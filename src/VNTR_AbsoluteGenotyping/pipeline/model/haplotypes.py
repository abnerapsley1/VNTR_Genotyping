"""
haplotypes.py — per-haplotype (allele) lengths from the combined prediction.

Read depth measures the SUM of a person's two alleles; the split has to come
from other information.  Two sources are combined per (locus, person):

  1. Allele catalog (known loci only).  The panel (vamos people minus held-out
     val/test people) gives each locus's observed haplotype alleles and their
     frequencies; a candidate allele pair (a, b) gets the Hardy-Weinberg prior
     p(a) p(b) (x2 if a != b).
  2. Allele-gap evidence.  A classifier predicts the bin of |a - b| from read
     signals that see the two alleles separately (spanning-pair insert-size
     spread / bimodality / 2-means split) plus locus context.  It is trained
     on training people, with cross-fitted intermediate predictions (plc_pred,
     transfer_pred) standing in for the in-sample final prediction.

  Candidate score (known loci):
     log prior(a, b) + log N(a + b; y_pred, sd) + w * log P(gap bin of b - a)
     sd = max(sd_floor, k * (y_hi - y_lo) / 3.29), from the conformal interval
  plus an equal-halves candidate (prior eps) so alleles absent from the
  catalog never force a catalogued pair.  Posterior over candidates gives the
  call, P(two alleles differ) and a 90% range per allele.

  Novel loci (no catalog may be used): total = round(y_pred), gap = the typical
  gap of the most likely gap bin (novel-safe classifier), het if
  P(gap = 0) < tau.

k, w, eps, tau are tuned on validation people; scored on test people.

Output: work/model/haplotype_predictions.parquet, work/model/haplotype_eval.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import time

import numpy as np
import polars as pl

from ..paths import load_config
from ..run.manifest import load_samples
from .assemble import _clean_vamos_col
from .calibration_transfer import novel_safe_features

GAP_EDGES = [0, 1, 2, 3, 5, 10, 20]           # gap bins: 0,1,2,3-4,5-9,10-19,>=20 units
N_BINS = len(GAP_EDGES)
TOP_ALLELES = 40                               # per locus; covers nearly all panel haplotypes
KNOWN_EXTRA = ["cov_h_exp", "cov_h_obs", "cov_n_alleles", "cov_n_distinct_hap_units",
               "cov_panel_core90_units", "cov_panel_cn_range_units", "cov_mean_combined_units",
               "cov_panel_median_units", "cov_panel_n_outliers"]
GRID = {"k": [0.5, 0.75, 1.0, 1.5], "w": [0.0, 0.5, 1.0], "eps": [1e-4, 1e-2]}


def gap_bin(g):
    return np.searchsorted(GAP_EDGES, np.asarray(g), side="right") - 1


# ---------------------------------------------------------------- data
def load(cfg):
    mdir = cfg.path("outputs.model_dir")
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg) if s.get("vamos_column")}
    sp = pl.read_parquet(mdir / "sample_split.parquet")
    held = [smap[s] for s in sp.filter(pl.col("fold") != "train")["sample_id"] if s in smap]
    gt = pl.read_parquet(cfg.path("outputs.ground_truth_dir") / "gt_long.parquet",
                         columns=["locus_id", "vamos_sample", "a1_units", "a2_units",
                                  "hap_lo_units", "hap_hi_units", "n_haps_called"])
    pred = pl.read_parquet(mdir / "predictions.parquet",
                           columns=["locus_id", "sample_id", "quadrant", "path", "y_pred", "y_lo", "y_hi",
                                    "plc_pred", "transfer_pred", "period"])
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    feats_novel = [c for c in novel_safe_features(tbl) if c != "period"] + ["transfer_pred"]
    feats_known = feats_novel + [c for c in KNOWN_EXTRA if c in tbl.columns] + ["plc_pred"]
    keep = sorted(set(feats_known) - {"transfer_pred", "plc_pred"})
    truth = (gt.filter(pl.col("n_haps_called") == 2)
             .join(pl.DataFrame({"sample_id": list(smap), "vamos_sample": list(smap.values())}), on="vamos_sample")
             .select("locus_id", "sample_id", "hap_lo_units", "hap_hi_units"))
    d = (pred.join(tbl.select(["locus_id", "sample_id"] + keep), on=["locus_id", "sample_id"], how="left")
         .join(truth, on=["locus_id", "sample_id"], how="left")
         .with_row_index("rid"))
    return d, gt.filter(~pl.col("vamos_sample").is_in(held)), feats_known, feats_novel


def catalog_pairs(panel: pl.DataFrame) -> pl.DataFrame:
    haps = pl.concat([panel.select("locus_id", u=pl.col("a1_units")),
                      panel.select("locus_id", u=pl.col("a2_units"))]).drop_nulls()
    freq = (haps.group_by("locus_id", "u").agg(c=pl.len())
            .with_columns(f=pl.col("c") / pl.col("c").sum().over("locus_id"))
            .sort(["locus_id", "f"], descending=[False, True])
            .with_columns(rk=pl.col("f").cum_count().over("locus_id"))
            .filter(pl.col("rk") <= TOP_ALLELES))
    return (freq.join(freq, on="locus_id", suffix="_b").filter(pl.col("u") <= pl.col("u_b"))
            .select("locus_id", a=pl.col("u"), b=pl.col("u_b"),
                    lp=(pl.col("f") * pl.col("f_b")
                        * pl.when(pl.col("u") != pl.col("u_b")).then(2.0).otherwise(1.0)).log())
            .with_columns(s=pl.col("a") + pl.col("b")))


# ---------------------------------------------------------------- gap model
def _classifier(seed):
    try:
        import xgboost as xgb
        return xgb.XGBClassifier(n_estimators=400, learning_rate=0.08, grow_policy="lossguide", max_leaves=63,
                                 max_depth=0, tree_method="hist", device="cuda", random_state=seed,
                                 objective="multi:softprob", eval_metric="mlogloss")
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        return HistGradientBoostingClassifier(max_iter=400, learning_rate=0.08, max_leaf_nodes=63,
                                              early_stopping=True, random_state=seed)


def gap_model(d: pl.DataFrame, feats: list[str], seed: int, n_train: int = 3_000_000):
    tr = d.filter((pl.col("quadrant") == "train") & pl.col("hap_lo_units").is_not_null())
    tr = tr.sample(min(n_train, tr.height), seed=seed)
    y = gap_bin(tr["hap_hi_units"].to_numpy() - tr["hap_lo_units"].to_numpy())
    m = _classifier(seed).fit(tr.select(feats).to_numpy().astype(np.float32), y)
    # typical gap per bin (novel path)
    g = tr["hap_hi_units"].to_numpy() - tr["hap_lo_units"].to_numpy()
    rep = np.array([float(np.median(g[y == b])) if (y == b).any() else GAP_EDGES[b] for b in range(N_BINS)])
    return m, rep


def gap_proba(m, d: pl.DataFrame, feats: list[str]) -> np.ndarray:
    out = np.empty((d.height, N_BINS), dtype=np.float32)
    for i in range(0, d.height, 2_000_000):
        out[i:i + 2_000_000] = m.predict_proba(d[i:i + 2_000_000].select(feats).to_numpy().astype(np.float32))
    return np.clip(out, 1e-6, 1)


# ---------------------------------------------------------------- known path
def score_known(rows: pl.DataFrame, pairs: pl.DataFrame, k: float, w: float, eps: float,
                sd_floor: float = 0.5, chunk: int = 300_000) -> pl.DataFrame:
    """rows: rid, locus_id, y_pred, y_lo, y_hi, gp0..gp6 -> call + posterior summaries per rid."""
    outs = []
    for c in rows.iter_slices(chunk):
        c = c.with_columns(sd=pl.max_horizontal(k * (pl.col("y_hi") - pl.col("y_lo")) / 3.29, pl.lit(sd_floor)))
        half_s = pl.col("y_pred").round().cast(pl.Int64)
        halves = c.select("rid", "locus_id", "y_pred", "sd", *[f"gp{b}" for b in range(N_BINS)],
                          a=half_s // 2, b=half_s - half_s // 2, lp=pl.lit(float(np.log(eps)))
                          ).with_columns(s=pl.col("a") + pl.col("b"))
        cand = (c.select("rid", "locus_id", "y_pred", "sd", *[f"gp{b}" for b in range(N_BINS)])
                .join(pairs, on="locus_id")
                .filter((pl.col("s") - pl.col("y_pred")).abs() <= 3 * pl.col("sd") + 0.5))
        cand = pl.concat([cand.select(halves.columns), halves])
        gb = pl.col("b") - pl.col("a")
        binx = pl.sum_horizontal([(gb >= e).cast(pl.Int32) for e in GAP_EDGES[1:]])
        cand = cand.with_columns(
            lgap=pl.concat_list([f"gp{b}" for b in range(N_BINS)]).list.get(binx).log())
        cand = cand.with_columns(
            score=pl.col("lp") - 0.5 * ((pl.col("s") - pl.col("y_pred")) / pl.col("sd")) ** 2 + w * pl.col("lgap"))
        cand = cand.with_columns(p=(pl.col("score") - pl.col("score").max().over("rid")).exp())
        cand = cand.with_columns(p=pl.col("p") / pl.col("p").sum().over("rid"))
        # collapse duplicate (a, b) (catalog pair == halves)
        cand = cand.group_by("rid", "a", "b").agg(pl.col("p").sum())
        top = (cand.sort("p", descending=True).group_by("rid").first()
               .select("rid", h_lo="a", h_hi="b", p_top="p"))
        het = cand.group_by("rid").agg(p_het=pl.col("p").filter(pl.col("a") != pl.col("b")).sum())
        q = []
        for col, nm in (("a", "lo"), ("b", "hi")):
            mm = (cand.group_by("rid", col).agg(pl.col("p").sum()).sort("rid", col)
                  .with_columns(cp=pl.col("p").cum_sum().over("rid")))
            q.append(mm.group_by("rid").agg(
                pl.col(col).filter(pl.col("cp") >= 0.05).first().alias(f"{nm}_q05"),
                pl.col(col).filter(pl.col("cp") >= 0.95).first().alias(f"{nm}_q95")))
        outs.append(top.join(het, on="rid").join(q[0], on="rid").join(q[1], on="rid"))
    return pl.concat(outs)


# ---------------------------------------------------------------- novel path
def call_novel(rows: pl.DataFrame, rep: np.ndarray, tau: float) -> pl.DataFrame:
    P = rows.select([f"gp{b}" for b in range(N_BINS)]).to_numpy()
    s = np.clip(np.round(rows["y_pred"].to_numpy()), 0, None)
    b = np.argmax(P[:, 1:], axis=1) + 1
    het = P[:, 0] < tau
    g = np.where(het, np.round(rep[b]), 0.0)
    g = np.minimum(g, s)
    lo = np.floor((s - g) / 2)
    hi = s - lo
    return rows.select("rid").with_columns(h_lo=pl.Series(lo), h_hi=pl.Series(hi), p_het=pl.Series(1 - P[:, 0]),
                                          p_top=pl.lit(None, dtype=pl.Float64),
                                          lo_q05=pl.lit(None, dtype=pl.Float64), lo_q95=pl.lit(None, dtype=pl.Float64),
                                          hi_q05=pl.lit(None, dtype=pl.Float64), hi_q95=pl.lit(None, dtype=pl.Float64))


# ---------------------------------------------------------------- scoring
def metrics(d: pl.DataFrame, lo: str = "h_lo", hi: str = "h_hi") -> dict:
    d = d.filter(pl.col("hap_lo_units").is_not_null() & pl.col(lo).is_not_null())
    if d.height == 0:
        return {"n": 0}
    tl, th = d["hap_lo_units"].to_numpy(), d["hap_hi_units"].to_numpy()
    cl, ch = d[lo].to_numpy(), d[hi].to_numpy()
    e = np.maximum(np.abs(cl - tl), np.abs(ch - th))
    het, chet = th != tl, ch != cl
    out = {"n": d.height, "frac_true_het": float(het.mean()),
           "both_exact": float(np.mean(e == 0)), "both_within_1": float(np.mean(e <= 1)),
           "both_within_2": float(np.mean(e <= 2)),
           "het_both_exact": float(np.mean(e[het] == 0)) if het.any() else None,
           "het_both_within_1": float(np.mean(e[het] <= 1)) if het.any() else None,
           "hom_both_exact": float(np.mean(e[~het] == 0)) if (~het).any() else None,
           "het_recall": float(np.mean(chet[het])) if het.any() else None,
           "het_precision": float(np.mean(het[chet])) if chet.any() else None,
           "mean_abs_err_lo": float(np.mean(np.abs(cl - tl))), "mean_abs_err_hi": float(np.mean(np.abs(ch - th))),
           "mean_abs_err_gap": float(np.mean(np.abs((ch - cl) - (th - tl))))}
    if "lo_q05" in d.columns and d["lo_q05"].is_not_null().any():
        k = d.filter(pl.col("lo_q05").is_not_null())
        cov = ((k["hap_lo_units"] >= k["lo_q05"]) & (k["hap_lo_units"] <= k["lo_q95"])
               & (k["hap_hi_units"] >= k["hi_q05"]) & (k["hap_hi_units"] <= k["hi_q95"]))
        out["both_in_90pct_range"] = float(cov.mean())
    if "p_het" in d.columns:
        ph = d["p_het"].to_numpy()
        for t in (0.5, 0.8, 0.95):
            sel = ph >= t
            out[f"p_het>={t}_share"] = float(sel.mean())
            out[f"p_het>={t}_truly_het"] = float(het[sel].mean()) if sel.any() else None
    return out


def build(cfg) -> None:
    t0 = time.time()
    mdir = cfg.path("outputs.model_dir")
    seed = int(cfg.get("model.random_seed", 20260831))
    d, panel, fk, fn = load(cfg)
    print(f"[{time.time()-t0:.0f}s] {d.height:,} person-sites; known features {len(fk)}, novel-safe {len(fn)}", flush=True)
    pairs = catalog_pairs(panel)
    print(f"[{time.time()-t0:.0f}s] catalog: {pairs.height:,} allele pairs", flush=True)

    known = d["path"] == "known"
    mk, _ = gap_model(d.filter(known), fk, seed)
    mn, rep = gap_model(d, fn, seed + 1)
    print(f"[{time.time()-t0:.0f}s] gap models fit; typical gap per bin {np.round(rep, 1).tolist()}", flush=True)
    gp = np.where(known.to_numpy()[:, None], gap_proba(mk, d, fk), gap_proba(mn, d, fn))
    d = d.with_columns([pl.Series(f"gp{b}", gp[:, b]) for b in range(N_BINS)])
    base_cols = ["rid", "locus_id", "y_pred", "y_lo", "y_hi"] + [f"gp{b}" for b in range(N_BINS)]

    # ---- tune known-path k, w, eps on validation people (subsample)
    val_k = d.filter((pl.col("quadrant") == "val") & known & pl.col("hap_lo_units").is_not_null())
    val_k = val_k.sample(min(400_000, val_k.height), seed=seed)
    best, tune_log = None, []
    for k, w, eps in itertools.product(GRID["k"], GRID["w"], GRID["eps"]):
        r = val_k.select("rid", "hap_lo_units", "hap_hi_units").join(
            score_known(val_k.select(base_cols), pairs, k, w, eps), on="rid")
        m = metrics(r)
        tune_log.append({"k": k, "w": w, "eps": eps, "both_exact": m["both_exact"], "het_both_exact": m["het_both_exact"]})
        print(f"  k={k} w={w} eps={eps}: both exact {m['both_exact']:.4f} (het {m['het_both_exact']:.3f})", flush=True)
        if best is None or m["both_exact"] > best[0]:
            best = (m["both_exact"], k, w, eps)
    _, K, Wt, EPS = best
    # ---- tune novel-path tau on validation people at novel loci
    val_n = d.filter((pl.col("quadrant") == "val") & ~known & pl.col("hap_lo_units").is_not_null())
    taus = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    res = {t: metrics(val_n.select("rid", "hap_lo_units", "hap_hi_units").join(call_novel(val_n, rep, t), on="rid"))["both_exact"]
           for t in taus} if val_n.height else {0.5: 0}
    TAU = max(res, key=res.get)
    print(f"[{time.time()-t0:.0f}s] tuned: k={K} w={Wt} eps={EPS}; novel tau={TAU} ({res})", flush=True)

    # ---- apply to everyone
    cols = ["rid", "h_lo", "h_hi", "p_top", "p_het", "lo_q05", "lo_q95", "hi_q05", "hi_q95"]
    calls = pl.concat([score_known(d.filter(known).select(base_cols), pairs, K, Wt, EPS).select(cols),
                       call_novel(d.filter(~known), rep, TAU).select(cols)], how="vertical_relaxed")
    out = d.select("rid", "locus_id", "sample_id", "quadrant", "path", "y_pred", "y_lo", "y_hi", "period",
                   "hap_lo_units", "hap_hi_units", *[f"gp{b}" for b in range(N_BINS)]).join(calls, on="rid", how="left").with_columns(
        half_lo=(pl.col("y_pred").round() / 2).floor(), half_hi=pl.col("y_pred").round() - (pl.col("y_pred").round() / 2).floor())
    out.drop("hap_lo_units", "hap_hi_units").write_parquet(mdir / "haplotype_predictions.parquet")

    report = {"params": {"k": K, "w": Wt, "eps": EPS, "tau": TAU, "top_alleles": TOP_ALLELES,
                         "gap_edges": GAP_EDGES, "typical_gap_per_bin": rep.tolist()},
              "tuning": tune_log, "quadrants": {}}
    for q in ["val", "test_novel_sample", "test_novel_locus", "test_novel_both"]:
        s = out.filter(pl.col("quadrant") == q)
        report["quadrants"][q] = {"model": metrics(s), "equal_halves": metrics(s, "half_lo", "half_hi")}
        a, b = report["quadrants"][q]["model"], report["quadrants"][q]["equal_halves"]
        if a.get("n"):
            print(f"{q:18s} both exact {a['both_exact']:.3f} (halves {b['both_exact']:.3f}) | "
                  f"het both exact {a['het_both_exact']:.3f} | het recall {a['het_recall']:.3f} "
                  f"precision {a['het_precision'] or float('nan'):.3f} | within 1 {a['both_within_1']:.3f} "
                  f"(halves {b['both_within_1']:.3f})", flush=True)
    (mdir / "haplotype_eval.json").write_text(json.dumps(report, indent=2, default=float))
    print(f"[{time.time()-t0:.0f}s] wrote {mdir/'haplotype_predictions.parquet'}, {mdir/'haplotype_eval.json'}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
