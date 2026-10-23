"""Does combining several short-read signals per repeat beat the single best signal?

Training people only (validation/test untouched).  ~5,000 repeats sampled across core-spread
groups.  5-fold cross-validation over people: every estimate comes from a model fitted
without that person's fold.

  single   current method: isotonic + Theil-Sen curve on the best single signal (|r|)
  ridge_iso  per-repeat ridge regression on all 11 signals -> combined score -> isotonic + Theil-Sen
  top3     curves on the 3 best signals, averaged with weights 1/resid_sd^2
  ridge    ridge regression alone (reference)

    python scripts/curve_combination_experiment.py [--n-loci 5000] [--jobs 24]
-> work/experiments/curve_combination/{per_row.parquet, summary.csv}
"""
import argparse
import hashlib
import sys
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from joblib import Parallel, delayed
from sklearn.linear_model import RidgeCV

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, "/data/GitHub_Repos/VNTR_Genotyping/src/vntrlen")
from vntrlen.curves import SIGNAL_SHORTLIST, best_fit, fit_one, predict_one  # noqa: E402

W = HERE / "work"
OUT = W / "experiments" / "curve_combination"
K, MIN_N = 5, 8
EDGES, LABELS = [1, 2, 5, 20, 100], ["0", "1", "2-4", "5-19", "20-99", ">=100"]


def fold(s: str) -> int:
    return int.from_bytes(hashlib.sha1(f"cc:{s}".encode()).digest()[:8], "big") % K


def ranked(sig, y):
    """Shortlist signals ranked by |train r| with their fitted curves."""
    out = []
    for c, x in sig.items():
        p = fit_one(x, y)
        if p is not None and np.isfinite(p["train_pearson"]):
            out.append((abs(p["train_pearson"]), c, p))
    return sorted(out, key=lambda t: -t[0])


def ridge_fit(X, y):
    med = np.nanmedian(X, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    Xf = np.where(np.isfinite(X), X, med)
    mu, sd = Xf.mean(0), Xf.std(0) + 1e-9
    m = RidgeCV(alphas=np.logspace(-2, 3, 12)).fit((Xf - mu) / sd, y)
    return lambda Z: m.predict((np.where(np.isfinite(Z), Z, med) - mu) / sd)


def locus_job(lid, X, y, folds):
    names = SIGNAL_SHORTLIST
    pred = {k: np.full(y.size, np.nan) for k in ("single", "ridge_iso", "top3", "ridge")}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for f in range(K):
            te, tr = folds == f, folds != f
            if te.sum() == 0 or np.isfinite(y[tr]).sum() < MIN_N:
                continue
            sig_tr = {c: X[tr, i] for i, c in enumerate(names)}
            sig_te = {c: X[te, i] for i, c in enumerate(names)}
            bf = best_fit(sig_tr, y[tr], names, MIN_N)
            if bf:
                pred["single"][te] = predict_one(bf[1], sig_te[bf[0]])
            rk = ranked(sig_tr, y[tr])[:3]
            if rk:
                P = np.array([predict_one(p, sig_te[c]) for _, c, p in rk])
                w = np.array([1 / max(p["resid_sd"], 1e-3) ** 2 for _, _, p in rk])[:, None] * np.isfinite(P)
                pred["top3"][te] = np.nansum(np.nan_to_num(P) * w, 0) / np.where(w.sum(0) > 0, w.sum(0), np.nan)
            if np.std(y[tr]) > 0:
                rf = ridge_fit(X[tr], y[tr])
                s_tr, s_te = rf(X[tr]), rf(X[te])
                pred["ridge"][te] = np.clip(s_te, 0, None)
                p = fit_one(s_tr, y[tr])
                pred["ridge_iso"][te] = predict_one(p, s_te) if p else pred["ridge"][te]
            else:                                  # truth constant in this fold's training people
                for k in ("ridge", "ridge_iso"):
                    pred[k][te] = y[tr][0]
    return lid, pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-loci", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=24)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["locus_id", "sample_id", "sample_fold", "combined_units", "cov_panel_core90_units"] + SIGNAL_SHORTLIST
    t = pl.read_parquet(W / "model" / "model_table.parquet", columns=cols).filter(pl.col("sample_fold") == "train")
    loci = t.group_by("locus_id").agg(core=pl.col("cov_panel_core90_units").first(), n=pl.len()) \
        .filter(pl.col("n") >= 50).with_columns(group=pl.col("core").cut(EDGES, labels=LABELS, left_closed=True))
    per = args.n_loci // len(LABELS)
    pick = pl.concat([g.sample(min(per, g.height), seed=20261023) for _, g in loci.group_by("group")])
    t = t.filter(pl.col("locus_id").is_in(pick["locus_id"].implode()))
    print(f"{pick.height:,} repeats; {t['sample_id'].n_unique()} training people", flush=True)
    jobs = []
    for (lid,), g in t.group_by(["locus_id"]):
        jobs.append(delayed(locus_job)(lid, g.select(SIGNAL_SHORTLIST).cast(pl.Float64).to_numpy(),
                                       g["combined_units"].to_numpy().astype(float),
                                       np.array([fold(s) for s in g["sample_id"]]),))
    res = dict(Parallel(n_jobs=args.jobs, batch_size=20)(jobs))
    rows = []
    for (lid,), g in t.group_by(["locus_id"]):
        p = res[lid]
        rows.append(g.select("locus_id", "sample_id", "combined_units").with_columns(
            [pl.Series(k, v) for k, v in p.items()]))
    R = pl.concat(rows).join(pick.select("locus_id", "group", "core"), on="locus_id")
    R.write_parquet(OUT / "per_row.parquet")

    methods = ["single", "ridge_iso", "top3", "ridge"]
    out = []
    for grp, g in [("all (sampled)", R)] + [(str(k), v) for (k,), v in R.group_by(["group"])]:
        for m in methods:
            e = (g[m] - g["combined_units"]).abs()
            d = g.with_columns(a=pl.col(m) - pl.col(m).mean().over("locus_id"),
                               b=pl.col("combined_units") - pl.col("combined_units").mean().over("locus_id")) \
                .filter(pl.col("a").is_finite())
            r = float(np.corrcoef(d["a"], d["b"])[0, 1]) if d.height > 2 and d["a"].std() > 0 else np.nan
            out.append(dict(group=grp, method=m, n_loci=g["locus_id"].n_unique(), mean_abs_err=float(e.mean()),
                            median_abs_err=float(e.median()), within_1=float((e <= 1).mean()),
                            within_locus_r=r, missing=float(g[m].is_null().mean() + g[m].is_nan().mean())))
    S = pl.DataFrame(out)
    S.write_csv(OUT / "summary.csv")
    order = {"all (sampled)": 0, **{l: i + 1 for i, l in enumerate(LABELS)}}
    S = S.with_columns(o=pl.col("group").replace_strict(order, return_dtype=pl.Int64)).sort("o", "method").drop("o")
    pl.Config.set_tbl_rows(40)
    pl.Config.set_tbl_width_chars(160)
    print(S.with_columns(pl.col(pl.Float64).round(4)))


if __name__ == "__main__":
    main()
