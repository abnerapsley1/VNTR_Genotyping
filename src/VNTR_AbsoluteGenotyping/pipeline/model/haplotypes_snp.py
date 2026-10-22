"""
haplotypes_snp.py — per-haplotype lengths using flanking phased SNP haplotypes.

A VNTR allele is usually inherited together with the SNPs around it, so a
person's two phased SNP haplotypes (1kGP 2022 phased panel) each predict the
repeat allele they carry.  This replaces the allele-frequency prior of
haplotypes.py with a haplotype-specific one, and adds phase: which SNP
haplotype carries the longer allele.

Per known locus (held-out loci are excluded by design: no panel there):
  SNPs     nearest N_FLANK SNPs on each side of the array (within 100 kb)
  match    for two haplotypes, the run of identical SNPs extending outward from
           the array on the left plus on the right (an IBD-length proxy)
  labels   training people's allele pair {a, b} is known but not which SNP
           haplotype carries which; EM over training people resolves it
           (leave-person-out kNN predictions of each haplotype's allele)
  predict  each held-out haplotype: alleles of its K best-matching training
           haplotypes (soft labels), smoothed with the panel allele frequencies
           (weight alpha)
  pair     prior(a, b) = P1(a) P2(b) + P1(b) P2(a); combined with the total
           likelihood and allele-gap evidence exactly as in haplotypes.py
           (k, w, eps taken from its tuning); K and alpha tuned on validation.
  reads    optional spanning-read term (features/spanning_alleles.py): read
           array lengths, shifted by a per-locus offset learned on training
           people (read length - vamos bp), scored against each candidate's
           allele lengths in bp; weight tuned on validation.

Reported variants (test people): equal halves, catalog (stage 1), SNP
haplotypes, catalog + spanning reads, SNP haplotypes + spanning reads.

Output: work/model/haplotype_snp_predictions.parquet, work/model/haplotype_snp_eval.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import polars as pl
from joblib import Parallel, delayed

from ..paths import load_config
from ..run.manifest import load_samples
from .assemble import _clean_vamos_col
from .haplotypes import GAP_EDGES, N_BINS, TOP_ALLELES, metrics

N_FLANK = 256
EM_ITERS = 4
BCF = "/data/miniconda3/envs/vntr_benchmark/bin/bcftools"


def load_snps(bcf: Path, samples: list[str]):
    """-> positions (n,), haplotypes (2 * len(samples), n) int8, in `samples` order."""
    hdr = subprocess.run([BCF, "query", "-l", str(bcf)], capture_output=True, text=True, check=True).stdout.split()
    idx = [hdr.index(s) for s in samples]
    # fast path: every genotype is phased diploid "a|b" -> fixed-width lines
    q = subprocess.run([BCF, "query", "-f", "%POS[\t%GT]\n", str(bcf)], capture_output=True, check=True).stdout
    lines = q.split(b"\n")[:-1]
    pos = np.array([int(ln.split(b"\t", 1)[0]) for ln in lines], dtype=np.int64)
    S = len(hdr)
    if all(len(ln) - ln.index(b"\t") == 4 * S for ln in lines):
        # each sample is "\tA|B": keep the leading tab so every field is 4 bytes
        body = np.frombuffer(b"".join(ln[ln.index(b"\t"):] for ln in lines), dtype=np.uint8).reshape(len(lines), S, 4)
        h1 = np.where(body[:, :, 1] == ord("1"), 1, np.where(body[:, :, 1] == ord("0"), 0, -1))
        h2 = np.where(body[:, :, 3] == ord("1"), 1, np.where(body[:, :, 3] == ord("0"), 0, -1))
        h = np.concatenate([h1[:, idx], h2[:, idx]], axis=1).astype(np.int8).T
    else:                                   # haploid calls present (chrX males): slow path
        rows = []
        for ln in lines:
            g = ln.decode().split("\t")[1:]
            rows.append([g[i][0] for i in idx] + [g[i][2] if len(g[i]) > 2 else g[i][0] for i in idx])
        a = np.array(rows, dtype="U1")
        h = np.where(a == "1", 1, np.where(a == "0", 0, -1)).astype(np.int8).T
    n = len(samples)
    order = np.empty(2 * n, dtype=int)
    order[0::2], order[1::2] = np.arange(n), np.arange(n) + n                  # person i -> rows 2i, 2i+1
    return np.asarray(pos), h[order]


def _match(X: np.ndarray, T: np.ndarray, nl: int) -> np.ndarray:
    """Outward identical-run length (left + right) between rows of X and rows of T."""
    eq = X[:, None, :] == T[None, :, :]
    left, right = eq[:, :, :nl], eq[:, :, nl:]
    return np.cumprod(left, axis=2).sum(2) + np.cumprod(right, axis=2).sum(2)


def _knn_dist(sim: np.ndarray, lab: np.ndarray, K: int) -> np.ndarray:
    """sim (x, t); lab (t, A) soft labels -> (x, A) allele distribution of the K best matches."""
    K = min(K, sim.shape[1])
    top = np.argpartition(-sim, K - 1, axis=1)[:, :K]
    return lab[top].mean(axis=1)


def read_ll(Lp: np.ndarray, A: np.ndarray, B: np.ndarray, eps_r: float = 0.05) -> np.ndarray:
    """Log-likelihood of offset-corrected spanning-read lengths Lp under each candidate pair (A, B) in bp:
    each read comes from either allele; exact length w.p. 1 - eps_r, else a two-sided geometric error."""
    def f(x):
        ax = np.abs(x)
        return np.where(ax == 0, 1 - eps_r, eps_r * 0.25 * 0.5 ** np.maximum(ax - 1, 0))
    L = Lp[None, :]
    return np.log(0.5 * f(L - A[:, None]) + 0.5 * f(L - B[:, None])).sum(1)


def locus_job(lid, H, nl, train_rows, train_alleles, held, freq_u, freq_f, Ks, alphas, k_sd, w, eps,
              wrs=(0.0,), offset=None, ubp=None):
    """H: (2S, m) haplotypes at this locus; train_rows: person indices with truth;
    train_alleles: (n_train, 2) unordered; held: list of dicts (person index, y_pred, y_lo, y_hi, gp*,
    optional 'lengths' = spanning-read array lengths).  offset: this locus's read-length offset
    (read length - vamos bp); ubp: unit allele -> bp.  Returns rows per held-out person x (K, alpha, wr)."""
    A = len(freq_u)
    ai = {u: i for i, u in enumerate(freq_u)}
    tr_h = np.concatenate([[2 * i, 2 * i + 1] for i in train_rows]) if len(train_rows) else np.array([], int)
    ta = np.array([[ai.get(a, -1), ai.get(b, -1)] for a, b in train_alleles])
    ok = (ta >= 0).all(1)
    tr_rows_ok = np.asarray(train_rows)[ok]
    ta = ta[ok]
    tr_h = np.concatenate([[2 * i, 2 * i + 1] for i in tr_rows_ok]) if len(tr_rows_ok) else np.array([], int)
    out = []
    if len(tr_h) < 20 or H.shape[1] < 4:
        return out
    T = H[tr_h]
    sim_tt = _match(T, T, nl)
    n_t = len(tr_rows_ok)
    own = np.repeat(np.arange(n_t), 2)
    sim_tt[own[:, None] == own[None, :]] = -1                           # leave-person-out
    q = np.full(n_t, 0.5)                                               # P(hap1 carries allele a)
    het = ta[:, 0] != ta[:, 1]
    K0 = Ks[len(Ks) // 2]
    for _ in range(EM_ITERS):
        lab = np.zeros((2 * n_t, A))
        lab[0::2][np.arange(n_t), ta[:, 0]] += q
        lab[0::2][np.arange(n_t), ta[:, 1]] += 1 - q
        lab[1::2][np.arange(n_t), ta[:, 0]] += 1 - q
        lab[1::2][np.arange(n_t), ta[:, 1]] += q
        P = 0.9 * _knn_dist(sim_tt, lab, K0) + 0.1 * freq_f
        p1a, p2b = P[0::2][np.arange(n_t), ta[:, 0]], P[1::2][np.arange(n_t), ta[:, 1]]
        p1b, p2a = P[0::2][np.arange(n_t), ta[:, 1]], P[1::2][np.arange(n_t), ta[:, 0]]
        q = np.where(het, p1a * p2b / np.clip(p1a * p2b + p1b * p2a, 1e-12, None), 0.5)
    lab = np.zeros((2 * n_t, A))
    lab[0::2][np.arange(n_t), ta[:, 0]] += q
    lab[0::2][np.arange(n_t), ta[:, 1]] += 1 - q
    lab[1::2][np.arange(n_t), ta[:, 0]] += 1 - q
    lab[1::2][np.arange(n_t), ta[:, 1]] += q
    phase_conf = float(np.mean(np.maximum(q, 1 - q)[het])) if het.any() else np.nan

    hp = np.array([r["person"] for r in held])
    X = H[np.concatenate([[2 * i, 2 * i + 1] for i in hp])]
    sim_xt = _match(X, T, nl)
    best_match = sim_xt.max(1)
    av, bv = np.triu_indices(A)                      # index pairs; alleles are in frequency order
    sv = freq_u[av] + freq_u[bv]
    gapb = np.searchsorted(GAP_EDGES, np.abs(freq_u[bv] - freq_u[av]), side="right") - 1
    if ubp:
        bpu = float(np.median([b / max(u, 1) for u, b in ubp.items()]))
    to_bp = (lambda us: np.array([ubp.get(int(u), u * bpu) for u in us], dtype=float)) if ubp else None
    for K, alpha, wr in itertools.product(Ks, alphas, wrs):
            D = _knn_dist(sim_xt, lab, K)
            P = (1 - alpha) * D + alpha * freq_f
            for j, r in enumerate(held):
                P1, P2 = P[2 * j], P[2 * j + 1]
                pr = np.where(av == bv, P1[av] * P2[bv], P1[av] * P2[bv] + P1[bv] * P2[av])
                sd = max(0.5, k_sd * (r["y_hi"] - r["y_lo"]) / 3.29)
                gp = np.array([r[f"gp{b}"] for b in range(N_BINS)])
                sc = np.log(np.clip(pr, 1e-300, None)) - 0.5 * ((sv - r["y_pred"]) / sd) ** 2 + w * np.log(gp[gapb])
                keep = np.abs(sv - r["y_pred"]) <= 3 * sd + 0.5
                hs = int(round(r["y_pred"]))
                ha, hb = hs // 2, hs - hs // 2
                sc_h = np.log(eps) - 0.5 * ((hs - r["y_pred"]) / sd) ** 2 + w * np.log(gp[min(hb - ha, N_BINS - 1)])
                cand_a = np.append(np.minimum(freq_u[av], freq_u[bv])[keep], ha)   # shorter allele
                cand_b = np.append(np.maximum(freq_u[av], freq_u[bv])[keep], hb)   # longer allele
                cs = np.append(sc[keep], sc_h)
                Lr = r.get("lengths")
                n_reads = 0 if Lr is None else len(Lr)
                if wr and n_reads and offset is not None and to_bp is not None:
                    Lp = np.asarray(Lr[:40], dtype=float) - offset
                    cs = cs + wr * read_ll(Lp, to_bp(cand_a), to_bp(cand_b))
                p = np.exp(cs - cs.max())
                p /= p.sum()
                t = int(np.argmax(p))
                lo_, hi_ = int(cand_a[t]), int(cand_b[t])
                # phase of the call: P(SNP hap 1 carries the longer allele)
                if lo_ != hi_ and lo_ in ai and hi_ in ai:
                    x1, x2 = P1[ai[hi_]] * P2[ai[lo_]], P1[ai[lo_]] * P2[ai[hi_]]
                    ph1 = x1 / max(x1 + x2, 1e-300)
                else:
                    ph1 = np.nan
                out.append((lid, r["sample_id"], K, alpha, wr, lo_, hi_, float(p[t]),
                            float(p[cand_a != cand_b].sum()), ph1,
                            float(min(best_match[2 * j], best_match[2 * j + 1])), phase_conf, n_reads))
    return out


def span_inputs(cfg, gt_bp: pl.DataFrame, smap: dict, train: list[str]):
    """Spanning-read lengths per (locus, person), per-locus read offset (read length - vamos bp,
    median over training people whose two alleles are equal), and per-locus unit -> bp map."""
    sdir = cfg.path("outputs.features_dir") / "spanning"
    files = {f.stem: f for f in sdir.glob("*.parquet")} if sdir.exists() else {}
    if not files:
        return {}, {}, {}
    sp = pl.concat([pl.read_parquet(f).with_columns(sample_id=pl.lit(s)) for s, f in files.items()])
    sp = sp.filter(pl.col("n_reads") > 0)
    lengths = {(r[0], r[1]): r[2] for r in sp.select("locus_id", "sample_id", "lengths").iter_rows()}
    tr = (sp.filter(pl.col("sample_id").is_in(train) & (pl.col("n_reads") >= 5))
          .join(pl.DataFrame({"sample_id": list(smap), "vamos_sample": list(smap.values())}), on="sample_id")
          .join(gt_bp.filter(pl.col("a1_bp") == pl.col("a2_bp")), on=["locus_id", "vamos_sample"])
          .with_columns(mode=pl.col("lengths").list.eval(pl.element().mode().first()).list.first())
          .with_columns(off=pl.col("mode") - pl.col("a1_bp")))
    offs = tr.group_by("locus_id").agg(off=pl.col("off").median(), n=pl.len()).filter(pl.col("n") >= 3)
    offset = dict(zip(offs["locus_id"], offs["off"].round().cast(pl.Int64)))
    ub = pl.concat([gt_bp.select("locus_id", u=pl.col("a1_units"), b=pl.col("a1_bp")),
                    gt_bp.select("locus_id", u=pl.col("a2_units"), b=pl.col("a2_bp"))]).drop_nulls()
    ub = ub.group_by("locus_id", "u").agg(pl.col("b").mode().first())
    ubp = {}
    for lid, u, b in ub.iter_rows():
        ubp.setdefault(lid, {})[int(u)] = int(b)
    print(f"  spanning reads: {len(files)} people; read offsets for {len(offset):,} loci "
          f"(offset 0 at {sum(v == 0 for v in offset.values()):,})", flush=True)
    return lengths, offset, ubp


def run(cfg, loci_filter: set | None, Ks, alphas, wrs, people: list[str]):
    mdir = cfg.path("outputs.model_dir")
    hp = json.loads((mdir / "haplotype_eval.json").read_text())["params"]
    k_sd, w, eps = hp["k"], hp["w"], hp["eps"]
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg) if s.get("vamos_column")}
    sp = pl.read_parquet(mdir / "sample_split.parquet")
    samples = sp["sample_id"].to_list()
    pidx = {s: i for i, s in enumerate(samples)}
    train = sp.filter(pl.col("fold") == "train")["sample_id"].to_list()
    held_vamos = [smap[s] for s in sp.filter(pl.col("fold") != "train")["sample_id"] if s in smap]
    gt = pl.read_parquet(cfg.path("outputs.ground_truth_dir") / "gt_long.parquet",
                         columns=["locus_id", "vamos_sample", "a1_units", "a2_units", "a1_bp", "a2_bp", "n_haps_called"])
    panel = gt.filter(~pl.col("vamos_sample").is_in(held_vamos))
    haps = pl.concat([panel.select("locus_id", u=pl.col("a1_units")), panel.select("locus_id", u=pl.col("a2_units"))]).drop_nulls()
    freq = (haps.group_by("locus_id", "u").agg(c=pl.len()).with_columns(f=pl.col("c") / pl.col("c").sum().over("locus_id"))
            .sort(["locus_id", "f"], descending=[False, True]).with_columns(rk=pl.col("f").cum_count().over("locus_id"))
            .filter(pl.col("rk") <= TOP_ALLELES).with_columns(f=pl.col("f") / pl.col("f").sum().over("locus_id")))
    v2s = {v: s for s, v in smap.items()}
    tr_truth = gt.filter((pl.col("n_haps_called") == 2)
                         & pl.col("vamos_sample").is_in([smap[s] for s in train if s in smap]))
    truth = (gt.filter(pl.col("n_haps_called") == 2)
             .join(pl.DataFrame({"sample_id": list(smap), "vamos_sample": list(smap.values())}), on="vamos_sample")
             .select("locus_id", "sample_id",
                     hap_lo_units=pl.min_horizontal("a1_units", "a2_units"),
                     hap_hi_units=pl.max_horizontal("a1_units", "a2_units")))
    lengths, offset, ubp = span_inputs(cfg, panel.select("locus_id", "vamos_sample", "a1_units", "a2_units", "a1_bp", "a2_bp")
                                       .filter(pl.col("vamos_sample").is_in([smap[s] for s in train if s in smap])
                                               | ~pl.col("vamos_sample").is_in(list(smap.values()))),
                                       smap, train) if any(wr > 0 for wr in wrs) else ({}, {}, {})
    hpred = pl.read_parquet(mdir / "haplotype_predictions.parquet").filter(
        (pl.col("path") == "known") & pl.col("sample_id").is_in(people)).join(
        truth, on=["locus_id", "sample_id"], how="left")
    idx = pl.read_csv(Path(cfg.repo_root) / "src" / "VNTR_Benchmark" / "work" / "danbing" / "tr_index.tsv", separator="\t")
    if loci_filter is not None:
        idx = idx.filter(pl.col("locus_id").is_in(list(loci_filter)))
    idx = idx.filter(pl.col("locus_id").is_in(hpred["locus_id"].unique().to_list()))
    fq = {lid: (g["u"].to_numpy(), g["f"].to_numpy()) for (lid,), g in freq.group_by(["locus_id"])}
    trd = {lid: g for (lid,), g in tr_truth.group_by(["locus_id"])}
    hd = {lid: g for (lid,), g in hpred.group_by(["locus_id"])}
    snpdir = Path(cfg.repo_root) / "data" / "1000Genomes" / "phased_snps"
    rows = []
    for (chrom,), g in idx.group_by(["chrom"], maintain_order=True):
        bcf = snpdir / f"{chrom}.bcf"
        if not Path(f"{bcf}.csi").exists():
            print(f"  {chrom}: no SNP file, skipped", flush=True)
            continue
        t0 = time.time()
        pos, H = load_snps(bcf, samples)
        jobs = []
        for r in g.iter_rows(named=True):
            lid = r["locus_id"]
            if lid not in fq or lid not in hd or lid not in trd:
                continue
            L = np.searchsorted(pos, r["array_start0"])
            R = np.searchsorted(pos, r["array_end0"], side="right")
            li = np.arange(max(0, L - N_FLANK), L)[::-1]
            li = li[pos[li] >= r["array_start0"] - 100_000]
            ri = np.arange(R, min(len(pos), R + N_FLANK))
            ri = ri[pos[ri] <= r["array_end0"] + 100_000]
            nl = len(li)
            Hl = np.concatenate([H[:, li], H[:, ri]], axis=1)
            t = trd[lid]
            tr_rows = [pidx[v2s[v]] for v in t["vamos_sample"]]
            tr_al = [tuple(sorted(x)) for x in zip(t["a1_units"], t["a2_units"])]
            held = [dict(rr, person=pidx[rr["sample_id"]], lengths=lengths.get((lid, rr["sample_id"])))
                    for rr in hd[lid].iter_rows(named=True)]
            u, f = fq[lid]
            jobs.append(delayed(locus_job)(lid, Hl, nl, tr_rows, tr_al, held, u, f, Ks, alphas, k_sd, w, eps,
                                           wrs, offset.get(lid), ubp.get(lid)))
        for res in Parallel(n_jobs=int(cfg.get("model.n_jobs", 16)), batch_size=16)(jobs):
            rows.extend(res)
        print(f"  {chrom}: {len(jobs):,} loci, {time.time()-t0:.0f}s", flush=True)
        del H
    cols = ["locus_id", "sample_id", "K", "alpha", "wr", "s_lo", "s_hi", "s_p_top", "s_p_het", "s_phase_hap1_long",
            "s_min_match", "s_train_phase_conf", "s_n_span_reads"]
    return pl.DataFrame(rows, schema=cols, orient="row"), hpred


def build(cfg, tune_loci: int = 6000, reads: bool = True) -> None:
    mdir = cfg.path("outputs.model_dir")
    sp = pl.read_parquet(mdir / "sample_split.parquet")
    val = sp.filter(pl.col("fold") == "val")["sample_id"].to_list()
    test = sp.filter(pl.col("fold") == "test")["sample_id"].to_list()
    Ks, alphas, wrs = [1, 2, 3, 5], [0.002, 0.005, 0.01, 0.02, 1.0], ([0.0, 0.5, 1.0] if reads else [0.0])
    loci = pl.read_parquet(mdir / "haplotype_predictions.parquet", columns=["locus_id", "path"]).filter(
        pl.col("path") == "known")["locus_id"].unique().sort()
    rng = np.random.default_rng(int(cfg.get("model.random_seed", 20260831)))
    tl = set(rng.choice(loci.to_numpy(), size=min(tune_loci, len(loci)), replace=False))
    print(f"tuning K x alpha x read weight on validation people at {len(tl):,} loci", flush=True)
    res, hpred = run(cfg, tl, Ks, alphas, wrs, val)
    tune = []
    for K, a, wr in itertools.product(Ks, alphas, wrs):
        r = hpred.join(res.filter((pl.col("K") == K) & (pl.col("alpha") == a) & (pl.col("wr") == wr)),
                       on=["locus_id", "sample_id"])
        m_s, m_c = metrics(r, "s_lo", "s_hi"), metrics(r)
        tune.append({"K": K, "alpha": a, "wr": wr, "both_exact": m_s["both_exact"], "catalog_both_exact": m_c["both_exact"],
                     "het_both_exact": m_s["het_both_exact"], "n": m_s["n"]})
        print(f"  K={K} alpha={a} wr={wr}: {m_s['both_exact']:.4f} (het {m_s['het_both_exact']:.3f}) vs stage 1 "
              f"{m_c['both_exact']:.4f}", flush=True)
    snp_only = max([t for t in tune if t["alpha"] < 1 and t["wr"] == 0], key=lambda x: x["both_exact"])
    K, a = snp_only["K"], snp_only["alpha"]
    wr_best = max([t for t in tune if t["K"] == K and t["alpha"] == a and t["wr"] > 0],
                  key=lambda x: x["both_exact"])["wr"] if reads else 0.0
    print(f"chosen K={K} alpha={a} read weight={wr_best}; running all known loci for val + test people", flush=True)
    res, hpred = run(cfg, None, [K], sorted({a, 1.0}), sorted({0.0, wr_best}), val + test)
    keys = ["locus_id", "sample_id"]
    pick = {"SNP haplotypes": (a, 0.0), "catalog + spanning reads": (1.0, wr_best),
            "SNP haplotypes + spanning reads": (a, wr_best)} if reads else {"SNP haplotypes": (a, 0.0)}
    out = hpred
    for name, (aa, ww) in pick.items():
        tag = {"SNP haplotypes": "snp", "catalog + spanning reads": "cat_span", "SNP haplotypes + spanning reads": "full"}[name]
        sub = res.filter((pl.col("alpha") == aa) & (pl.col("wr") == ww)).drop("K", "alpha", "wr")
        sub = sub.rename({c: f"{tag}_{c[2:]}" for c in sub.columns if c.startswith("s_")})
        out = out.join(sub, on=keys, how="left")
    out.drop("hap_lo_units", "hap_hi_units").write_parquet(mdir / "haplotype_snp_predictions.parquet")
    report = {"params": {"K": K, "alpha": a, "read_weight": wr_best, "n_flank": N_FLANK, "em_iters": EM_ITERS},
              "tuning": tune, "quadrants": {}}
    for q in ["val", "test_novel_sample"]:
        last = "full" if reads else "snp"
        s = out.filter((pl.col("quadrant") == q) & pl.col(f"{last}_lo").is_not_null())
        rep = {"equal halves": metrics(s, "half_lo", "half_hi"), "catalog (stage 1)": metrics(s)}
        for name, tag in (("SNP haplotypes", "snp"), ("catalog + spanning reads", "cat_span"),
                          ("SNP haplotypes + spanning reads", "full")):
            if f"{tag}_lo" in s.columns:
                rep[name] = metrics(s, f"{tag}_lo", f"{tag}_hi")
        if reads:
            span = s.filter(pl.col("full_n_span_reads") > 0)
            rep["_at loci with spanning reads"] = {n: metrics(span, *c) for n, c in (
                ("catalog (stage 1)", ("h_lo", "h_hi")), ("catalog + spanning reads", ("cat_span_lo", "cat_span_hi")),
                ("SNP haplotypes + spanning reads", ("full_lo", "full_hi")))}
        report["quadrants"][q] = rep
        print(f"{q}: " + " | ".join(f"{n}: {m['both_exact']:.3f} (het {m['het_both_exact']:.3f})"
                                  for n, m in rep.items() if not n.startswith("_")), flush=True)
    report["params"]["reads"] = reads
    (mdir / "haplotype_snp_eval.json").write_text(json.dumps(report, indent=2, default=float))
    print(f"wrote {mdir/'haplotype_snp_predictions.parquet'}, {mdir/'haplotype_snp_eval.json'}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--tune-loci", type=int, default=6000)
    ap.add_argument("--no-reads", action="store_true", help="SNP haplotypes only (no spanning-read term)")
    args = ap.parse_args(argv)
    build(load_config(args.config), args.tune_loci, reads=not args.no_reads)


if __name__ == "__main__":
    main()
