"""Is our tool significantly better than danbing-tk?  Paired tests on the test people.

Both tools are scored on exactly the same (person, locus) calls: sites where
both made a call and the truth is known (compare_tools.py calls.parquet).
Calls from one person, or from one locus, are not independent, so tests use
people or loci as units:

  per person  each of the 52 test people gets one value per measure and tool;
              Wilcoxon signed-rank test on the 52 paired differences
              (+ how many people favour each tool)
  per locus   each locus gets each tool's mean absolute error over people;
              Wilcoxon signed-rank on the paired locus differences
  bootstrap   95% CI for the difference (ours - danbing), resampling people
              (500 draws; 200 per stratum) and, for the overall row, loci (500 draws)

Measures: median relative error, share within 5% of true length, mean
absolute error (bp), within-locus r (person-to-person tracking), and the share
of rare long/short alleles called in the right direction.  Overall and per
stratum (core spread bp, site type, reference array length, segdup); p-values
Holm-corrected across all tests in the table.

    python src/VNTR_Benchmark/stats_ours_vs_danbing.py
-> work/compare/stats_ours_vs_danbing.csv
"""
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[2]
CMP = ROOT / "src" / "VNTR_Benchmark" / "work" / "compare"
N_BOOT = 500
RNG = np.random.default_rng(20260831)


def prep() -> pl.DataFrame:
    d = pl.read_parquet(CMP / "calls.parquet").filter(
        pl.col("ours_bp").is_finite() & pl.col("danbing_bp").is_finite() & (pl.col("combined_bp") > 0))
    x = pl.col("combined_bp").cast(pl.Float64)
    for t in ("ours", "danbing"):
        d = d.with_columns((pl.col(f"{t}_bp") - x).abs().alias(f"{t}_ae"))
        d = d.with_columns((pl.col(f"{t}_ae") / x).alias(f"{t}_re"))
    # rare-carrier direction (as in compare_tools.py): call beyond the panel core edge on the right side
    bpr = x / pl.max_horizontal(pl.col("combined_units"), pl.lit(1))
    edge = pl.when(pl.col("rare_up")).then(pl.col("cov_panel_p95_units")).otherwise(pl.col("cov_panel_p5_units")) * bpr
    for t in ("ours", "danbing"):
        d = d.with_columns(pl.when(pl.col("rare_carrier")).then(
            pl.when(pl.col("rare_up")).then(pl.col(f"{t}_bp") > edge).otherwise(pl.col(f"{t}_bp") < edge)
        ).alias(f"{t}_rare_ok"))
    return d


def within_r(locus: np.ndarray, y: np.ndarray, x: np.ndarray) -> float:
    """Pooled within-locus correlation: centre call and truth on each locus mean."""
    _, inv = np.unique(locus, return_inverse=True)
    cnt = np.bincount(inv)
    yc = y - (np.bincount(inv, y) / cnt)[inv]
    xc = x - (np.bincount(inv, x) / cnt)[inv]
    den = np.sqrt((yc ** 2).sum() * (xc ** 2).sum())
    return float((yc * xc).sum() / den) if den > 0 else np.nan


def measures(g: pl.DataFrame, t: str) -> dict:
    out = {"median_rel_err": float(g[f"{t}_re"].median()),
           "within_5pct": float((g[f"{t}_re"] <= 0.05).mean()),
           "mean_abs_err_bp": float(g[f"{t}_ae"].mean())}
    out["within_locus_r"] = within_r(g["locus_id"].to_physical().to_numpy() if g["locus_id"].dtype == pl.Categorical
                                     else g["locus_id"].cast(pl.Categorical).to_physical().to_numpy(),
                                     g[f"{t}_bp"].to_numpy().astype(float), g["combined_bp"].to_numpy().astype(float))
    r = g.filter(pl.col("rare_carrier"))
    out["rare_right_direction"] = float(r[f"{t}_rare_ok"].mean()) if r.height >= 20 else np.nan
    return out


MEASURES = ["median_rel_err", "within_5pct", "mean_abs_err_bp", "within_locus_r", "rare_right_direction"]
LOWER_BETTER = {"median_rel_err", "mean_abs_err_bp"}


def per_unit(g: pl.DataFrame, unit: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Per person (or locus) values of each measure for both tools."""
    vals = {m: ([], []) for m in MEASURES}
    for _, h in g.group_by(unit):
        if h.height < 30:
            continue
        a, b = measures(h, "ours"), measures(h, "danbing")
        for m in MEASURES:
            if np.isfinite(a[m]) and np.isfinite(b[m]):
                vals[m][0].append(a[m])
                vals[m][1].append(b[m])
    return {m: (np.array(v[0]), np.array(v[1])) for m, v in vals.items()}


def boot(g: pl.DataFrame, unit: str, n: int) -> dict[str, tuple[float, float]]:
    """95% CI of (ours - danbing) resampling whole `unit`s."""
    ids = g[unit].unique().to_numpy()
    parts = {k: h for (k,), h in g.group_by([unit])}
    draws = {m: [] for m in MEASURES}
    for _ in range(n):
        pick = RNG.choice(ids, size=len(ids), replace=True)
        s = pl.concat([parts[k] for k in pick])
        if unit == "locus_id":            # keep loci distinct for the within-locus centring
            s = s.with_columns(locus_id=pl.col("locus_id") + "_" + pl.int_range(pl.len()).over("locus_id").cast(pl.Utf8))
        a, b = measures(s, "ours"), measures(s, "danbing")
        for m in MEASURES:
            draws[m].append(a[m] - b[m])
    return {m: tuple(np.nanpercentile(v, [2.5, 97.5])) for m, v in draws.items()}


def holm(p: np.ndarray) -> np.ndarray:
    o = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for i, k in enumerate(o):
        run = max(run, (len(p) - i) * p[k])
        adj[k] = min(1.0, run)
    return adj


def main():
    d = prep()
    print(f"{d.height:,} paired calls; {d['sample_id'].n_unique()} people; {d['locus_id'].n_unique():,} loci")
    strata = [("all", None), ("core_bp_bin", "core_bp_bin"), ("site_type", "site_type"),
              ("span_bin", "span_bin"), ("segdup", "segdup")]
    rows = []
    for name, col in strata:
        groups = [("all", d)] if col is None else [(str(k), g) for (k,), g in d.group_by([col], maintain_order=True)]
        for lab, g in groups:
            if lab in ("None", "null") or g.height < 1000:
                continue
            est_o, est_d = measures(g, "ours"), measures(g, "danbing")
            pp = per_unit(g, "sample_id")
            n_boot = N_BOOT if name == "all" else 200
            ci_p = boot(g, "sample_id", n_boot)
            # per-locus test on mean absolute error only (other measures need many people per locus)
            lm = g.group_by("locus_id").agg(o=pl.col("ours_ae").mean(), b=pl.col("danbing_ae").mean())
            p_loc = wilcoxon(lm["o"], lm["b"]).pvalue if lm.height > 10 else np.nan
            for m in MEASURES:
                a, b = pp[m]
                better = (a < b) if m in LOWER_BETTER else (a > b)
                p = wilcoxon(a, b).pvalue if len(a) >= 10 and np.any(a != b) else np.nan
                rows.append({"stratum": name, "group": lab, "measure": m, "n_calls": g.height,
                             "n_loci": g["locus_id"].n_unique(), "ours": est_o[m], "danbing": est_d[m],
                             "diff_ours_minus_danbing": est_o[m] - est_d[m],
                             "ci95_people_lo": ci_p[m][0], "ci95_people_hi": ci_p[m][1],
                             "people_favouring_ours": int(better.sum()), "people_tested": len(a),
                             "p_wilcoxon_people": p,
                             "p_wilcoxon_loci_mean_abs_err": p_loc if m == "mean_abs_err_bp" else np.nan})
            print(f"  {name}/{lab}: done", flush=True)
    res = pl.DataFrame(rows)
    p = res["p_wilcoxon_people"].to_numpy()
    ok = np.isfinite(p)
    adj = np.full_like(p, np.nan)
    adj[ok] = holm(p[ok])
    res = res.with_columns(p_holm=pl.Series(adj))
    # locus-resampled CI for the overall row
    ci_l = boot(d, "locus_id", 500)
    res = res.with_columns(
        ci95_loci_lo=pl.when(pl.col("stratum") == "all").then(pl.col("measure").replace_strict({m: ci_l[m][0] for m in MEASURES}))
        .otherwise(None),
        ci95_loci_hi=pl.when(pl.col("stratum") == "all").then(pl.col("measure").replace_strict({m: ci_l[m][1] for m in MEASURES}))
        .otherwise(None))
    res.write_csv(CMP / "stats_ours_vs_danbing.csv")
    pl.Config.set_tbl_rows(100)
    pl.Config.set_tbl_width_chars(250)
    print(res.filter(pl.col("stratum") == "all").select(
        "measure", "ours", "danbing", "ci95_people_lo", "ci95_people_hi", "ci95_loci_lo", "ci95_loci_hi",
        "people_favouring_ours", "people_tested", "p_holm"))
    print(f"wrote {CMP / 'stats_ours_vs_danbing.csv'}")


if __name__ == "__main__":
    main()
