"""Supplementary paired tests (ours vs danbing-tk) where per-person tests do not apply.

  within-locus tracking  per locus with true variation among the 52 test people:
                         Pearson r of call vs truth for each tool; Wilcoxon
                         signed-rank on the paired per-locus r (overall and by stratum)
  rare alleles           McNemar exact test on paired right/wrong direction calls
                         for rare long/short allele carriers (overall and by stratum)

    python src/VNTR_Benchmark/stats_supplement.py -> work/compare/stats_supplement.csv
"""
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats_ours_vs_danbing import CMP, prep  # noqa: E402


def locus_r(g: pl.DataFrame) -> pl.DataFrame:
    return (g.group_by("locus_id").agg(
        n=pl.len(), sd_t=pl.col("combined_bp").cast(pl.Float64).std(),
        r_o=pl.corr("ours_bp", pl.col("combined_bp").cast(pl.Float64)),
        r_d=pl.corr("danbing_bp", pl.col("combined_bp").cast(pl.Float64)))
        .filter((pl.col("n") >= 20) & (pl.col("sd_t") > 0))
        .with_columns(pl.col("r_o", "r_d").fill_nan(0.0).fill_null(0.0)))   # constant call -> r = 0


def main():
    d = prep()
    rows = []
    for name, col in [("all", None), ("core_bp_bin", "core_bp_bin"), ("site_type", "site_type"),
                      ("span_bin", "span_bin"), ("segdup", "segdup")]:
        groups = [("all", d)] if col is None else [(str(k), g) for (k,), g in d.group_by([col], maintain_order=True)]
        for lab, g in groups:
            if lab in ("None", "null"):
                continue
            lr = locus_r(g)
            if lr.height >= 10:
                rows.append({"stratum": name, "group": lab, "test": "per-locus r (Wilcoxon)", "n": lr.height,
                             "ours": float(lr["r_o"].median()), "danbing": float(lr["r_d"].median()),
                             "favouring_ours": int((lr["r_o"] > lr["r_d"]).sum()),
                             "favouring_danbing": int((lr["r_d"] > lr["r_o"]).sum()),
                             "p": float(wilcoxon(lr["r_o"], lr["r_d"]).pvalue)})
            rc = g.filter(pl.col("rare_carrier"))
            if rc.height >= 10:
                o, b = rc["ours_rare_ok"].to_numpy(), rc["danbing_rare_ok"].to_numpy()
                n01, n10 = int((~o & b).sum()), int((o & ~b).sum())
                p = binomtest(n10, n10 + n01).pvalue if n10 + n01 else np.nan
                rows.append({"stratum": name, "group": lab, "test": "rare alleles right direction (McNemar)",
                             "n": rc.height, "ours": float(o.mean()), "danbing": float(b.mean()),
                             "favouring_ours": n10, "favouring_danbing": n01, "p": p})
    res = pl.DataFrame(rows)
    p = res["p"].to_numpy()
    o = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for i, k in enumerate(o):
        run = max(run, (len(p) - i) * p[k])
        adj[k] = min(1.0, run)
    res = res.with_columns(p_holm=pl.Series(adj))
    res.write_csv(CMP / "stats_supplement.csv")
    pl.Config.set_tbl_rows(100)
    pl.Config.set_tbl_width_chars(220)
    print(res.with_columns(pl.col("ours", "danbing").round(3)))


if __name__ == "__main__":
    main()
