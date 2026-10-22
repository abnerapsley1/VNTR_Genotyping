"""
panel_stats.py — per-locus long-read panel statistics, EXCLUDING held-out people.

The covariates built from gt_locus_summary use all 416 vamos individuals, which
include the short-read val/test samples.  Panel min/max (and to a lesser
degree mean/var/allele counts) then carry a held-out person's own answer when
they hold a rare allele.  In deployment the panel never contains the person
being genotyped, so for modelling the panel is recomputed here from gt_long
without the val/test samples (the training samples stay in, as they would in
the shipped panel).

Besides the original summaries it adds a robust split of "how variable is this
locus" into two parts:
  panel_core90_units   p95 - p5 of combined length: spread of the bulk of people
                       (a handful of extreme individuals cannot move it)
  panel_n_outliers     people beyond the core by > 2 x max(core, 1) units
  panel_outlier_frac   the same as a fraction of called people
"""

from __future__ import annotations

import polars as pl

# panel-derived covariate columns this module (re)defines, without the cov_ prefix
PANEL_COLUMNS = [
    "n_called", "n_alleles", "h_exp", "h_obs",
    "mean_combined_units", "var_combined_units", "min_combined_units", "max_combined_units",
    "mean_combined_bp", "var_combined_bp", "n_distinct_hap_units", "panel_cn_range_units",
    "panel_median_units", "panel_p5_units", "panel_p95_units", "panel_core90_units",
    "panel_n_outliers", "panel_outlier_frac",
]


def panel_stats(gt_long: pl.DataFrame, exclude_vamos_samples: set[str]) -> pl.DataFrame:
    """One row per locus_id with PANEL_COLUMNS, from gt_long minus excluded samples."""
    g = gt_long.filter(~pl.col("vamos_sample").is_in(list(exclude_vamos_samples)))
    comb = pl.col("combined_units")
    called = g.filter(comb.is_not_null())

    core = called.group_by("locus_id").agg(
        n_called=pl.len(),
        mean_combined_units=comb.mean(),
        var_combined_units=comb.var(),
        min_combined_units=comb.min(),
        max_combined_units=comb.max(),
        mean_combined_bp=pl.col("combined_bp").mean(),
        var_combined_bp=pl.col("combined_bp").var(),
        h_obs=(pl.col("a1_units") != pl.col("a2_units")).mean(),
        panel_median_units=comb.median(),
        panel_p5_units=comb.quantile(0.05, "nearest"),
        panel_p95_units=comb.quantile(0.95, "nearest"),
    ).with_columns(
        panel_cn_range_units=pl.col("max_combined_units") - pl.col("min_combined_units"),
        panel_core90_units=pl.col("panel_p95_units") - pl.col("panel_p5_units"),
    )

    w = 2 * pl.max_horizontal(pl.col("panel_core90_units"), pl.lit(1))
    out = (called.join(core.select("locus_id", "panel_p5_units", "panel_p95_units",
                                   "panel_core90_units"), on="locus_id")
           .group_by("locus_id")
           .agg(panel_n_outliers=((comb > pl.col("panel_p95_units") + w)
                                  | (comb < pl.col("panel_p5_units") - w)).sum()))

    # haplotype-level allele stats: expected heterozygosity over allele lengths,
    # and the number of distinct ALTANNO alleles actually carried by kept people
    haps = pl.concat([
        g.select("locus_id", pl.col("a1_units").alias("u"),
                 pl.col("gt_raw").str.split("/").list.get(0, null_on_oob=True).alias("allele")),
        g.select("locus_id", pl.col("a2_units").alias("u"),
                 pl.col("gt_raw").str.split("/").list.get(1, null_on_oob=True).alias("allele")),
    ]).filter(pl.col("u").is_not_null())
    h_exp = (haps.group_by("locus_id", "u").agg(c=pl.len())
             .with_columns(p=pl.col("c") / pl.col("c").sum().over("locus_id"))
             .group_by("locus_id").agg(h_exp=1 - (pl.col("p") ** 2).sum()))
    n_all = (haps.filter(pl.col("allele").is_not_null() & (pl.col("allele") != "."))
             .group_by("locus_id").agg(n_alleles=pl.col("allele").n_unique()))

    # over every person with >= 1 haplotype called (as parse_vamos defines it)
    n_hap = g.group_by("locus_id").agg(n_distinct_hap_units=pl.col("hap_hi_units").n_unique())

    res = (core.join(out, on="locus_id", how="left")
           .join(n_hap, on="locus_id", how="left")
           .join(h_exp, on="locus_id", how="left")
           .join(n_all, on="locus_id", how="left")
           .with_columns(panel_outlier_frac=pl.col("panel_n_outliers") / pl.col("n_called")))
    return res.select(["locus_id"] + PANEL_COLUMNS)
