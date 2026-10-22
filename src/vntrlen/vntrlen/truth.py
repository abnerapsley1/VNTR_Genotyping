"""Ground-truth input for model building, and the panel statistics derived from it.

Truth file (tab-separated, header):
  sample_id                      matches the sample sheet
  locus_id  or  chrom/start/end  start 0-based, as in the BED
  length                         combined length of both alleles
  unit                           'units' (repeat units) or 'bp'
  allele1, allele2               optional per-allele lengths in the same unit; enable
                                 the allele-count and heterozygosity statistics

units <-> bp uses the repeat's mean motif length (an approximation when an allele
mixes motifs of different lengths).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl

PANEL_COLUMNS = [
    "cov_n_called", "cov_n_alleles", "cov_h_exp", "cov_h_obs",
    "cov_mean_combined_units", "cov_var_combined_units", "cov_min_combined_units", "cov_max_combined_units",
    "cov_mean_combined_bp", "cov_var_combined_bp", "cov_n_distinct_hap_units", "cov_panel_cn_range_units",
    "cov_panel_median_units", "cov_panel_p5_units", "cov_panel_p95_units", "cov_panel_core90_units",
    "cov_panel_n_outliers", "cov_panel_outlier_frac",
]


def read_truth(path: str | Path, bed: pl.DataFrame) -> pl.DataFrame:
    t = pl.read_csv(path, separator="\t", infer_schema_length=0, comment_prefix="#")
    need = {"sample_id", "length", "unit"}
    if not need <= set(t.columns):
        raise ValueError(f"{path}: needs columns {sorted(need)} plus locus_id or chrom/start/end")
    if "locus_id" not in t.columns:
        t = t.with_columns(locus_id=pl.col("chrom") + ":" + (pl.col("start").cast(pl.Int64) + 1)
                           .cast(pl.Utf8) + "-" + pl.col("end"))
    bad = set(t["unit"].str.to_lowercase().unique()) - {"units", "bp"}
    if bad:
        raise ValueError(f"{path}: unit must be 'units' or 'bp'; found {sorted(bad)}")
    mlen = bed.select("locus_id", bpu=pl.col("motifs").str.split(",")
                      .list.eval(pl.element().str.len_chars()).list.mean())
    t = t.join(mlen, on="locus_id", how="inner")
    is_bp = pl.col("unit").str.to_lowercase() == "bp"
    num = lambda c: pl.col(c).cast(pl.Float64, strict=False)  # noqa: E731
    cols = [pl.when(is_bp).then(num("length") / pl.col("bpu")).otherwise(num("length")).alias("combined_units"),
            pl.when(is_bp).then(num("length")).otherwise(num("length") * pl.col("bpu")).alias("combined_bp")]
    for a in ("allele1", "allele2"):
        if a in t.columns:
            cols.append(pl.when(is_bp).then(num(a) / pl.col("bpu")).otherwise(num(a)).round(0).alias(f"{a}_units"))
    t = t.with_columns(cols).filter(pl.col("combined_units").is_finite())
    keep = ["sample_id", "locus_id", "combined_units", "combined_bp"] + [c for c in t.columns if c.endswith("_units")
                                                                        and c.startswith("allele")]
    return t.select(keep)


def panel_stats(truth: pl.DataFrame, exclude: set[str] | None = None) -> pl.DataFrame:
    """Per-repeat panel statistics from training truth (optionally leaving people out)."""
    g = truth.filter(~pl.col("sample_id").is_in(list(exclude or set())))
    comb = pl.col("combined_units")
    core = g.group_by("locus_id").agg(
        cov_n_called=pl.len().cast(pl.Float64),
        cov_mean_combined_units=comb.mean(), cov_var_combined_units=comb.var(),
        cov_min_combined_units=comb.min(), cov_max_combined_units=comb.max(),
        cov_mean_combined_bp=pl.col("combined_bp").mean(), cov_var_combined_bp=pl.col("combined_bp").var(),
        cov_panel_median_units=comb.median(),
        cov_panel_p5_units=comb.quantile(0.05, "nearest"), cov_panel_p95_units=comb.quantile(0.95, "nearest"),
    ).with_columns(cov_panel_cn_range_units=pl.col("cov_max_combined_units") - pl.col("cov_min_combined_units"),
                   cov_panel_core90_units=pl.col("cov_panel_p95_units") - pl.col("cov_panel_p5_units"))
    w = 2 * pl.max_horizontal(pl.col("cov_panel_core90_units"), pl.lit(1))
    out = (g.join(core.select("locus_id", "cov_panel_p5_units", "cov_panel_p95_units", "cov_panel_core90_units"),
                  on="locus_id")
           .group_by("locus_id").agg(cov_panel_n_outliers=((comb > pl.col("cov_panel_p95_units") + w)
                                                          | (comb < pl.col("cov_panel_p5_units") - w)).sum()
                                     .cast(pl.Float64)))
    core = core.join(out, on="locus_id").with_columns(
        cov_panel_outlier_frac=pl.col("cov_panel_n_outliers") / pl.col("cov_n_called"))
    if {"allele1_units", "allele2_units"} <= set(g.columns):
        haps = pl.concat([g.select("locus_id", u=pl.col("allele1_units")),
                          g.select("locus_id", u=pl.col("allele2_units"))]).drop_nulls()
        hx = (haps.group_by("locus_id", "u").agg(c=pl.len())
              .with_columns(p=pl.col("c") / pl.col("c").sum().over("locus_id"))
              .group_by("locus_id").agg(cov_h_exp=1 - (pl.col("p") ** 2).sum(),
                                        cov_n_distinct_hap_units=pl.len().cast(pl.Float64),
                                        cov_n_alleles=pl.len().cast(pl.Float64)))
        ho = g.drop_nulls(["allele1_units", "allele2_units"]).group_by("locus_id").agg(
            cov_h_obs=(pl.col("allele1_units") != pl.col("allele2_units")).mean())
        core = core.join(hx, on="locus_id", how="left").join(ho, on="locus_id", how="left")
    for c in PANEL_COLUMNS:
        if c not in core.columns:
            core = core.with_columns(pl.lit(None, dtype=pl.Float64).alias(c))
    return core.select(["locus_id"] + PANEL_COLUMNS)
