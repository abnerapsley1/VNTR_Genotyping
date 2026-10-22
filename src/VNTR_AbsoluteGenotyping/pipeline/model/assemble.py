"""
assemble.py — build the (locus, sample) modelling table.

Joins, on (locus_id[, sample_id]):
  - features/<sample>.parquet    (per-locus signals, all samples stacked)
  - features/<sample>.nuisance.json  -> sample-level columns (prefix nuis_)
  - work/ground_truth/gt_long.parquet -> target (combined_units, combined_bp)
  - work/covariates/locus_covariates.parquet -> locus covariates (prefix cov_);
    panel-derived ones recomputed without held-out people (panel_stats.py)
  - work/model/locus_split.parquet / sample_split.parquet -> folds + quadrant

Keeps only rows with extraction_ok AND n_haps_called == 2 (combined length defined).
Writes work/model/model_table.parquet and prints a quadrant x target summary.

vamos_column in the manifest maps a sample_id to its vamos_sample; ground truth
is looked up by that.
"""

from __future__ import annotations

import argparse
import os
import json

import numpy as np
import polars as pl

from ..paths import load_config
from ..run.manifest import load_samples
from .panel_stats import PANEL_COLUMNS, panel_stats
from .splits import quadrant_expr

_FEATURE_DROP = {"array_start0", "array_end0", "extraction_ok", "note",
                 "n_reads_fetched"}
_COV_DROP = ["chrom", "start", "end", "gene", "locus_key", "match_type"]


def _clean_vamos_col(col: str) -> str:
    left = col.split("/")[0]
    for suf in ("_hap1", "_hap2"):
        if left.endswith(suf):
            left = left[: -len(suf)]
    return left.rstrip("_")


def build(cfg) -> str:
    gt_dir = cfg.path("outputs.ground_truth_dir")
    feat_dir = cfg.path("outputs.features_dir")
    mdir = cfg.ensure_dir("outputs.model_dir")

    # ---- ground truth (locus, vamos_sample) -> combined length ----
    gt = (pl.read_parquet(gt_dir / "gt_long.parquet")
          .filter(pl.col("n_haps_called") == 2)
          .select(["locus_id", "vamos_sample", "combined_units", "combined_bp",
                   "hap_hi_units", "hap_lo_units"]))

    # ---- sample_id <-> vamos_sample map from the manifest ----
    smap = {}
    for s in load_samples(cfg):
        if s.get("vamos_column"):
            smap[s["sample_id"]] = _clean_vamos_col(s["vamos_column"])
    if not smap:
        raise SystemExit("no vamos_column entries in samples.tsv — cannot join ground truth")

    # ---- features + nuisance, per sample ----
    frames = []
    for sid, vs in smap.items():
        fp = feat_dir / f"{sid}.parquet"
        if not fp.exists():
            print(f"  (skip {sid}: no features parquet)")
            continue
        f = pl.read_parquet(fp)
        f = f.filter(pl.col("extraction_ok") == 1.0)
        np_ = feat_dir / f"{sid}.nuisance.json"
        nuis = json.loads(np_.read_text()) if np_.exists() else {}
        for k, v in nuis.items():
            if k in ("sample_id", "gc_depth_curve"):
                continue
            f = f.with_columns(pl.lit(v).alias(f"nuis_{k}"))
        f = f.with_columns(pl.lit(sid).alias("sample_id"),
                           pl.lit(vs).alias("vamos_sample"))
        frames.append(f)
    if not frames:
        raise SystemExit("no per-sample feature parquets found — run extraction first")
    feats = pl.concat(frames, how="diagonal_relaxed")

    # ---- splits ----
    lsplit = pl.read_parquet(mdir / "locus_split.parquet").select(
        pl.col("locus_id"), pl.col("fold").alias("locus_fold"))
    ssplit = pl.read_parquet(mdir / "sample_split.parquet").select(
        pl.col("sample_id"), pl.col("fold").alias("sample_fold"))

    # ---- covariates ----
    # panel-derived columns are recomputed WITHOUT the held-out (val/test)
    # people, so a held-out person's own long-read call never reaches a feature
    cov = pl.read_parquet(cfg.path("outputs.covariates_dir") / "locus_covariates.parquet")
    cov = cov.drop([c for c in _COV_DROP + PANEL_COLUMNS if c in cov.columns])
    held_out = {smap[s] for s in ssplit.filter(pl.col("sample_fold") != "train")["sample_id"]
                if s in smap}
    panel = panel_stats(pl.read_parquet(gt_dir / "gt_long.parquet"), held_out)
    print(f"panel statistics: {len(held_out)} held-out people excluded from the "
          f"long-read panel")
    cov = cov.join(panel, on="locus_id", how="left")
    cov = cov.rename({c: f"cov_{c}" for c in cov.columns if c != "locus_id"})

    tbl = (feats
           .join(gt, on=["locus_id", "vamos_sample"], how="inner")
           .join(cov, on="locus_id", how="left")
           .join(lsplit, on="locus_id", how="left")
           .join(ssplit, on="sample_id", how="left")
           .with_columns(locus_fold=pl.col("locus_fold").fill_null("train"),
                         sample_fold=pl.col("sample_fold").fill_null("train"))
           .with_columns(quadrant=quadrant_expr())
           .drop([c for c in _FEATURE_DROP if c in feats.columns]))

    out = mdir / "model_table.parquet"
    tbl.write_parquet(out)
    print(f"wrote {out}  ({tbl.height:,} rows x {tbl.width} cols)")
    print(tbl.group_by("quadrant").agg(
        n=pl.len(),
        n_loci=pl.col("locus_id").n_unique(),
        n_samples=pl.col("sample_id").n_unique(),
        median_combined_units=pl.col("combined_units").median(),
    ).sort("quadrant"))
    return str(out)


# Population summaries precomputed by external long-read studies (ANOVA tables and the
# population columns of the mutation-rate table).  They were computed over panels that
# likely include our validation/test people, so they are never used as features.
# Dropping them changed test accuracy by <= 0.02 repeats (work/model_withext kept them).
EXTERNAL_POPULATION_FEATURES = {
    "cov_anno_len_between_var",
    "cov_anno_len_effect_size",
    "cov_anno_len_f",
    "cov_anno_len_grand_mean",
    "cov_anno_len_p",
    "cov_anno_len_within_var",
    "cov_ed2_major_between_var",
    "cov_ed2_major_effect_size",
    "cov_ed2_major_f",
    "cov_ed2_major_grand_mean",
    "cov_ed2_major_p",
    "cov_ed2_major_within_var",
    "cov_mr_all_gc",
    "cov_mr_ave_block",
    "cov_mr_hete_comp",
    "cov_mr_hete_len",
    "cov_mr_homo_comp",
    "cov_mr_homo_len",
    "cov_mr_motif_entropy",
    "cov_mr_mut_rate",
    "cov_mr_mut_rate_normlen",
    "cov_mr_num_pattern",
    "cov_mr_snp_density",
    "cov_mr_std_len_units",
}


def feature_columns(tbl: pl.DataFrame) -> list[str]:
    """Numeric predictor columns: signals + nuis_* + cov_* (excludes ids/target/folds).

    EXTERNAL_POPULATION_FEATURES are always dropped; VNTR_EXCLUDE_FEATURES
    (comma-separated) drops further columns, for ablation runs."""
    exclude = {"sample_id", "locus_id", "vamos_sample", "chrom", "note",
               "combined_units", "combined_bp", "hap_hi_units", "hap_lo_units",
               "locus_fold", "sample_fold", "quadrant"}
    exclude |= EXTERNAL_POPULATION_FEATURES
    exclude |= {c for c in os.environ.get("VNTR_EXCLUDE_FEATURES", "").split(",") if c}
    cols = []
    for c, dt in zip(tbl.columns, tbl.dtypes):
        if c in exclude:
            continue
        if dt.is_numeric():
            cols.append(c)
    return cols


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
