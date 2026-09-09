"""
splits.py — frozen, deterministic train/val/test assignment.

Loci and samples are split INDEPENDENTLY (hash of id + seed), so evaluation can
separate two kinds of generalisation:

    quadrant            loci          samples      what it tests
    ----------------    -----------   ----------   -----------------------------
    train               train         train        (fitting)
    test_novel_sample   train         test         new individual, known locus
    test_novel_locus    test          train        new locus (arbitrary user BED)
    test_novel_both     test          test         both
    val_*               val fold on either axis -> hyper-parameter tuning

Outputs: work/model/locus_split.parquet, work/model/sample_split.parquet
"""

from __future__ import annotations

import argparse
import hashlib

import polars as pl

from ..paths import load_config


def _bucket(key: str, seed: int) -> float:
    h = hashlib.sha1(f"{seed}:{key}".encode()).digest()
    return int.from_bytes(h[:8], "big") / 2**64


def _assign(ids: list[str], seed: int, test_frac: float, val_frac: float) -> pl.DataFrame:
    rows = []
    for x in ids:
        u = _bucket(x, seed)
        fold = "test" if u < test_frac else ("val" if u < test_frac + val_frac else "train")
        rows.append(dict(id=x, u=u, fold=fold))
    return pl.DataFrame(rows)


def build(cfg) -> tuple:
    seed = int(cfg.get("model.random_seed", 20260831))
    test_l = float(cfg.get("model.holdout_loci_frac", 0.20))
    test_s = float(cfg.get("model.holdout_samples_frac", 0.20))
    val_l = float(cfg.get("model.val_loci_frac", 0.10))
    val_s = float(cfg.get("model.val_samples_frac", 0.10))

    cov = pl.read_parquet(cfg.path("outputs.covariates_dir") / "locus_covariates.parquet")
    loci = cov["locus_id"].to_list()

    feat_dir = cfg.path("outputs.features_dir")
    sample_ids = sorted(p.stem for p in feat_dir.glob("*.parquet")
                        if not p.stem.startswith("_"))
    if not sample_ids:
        print("WARNING: no features/*.parquet yet — sample split empty until extraction runs")

    ldf = _assign(loci, seed, test_l, val_l).rename({"id": "locus_id"})
    sdf = _assign(sample_ids, seed + 1, test_s, val_s).rename({"id": "sample_id"})

    mdir = cfg.ensure_dir("outputs.model_dir")
    ldf.write_parquet(mdir / "locus_split.parquet")
    sdf.write_parquet(mdir / "sample_split.parquet")
    print("loci:", ldf["fold"].value_counts().sort("fold"))
    print("samples:", sdf["fold"].value_counts().sort("fold") if sdf.height else "(none)")
    print(f"wrote {mdir/'locus_split.parquet'} and {mdir/'sample_split.parquet'}")
    return ldf, sdf


def quadrant_expr() -> pl.Expr:
    """polars expr: given columns locus_fold, sample_fold -> quadrant label."""
    lf, sf = pl.col("locus_fold"), pl.col("sample_fold")
    return (
        pl.when((lf == "train") & (sf == "train")).then(pl.lit("train"))
        .when((lf == "val") | (sf == "val")).then(pl.lit("val"))
        .when((lf == "train") & (sf == "test")).then(pl.lit("test_novel_sample"))
        .when((lf == "test") & (sf == "train")).then(pl.lit("test_novel_locus"))
        .when((lf == "test") & (sf == "test")).then(pl.lit("test_novel_both"))
        .otherwise(pl.lit("other"))
    )


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
