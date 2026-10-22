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

Individuals are taken from the manifest and assigned per family
(config/sample_families.tsv), so relatives always share a fold.

Outputs: work/model/locus_split.parquet, work/model/sample_split.parquet
"""

from __future__ import annotations

import argparse
import hashlib

import polars as pl

from ..paths import load_config
from ..run.manifest import load_samples


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

    # Individuals come from the manifest, so the split is fixed before extraction.
    # Folds are drawn per FAMILY (config/sample_families.tsv): relatives share
    # alleles, and a parent in training with a child in test would inflate test
    # scores.  Unrelated individuals keep their per-sample fold.
    sample_ids = sorted(s["sample_id"] for s in load_samples(cfg))
    fam_path = cfg.path("inputs.samples_tsv").with_name("sample_families.tsv")
    family = {}
    if fam_path.exists():
        for line in fam_path.read_text().splitlines():
            if line.startswith("#") or line.startswith("sample_id\t") or not line.strip():
                continue
            sid, fid = line.split("\t")[:2]
            family[sid] = fid
    # 1kGP gives every sample a family code even when no relative is in the
    # panel; only group by it when >= 2 panel members share it, so unrelated
    # individuals keep the per-sample fold they had before families were used.
    n_in_panel = {}
    for s_ in sample_ids:
        n_in_panel[family.get(s_, s_)] = n_in_panel.get(family.get(s_, s_), 0) + 1
    key = {s_: (family[s_] if s_ in family and n_in_panel[family[s_]] > 1 else s_) for s_ in sample_ids}
    fams = _assign(sorted(set(key.values())), seed + 1, test_s, val_s)
    fold_of = dict(zip(fams["id"], fams["fold"]))
    sdf = pl.DataFrame({"sample_id": sample_ids, "family_id": [key[s_] for s_ in sample_ids],
                        "fold": [fold_of[key[s_]] for s_ in sample_ids]})

    ldf = _assign(loci, seed, test_l, val_l).rename({"id": "locus_id"})

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
