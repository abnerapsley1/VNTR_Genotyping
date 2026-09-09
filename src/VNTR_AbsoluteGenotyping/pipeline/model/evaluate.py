"""
evaluate.py — score work/model/predictions.parquet against ground truth.

Per quadrant (train / test_novel_sample / test_novel_locus / test_novel_both):
  pooled Pearson & Spearman
  median per-locus Pearson & Spearman   (the number that matters for genotyping)
  MAE (motif units), median relative error
  calibration slope  (OLS y_true ~ y_pred; ideal 1.0)
  |err| <= 1 unit, <= 2 units fractions
  90% interval coverage
  fraction of loci with per-locus Pearson > 0.5 / > 0.7
Also a PASS-tier-only view, and a stratification by segdup / size.

Writes work/model/evaluation.json and prints a table.  Compares against the
naive baseline (density_ratio_d20_w2000 * 2 * ref_units) where available.
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import polars as pl
from scipy.stats import pearsonr, spearmanr

from ..paths import load_config

_QUADS = ["train", "val", "test_novel_sample", "test_novel_locus", "test_novel_both"]


def _per_locus_corr(df: pl.DataFrame, col_pred: str, min_n: int = 5):
    prs, sps = [], []
    for (_l,), g in df.group_by(["locus_id"]):
        if g.height < min_n:
            continue
        yt = g["combined_units"].to_numpy()
        yp = g[col_pred].to_numpy()
        if np.std(yt) < 1e-9 or np.std(yp) < 1e-9:
            continue
        prs.append(pearsonr(yp, yt).statistic)
        sps.append(spearmanr(yp, yt).statistic)
    return np.array(prs), np.array(sps)


def _block(df: pl.DataFrame, pred_col: str) -> dict:
    d = df.filter(pl.col(pred_col).is_not_null() & pl.col("combined_units").is_not_null())
    if d.height < 5:
        return {"n": d.height}
    yt = d["combined_units"].to_numpy()
    yp = d[pred_col].to_numpy()
    err = yp - yt
    out = {
        "n": int(d.height),
        "n_loci": int(d["locus_id"].n_unique()),
        "pooled_pearson": float(pearsonr(yp, yt).statistic),
        "pooled_spearman": float(spearmanr(yp, yt).statistic),
        "mae_units": float(np.median(np.abs(err))),
        "median_rel_err": float(np.median(np.abs(err) / np.clip(yt, 1, None))),
        "calib_slope": float(np.polyfit(yp, yt, 1)[0]),
        "frac_within_1u": float(np.mean(np.abs(err) <= 1)),
        "frac_within_2u": float(np.mean(np.abs(err) <= 2)),
    }
    pr, sp = _per_locus_corr(d, pred_col)
    if pr.size:
        out.update(
            median_perlocus_pearson=float(np.median(pr)),
            median_perlocus_spearman=float(np.median(sp)),
            frac_loci_pearson_gt_0p5=float(np.mean(pr > 0.5)),
            frac_loci_pearson_gt_0p7=float(np.mean(pr > 0.7)),
            n_loci_scored=int(pr.size),
        )
    if {"y_lo", "y_hi"}.issubset(d.columns):
        cov = np.mean((yt >= d["y_lo"].to_numpy()) & (yt <= d["y_hi"].to_numpy()))
        out["interval_coverage_90"] = float(cov)
    return out


def build(cfg) -> None:
    mdir = cfg.path("outputs.model_dir")
    pred = pl.read_parquet(mdir / "predictions.parquet")

    # naive baseline for context
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    if {"depth_ratio_d20_w2000", "cov_ref_units_est"}.issubset(tbl.columns):
        base = tbl.select("locus_id", "sample_id", "combined_units").with_columns(
            naive_pred=(tbl["depth_ratio_d20_w2000"] * 2.0 * tbl["cov_ref_units_est"]))
        pred = pred.join(base.select("locus_id", "sample_id", "naive_pred"),
                         on=["locus_id", "sample_id"], how="left")

    report: dict = {"quadrants": {}}
    for q in _QUADS:
        dq = pred.filter(pl.col("quadrant") == q)
        if dq.height < 5:
            continue
        entry = {"model": _block(dq, "y_pred")}
        if "naive_pred" in dq.columns:
            entry["naive_baseline"] = _block(dq, "naive_pred")
        if "transfer_pred" in dq.columns:
            entry["transfer_only"] = _block(dq, "transfer_pred")
        if "plc_pred" in dq.columns:
            entry["plc_only"] = _block(dq.filter(pl.col("plc_pred").is_not_null()), "plc_pred")
        pass_only = dq.filter(pl.col("tier") == "PASS")
        if pass_only.height >= 5:
            entry["model_PASS_tier"] = _block(pass_only, "y_pred")
            entry["PASS_fraction"] = float(pass_only.height / dq.height)
        report["quadrants"][q] = entry

    # segdup / size stratification on the hardest quadrant
    hard = pred.filter(pl.col("quadrant") == "test_novel_locus")
    if hard.height >= 20 and "cov_is_segdup" in hard.columns:
        report["test_novel_locus_by_segdup"] = {
            str(k): _block(hard.filter(pl.col("cov_is_segdup") == k), "y_pred")
            for k in (0, 1) if hard.filter(pl.col("cov_is_segdup") == k).height >= 5
        }

    (mdir / "evaluation.json").write_text(json.dumps(report, indent=2, default=float))
    _print(report)
    print(f"\nwrote {mdir/'evaluation.json'}")


def _print(report: dict) -> None:
    for q, entry in report["quadrants"].items():
        print(f"\n=== {q} ===")
        for name, m in entry.items():
            if not isinstance(m, dict):
                print(f"  {name}: {m:.3f}")
                continue
            if "pooled_pearson" not in m:
                print(f"  {name}: n={m.get('n')}")
                continue
            print(f"  {name:18s} n={m['n']:>7,} loci={m.get('n_loci','?'):>5} "
                  f"MAE={m['mae_units']:6.2f}u relErr={m['median_rel_err']:5.2f} "
                  f"slope={m['calib_slope']:+.2f} "
                  f"perLocus_r={m.get('median_perlocus_pearson', float('nan')):+.3f} "
                  f"r>0.5={m.get('frac_loci_pearson_gt_0p5', float('nan')):.2f} "
                  f"cov90={m.get('interval_coverage_90', float('nan')):.2f}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
