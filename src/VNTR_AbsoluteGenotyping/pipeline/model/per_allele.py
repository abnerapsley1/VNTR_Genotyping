"""
per_allele.py — evaluate / gate the opt-in per-allele estimate.

Signal extraction already emits a spanning-pair mixture estimate
(pa_units_a1 <= pa_units_a2 after ordering, pa_confidence, pa_mix_separation).
This module:
  * orders the estimate against ground-truth hap_lo/hap_hi,
  * picks a confidence threshold on the train quadrant that maximises
    "both alleles within tolerance" while keeping >= min call rate,
  * reports call rate + accuracy per quadrant at that threshold.

Per-allele is only expected to fire for arrays shorter than ~1 fragment with
two distinguishable alleles; low call rates are normal and acceptable.

Output: work/model/per_allele_eval.json
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import polars as pl

from ..paths import load_config


def build(cfg, tol_units: float = 2.0, min_call_rate: float = 0.02) -> None:
    mdir = cfg.path("outputs.model_dir")
    tbl = pl.read_parquet(mdir / "model_table.parquet")
    need = {"pa_units_a1", "pa_units_a2", "pa_confidence",
            "hap_lo_units", "hap_hi_units", "quadrant"}
    if not need.issubset(tbl.columns):
        raise SystemExit(f"model_table missing per-allele columns: {need - set(tbl.columns)}")

    d = tbl.filter(pl.col("pa_units_a1").is_not_null()
                   & pl.col("hap_lo_units").is_not_null()).with_columns(
        pa_lo=pl.min_horizontal("pa_units_a1", "pa_units_a2"),
        pa_hi=pl.max_horizontal("pa_units_a1", "pa_units_a2"))
    if d.height < 20:
        (mdir / "per_allele_eval.json").write_text(json.dumps(
            {"note": "insufficient per-allele observations", "n": d.height}, indent=2))
        print("per-allele: insufficient observations")
        return

    d = d.with_columns(
        err_lo=(pl.col("pa_lo") - pl.col("hap_lo_units")).abs(),
        err_hi=(pl.col("pa_hi") - pl.col("hap_hi_units")).abs())
    d = d.with_columns(both_ok=((pl.col("err_lo") <= tol_units)
                                & (pl.col("err_hi") <= tol_units)).cast(pl.Int8))

    tr = d.filter(pl.col("quadrant") == "train")
    thr = _pick_threshold(tr, min_call_rate) if tr.height >= 20 else 0.3

    report = {"tol_units": tol_units, "confidence_threshold": thr, "quadrants": {}}
    for q in ["train", "test_novel_sample", "test_novel_locus", "test_novel_both"]:
        s = d.filter(pl.col("quadrant") == q)
        if s.height < 10:
            continue
        called = s.filter(pl.col("pa_confidence") >= thr)
        report["quadrants"][q] = {
            "n_with_estimate": int(s.height),
            "call_rate": float(called.height / s.height),
            "both_within_tol": float(called["both_ok"].mean()) if called.height else None,
            "median_err_lo": float(called["err_lo"].median()) if called.height else None,
            "median_err_hi": float(called["err_hi"].median()) if called.height else None,
        }

    (mdir / "per_allele_eval.json").write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=float))
    print(f"wrote {mdir/'per_allele_eval.json'}")


def _pick_threshold(tr: pl.DataFrame, min_call_rate: float) -> float:
    best = (0.0, -1.0)
    n = tr.height
    for thr in np.linspace(0.05, 0.95, 19):
        c = tr.filter(pl.col("pa_confidence") >= thr)
        if c.height / n < min_call_rate or c.height < 10:
            continue
        acc = float(c["both_ok"].mean())
        if acc > best[1]:
            best = (float(thr), acc)
    return best[0] or 0.3


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--tol-units", type=float, default=2.0)
    args = ap.parse_args(argv)
    build(load_config(args.config), tol_units=args.tol_units)


if __name__ == "__main__":
    main()
