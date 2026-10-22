"""
per_allele.py — opt-in per-allele estimate, built on the combined prediction.

Why not the raw spanning-pair split from extraction (pa_units_a1/a2)?  It
always 2-means-splits the spanning inserts, so a homozygous locus gets two
alleles ~1.6 insert-SDs apart (~80 bp at 30x/150bp), and its confidence
score does not separate het from hom.  It also mapped the insert shift to an
allele with the wrong sign: TLEN is in reference coordinates, so an allele
D bp longer than the reference shortens the observed insert by D.  On the
real panel ~85% of (locus, sample) pairs are homozygous and the median het
allele difference is ~22 bp — far below the insert-size spread — so short
reads can only resolve large allele differences.

Estimator (sign-free; uses only the 2-means gap):
  sigma  = 1.4826 * sample insert MAD
  z      = |shift_a2 - shift_a1| / sigma        (~1.6 for a single allele)
  het    if z > z_het (tuned on train) and >= min_pairs spanning pairs
  d_bp   = sqrt(max(0, gap^2 - (1.596 sigma)^2))  (deconvolved allele gap)
  alleles = (y_pred -/+ d_units) / 2   (hom: y_pred / 2 each)
centred on the combined y_pred from global_residual.

Output: work/model/per_allele_eval.json, work/model/per_allele_predictions.parquet
"""

from __future__ import annotations

import argparse
import json

import numpy as np
import polars as pl

from ..paths import load_config

_HOM_GAP_SD = 2 * np.sqrt(2 / np.pi)      # E|m2 - m1| / sigma for a 2-means split of N(0, sigma)


def _alleles(d: pl.DataFrame, z_het: float, min_pairs: int) -> pl.DataFrame:
    sigma = 1.4826 * d["nuis_insert_mad"].to_numpy().astype(float)
    gap = np.abs(d["pa_bp_shift_a2"].to_numpy() - d["pa_bp_shift_a1"].to_numpy())
    z = gap / np.clip(sigma, 1e-6, None)
    n = d["pa_spanning_n"].to_numpy()
    het = np.isfinite(z) & (z > z_het) & (n >= min_pairs)
    d_bp = np.sqrt(np.clip(gap ** 2 - (_HOM_GAP_SD * sigma) ** 2, 0, None))
    period = np.clip(d["period"].to_numpy().astype(float), 1, None)
    d_units = np.where(het, d_bp / period, 0.0)
    total = d["y_pred"].to_numpy()
    lo = np.clip((total - d_units) / 2, 0, None)
    hi = np.clip((total + d_units) / 2, 0, None)
    return d.with_columns(pa_het_call=pl.Series(het), pa_gap_z=pl.Series(z),
                          pa_lo=pl.Series(lo), pa_hi=pl.Series(hi))


def _score(d: pl.DataFrame, tol: float) -> dict:
    tl, th = d["hap_lo_units"].to_numpy(), d["hap_hi_units"].to_numpy()
    lo, hi = d["pa_lo"].to_numpy(), d["pa_hi"].to_numpy()
    ok = (np.abs(lo - tl) <= tol) & (np.abs(hi - th) <= tol)
    true_het = th != tl
    call_het = d["pa_het_call"].to_numpy()
    return {
        "n": int(d.height),
        "both_within_tol": float(ok.mean()),
        "median_err_lo": float(np.median(np.abs(lo - tl))),
        "median_err_hi": float(np.median(np.abs(hi - th))),
        "frac_true_het": float(true_het.mean()),
        "het_call_rate": float(call_het.mean()),
        "het_precision": float(true_het[call_het].mean()) if call_het.any() else None,
        "het_recall": float(call_het[true_het].mean()) if true_het.any() else None,
        "both_within_tol_true_het": float(ok[true_het].mean()) if true_het.any() else None,
    }


def build(cfg, tol_units: float = 2.0, min_pairs: int = 10) -> None:
    mdir = cfg.path("outputs.model_dir")
    pred = pl.read_parquet(mdir / "predictions.parquet").select(
        "locus_id", "sample_id", "quadrant", "path", "tier", "y_pred", "period")
    tbl = pl.read_parquet(mdir / "model_table.parquet", columns=[
        "locus_id", "sample_id", "pa_bp_shift_a1", "pa_bp_shift_a2", "pa_spanning_n",
        "nuis_insert_mad", "hap_lo_units", "hap_hi_units"])
    d = pred.join(tbl, on=["locus_id", "sample_id"], how="inner").with_columns(
        pa_spanning_n=pl.col("pa_spanning_n").fill_null(0),
        pa_bp_shift_a1=pl.col("pa_bp_shift_a1").fill_nan(None).fill_null(0.0),
        pa_bp_shift_a2=pl.col("pa_bp_shift_a2").fill_nan(None).fill_null(0.0))
    d = d.filter(pl.col("hap_lo_units").is_not_null() & pl.col("y_pred").is_finite())

    # tune the het threshold on train rows; z_het = inf means "always homozygous"
    tr = d.filter(pl.col("quadrant") == "train")
    grid = [np.inf] + [float(x) for x in np.arange(2.0, 6.01, 0.25)]
    scores = {z: _score(_alleles(tr, z, min_pairs), tol_units)["both_within_tol"] for z in grid}
    z_het = max(scores, key=scores.get)
    print("per-allele: train both-within-tol by z_het: " +
          ", ".join(f"{k:g}:{v:.3f}" for k, v in scores.items()))
    print(f"per-allele: chosen z_het = {z_het:g}")

    d = _alleles(d, z_het, min_pairs)
    report = {"tol_units": tol_units, "z_het": z_het, "min_pairs": min_pairs,
              "quadrants": {}}
    for q in ["train", "val", "test_novel_sample", "test_novel_locus", "test_novel_both"]:
        s = d.filter(pl.col("quadrant") == q)
        if s.height < 10:
            continue
        entry = {"all": _score(s, tol_units)}
        p = s.filter(pl.col("tier") == "PASS")
        if p.height >= 10:
            entry["PASS_tier"] = _score(p, tol_units)
        report["quadrants"][q] = entry

    d.select("locus_id", "sample_id", "quadrant", "path", "tier", "pa_het_call",
             "pa_gap_z", "pa_lo", "pa_hi").write_parquet(mdir / "per_allele_predictions.parquet")
    (mdir / "per_allele_eval.json").write_text(json.dumps(report, indent=2, default=float))
    for q, e in report["quadrants"].items():
        a = e["all"]
        print(f"  {q:18s} both±{tol_units:g}u={a['both_within_tol']:.3f}  "
              f"true-het={a['frac_true_het']:.3f}  het-call={a['het_call_rate']:.3f}  "
              f"het-prec={a['het_precision'] if a['het_precision'] is not None else float('nan'):.3f}  "
              f"het-recall={a['het_recall'] if a['het_recall'] is not None else float('nan'):.3f}")
    print(f"wrote {mdir/'per_allele_eval.json'}, {mdir/'per_allele_predictions.parquet'}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--tol-units", type=float, default=2.0)
    args = ap.parse_args(argv)
    build(load_config(args.config), tol_units=args.tol_units)


if __name__ == "__main__":
    main()
