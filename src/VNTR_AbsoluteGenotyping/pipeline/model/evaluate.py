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
    """Per-locus Pearson/Spearman over loci whose TRUTH varies.

    A prediction that is constant within such a locus carries no ranking
    information and scores r = 0 (rather than being dropped, which would let
    methods be scored on different, easier locus sets).
    """
    prs, sps = [], []
    for (_l,), g in df.group_by(["locus_id"]):
        if g.height < min_n:
            continue
        yt = g["combined_units"].to_numpy()
        yp = g[col_pred].to_numpy()
        if np.std(yt) < 1e-9:
            continue
        if np.std(yp) < 1e-9:
            prs.append(0.0)
            sps.append(0.0)
            continue
        prs.append(pearsonr(yp, yt).statistic)
        sps.append(spearmanr(yp, yt).statistic)
    return np.array(prs), np.array(sps)


def _within_locus_r(d: pl.DataFrame, col_pred: str) -> float:
    """Pearson r of locus-demeaned prediction vs locus-demeaned truth, pooled
    over all rows.  Stable with few samples per locus, unlike per-locus r."""
    w = d.with_columns(
        _yt=pl.col("combined_units") - pl.col("combined_units").mean().over("locus_id"),
        _yp=pl.col(col_pred) - pl.col(col_pred).mean().over("locus_id"))
    a, b = w["_yp"].to_numpy(), w["_yt"].to_numpy()
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    return float(pearsonr(a, b).statistic)


def _block(df: pl.DataFrame, pred_col: str) -> dict:
    d = df.filter(pl.col(pred_col).is_finite() & pl.col("combined_units").is_not_null())
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
        "mean_abs_err_units": float(np.mean(np.abs(err))),
        "median_rel_err": float(np.median(np.abs(err) / np.clip(yt, 1, None))),
        "calib_slope": float(np.polyfit(yp, yt, 1)[0]),
        "frac_within_1u": float(np.mean(np.abs(err) <= 1)),
        "frac_within_2u": float(np.mean(np.abs(err) <= 2)),
    }
    out["within_locus_r"] = _within_locus_r(d, pred_col)
    pr, sp = _per_locus_corr(d, pred_col)
    if pr.size:
        out.update(
            median_perlocus_pearson=float(np.median(pr)),
            median_perlocus_spearman=float(np.median(sp)),
            frac_loci_pearson_gt_0p5=float(np.mean(pr > 0.5)),
            frac_loci_pearson_gt_0p7=float(np.mean(pr > 0.7)),
            n_loci_scored=int(pr.size),
        )
    if pred_col == "y_pred" and {"y_lo", "y_hi"}.issubset(d.columns):
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
    # panel-mean baseline: the locus's long-read panel mean, no reads at all.
    # Only defined for KNOWN loci (a novel locus has no panel). It has no
    # per-locus correlation (constant within locus); compare it on error.
    if "cov_mean_combined_units" in pred.columns:
        pred = pred.with_columns(panel_mean_pred=pl.when(pl.col("path") == "known")
                                 .then(pl.col("cov_mean_combined_units")))

    report: dict = {"quadrants": {}}
    for q in _QUADS:
        dq = pred.filter(pl.col("quadrant") == q)
        if dq.height < 5:
            continue
        entry = {"model": _block(dq, "y_pred")}
        if "naive_pred" in dq.columns:
            entry["naive_baseline"] = _block(dq, "naive_pred")
        if "panel_mean_pred" in dq.columns and dq["panel_mean_pred"].is_not_null().any():
            entry["panel_mean_baseline"] = _block(dq, "panel_mean_pred")
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

    report["by_site_variability"] = _by_variability(pred)
    report["outlier_people"] = _outlier_people(pred)

    (mdir / "evaluation.json").write_text(json.dumps(report, indent=2, default=float))
    _print(report)
    _print_variability(report)
    print(f"\nwrote {mdir/'evaluation.json'}")


def _short_block(d: pl.DataFrame, baseline: str | None) -> dict:
    yt = d["combined_units"].to_numpy()
    yp = d["y_pred"].to_numpy()
    out = {"n": int(d.height), "n_loci": int(d["locus_id"].n_unique()),
           "mean_abs_err": float(np.mean(np.abs(yp - yt))),
           "frac_within_1u": float(np.mean(np.abs(yp - yt) <= 1)),
           "median_rel_err": float(np.median(np.abs(yp - yt) / np.clip(yt, 1, None))),
           "within_locus_r": _within_locus_r(d, "y_pred"),
           "interval_coverage_90": float(np.mean((yt >= d["y_lo"].to_numpy())
                                                 & (yt <= d["y_hi"].to_numpy()))),
           "median_interval_width": float(np.median(d["y_hi"].to_numpy() - d["y_lo"].to_numpy()))}
    pr, _ = _per_locus_corr(d, "y_pred")
    out["median_perlocus_pearson"] = float(np.median(pr)) if pr.size else None
    if baseline is not None:
        b = d[baseline].to_numpy()
        ok = np.isfinite(b)
        out["baseline"] = baseline
        out["baseline_mean_abs_err"] = float(np.mean(np.abs(b[ok] - yt[ok]))) if ok.any() else None
    return out


def _by_variability(pred: pl.DataFrame) -> dict:
    """Test performance by locus core spread (p95 - p5 of the long-read panel
    without held-out people) and by outlier status.  Panel stats are used here
    only to stratify the evaluation, for both paths."""
    if "cov_panel_core90_units" not in pred.columns:
        return {}
    edges, labels = [1, 2, 5, 20, 100], ["0", "1", "2-4", "5-19", "20-99", ">=100"]
    p = pred.with_columns(
        core_bin=pl.col("cov_panel_core90_units").cut(edges, labels=labels, left_closed=True),
        site_type=pl.when(pl.col("cov_panel_core90_units") <= 1)
        .then(pl.when(pl.col("cov_panel_n_outliers") > 0).then(pl.lit("stable core + outliers"))
              .otherwise(pl.lit("stable")))
        .otherwise(pl.when(pl.col("cov_panel_n_outliers") > 0).then(pl.lit("broad spread + outliers"))
                   .otherwise(pl.lit("broad spread"))))
    rep: dict = {}
    for q, base in (("test_novel_sample", "cov_mean_combined_units"),
                    ("test_novel_locus", "naive_pred")):
        d = p.filter(pl.col("quadrant") == q)
        if d.height < 50:
            continue
        b = base if base in d.columns else None
        rep[q] = {
            "by_core_spread": {lab: _short_block(d.filter(pl.col("core_bin") == lab), b)
                               for lab in labels
                               if d.filter(pl.col("core_bin") == lab).height >= 50},
            "by_site_type": {t: _short_block(d.filter(pl.col("site_type") == t), b)
                             for t in ("stable", "broad spread", "stable core + outliers",
                                       "broad spread + outliers")
                             if d.filter(pl.col("site_type") == t).height >= 50},
        }
    return rep


def _outlier_people(pred: pl.DataFrame) -> dict:
    """Can the tool tell when a held-out person carries a rare extreme allele?

    An outlier person is beyond the panel core (p5..p95, panel excludes them)
    by > 2 x max(core, 1) units.  For them:
      right_direction  rounded prediction is outside the core on the true side
      halfway          prediction moved >= half the distance from the core edge
                       to the true value
      called_outlier   prediction is itself beyond the outlier margin
    false_alarm_rate_non_outliers: non-outlier people whose prediction is
    beyond the outlier margin (same threshold as called_outlier).
    """
    need = {"cov_panel_p5_units", "cov_panel_p95_units", "cov_panel_core90_units"}
    if not need.issubset(pred.columns):
        return {}
    w = 2 * pl.max_horizontal(pl.col("cov_panel_core90_units"), pl.lit(1))
    p = pred.with_columns(
        hi_out=pl.col("combined_units") > pl.col("cov_panel_p95_units") + w,
        lo_out=pl.col("combined_units") < pl.col("cov_panel_p5_units") - w,
        pred_hi=pl.col("y_pred").round() > pl.col("cov_panel_p95_units"),
        pred_lo=pl.col("y_pred").round() < pl.col("cov_panel_p5_units"),
        half_hi=pl.col("y_pred") >= pl.col("cov_panel_p95_units")
        + 0.5 * (pl.col("combined_units") - pl.col("cov_panel_p95_units")),
        half_lo=pl.col("y_pred") <= pl.col("cov_panel_p5_units")
        - 0.5 * (pl.col("cov_panel_p5_units") - pl.col("combined_units")),
        pred_far_hi=pl.col("y_pred") > pl.col("cov_panel_p95_units") + w,
        pred_far_lo=pl.col("y_pred") < pl.col("cov_panel_p5_units") - w)
    p = p.with_columns(pred_far=pl.col("pred_far_hi") | pl.col("pred_far_lo"))
    rep = {}
    for q in ("test_novel_sample", "test_novel_locus"):
        d = p.filter(pl.col("quadrant") == q)
        o = d.filter(pl.col("hi_out") | pl.col("lo_out"))
        n = d.filter(~(pl.col("hi_out") | pl.col("lo_out")))
        if o.height < 10:
            continue
        side = lambda hi, lo: ((o["hi_out"] & o[hi]) | (o["lo_out"] & o[lo])).to_numpy()
        yt, yp = o["combined_units"].to_numpy(), o["y_pred"].to_numpy()
        rep[q] = {
            "n_outlier_person_sites": int(o.height),
            "frac_of_rows": float(o.height / d.height),
            "right_direction": float(side("pred_hi", "pred_lo").mean()),
            "halfway": float(side("half_hi", "half_lo").mean()),
            "called_outlier": float(side("pred_far_hi", "pred_far_lo").mean()),
            "mean_abs_err": float(np.mean(np.abs(yp - yt))),
            "interval_coverage_90": float(np.mean((yt >= o["y_lo"].to_numpy())
                                                  & (yt <= o["y_hi"].to_numpy()))),
            "false_alarm_rate_non_outliers": float(n["pred_far"].mean()),
        }
        if "cov_mean_combined_units" in o.columns and o["cov_mean_combined_units"].is_not_null().any():
            rep[q]["panel_mean_abs_err"] = float(np.nanmean(
                np.abs(o["cov_mean_combined_units"].to_numpy().astype(float) - yt)))
    return rep


def _print_variability(report: dict) -> None:
    for q, e in report.get("by_site_variability", {}).items():
        for key, title in (("by_core_spread", "core spread (units)"), ("by_site_type", "site type")):
            print(f"\n=== {q}: by {title} ===")
            for lab, m in e[key].items():
                b = m.get("baseline_mean_abs_err")
                print(f"  {lab:24s} loci={m['n_loci']:>6,} meanAE={m['mean_abs_err']:7.2f} "
                      f"baseline={b if b is not None else float('nan'):7.2f} "
                      f"±1u={m['frac_within_1u']:.2f} relErr={m['median_rel_err']:.3f} "
                      f"siteR={m['median_perlocus_pearson'] if m['median_perlocus_pearson'] is not None else float('nan'):+.2f} "
                      f"cov90={m['interval_coverage_90']:.2f} width={m['median_interval_width']:.0f}")
    for q, m in report.get("outlier_people", {}).items():
        print(f"\n=== {q}: outlier people (rare extreme alleles) ===")
        print("  " + "  ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                               for k, v in m.items()))


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
                  f"medAE={m['mae_units']:5.2f}u meanAE={m['mean_abs_err_units']:6.2f}u "
                  f"relErr={m['median_rel_err']:5.2f} "
                  f"slope={m['calib_slope']:+.2f} "
                  f"perLocus_r={m.get('median_perlocus_pearson', float('nan')):+.3f} "
                  f"withinLocus_r={m.get('within_locus_r', float('nan')):+.3f} "
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
