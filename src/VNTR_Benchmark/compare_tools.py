"""Head-to-head comparison of our tool, GangSTR and danbing-tk on the test
individuals of the current sample split, against vamos long-read truth.

Common scale: combined (both-allele) array length in bp.
  ours       y_pred x bp-per-repeat (known loci: panel mean bp / mean units,
             held-out individuals excluded; novel loci: BED motif period)
  GangSTR    2 x reference span + (a1 + a2 - 2 x ref copies) x period
             (genotypes read as change from the reference; whole-copy counts
             times period would drop partial/impure units)
  danbing-tk per-locus regression of length on depth-normalised k-mer counts,
             trained on the training individuals (fit_danbing_length.py)

A tool with no predictions yet (e.g. danbing-tk before its length model is fit)
is left out; rerun once its file exists and it is added everywhere.

Writes work/compare/tool_metrics.csv (one row per tool x stratum x scope) and
work/compare/calls.parquet (one row per test person x locus: truth, each tool's
call, strata) for make_compare_figures.py.

    python src/VNTR_Benchmark/compare_tools.py [--danbing PATH] [--out DIR]
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))
from pipeline.model.assemble import _clean_vamos_col  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402

PIPE = ROOT / "src" / "VNTR_AbsoluteGenotyping" / "work"
BENCH = ROOT / "src" / "VNTR_Benchmark" / "work"
OUT = BENCH / "compare"
TOOLS = ["ours", "gangstr", "danbing"]


def truth_and_ours(test: list[str]) -> pl.DataFrame:
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(load_config(None))}
    gt = (pl.read_parquet(PIPE / "ground_truth" / "gt_long.parquet",
                          columns=["locus_id", "vamos_sample", "combined_bp", "combined_units", "n_haps_called"])
          .filter(pl.col("vamos_sample").is_in([smap[s] for s in test]) & (pl.col("n_haps_called") == 2))
          .join(pl.DataFrame({"sample_id": test, "vamos_sample": [smap[s] for s in test]}), on="vamos_sample"))
    pred = pl.read_parquet(PIPE / "model" / "predictions.parquet").filter(pl.col("sample_id").is_in(test))
    tbl = pl.read_parquet(PIPE / "model" / "model_table.parquet",
                          columns=["locus_id", "sample_id", "cov_mean_combined_bp", "cov_mean_combined_units",
                                   "cov_panel_core90_units", "cov_panel_n_outliers", "cov_panel_p5_units",
                                   "cov_panel_p95_units", "cov_is_segdup"]).filter(pl.col("sample_id").is_in(test))
    d = (gt.join(pred.select("locus_id", "sample_id", "y_pred", "period", "path"), on=["locus_id", "sample_id"], how="left")
         .join(tbl, on=["locus_id", "sample_id"], how="left"))
    bpr = (pl.when((pl.col("path") == "known") & (pl.col("cov_mean_combined_units") > 0))
           .then(pl.col("cov_mean_combined_bp") / pl.col("cov_mean_combined_units"))
           .otherwise(pl.col("period").cast(pl.Float64)))
    return d.with_columns(ours_bp=pl.col("y_pred") * bpr)


def add_gangstr(d: pl.DataFrame) -> pl.DataFrame:
    idx = pl.read_csv(BENCH / "gangstr" / "regions_index.tsv", separator="\t").with_columns(
        span=pl.col("end") - pl.col("start") + 1)
    parts = []
    for f in sorted((BENCH / "gangstr" / "calls").glob("*.tsv")):
        g = pl.read_csv(f, separator="\t", infer_schema_length=0).filter(pl.col("a1_copies") != ".")
        parts.append(g.with_columns(sample_id=pl.lit(f.stem)))
    if not parts:
        return d.with_columns(gangstr_bp=pl.lit(None, dtype=pl.Float64), gangstr_q=pl.lit(None, dtype=pl.Float64))
    g = pl.concat(parts).with_columns(
        pl.col("a1_copies", "a2_copies", "period", "ref_copies", "q").cast(pl.Float64, strict=False))
    # haploid calls (male chrX/chrY) carry one allele
    allele = lambda c: pl.col("span") + (pl.col(c) - pl.col("ref_copies")) * pl.col("period")  # noqa: E731
    g = g.join(idx.select("locus_id", "span"), on="locus_id").with_columns(
        gangstr_bp=allele("a1_copies") + pl.when(pl.col("a2_copies").is_not_null())
        .then(allele("a2_copies")).otherwise(0), gangstr_q=pl.col("q"))
    return d.join(g.select("locus_id", "sample_id", "gangstr_bp", "gangstr_q"), on=["locus_id", "sample_id"], how="left")


def add_danbing(d: pl.DataFrame, f: Path) -> pl.DataFrame:
    if not f.exists():
        return d.with_columns(danbing_bp=pl.lit(None, dtype=pl.Float64))
    p = pl.read_parquet(f).select("locus_id", "sample_id", "danbing_bp")
    return d.join(p, on=["locus_id", "sample_id"], how="left")


def add_strata(d: pl.DataFrame) -> pl.DataFrame:
    span = pl.read_csv(BENCH / "gangstr" / "regions_index.tsv", separator="\t").select(
        "locus_id", ref_span=pl.col("end") - pl.col("start") + 1)
    core_bp = pl.read_parquet(PIPE / "model" / "panel_core_bp.parquet").select("locus_id", "panel_core90_bp")
    d = d.join(span, on="locus_id", how="left").join(core_bp, on="locus_id", how="left")
    w = 2 * pl.max_horizontal(pl.col("cov_panel_core90_units"), pl.lit(1))
    return d.with_columns(
        core_units_bin=pl.col("cov_panel_core90_units").cut([1, 2, 5, 20, 100], labels=["0", "1", "2-4", "5-19", "20-99", ">=100"], left_closed=True),
        core_bp_bin=pl.col("panel_core90_bp").cut([1, 25, 100, 300, 1000], labels=["0", "1-24", "25-99", "100-299", "300-999", ">=1000"], left_closed=True),
        motif_bin=pl.col("period").cut([11, 21, 51], labels=["7-10", "11-20", "21-50", ">50"], left_closed=True),
        span_bin=pl.col("ref_span").cut([150, 450, 1000], labels=["<150 (read)", "150-449 (fragment)", "450-999", ">=1000"], left_closed=True),
        segdup=pl.when(pl.col("cov_is_segdup") == 1).then(pl.lit("segdup")).otherwise(pl.lit("unique")),
        site_type=pl.when(pl.col("cov_panel_core90_units") <= 1)
        .then(pl.when(pl.col("cov_panel_n_outliers") > 0).then(pl.lit("stable core + outliers")).otherwise(pl.lit("stable")))
        .otherwise(pl.when(pl.col("cov_panel_n_outliers") > 0).then(pl.lit("broad + outliers")).otherwise(pl.lit("broad spread"))),
        rare_carrier=(pl.col("combined_units") > pl.col("cov_panel_p95_units") + w)
        | (pl.col("combined_units") < pl.col("cov_panel_p5_units") - w),
        rare_up=pl.col("combined_units") > pl.col("cov_panel_p95_units") + w)


def within_r(g: pl.DataFrame, c: str) -> float:
    w = g.with_columns(a=pl.col(c) - pl.col(c).mean().over("locus_id"),
                       b=pl.col("combined_bp") - pl.col("combined_bp").mean().over("locus_id"))
    a, b = w["a"].to_numpy(), w["b"].to_numpy()
    return float(pearsonr(a, b).statistic) if len(a) > 2 and a.std() > 0 and b.std() > 0 else np.nan


def metrics(g: pl.DataFrame, tool: str, n_total: int) -> dict:
    c = f"{tool}_bp"
    h = g.filter(pl.col(c).is_finite())
    if h.height < 20:
        return {"n_calls": h.height}
    y, x = h[c].to_numpy(), h["combined_bp"].to_numpy().astype(float)
    err = np.abs(y - x)
    rc = h.filter(pl.col("rare_carrier"))
    out = {"n_calls": h.height, "call_rate": h.height / max(n_total, 1),
           "mean_abs_err_bp": float(err.mean()), "median_abs_err_bp": float(np.median(err)),
           "median_rel_err": float(np.median(err / np.clip(x, 1, None))),
           "frac_within_5pct": float(np.mean(err / np.clip(x, 1, None) <= 0.05)),
           "within_locus_r": within_r(h, c)}
    if rc.height >= 10:
        up = rc["rare_up"].to_numpy()
        y0 = rc[c].to_numpy()
        ref = np.where(up, rc["cov_panel_p95_units"].to_numpy(), rc["cov_panel_p5_units"].to_numpy())
        bpr = rc["combined_bp"].to_numpy() / np.clip(rc["combined_units"].to_numpy(), 1, None)
        edge_bp = ref * bpr
        out["rare_n"] = rc.height
        out["rare_right_direction"] = float(np.mean(np.where(up, y0 > edge_bp, y0 < edge_bp)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--danbing", type=Path, default=BENCH / "danbing" / "length_predictions.parquet")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    sp = pl.read_parquet(PIPE / "model" / "sample_split.parquet")
    test = sorted(sp.filter(pl.col("fold") == "test")["sample_id"])
    d = add_strata(add_danbing(add_gangstr(truth_and_ours(test)), args.danbing))
    have = {t: d[f"{t}_bp"].is_finite().any() for t in TOOLS}
    samples_with = {t: d.filter(pl.col(f"{t}_bp").is_finite())["sample_id"].n_unique() for t in TOOLS}
    print("test individuals with calls:", samples_with)
    present = [t for t in TOOLS if have[t]]
    # restrict to individuals every present tool has run on, for a like-for-like table
    common_samples = set(test)
    for t in present:
        common_samples &= set(d.filter(pl.col(f"{t}_bp").is_finite())["sample_id"].unique())
    d = d.filter(pl.col("sample_id").is_in(list(common_samples)))
    all_called = pl.all_horizontal([pl.col(f"{t}_bp").is_finite() for t in present])
    rows = []
    strata = [("all", None), ("core_units_bin", "core_units_bin"), ("core_bp_bin", "core_bp_bin"),
              ("site_type", "site_type"), ("motif_bin", "motif_bin"), ("span_bin", "span_bin"),
              ("segdup", "segdup"), ("path", "path")]
    for name, col in strata:
        groups = [("all", d)] if col is None else [(str(k[0]), g) for k, g in d.group_by([col], maintain_order=True)]
        for lab, g in groups:
            if lab in ("None", "null"):
                continue
            for scope, gg in (("own_calls", g), ("all_tools_called", g.filter(all_called))):
                for t in present:
                    m = metrics(gg, t, g.height)
                    rows.append({"stratum": name, "group": lab, "scope": scope, "tool": t,
                                 "n_sites": gg["locus_id"].n_unique(), "n_people": len(common_samples), **m})
    args.out.mkdir(parents=True, exist_ok=True)
    res = pl.DataFrame(rows, infer_schema_length=None)
    res.write_csv(args.out / "tool_metrics.csv")
    d.select("locus_id", "sample_id", "combined_bp", "combined_units", "path", "period", "ref_span",
             *[f"{t}_bp" for t in present], "core_units_bin", "core_bp_bin", "motif_bin", "span_bin",
             "segdup", "site_type", "rare_carrier", "rare_up", "cov_panel_p5_units", "cov_panel_p95_units"
             ).write_parquet(args.out / "calls.parquet")
    show = res.filter((pl.col("stratum") == "all"))
    print(show.select("scope", "tool", "n_calls", "call_rate", "mean_abs_err_bp", "median_rel_err",
                      "frac_within_5pct", "within_locus_r", "rare_right_direction"))
    print(f"wrote {args.out / 'tool_metrics.csv'} ({res.height} rows; {len(common_samples)} individuals; "
          f"tools: {', '.join(present)})")


if __name__ == "__main__":
    main()
