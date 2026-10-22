"""Performance reports for built model sets (per repeat TSV + one-page HTML summary)."""
from __future__ import annotations

import html as _html
from pathlib import Path

import numpy as np
import polars as pl


def per_repeat(rows: pl.DataFrame) -> pl.DataFrame:
    """Per-repeat metrics from cross-validated estimates (est_units, ci90_*, combined_units)."""
    err = (pl.col("est_units") - pl.col("combined_units")).abs()
    cov = (pl.col("combined_units") >= pl.col("ci90_lo_units")) & (pl.col("combined_units") <= pl.col("ci90_hi_units"))
    g = rows.filter(pl.col("est_units").is_finite()).group_by("locus_id").agg(
        perf_n_people=pl.len().cast(pl.Float64),
        perf_mean_abs_err_units=err.mean(),
        perf_median_abs_err_units=err.median(),
        perf_within_1_unit=(err <= 1).mean(),
        perf_ci90_coverage=cov.mean(),
        perf_mean_ci90_width_units=(pl.col("ci90_hi_units") - pl.col("ci90_lo_units")).mean(),
        perf_truth_sd_units=pl.col("combined_units").std(),
        perf_pearson_r=pl.corr("est_units", "combined_units"),
    )
    # r is undefined where the truth (or the estimate) does not vary across people
    return g.with_columns(perf_pearson_r=pl.when(pl.col("perf_truth_sd_units") > 0)
                          .then(pl.col("perf_pearson_r").fill_nan(0.0)))


def html(perf: pl.DataFrame, rows: pl.DataFrame, path: Path, name: str, retrained: bool) -> None:
    err = np.abs(rows["est_units"].to_numpy() - rows["combined_units"].to_numpy())
    cov = ((rows["combined_units"] >= rows["ci90_lo_units"]) & (rows["combined_units"] <= rows["ci90_hi_units"])).mean()
    q = perf["perf_mean_abs_err_units"].drop_nulls()
    worst = perf.sort("perf_mean_abs_err_units", descending=True, nulls_last=True).head(25)
    def tr(r):
        return "<tr>" + "".join(f"<td>{_html.escape(str(r[c] if not isinstance(r[c], float) else round(r[c], 3)))}</td>"
                                for c in ("locus_id", "n_train_people", "curve_signal", "perf_mean_abs_err_units",
                                          "perf_within_1_unit", "perf_pearson_r", "perf_ci90_coverage")) + "</tr>"
    body = f"""<!doctype html><html><head><meta charset="utf-8"><title>vntrlen build: {_html.escape(name)}</title>
<style>body{{font-family:system-ui,sans-serif;max-width:960px;margin:24px auto;padding:0 16px;color:#1b1b1b;background:#fcfcfb}}
table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{border-bottom:1px solid #e4e3df;padding:4px 6px;text-align:left}}
.k{{display:inline-block;margin:0 24px 12px 0}}.k b{{display:block;font-size:22px}}</style></head><body>
<h1>Model set: {_html.escape(name)}</h1>
<p>All numbers are 5-fold cross-validated over people (families kept together): each person is scored by
models and interval calibration that did not see them. Correction models were
{'retrained on your data' if retrained else 'reused from the shipped models; intervals recalibrated on your data'}.</p>
<div class="k"><b>{perf.height:,}</b>repeats</div>
<div class="k"><b>{rows['sample_id'].n_unique():,}</b>people</div>
<div class="k"><b>{err.mean():.2f}</b>mean abs. error (units)</div>
<div class="k"><b>{np.mean(err <= 1):.1%}</b>within 1 unit</div>
<div class="k"><b>{cov:.1%}</b>90% interval coverage</div>
<p>Per-repeat mean absolute error: median {q.median():.2f}, 90th percentile {q.quantile(0.9):.2f} units.</p>
<h2>25 repeats with the largest error</h2>
<table><tr><th>repeat</th><th>people</th><th>curve signal</th><th>mean abs. error</th><th>within 1</th>
<th>r</th><th>90% coverage</th></tr>{''.join(tr(r) for r in worst.iter_rows(named=True))}</table>
<p>Full table: performance_per_repeat.tsv</p></body></html>"""
    path.write_text(body)
