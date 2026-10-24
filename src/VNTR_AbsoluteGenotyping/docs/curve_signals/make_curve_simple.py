"""Simplified calibration-curve figures for one repeat: the whole curve, and a zoom on the
panel's core spread (5th-95th percentile of combined length).

    python docs/curve_signals/make_curve_simple.py [locus_id]
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import polars as pl

sys.path.insert(0, "/data/GitHub_Repos/VNTR_Genotyping/src/vntrlen")
from vntrlen import curves  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

D = Path(__file__).resolve().parent
W = D.parents[1] / "work"
BUNDLE = Path("/data/GitHub_Repos/VNTR_Genotyping/data/vntrlen_bundle/52ec291b8f70/shipped_GRCh38_v1/loci.parquet")
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
XLABEL = {"ridge": "Combined short-read signal score","depth_ratio_d0_w2000": "Depth ratio (all reads)", "depth_ratio_d20_w5000": "Depth ratio (MAPQ ≥ 20)",
          "irr_count_norm": "In-repeat reads (normalised)"}


def pretty(lid: str) -> str:
    c, r = lid.split(":")
    a, b = r.split("-")
    return f"{c}:{int(a):,}–{int(b):,}"


def draw(ax, tr, te, sig, p, xlim, ylim, legend):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=10)
    x0, x1 = p["x_min"], p["x_max"]
    lin = lambda x: p["lin_slope"] * x + p["lin_intercept"]  # noqa: E731
    ax.scatter(tr[sig], tr["combined_units"], s=20, color="#a9a8a3", lw=0, zorder=2, label="Training")
    ax.scatter(te[sig], te["combined_units"], s=36, facecolor="none", edgecolor=AQUA, lw=1.5, zorder=5, label="Test")
    xg = np.linspace(max(x0, xlim[0]), min(x1, xlim[1]), 2000)
    ax.plot(xg, curves.predict_one(p, xg), color=BLUE, lw=2.2, zorder=4, label="Isotonic")
    for a, b in ((xlim[0], x0), (x1, xlim[1])):
        if b > a:
            xs = np.array([a, b])
            ax.plot(xs, lin(xs), color=ORANGE, lw=2.2, ls="--", zorder=4,
                    label="Theil–Sen" if a == xlim[0] or x0 <= xlim[0] else None)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(XLABEL.get(sig, sig), fontsize=11, color=INK2)
    ax.set_ylabel("Combined length (repeat units)", fontsize=11, color=INK2)
    if legend:
        h, l = ax.get_legend_handles_labels()
        seen = {}
        for hh, ll in zip(h, l):
            seen.setdefault(ll, hh)
        ax.legend(seen.values(), seen.keys(), frameon=False, fontsize=10, loc="upper left")


def main(lid="chr19:56744829-56745871"):
    L = pl.read_parquet(BUNDLE).filter(pl.col("locus_id") == lid).row(0, named=True)
    sig, p = json.loads(L["curves"])["pe"]
    sigcols = p["signals"] if sig == "ridge" else [sig]
    t = pl.read_parquet(W / "model" / "model_table.parquet", columns=["locus_id", "sample_fold", "combined_units", *sigcols]) \
        .filter(pl.col("locus_id") == lid)
    if sig == "ridge":                      # x axis: the per-repeat combined signal score
        X = t.select(p["signals"]).cast(pl.Float64).to_numpy()
        t = t.with_columns(ridge=pl.Series(curves._score(p, X)))
    tr, te = t.filter(pl.col("sample_fold") == "train"), t.filter(pl.col("sample_fold") == "test")
    stem = lid.replace(":", "_").replace("-", "_")

    # whole curve
    x0, x1 = p["x_min"], p["x_max"]
    pad = 0.25 * (x1 - x0)
    fig, ax = plt.subplots(figsize=(7, 5.2), facecolor=SURF)
    ymax = max(tr["combined_units"].max(), te["combined_units"].max())
    lin_hi = p["lin_slope"] * (x1 + pad) + p["lin_intercept"]
    draw(ax, tr, te, sig, p, (max(0.0, x0 - pad), x1 + pad), (0, max(ymax, lin_hi) * 1.04), legend=True)
    ax.set_title(pretty(lid), fontsize=12, color=INK, loc="left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(D / f"curve_simple_full_{stem}.{ext}", dpi=200, facecolor=SURF)
    plt.close(fig)

    # zoom on the core spread (panel 5th-95th percentile of combined length)
    lo, hi = L["cov_panel_p5_units"], L["cov_panel_p95_units"]
    m = (lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    inside = pl.concat([tr, te]).filter(pl.col("combined_units").is_between(*m))
    xs = inside[sig].drop_nulls()
    xm = (xs.min() - 0.05 * (xs.max() - xs.min()), xs.max() + 0.05 * (xs.max() - xs.min()))
    fig, ax = plt.subplots(figsize=(7, 5.2), facecolor=SURF)
    draw(ax, tr, te, sig, p, xm, m, legend=True)
    ax.axhspan(lo, hi, color="#eef4fc", zorder=0)
    ax.set_title(f"{pretty(lid)}: core spread", fontsize=12, color=INK, loc="left")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(D / f"curve_simple_core_{stem}.{ext}", dpi=200, facecolor=SURF)
    plt.close(fig)
    print(f"wrote curve_simple_full_{stem} and curve_simple_core_{stem} (core {lo:.0f}-{hi:.0f} units)")


if __name__ == "__main__":
    main(*sys.argv[1:])
