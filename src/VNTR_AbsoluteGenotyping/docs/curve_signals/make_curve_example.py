"""Example calibration curve for one highly variable, well-fitting repeat: training people,
isotonic fit within their signal range, Theil-Sen line beyond both edges, test people.

    python docs/curve_signals/make_curve_example.py [locus_id]
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
LABEL = {"depth_ratio_d20_w5000": "depth ratio (MAPQ >= 20, 5 kb flanks)",
         "depth_ratio_d0_w2000": "depth ratio (all reads, 2 kb flanks)",
         "irr_count_norm": "in-repeat reads per 1x flank depth"}


def main(lid="chr5:135778989-135785061"):
    L = pl.read_parquet(BUNDLE).filter(pl.col("locus_id") == lid).row(0, named=True)
    sig, p = json.loads(L["curves"])["pe"]
    t = pl.read_parquet(W / "model" / "model_table.parquet", columns=["locus_id", "sample_id", "sample_fold",
                                                                     "combined_units", sig]).filter(pl.col("locus_id") == lid)
    tr, te = t.filter(pl.col("sample_fold") == "train"), t.filter(pl.col("sample_fold") == "test")
    xi, yi = np.array(p["iso_x"]), np.array(p["iso_y"])
    x0, x1 = p["x_min"], p["x_max"]
    pad = 0.35 * (x1 - x0)
    lpad = min(pad, x0) if x0 > 0 else pad          # signals here are ratios/counts: no negative axis
    lin = lambda x: p["lin_slope"] * x + p["lin_intercept"]  # noqa: E731

    fig, ax = plt.subplots(figsize=(9.5, 6.2), facecolor=SURF)
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.grid(color=GRID, lw=0.6)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.axvspan(x0 - lpad, x0, color="#f3f2ee", zorder=0)
    ax.axvspan(x1, x1 + pad, color="#f3f2ee", zorder=0)
    ax.scatter(tr[sig], tr["combined_units"], s=22, color="#9b9a95", alpha=0.75, lw=0, zorder=2,
               label=f"training people (n={tr.height})")
    xg = np.linspace(x0, x1, 2000)               # the tool's own prediction function (linear between knots)
    ax.plot(xg, curves.predict_one(p, xg), color=BLUE, lw=2.2, zorder=4, label="isotonic fit (within training range)")
    xl, xr = np.array([x0 - lpad, x0]), np.array([x1, x1 + pad])
    ax.plot(xl, lin(xl), color=ORANGE, lw=2.2, ls="--", zorder=4, label="Theil-Sen line (beyond the edges)")
    ax.plot(xr, lin(xr), color=ORANGE, lw=2.2, ls="--", zorder=4)
    xs = np.linspace(x0, x1, 50)
    ax.plot(xs, lin(xs), color=ORANGE, lw=1, ls=":", alpha=0.7, zorder=3, label="Theil-Sen line within range (not used)")
    ax.scatter(te[sig], te["combined_units"], s=40, facecolor="none", edgecolor=AQUA, lw=1.6, zorder=5,
               label=f"held-out test people, true length (n={te.height})")
    for x in (x0, x1):
        ax.axvline(x, color=INK2, lw=0.8, ls="-", alpha=0.5)
    ylo, yhi = ax.get_ylim()
    ax.text(x0 - lpad / 2, yhi, "beyond\ntraining\nrange", ha="center", va="top", fontsize=8.5, color=INK2)
    ax.text(x1 + pad / 2, yhi, "beyond training range:\nTheil-Sen line", ha="center", va="top", fontsize=8.5, color=INK2)
    ax.set_xlim(x0 - lpad, x1 + pad)
    ax.set_xlabel(f"short-read signal: {LABEL.get(sig, sig)}", fontsize=10, color=INK2)
    ax.set_ylabel("true combined length, both alleles (repeat units, long reads)", fontsize=10, color=INK2)
    ax.set_title(f"Calibration curve for {lid} (motif {L['period']} bp; panel core spread "
                 f"{L['cov_panel_core90_units']:.0f} repeats)\ntraining r = {p['train_pearson']:.3f}; test people "
                 f"r = {L['perf_pearson_r']:.4f}, mean abs. error {L['perf_mean_abs_err_units']:.2f} repeats",
                 fontsize=11, color=INK, loc="left")
    ax.legend(frameon=False, fontsize=8.8, loc="lower right")
    fig.tight_layout()
    stem = "curve_example_" + lid.replace(":", "_").replace("-", "_")
    for ext in ("png", "pdf"):
        fig.savefig(D / f"{stem}.{ext}", dpi=180, facecolor=SURF)
    print(f"wrote {D / stem}.png; slope {p['lin_slope']:.2f} repeats per unit signal; "
          f"line value at upper edge {lin(x1):.1f} vs isotonic {yi[-1]:.1f}")


if __name__ == "__main__":
    main(*sys.argv[1:])
