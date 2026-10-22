"""Figures for the tool comparison (reads compare_tools.py outputs).

Tools are drawn from whatever compare_tools.py found: once danbing-tk's length
predictions exist, rerun compare_tools.py and this script and it appears in
every figure, in its fixed colour.

  compare_headline.png       call rate / accuracy / person-to-person tracking /
                             rare-allele direction, all sites and sites every tool called
  compare_by_<stratum>.png   the same measures across core spread (bp, repeats),
                             reference array length, motif length, site type
  compare_scatter.png        truth vs call (bp), one panel per tool, density

    python src/VNTR_Benchmark/make_compare_figures.py [--dir work/compare]
"""
import argparse
import warnings
from pathlib import Path

import matplotlib
import numpy as np
import polars as pl

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
SURF, INK, INK2, GRID, REF = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#8a8984"
# fixed categorical order (validated palette): colour follows the tool, never its rank
TOOL_STYLE = {"ours": ("Our tool", "#2a78d6"), "gangstr": ("GangSTR", "#eb6834"),
              "danbing": ("danbing-tk", "#1baf7a")}
MEASURES = [("call_rate", "Sites called", "pct", "higher is better"),
            ("median_rel_err", "Median error\n(% of true length)", "pct", "lower is better"),
            ("frac_within_5pct", "Calls within 5%\nof true length", "pct", "higher is better"),
            ("within_locus_r", "Tracks person-to-person\ndifferences (r)", "r", "higher is better"),
            ("rare_right_direction", "Rare long/short alleles:\ncalled in right direction", "pct", "higher is better")]
STRATA = {"core_bp_bin": ("core spread across the long-read panel (bp)", None),
          "core_units_bin": ("core spread across the long-read panel (repeats)", None),
          "span_bin": ("reference array length (bp); short-read fragments ~450 bp", None),
          "motif_bin": ("motif length (bp)", None),
          "site_type": ("site type", ["stable", "stable core + outliers", "broad spread", "broad + outliers"])}


def style(ax):
    ax.set_facecolor(SURF)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def fmt(v, kind):
    if v is None or not np.isfinite(v):
        return "n/a"
    return f"{v * 100:.1f}%" if kind == "pct" else f"{v:.2f}"


def headline(m: pl.DataFrame, tools: list[str], out: Path, n_people: int):
    fig, axes = plt.subplots(2, len(MEASURES), figsize=(3.0 * len(MEASURES), 6.4), facecolor=SURF)
    for row, (scope, title) in enumerate((("own_calls", "Every site each tool called"),
                                          ("all_tools_called", "Only sites every tool called"))):
        s = m.filter((pl.col("stratum") == "all") & (pl.col("scope") == scope))
        for j, (col, lab, kind, better) in enumerate(MEASURES):
            ax = axes[row, j]
            style(ax)
            vals = [s.filter(pl.col("tool") == t)[col].item() if col in s.columns and
                    s.filter(pl.col("tool") == t).height else np.nan for t in tools]
            vals = [np.nan if v is None else v for v in vals]
            x = np.arange(len(tools))
            ax.bar(x, np.nan_to_num(vals), width=0.62, color=[TOOL_STYLE[t][1] for t in tools],
                   edgecolor=SURF, linewidth=2)
            top = max([v for v in vals if np.isfinite(v)] + [0.01])
            for xi, v in zip(x, vals):
                ax.text(xi, (v if np.isfinite(v) else 0) + top * 0.03, fmt(v, kind), ha="center",
                        va="bottom", fontsize=8.5, color=INK)
            ax.set_xticks(x, [TOOL_STYLE[t][0] for t in tools], fontsize=8, color=INK)
            ax.set_ylim(0, top * 1.22 if kind == "pct" else max(1.0, top * 1.22))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(
                (lambda v, _: f"{v * 100:.0f}%") if kind == "pct" else (lambda v, _: f"{v:.1f}")))
            if row == 0:
                ax.set_title(f"{lab}\n", fontsize=9.5, color=INK, loc="left")
                ax.text(0, 1.02, better, transform=ax.transAxes, fontsize=7.5, color=INK2)
            if j == 0:
                ax.set_ylabel(title, fontsize=9, color=INK)
    fig.suptitle(f"Tool comparison on {n_people} held-out test people (truth: vamos long-read calls)",
                 fontsize=12, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=170, facecolor=SURF)
    plt.close(fig)


def by_stratum(m: pl.DataFrame, tools: list[str], stratum: str, out: Path):
    desc, order = STRATA[stratum]
    s = m.filter((pl.col("stratum") == stratum) & (pl.col("scope") == "own_calls"))
    if s.height == 0:
        return
    groups = order or _natural(s["group"].unique().to_list())
    groups = [g for g in groups if g in set(s["group"])]
    nsites = {g: s.filter(pl.col("group") == g)["n_sites"].max() for g in groups}
    meas = [mm for mm in MEASURES if mm[0] != "rare_right_direction"]
    fig, axes = plt.subplots(2, 2, figsize=(max(11, 1.5 * len(groups) * 2 + 2), 8.2), facecolor=SURF)
    x = np.arange(len(groups))
    w = 0.8 / len(tools)
    for ax, (col, lab, kind, better) in zip(axes.flat, meas):
        style(ax)
        for i, t in enumerate(tools):
            v = [s.filter((pl.col("group") == g) & (pl.col("tool") == t))[col].item()
                 if s.filter((pl.col("group") == g) & (pl.col("tool") == t)).height else np.nan
                 for g in groups]
            v = np.array([np.nan if a is None else a for a in v], dtype=float)
            ax.bar(x + (i - (len(tools) - 1) / 2) * w, np.nan_to_num(v), width=w, color=TOOL_STYLE[t][1],
                   edgecolor=SURF, linewidth=1.5, label=TOOL_STYLE[t][0])
        ax.set_xticks(x, [f"{g}\n{nsites[g]:,} sites" for g in groups], fontsize=7.5, color=INK2)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(
            (lambda v, _: f"{v * 100:.0f}%") if kind == "pct" else (lambda v, _: f"{v:.1f}")))
        if kind == "r":
            ax.set_ylim(min(0, ax.get_ylim()[0]), 1)
            ax.axhline(0, color=REF, lw=0.8)
        ax.set_title(f"{lab.replace(chr(10), ' ')}  ({better})", fontsize=9.5, color=INK, loc="left")
    h, l = axes.flat[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=9, loc="upper right", ncol=len(tools), bbox_to_anchor=(0.99, 0.985))
    fig.suptitle(f"By {desc}", fontsize=11.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=170, facecolor=SURF)
    plt.close(fig)


def _natural(labels):
    def key(s):
        num = "".join(ch if ch.isdigit() else " " for ch in s).split()
        return (s.startswith(">"), int(num[0]) if num else -1)
    return sorted(labels, key=key)


def scatter(calls: pl.DataFrame, tools: list[str], out: Path):
    fig, axes = plt.subplots(1, len(tools), figsize=(4.3 * len(tools), 4.5), facecolor=SURF, squeeze=False)
    lo, hi = 10, 3e4
    for ax, t in zip(axes[0], tools):
        style(ax)
        ax.grid(color=GRID, lw=0.6)
        c = calls.filter(pl.col(f"{t}_bp").is_finite() & (pl.col(f"{t}_bp") > 0) & (pl.col("combined_bp") > 0))
        x = np.clip(c["combined_bp"].to_numpy().astype(float), lo, hi)
        y = np.clip(c[f"{t}_bp"].to_numpy().astype(float), lo, hi)
        hb = ax.hexbin(x, y, xscale="log", yscale="log", gridsize=70, mincnt=1, norm=LogNorm(),
                       cmap=_ramp(TOOL_STYLE[t][1]), extent=(np.log10(lo), np.log10(hi), np.log10(lo), np.log10(hi)),
                       linewidths=0)
        ax.plot([lo, hi], [lo, hi], color=REF, lw=1, ls="--")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        rel = np.median(np.abs(y - x) / x)
        ax.set_title(f"{TOOL_STYLE[t][0]}\n{len(x):,} calls · median error {rel * 100:.1f}% of length",
                     fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("true array length, both alleles (bp, vamos)", fontsize=8, color=INK2)
        cb = fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=7, colors=INK2)
        cb.set_label("person-site calls per cell", fontsize=7, color=INK2)
    axes[0][0].set_ylabel("tool's call (bp)", fontsize=8, color=INK2)
    fig.suptitle("Called vs true length, held-out test people (dashed line = perfect)",
                 fontsize=11.5, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=170, facecolor=SURF)
    plt.close(fig)


def _ramp(hex_color):
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    r, g, b = to_rgb(hex_color)
    return LinearSegmentedColormap.from_list("t", [(0.86 + 0.14 * r, 0.86 + 0.14 * g, 0.86 + 0.14 * b),
                                                   (r, g, b), (r * 0.45, g * 0.45, b * 0.45)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, default=ROOT / "src" / "VNTR_Benchmark" / "work" / "compare")
    args = ap.parse_args()
    m = pl.read_csv(args.dir / "tool_metrics.csv", infer_schema_length=None)
    calls = pl.read_parquet(args.dir / "calls.parquet")
    tools = [t for t in TOOL_STYLE if t in set(m["tool"])]
    n_people = int(m["n_people"].max())
    figdir = args.dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    headline(m, tools, figdir / "compare_headline.png", n_people)
    for st in STRATA:
        by_stratum(m, tools, st, figdir / f"compare_by_{st.replace('_bin', '')}.png")
    scatter(calls, tools, figdir / "compare_scatter.png")
    print(f"wrote figures for {', '.join(tools)} -> {figdir}")


if __name__ == "__main__":
    main()
