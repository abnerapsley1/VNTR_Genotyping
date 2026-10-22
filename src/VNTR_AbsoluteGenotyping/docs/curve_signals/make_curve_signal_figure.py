"""Heatmap: which short-read signal each repeat's calibration curve uses, by core spread.
    python docs/curve_signals/make_curve_signal_figure.py
"""
from pathlib import Path

import matplotlib
import numpy as np
import polars as pl

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

D = Path(__file__).resolve().parent
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
LABEL = {
    "depth_ratio_d20_w2000": "Depth ratio, MAPQ>=20, 2 kb flanks",
    "depth_ratio_d20_w1000": "Depth ratio, MAPQ>=20, 1 kb flanks",
    "depth_ratio_d20_w5000": "Depth ratio, MAPQ>=20, 5 kb flanks",
    "depth_ratio_d1_w2000": "Depth ratio, MAPQ>=1, 2 kb flanks",
    "depth_ratio_d0_w2000": "Depth ratio, all reads (incl. MAPQ 0), 2 kb",
    "irr_count_norm": "In-repeat reads (normalised)",
    "partial_repeat_frac": "Reads straddling a repeat edge (fraction)",
    "anchored_irr_count_norm": "In-repeat reads, mate unmapped or on another chr.",
    "n_spanning_norm": "Read pairs spanning the repeat (normalised)",
    "spanning_tlen_shift": "Fragment-length shift of spanning pairs",
    "(no curve)": "(no curve fitted)",
}
ORDER = list(LABEL)
RAMP = LinearSegmentedColormap.from_list("blue", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])


def main():
    t = pl.read_csv(D / "curve_signal_by_core_spread.csv").filter(pl.col("curve") == "paired-end curve")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2), facecolor=SURF, gridspec_kw={"wspace": 0.05})
    for ax, (grouping, groups, unit) in zip(axes, (
            ("core spread (repeat units)", ["0", "1", "2-4", "5-19", "20-99", ">=100"], "repeats"),
            ("core spread (bp)", ["0", "1-24", "25-99", "100-299", "300-999", ">=1000"], "bp"))):
        s = t.filter(pl.col("grouping") == grouping)
        n = s.group_by("group").agg(pl.col("n").sum())
        nmap = dict(n.iter_rows())
        M = np.zeros((len(ORDER), len(groups)))
        for r in s.iter_rows(named=True):
            if r["signal"] in ORDER and r["group"] in groups:
                M[ORDER.index(r["signal"]), groups.index(r["group"])] = r["share"]
        ax.imshow(M, cmap=RAMP, vmin=0, vmax=0.5, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if M[i, j] >= 0.005:
                    ax.text(j, i, f"{M[i, j] * 100:.0f}%", ha="center", va="center", fontsize=8.5,
                            color="white" if M[i, j] > 0.28 else INK)
        ax.set_xticks(range(len(groups)), [f"{g}\n(n={nmap.get(g, 0):,})" for g in groups], fontsize=8,
                      color=INK2)
        ax.set_yticks(range(len(ORDER)), [LABEL[o] for o in ORDER] if ax is axes[0] else [], fontsize=9,
                      color=INK)
        ax.set_xlabel(f"core spread across the long-read panel ({unit})", fontsize=9.5, color=INK2)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"Grouped by core spread in {unit}", fontsize=10.5, color=INK, loc="left")
    fig.suptitle("Which short-read signal each repeat's calibration curve uses (share of repeats in each column)",
                 fontsize=12.5, color=INK, x=0.01, ha="left", y=0.99)
    fig.text(0.01, 0.01, "Shipped models: one curve per repeat, fitted on 212 training people; each curve uses the "
             "signal that best tracks the long-read truth at that repeat. Columns sum to 100%.\n"
             "(The 11th shortlist signal, repeat depth / flank depth at MAPQ>=20, equals the 2 kb depth ratio and is never picked.)",
             fontsize=8.5, color=INK2)
    fig.subplots_adjust(left=0.24, right=0.99, top=0.9, bottom=0.16)
    for ext in ("png", "pdf"):
        fig.savefig(D / f"curve_signal_by_core_spread.{ext}", dpi=180, facecolor=SURF)


if __name__ == "__main__":
    main()
