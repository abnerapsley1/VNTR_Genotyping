"""Schematic figures of the combined-length tool's statistical methods.

    python docs/figures/make_methods_figures.py
-> docs/figures/methods_overview.{png,pdf}, methods_validation_interval.{png,pdf}
"""
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

OUT = Path(__file__).resolve().parent
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"           # known path, new path, output
TINT = {BLUE: "#e3eefb", ORANGE: "#fde9e0", AQUA: "#dff4ec", "grey": "#efeeea"}


def box(ax, x, y, w, h, title, body="", color="grey", title_size=10.5, body_size=8.6):
    edge = INK2 if color == "grey" else color
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.012",
                                fc=TINT[color], ec=edge, lw=1.3))
    ax.text(x + 0.012, y + h - 0.016, title, ha="left", va="top", fontsize=title_size,
            fontweight="bold", color=INK)
    if body:
        ax.text(x + 0.012, y + h - 0.055, body, ha="left", va="top", fontsize=body_size,
                color=INK2, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1, color=INK2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
                                 lw=1.4, color=color, shrinkA=2, shrinkB=2))


def overview():
    fig = plt.figure(figsize=(13, 9.2), facecolor=SURF)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.975, "How the tool estimates a repeat's combined length (both alleles) from short reads",
            fontsize=14, fontweight="bold", color=INK, va="top")

    # inputs
    box(ax, 0.02, 0.79, 0.30, 0.13, "Input: short-read BAM/CRAM",
        "30x Illumina, aligned to GRCh38\n(training data: 1000 Genomes NYGC CRAMs)")
    box(ax, 0.35, 0.79, 0.30, 0.13, "Repeat annotation (per locus)",
        "motif, reference length, GC, segmental dup.\nknown repeats only: long-read panel\nstatistics (mean, core spread, outliers)")
    box(ax, 0.68, 0.79, 0.30, 0.13, "Training truth (not needed to run)",
        "vamos long-read genotypes: combined\nrepeat units per person, 290 people\n(training, calibration and testing)")

    # features
    box(ax, 0.02, 0.57, 0.63, 0.17, "1  Read-level signals, extracted for every repeat",
        "• read depth over the repeat vs flanking sequence, at 4 mapping-quality levels\n"
        "• reads made of the repeat's motif (in-repeat reads, anchored pairs)\n"
        "• read pairs spanning the repeat: fragment-length shift and spread\n"
        "• divided by the person's genome-wide depth, fragment size and read length")
    arrow(ax, 0.17, 0.79, 0.17, 0.74)
    arrow(ax, 0.50, 0.79, 0.50, 0.74)

    # two paths
    box(ax, 0.02, 0.225, 0.40, 0.275, "2a  Known repeats (in the long-read catalog)",
        "Per-repeat calibration curve: monotone (isotonic) fit\n"
        "of truth on the best signal, robust-linear (Theil-Sen)\n"
        "beyond the training range; fitted per repeat\n\n"
        "Anchor = this curve's estimate\n"
        "Features = signals + annotation + panel statistics\n"
        "+ curve and transfer-model estimates", color=BLUE)
    box(ax, 0.45, 0.225, 0.33, 0.275, "2b  New repeats (not in catalog)",
        "No per-repeat curve or panel statistics\n\n"
        "Anchor = 2 x reference copy number\n"
        "Features = signals + reference-derived\n"
        "annotation + transfer-model estimate\n"
        "(gradient-boosted trees trained across\nall repeats)", color=ORANGE)
    arrow(ax, 0.20, 0.57, 0.20, 0.50, BLUE)
    arrow(ax, 0.50, 0.57, 0.60, 0.50, ORANGE)
    ax.text(0.33, 0.535, "is the repeat in the catalog?", fontsize=9, color=INK2, ha="center", style="italic")

    # model
    box(ax, 0.02, 0.02, 0.76, 0.175, "3  Gradient-boosted trees model the correction to the anchor",
        "target = log(true + 1) - log(anchor + 1):  one model for the point estimate (squared error),\n"
        "two for the 5th and 95th percentiles (quantile loss)\n\n"
        "4  Split-conformal calibration on held-out data (validation people for known repeats, validation repeats\n"
        "for new ones): the 5-95% range is widened or narrowed per repeat group so 90% of true values fall inside,\n"
        "then rounded out to whole repeats", body_size=8.6)
    arrow(ax, 0.20, 0.225, 0.20, 0.195, BLUE)
    arrow(ax, 0.60, 0.225, 0.60, 0.195, ORANGE)

    # outputs
    box(ax, 0.81, 0.02, 0.17, 0.72, "Output per repeat", "", color=AQUA)
    items = [("Point estimate", "combined repeat units\n(both alleles), decimal"),
             ("90% interval", "whole repeat units;\ncontains the truth for\n91-92% of test calls"),
             ("Tier", "PASS: interval narrower\nthan the estimate and\nrepeat not unusual\nLOW_CONF: otherwise\nNA: no usable reads")]
    yy = 0.66
    for t, b in items:
        ax.text(0.825, yy, t, fontsize=10, fontweight="bold", color=INK, va="top")
        ax.text(0.825, yy - 0.032, b, fontsize=8.6, color=INK2, va="top", linespacing=1.35)
        yy -= 0.19
    arrow(ax, 0.78, 0.08, 0.81, 0.08, AQUA)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"methods_overview.{ext}", dpi=200, facecolor=SURF)
    plt.close(fig)


def validation_and_interval():
    fig = plt.figure(figsize=(14, 6.8), facecolor=SURF)
    fig.text(0.02, 0.965, "How the model was trained and checked, and how the interval and tier are built",
             fontsize=14, fontweight="bold", color=INK, va="top")

    # left: people x repeats grid (widths schematic, not to scale)
    ax = fig.add_axes([0.13, 0.20, 0.36, 0.62])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    cols = [(0.0, 0.46, "212 training\npeople"), (0.46, 0.70, "26 validation\npeople"), (0.70, 1.0, "52 test\npeople")]
    rows = [(0.40, 1.0, "92% of repeats\n(trained on)"), (0.20, 0.40, "3% validation\nrepeats"),
            (0.0, 0.20, "5% test\nrepeats")]
    GREY = "#c3c2b7"
    fill = {(0, 0): (BLUE, "fit all\nmodels"), (1, 0): ("#86b6ef", "calibrate\n90% intervals\n(known path)"),
            (2, 0): (AQUA, "TEST\nknown repeats,\nnew people"),
            (0, 1): ("#f3a988", "calibrate 90% intervals\n(new-repeat path)"), (1, 1): (GREY, ""), (2, 1): (GREY, ""),
            (0, 2): (ORANGE, "TEST\nnew repeats"), (1, 2): (GREY, ""), (2, 2): (AQUA, "TEST\nboth new")}
    for ri, (y0, y1, rl) in enumerate(rows):
        for ci, (x0, x1, cl) in enumerate(cols):
            c, lab = fill[(ci, ri)]
            ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=c, ec=SURF, lw=3))
            ax.text((x0 + x1) / 2, (y0 + y1) / 2, lab, ha="center", va="center", fontsize=9,
                    color="white", fontweight="bold", linespacing=1.3)
        ax.text(-0.03, (y0 + y1) / 2, rl, ha="right", va="center", fontsize=9.5, color=INK2)
    for x0, x1, cl in cols:
        ax.text((x0 + x1) / 2, 1.03, cl, ha="center", va="bottom", fontsize=9.5, color=INK2)
    ax.text(0.5, -0.06, "Box sizes are schematic, not to scale; grey = not used.\n"
            "People are split by family (relatives never cross groups); repeats by a fixed hash.\n"
            "Model estimates for training people used inside the model are cross-fitted (5 folds).",
            ha="center", va="top", fontsize=8.6, color=INK2, linespacing=1.4)

    # right: interval construction for one call
    ax = fig.add_axes([0.57, 0.30, 0.40, 0.50])
    ax.set_facecolor(SURF)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_yticks([])
    ax.set_xlim(23, 37)
    ax.set_ylim(0.4, 4.0)
    ax.tick_params(colors=INK2, labelsize=8.5)
    ax.set_xlabel("combined repeat units (both alleles)", fontsize=9, color=INK2)
    est = 29.6
    steps = [(3.3, 26.9, 32.4, "1  raw 5-95% range from the quantile models", "#86b6ef"),
             (2.2, 26.2, 33.3, "2  shifted by the conformal correction for this repeat's group", BLUE),
             (1.1, 26, 34, "3  rounded out to whole repeats = reported 90% interval (26-34)", AQUA)]
    for y, lo, hi, lab, c in steps:
        ax.plot([lo, hi], [y, y], color=c, lw=8, solid_capstyle="round")
        ax.text(23.2, y + 0.38, lab, fontsize=9, color=INK, va="center",
                bbox=dict(fc=SURF, ec="none", pad=1.5), zorder=3)
    ax.plot([est] * 2, [0.7, 3.95], color=INK, lw=1.2, ls="--")
    ax.text(est + 0.15, 3.85, f"point estimate {est}", ha="left", fontsize=8.6, color=INK)
    ax.set_title("Example: building one call's 90% interval", fontsize=10.5, color=INK, loc="left",
                 fontweight="bold", pad=12)
    fig.text(0.57, 0.12, "Tier for this call:  PASS, because the interval width (34 - 26 = 8) is no larger than\n"
             "the estimate (29.6) and the repeat is not in the top 10% of the novelty score.\n"
             "Otherwise LOW_CONF;  NA when the repeat could not be read.", fontsize=9, color=INK2, va="top",
             linespacing=1.4)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"methods_validation_interval.{ext}", dpi=200, facecolor=SURF)
    plt.close(fig)


if __name__ == "__main__":
    overview()
    validation_and_interval()
    print(f"wrote figures to {OUT}")
