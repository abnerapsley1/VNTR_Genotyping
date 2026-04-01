"""
build_flow_diagram.py — Generate a flow diagram of the VNTR Genotyping Tool pipeline.

Output: VNTR_GenotypingTool_Flow.png in the same directory as this script.

Usage:
    python build_flow_diagram.py
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------------------
# Output path
# ---------------------------------------------------------------------------

OUT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUT_FILE = os.path.join(OUT_DIR, "VNTR_GenotypingTool_Flow.png")

# ---------------------------------------------------------------------------
# Coordinate system
#   x : 0 → 10   (data units)
#   y : 80 → 200 (data units, top = 200, bottom = 80)
# ---------------------------------------------------------------------------

FIG_W = 20
FIG_H = 28

# Colours
C_TERM    = "#2C3E50"
C_SETUP   = "#1A5276"
C_NORM    = "#1E8449"
C_IO      = "#7D3C98"
C_COMPUTE = "#9A7D0A"
C_DECIDE  = "#CB4335"
C_OUTPUT  = "#1F618D"
C_TEXT    = "white"
C_ARROW   = "#2C2C2C"
C_BORDER  = "#111111"
C_BG      = "#F4F6F7"

FONT   = 10.5
FONT_S = 9.0
FONT_XS = 8.0


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def box(ax, cx, cy, w, h, label, color, sublabel=None, fontsize=FONT):
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.10",
        linewidth=1.2, edgecolor=C_BORDER, facecolor=color, zorder=3,
    )
    ax.add_patch(rect)
    if sublabel:
        ax.text(cx, cy + h * 0.15, label,
                ha="center", va="center", fontsize=fontsize, color=C_TEXT,
                fontweight="bold", zorder=4, multialignment="center")
        ax.text(cx, cy - h * 0.28, sublabel,
                ha="center", va="center", fontsize=fontsize - 1.5, color=C_TEXT,
                style="italic", zorder=4, multialignment="center")
    else:
        ax.text(cx, cy, label,
                ha="center", va="center", fontsize=fontsize, color=C_TEXT,
                fontweight="bold", zorder=4, multialignment="center")


def diamond(ax, cx, cy, w, h, label, color, fontsize=FONT):
    xs = [cx, cx + w/2, cx, cx - w/2, cx]
    ys = [cy + h/2, cy, cy - h/2, cy, cy + h/2]
    ax.fill(xs, ys, color=color, zorder=3)
    ax.plot(xs, ys, color=C_BORDER, lw=1.2, zorder=4)
    ax.text(cx, cy, label,
            ha="center", va="center", fontsize=fontsize, color=C_TEXT,
            fontweight="bold", zorder=5, multialignment="center")


def arr(ax, x0, y0, x1, y1, label=None, lside="right", dashed=False):
    ap = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle="-|>",
        color="#AAAAAA" if dashed else C_ARROW,
        lw=1.4 if dashed else 1.8,
        mutation_scale=13 if dashed else 16,
        connectionstyle="arc3,rad=0.0",
        linestyle="dashed" if dashed else "solid",
        zorder=2,
    )
    ax.add_patch(ap)
    if label:
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        dx = 0.20 if lside == "right" else -0.20
        ax.text(mx + dx, my, label,
                ha="left" if lside == "right" else "right",
                va="center", fontsize=FONT_XS, color=C_ARROW, style="italic", zorder=5)


def hline(ax, x0, y, x1):
    ax.plot([x0, x1], [y, y], color=C_ARROW, lw=1.8, zorder=2)


def vline(ax, x, y0, y1, dashed=False):
    ax.plot([x, x], [y0, y1],
            color="#AAAAAA" if dashed else C_ARROW,
            lw=1.4 if dashed else 1.8,
            ls="--" if dashed else "-", zorder=2)


def section(ax, y, text):
    ax.axhline(y, color="#CCCCCC", lw=0.9, zorder=1)
    ax.text(0.05, y + 0.30, text,
            ha="left", va="bottom", fontsize=FONT_XS, color="#777777",
            style="italic", zorder=2)


def branch_label(ax, x, y, text, ha="left"):
    ax.text(x, y, text, ha=ha, va="bottom",
            fontsize=FONT_XS, color=C_ARROW, style="italic", zorder=5)


# ---------------------------------------------------------------------------
# Main diagram
# ---------------------------------------------------------------------------

def build():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.set_xlim(0, 10)
    ax.set_ylim(116, 200)
    ax.axis("off")
    fig.patch.set_facecolor(C_BG)
    ax.set_facecolor(C_BG)

    # Layout constants
    CX  = 5.0    # main spine x
    LX  = 1.9    # left branch x
    RX  = 8.1    # right branch x
    BW  = 4.0    # standard box width
    SBW = 2.6    # side-branch box width
    BH  = 1.6    # standard box height
    DBH = 2.0    # diamond height
    DBW = 3.0    # diamond width

    # ── Title ──────────────────────────────────────────────────────────
    ax.text(5, 199.0, "VNTR Genotyping Tool — Processing Flow",
            ha="center", va="top", fontsize=16, fontweight="bold", color=C_TERM)
    ax.text(5, 197.2, "count_vntrs()  ·  _counting.py",
            ha="center", va="top", fontsize=11, color="#555555", style="italic")

    # ── Entry point ────────────────────────────────────────────────────
    section(ax, 195.8, "Entry point")
    box(ax, CX, 194.5, BW, BH, "START  —  count_vntrs() called", C_TERM)
    arr(ax, CX, 194.5 - BH/2, CX, 192.5 + BH/2)

    box(ax, CX, 192.5, BW, BH,
        "Print run header\n(samples, norm method, output path)", C_SETUP)
    arr(ax, CX, 192.5 - BH/2, CX, 190.2 + BH/2)

    # ── Step 1 ─────────────────────────────────────────────────────────
    section(ax, 189.8, "Step 1 · Region Loading")
    box(ax, CX, 190.2, BW, BH,
        "build_regions()\nLoad & filter VNTR BED file", C_SETUP,
        sublabel="--default / --gene / --vntr / --regions")
    arr(ax, CX, 190.2 - BH/2, CX, 187.5 + DBH/2)

    diamond(ax, CX, 187.5, DBW, DBH, "--chrom\nfilter?", C_DECIDE)
    # Yes → left
    hline(ax, CX - DBW/2, 187.5, LX + SBW/2)
    branch_label(ax, LX + SBW/2 + 0.08, 187.5 + 0.3, "Yes")
    box(ax, LX, 187.5, SBW, BH, "Filter to requested\nchromosome(s)", C_SETUP)
    arr(ax, LX, 187.5 - BH/2, CX, 185.2 + BH/2)
    # No → down
    arr(ax, CX, 187.5 - DBH/2, CX, 185.2 + BH/2, label="No")

    box(ax, CX, 185.2, BW, BH, "N VNTRs ready for genotyping", C_SETUP)
    arr(ax, CX, 185.2 - BH/2, CX, 182.5 + DBH/2)

    # ── Step 2 ─────────────────────────────────────────────────────────
    section(ax, 182.1, "Step 2 · Normalization Setup")
    diamond(ax, CX, 182.5, DBW, DBH, "norm_method?", C_DECIDE)

    # — gene branch (left) —
    hline(ax, CX - DBW/2, 182.5, LX + SBW/2)
    branch_label(ax, LX + SBW/2 + 0.08, 182.5 + 0.3, "gene")
    arr(ax, LX + SBW/2, 182.5, LX, 180.0 + BH/2)

    box(ax, LX, 180.0, SBW, BH,
        "parse_gtf()\nLoad gene annotation", C_NORM,
        sublabel="bundled GENCODE v38 or custom GTF")
    arr(ax, LX, 180.0 - BH/2, LX, 177.5 + BH/2)

    box(ax, LX, 177.5, SBW, BH,
        "assign_nearest_genes()\nBisect-index gene lookup", C_NORM,
        sublabel="O(N·K) vs old O(N×G)")
    arr(ax, LX, 177.5 - BH/2, LX, 175.0 + BH/2)

    box(ax, LX, 175.0, SBW, BH,
        "Build normalization backgrounds\nEffective gene lengths", C_NORM)
    arr(ax, LX, 175.0 - BH/2, LX, 172.5 + BH/2)

    box(ax, LX, 172.5, SBW, BH,
        "parse_psl()  [optional]\nAlt-contig gene regions", C_NORM,
        sublabel="bundled GRCh38 PSL or custom")
    arr(ax, LX, 172.5 - BH/2, LX, 169.5)
    hline(ax, LX, 169.5, CX)

    # — local branch (right) —
    hline(ax, CX + DBW/2, 182.5, RX - SBW/2)
    branch_label(ax, RX - SBW/2 - 0.08, 182.5 + 0.3, "local", ha="right")
    arr(ax, RX - SBW/2, 182.5, RX, 180.0 + BH/2)

    box(ax, RX, 180.0, SBW, BH,
        "Local window setup\nnorm_window bp flanks", C_NORM,
        sublabel="no GTF required")
    arr(ax, RX, 180.0 - BH/2, RX, 169.5)
    hline(ax, CX, 169.5, RX)

    # — no norm bypass (dashed) —
    vline(ax, 9.5, 182.5 - DBH/2, 169.5, dashed=True)
    ax.text(9.56, (182.5 - DBH/2 + 169.5) / 2, "no norm\n(raw counts)",
            ha="left", va="center", fontsize=FONT_XS, color="#888888", style="italic")
    hline(ax, CX, 169.5, 9.5)

    arr(ax, CX, 169.5, CX, 167.8 + BH/2)

    # ── Step 3 ─────────────────────────────────────────────────────────
    section(ax, 167.4, "Step 3 · Per-Sample Processing")
    box(ax, CX, 167.8, BW, BH,
        "Determine output metric\npredicted_copies / density_ratio / read_count", C_COMPUTE)
    arr(ax, CX, 167.8 - BH/2, CX, 165.0 + DBH/2)

    diamond(ax, CX, 165.0, DBW, DBH, "workers > 1\n& samples > 1?", C_DECIDE)

    # Yes → parallel (left)
    hline(ax, CX - DBW/2, 165.0, LX + SBW/2)
    branch_label(ax, LX + SBW/2 + 0.08, 165.0 + 0.3, "Yes")
    arr(ax, LX + SBW/2, 165.0, LX, 162.5 + BH/2)
    box(ax, LX, 162.5, SBW, BH,
        "ProcessPoolExecutor\nN worker processes", C_IO,
        sublabel="shared data via initializer")
    arr(ax, LX, 162.5 - BH/2, LX, 159.5)
    hline(ax, LX, 159.5, CX)

    # No → serial (straight down)
    arr(ax, CX, 165.0 - DBH/2, CX, 162.5 + BH/2, label="No")
    box(ax, CX, 162.5, BW, BH,
        "Serial sample loop\n(vntr_threads if single sample)", C_IO)
    arr(ax, CX, 162.5 - BH/2, CX, 159.5)

    arr(ax, CX, 159.5, CX, 157.5 + BH/2)
    box(ax, CX, 157.5, BW, BH,
        "_process_sample_task()\nper-sample worker", C_IO,
        sublabel="[i/N] sample_name")
    arr(ax, CX, 157.5 - BH/2, CX, 154.8 + DBH/2)

    # ── Step 3a ────────────────────────────────────────────────────────
    section(ax, 154.4, "Step 3a · Inside each sample")
    diamond(ax, CX, 154.8, DBW, DBH, "normalize?", C_DECIDE)

    # No → raw counts (right)
    hline(ax, CX + DBW/2, 154.8, RX - SBW/2)
    branch_label(ax, RX - SBW/2 - 0.08, 154.8 + 0.3, "No", ha="right")
    arr(ax, RX - SBW/2, 154.8, RX, 152.0 + BH/2)
    box(ax, RX, 152.0, SBW, BH,
        "Fetch raw read counts\nget_read_names() × N_vntrs", C_IO)
    arr(ax, RX, 152.0 - BH/2, RX, 143.5)
    hline(ax, CX, 143.5, RX)

    # Yes → fetch VNTR reads (down)
    arr(ax, CX, 154.8 - DBH/2, CX, 152.0 + BH/2, label="Yes")
    box(ax, CX, 152.0, BW, BH,
        "Fetch VNTR reads\nget_read_names() × N_vntrs", C_IO,
        sublabel="ThreadPoolExecutor if vntr_workers > 1")
    arr(ax, CX, 152.0 - BH/2, CX, 149.3 + DBH/2)

    diamond(ax, CX, 149.3, DBW, DBH, "norm_method?", C_DECIDE)

    # gene → left
    hline(ax, CX - DBW/2, 149.3, LX + SBW/2)
    branch_label(ax, LX + SBW/2 + 0.08, 149.3 + 0.3, "gene")
    arr(ax, LX + SBW/2, 149.3, LX, 146.8 + BH/2)
    box(ax, LX, 146.8, SBW, BH,
        "Fetch gene background reads\nexclude overlapping VNTR reads", C_IO)
    arr(ax, LX, 146.8 - BH/2, LX, 143.5)
    hline(ax, LX, 143.5, CX)

    # local → right
    hline(ax, CX + DBW/2, 149.3, RX - SBW/2)
    branch_label(ax, RX - SBW/2 - 0.08, 149.3 + 0.3, "local", ha="right")
    arr(ax, RX - SBW/2, 149.3, RX, 146.8 + BH/2)
    box(ax, RX, 146.8, SBW, BH,
        "Fetch window reads\n(norm_window bp flanks)", C_IO)
    arr(ax, RX, 146.8 - BH/2, RX, 143.5)
    hline(ax, CX, 143.5, RX)

    # converge → compute
    arr(ax, CX, 143.5, CX, 141.5 + BH/2)
    box(ax, CX, 141.5, BW, BH,
        "Compute per-VNTR metric\ndensity_ratio() → predicted_copies", C_COMPUTE,
        sublabel="all from read cache — no further I/O")
    arr(ax, CX, 141.5 - BH/2, CX, 139.0 + BH/2)

    box(ax, CX, 139.0, BW, BH,
        "Append row to results\n{sample, metric, vntr_1 … vntr_N}", C_COMPUTE)
    arr(ax, CX, 139.0 - BH/2, CX, 136.3 + DBH/2)

    diamond(ax, CX, 136.3, DBW, DBH, "More\nsamples?", C_DECIDE)

    # Loop back arc (right side)
    ax.annotate("",
                xy=(CX + BW/2 + 0.12, 157.5),
                xytext=(CX + DBW/2, 136.3),
                arrowprops=dict(
                    arrowstyle="-|>", color=C_ARROW, lw=1.6,
                    connectionstyle="arc3,rad=-0.35",
                    mutation_scale=15,
                ))
    ax.text(8.85, 147.0, "Yes\n(next sample)",
            ha="left", va="center", fontsize=FONT_XS, color=C_ARROW, style="italic")

    arr(ax, CX, 136.3 - DBH/2, CX, 133.8 + BH/2, label="No")

    # ── Step 4 ─────────────────────────────────────────────────────────
    section(ax, 133.4, "Step 4 · Output")
    box(ax, CX, 133.8, BW, BH,
        "Build pandas DataFrame\ncolumns: sample, metric, <vntr_1> … <vntr_N>", C_OUTPUT)
    arr(ax, CX, 133.8 - BH/2, CX, 131.1 + DBH/2)

    diamond(ax, CX, 131.1, DBW, DBH, "output_csv\nspecified?", C_DECIDE)

    # Yes → left
    hline(ax, CX - DBW/2, 131.1, LX + SBW/2)
    branch_label(ax, LX + SBW/2 + 0.08, 131.1 + 0.3, "Yes")
    arr(ax, LX + SBW/2, 131.1, LX, 128.5 + BH/2)
    box(ax, LX, 128.5, SBW, BH,
        "Write CSV\ndf.to_csv(output_csv)", C_OUTPUT)
    arr(ax, LX, 128.5 - BH/2, LX, 126.0)
    hline(ax, LX, 126.0, CX)

    # No → down
    arr(ax, CX, 131.1 - DBH/2, CX, 126.0, label="No")
    arr(ax, CX, 126.0, CX, 124.0 + BH/2)

    box(ax, CX, 124.0, BW, BH,
        "Print completion summary\nReturn DataFrame to caller", C_OUTPUT)
    arr(ax, CX, 124.0 - BH/2, CX, 121.5 + BH/2)

    box(ax, CX, 121.5, BW, BH, "END", C_TERM)

    # ── Legend ─────────────────────────────────────────────────────────
    legend_items = [
        (C_TERM,    "Start / End"),
        (C_SETUP,   "Setup / Loading"),
        (C_NORM,    "Normalization setup"),
        (C_IO,      "I/O (file reads)"),
        (C_COMPUTE, "Computation"),
        (C_DECIDE,  "Decision"),
        (C_OUTPUT,  "Output"),
    ]
    lx, ly = 0.3, 118.5
    ax.text(lx, ly + 1.2, "Legend:", fontsize=FONT_S, color="#444444", fontweight="bold")
    for color, label in legend_items:
        rect = FancyBboxPatch((lx, ly - 0.5), 1.0, 1.0,
                              boxstyle="round,pad=0.06",
                              facecolor=color, edgecolor=C_BORDER, lw=0.9)
        ax.add_patch(rect)
        ax.text(lx + 1.2, ly + 0.0, label,
                va="center", fontsize=FONT_S, color="#333333")
        lx += 1.35 + len(label) * 0.09

    plt.tight_layout(pad=0.5)
    plt.savefig(OUT_FILE, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"Flow diagram saved to: {OUT_FILE}")


if __name__ == "__main__":
    build()
