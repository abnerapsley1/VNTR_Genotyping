"""Heatmap: how much weight each short-read signal gets in the ridge-combined calibration
curves, by core spread.  Weight share = |standardised ridge coefficient| / sum over signals,
averaged over the repeats in each group (repeats where everyone has the same length excluded).

    python docs/curve_signals/make_curve_weight_figure.py
"""
import json
from pathlib import Path

import matplotlib
import numpy as np
import polars as pl

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

D = Path(__file__).resolve().parent
W = D.parents[1] / "work"
BUNDLE = Path("/data/GitHub_Repos/VNTR_Genotyping/data/vntrlen_bundle/52ec291b8f70/shipped_GRCh38_v1/loci.parquet")
SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
LABEL = {
    "depth_ratio_d20_w2000": "Depth ratio, MAPQ>=20, 2 kb flanks",
    "depth_ratio_d20_w1000": "Depth ratio, MAPQ>=20, 1 kb flanks",
    "depth_ratio_d20_w5000": "Depth ratio, MAPQ>=20, 5 kb flanks",
    "depth_ratio_d1_w2000": "Depth ratio, MAPQ>=1, 2 kb flanks",
    "depth_ratio_d0_w2000": "Depth ratio, all reads (incl. MAPQ 0), 2 kb",
    "array_depth_d20_norm": "Repeat depth / flank depth (MAPQ>=20)",
    "irr_count_norm": "In-repeat reads (normalised)",
    "partial_repeat_frac": "Reads straddling a repeat edge (fraction)",
    "anchored_irr_count_norm": "In-repeat reads, mate unmapped or on another chr.",
    "n_spanning_norm": "Read pairs spanning the repeat (normalised)",
    "spanning_tlen_shift": "Fragment-length shift of spanning pairs",
}
RAMP = LinearSegmentedColormap.from_list("blue", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
SCHEMES = (("repeats", "cov_panel_core90_units", [1, 2, 5, 20, 100], ["0", "1", "2-4", "5-19", "20-99", ">=100"]),
           ("bp", "panel_core90_bp", [1, 25, 100, 300, 1000], ["0", "1-24", "25-99", "100-299", "300-999", ">=1000"]))


def shares(cj):
    if not cj:
        return None
    e = json.loads(cj).get("pe")
    if not e or e[0] != "ridge":
        return None
    c = np.abs(np.asarray(e[1]["coef"], dtype=float))
    if c.sum() <= 0:
        return None
    return dict(zip(e[1]["signals"], (c / c.sum()).tolist()))


def main():
    L = pl.read_parquet(BUNDLE).join(pl.read_parquet(W / "model" / "panel_core_bp.parquet")
                                     .select("locus_id", "panel_core90_bp"), on="locus_id", how="left")
    S = [shares(c) for c in L["curves"].to_list()]
    keep = [i for i, s in enumerate(S) if s]
    sig = list(LABEL)
    M = np.array([[S[i].get(c, 0.0) for c in sig] for i in keep])
    L = L[keep]
    rows, tables = [], {}
    for unit, col, edges, labels in SCHEMES:
        g = np.digitize(np.nan_to_num(L[col].to_numpy().astype(float), nan=0.0), edges)
        T = np.array([M[g == k].mean(0) if (g == k).any() else np.zeros(len(sig)) for k in range(len(labels))]).T
        n = [int((g == k).sum()) for k in range(len(labels))]
        tables[unit] = (T, labels, n)
        for k, lab in enumerate(labels):
            for j, s in enumerate(sig):
                rows.append(dict(grouping=f"core spread ({unit})", group=lab, n_repeats=n[k], signal=s,
                                 mean_weight_share=round(float(T[j, k]), 4)))
    pl.DataFrame(rows).write_csv(D / "curve_weight_by_core_spread.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2), facecolor=SURF, gridspec_kw={"wspace": 0.05})
    for ax, (unit, (T, labels, n)) in zip(axes, tables.items()):
        ax.imshow(T, cmap=RAMP, vmin=0, vmax=0.35, aspect="auto")
        for i in range(T.shape[0]):
            for j in range(T.shape[1]):
                ax.text(j, i, f"{T[i, j] * 100:.0f}%", ha="center", va="center", fontsize=8.5,
                        color="white" if T[i, j] > 0.2 else INK)
        ax.set_xticks(range(len(labels)), [f"{l}\n(n={k:,})" for l, k in zip(labels, n)], fontsize=8, color=INK2)
        ax.set_yticks(range(len(sig)), [LABEL[s] for s in sig] if ax is axes[0] else [], fontsize=9, color=INK)
        ax.set_xlabel(f"core spread across the long-read panel ({unit})", fontsize=9.5, color=INK2)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"Grouped by core spread in {unit}", fontsize=10.5, color=INK, loc="left")
    fig.suptitle("Weight each short-read signal gets in the combined calibration curves (average share per repeat)",
                 fontsize=12.5, color=INK, x=0.01, ha="left", y=0.99)
    fig.text(0.01, 0.01, "Shipped models: per repeat, the 11 signals are combined by ridge regression (standardised) "
             "before the monotone curve is fitted. Share = |weight| / sum of |weights|; columns sum to 100%.\n"
             "Repeats where all training people have the same length are excluded. 'Repeat depth / flank depth "
             "(MAPQ>=20)' equals the 2 kb MAPQ>=20 depth ratio, so the two share their weight.",
             fontsize=8.5, color=INK2)
    fig.subplots_adjust(left=0.25, right=0.99, top=0.9, bottom=0.16)
    for ext in ("png", "pdf"):
        fig.savefig(D / f"curve_weight_by_core_spread.{ext}", dpi=180, facecolor=SURF)
    print(f"{len(keep):,} repeats with combined curves")


if __name__ == "__main__":
    main()
