import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

xlsx_path = os.path.join(os.path.dirname(__file__), "copy_number_variable_vntrs_summary.xlsx")
out_dir   = os.path.dirname(__file__)

print(f"Loading {xlsx_path}...")
df = pd.read_excel(xlsx_path)

for col in ["cons_len", "mean_hap1", "mean_hap2", "var_hap1", "var_hap2"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["cons_x_hap1"] = df["cons_len"] * df["mean_hap1"]
df["cons_x_hap2"] = df["cons_len"] * df["mean_hap2"]

# (x_col, y_col, x_label, y_label)
plots = [
    ("cons_len",    "mean_hap1", "Consensus Sequence Length (bp)",              "Mean Haplotype 1 Copy Number"),
    ("cons_len",    "mean_hap2", "Consensus Sequence Length (bp)",              "Mean Haplotype 2 Copy Number"),
    ("cons_len",    "var_hap1",  "Consensus Sequence Length (bp)",              "Variance of Haplotype 1 Copy Number"),
    ("cons_len",    "var_hap2",  "Consensus Sequence Length (bp)",              "Variance of Haplotype 2 Copy Number"),
    ("mean_hap1",  "var_hap1",  "Mean Haplotype 1 Copy Number",                "Variance of Haplotype 1 Copy Number"),
    ("mean_hap2",  "var_hap2",  "Mean Haplotype 2 Copy Number",                "Variance of Haplotype 2 Copy Number"),
    ("cons_x_hap1", "var_hap1", "Consensus Length × Mean Hap1 Copy Number",   "Variance of Haplotype 1 Copy Number"),
    ("cons_x_hap2", "var_hap2", "Consensus Length × Mean Hap2 Copy Number",   "Variance of Haplotype 2 Copy Number"),
]

fig, axes = plt.subplots(2, 4, figsize=(26, 12))
axes = axes.flatten()

for ax, (xcol, ycol, xlabel, ylabel) in zip(axes, plots):
    plot_data = df[[xcol, ycol]].dropna()
    x = plot_data[xcol].values
    y = plot_data[ycol].values

    ax.scatter(x, y, s=4, alpha=0.3, color="steelblue", linewidths=0)

    # Line of best fit + slope p-value
    slope, intercept, r_value, p_value, _ = linregress(x, y)
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="crimson", linewidth=1.5,
            label=f"slope = {slope:.4f}\np = {p_value:.2e}")
    ax.legend(fontsize=9, loc="upper right")

    # Zoom to 1st–99th percentile to suppress outlier distortion
    ax.set_xlim(np.percentile(x, 1), np.percentile(x, 99))
    ax.set_ylim(np.percentile(y, 1), np.percentile(y, 99))

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(f"{xlabel.split(' (')[0]} vs. {ylabel}", fontsize=10)

fig.suptitle(
    "Copy-Number-Variable VNTRs\n(vamosExpanded v3.0 GRCh38 oriMotifs)",
    fontsize=13,
)
plt.tight_layout()

scatter_path = os.path.join(out_dir, "cn_variable_scatter_plots.png")
plt.savefig(scatter_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Scatter plots saved to: {scatter_path}")

