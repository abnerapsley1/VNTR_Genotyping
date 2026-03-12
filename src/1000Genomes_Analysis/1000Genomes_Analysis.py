import os
import csv
import statistics
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# File paths
data_dir = os.path.join(os.path.dirname(__file__), "../../data/1000Genomes")
vcf_path = os.path.join(data_dir, "vamosExpanded_v3.0_GRCh38_oriMotifs.vcf")
tsv_path = os.path.join(data_dir, "vamosExpanded_v3.0_GRCh38_oriMotifs.tsv")
out_dir  = os.path.dirname(__file__)

# Standard chromosome order
chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# ---------------------------------------------------------------------------
# 1. VNTR counts per chromosome
# ---------------------------------------------------------------------------
chrom_counts = Counter()
with open(vcf_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        chrom = line.split("\t", 1)[0]
        chrom_counts[chrom] += 1

print(f"{'Chromosome':<12} {'VNTR Entries':>12}")
print("-" * 26)
total = 0
for chrom in chrom_order:
    if chrom in chrom_counts:
        count = chrom_counts[chrom]
        print(f"{chrom:<12} {count:>12,}")
        total += count
for chrom, count in sorted(chrom_counts.items()):
    if chrom not in chrom_order:
        print(f"{chrom:<12} {count:>12,}")
        total += count
print("-" * 26)
print(f"{'Total':<12} {total:>12,}")

# ---------------------------------------------------------------------------
# 2. Load TSV and filter to VNTRs >= 6 bp
# ---------------------------------------------------------------------------
df = pd.read_csv(
    tsv_path, sep="\t", header=None, usecols=[0, 1, 6, 13, 14, 15],
    names=["chrom", "pos", "cons_len", "bio_context_1", "bio_context_2", "bio_context_3"]
)

# Restrict to standard chromosomes and apply categorical ordering
df = df[df["chrom"].isin(chrom_order)].copy()
df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
df.sort_values("chrom", inplace=True)

df_filtered = df[df["cons_len"] > 6].copy()

# ---------------------------------------------------------------------------
# 3. Statistics before and after filtering
# ---------------------------------------------------------------------------
def print_stats(label, data):
    s = data["cons_len"]
    print(f"\n{label}")
    print("-" * 40)
    print(f"  Count   : {len(s):>10,}")
    print(f"  Mean    : {s.mean():>10.2f} bp")
    print(f"  Median  : {s.median():>10.2f} bp")
    print(f"  Std dev : {s.std():>10.2f} bp")
    print(f"  Min     : {s.min():>10,} bp")
    print(f"  Max     : {s.max():>10,} bp")
    print(f"  Q1      : {s.quantile(0.25):>10.2f} bp")
    print(f"  Q3      : {s.quantile(0.75):>10.2f} bp")


print_stats("VNTR Statistics — Before Filtering", df)
print_stats("VNTR Statistics — After Filtering (> 6 bp)", df_filtered)
print(f"\n  Removed : {len(df) - len(df_filtered):>10,} VNTRs")

# Per-chromosome counts after filtering
present_chroms = [c for c in chrom_order if c in df_filtered["chrom"].values]
filtered_counts = df_filtered["chrom"].value_counts()

print(f"\n{'Chromosome':<12} {'VNTRs (> 6 bp)':>16}")
print("-" * 30)
filtered_total = 0
for chrom in present_chroms:
    count = filtered_counts.get(chrom, 0)
    print(f"{chrom:<12} {count:>16,}")
    filtered_total += count
print("-" * 30)
print(f"{'Total':<12} {filtered_total:>16,}")

# ---------------------------------------------------------------------------
# 4. Unique biological context values (columns 14-16) for filtered VNTRs
# ---------------------------------------------------------------------------
bio_cols = {
    "Column 14 (bio_context_1)": "bio_context_1",
    "Column 15 (bio_context_2)": "bio_context_2",
    "Column 16 (bio_context_3)": "bio_context_3",
}

print("\nUnique biological context values (filtered VNTRs > 6 bp):")
for label, col in bio_cols.items():
    unique_vals = sorted(df_filtered[col].dropna().unique())
    print(f"\n  {label}:")
    for val in unique_vals:
        count = (df_filtered[col] == val).sum()
        print(f"    {val:<20} {count:>10,}")

# ---------------------------------------------------------------------------
# 5. Genotype variance across individuals (filtered VNTRs)
# ---------------------------------------------------------------------------
# Build a set of (chrom, pos) keys from the filtered set for fast lookup
filtered_keys = set(zip(df_filtered["chrom"].astype(str), df_filtered["pos"].astype(str)))

# Stream through the VCF; for each filtered VNTR determine if genotype varies
print("\nScanning VCF for genotype variance (this may take a few minutes)...")
variance_map = {}  # (chrom, pos) -> True (variable) / False (invariant)

with open(vcf_path, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        key = (fields[0], fields[1])
        if key not in filtered_keys:
            continue
        # Collect unique non-missing allele indices across all samples
        alleles_seen = set()
        for gt in fields[9:]:
            for allele in gt.split("/"):
                if allele != ".":
                    alleles_seen.add(allele)
                    if len(alleles_seen) > 1:
                        break  # variance confirmed; no need to scan further
            if len(alleles_seen) > 1:
                break
        variance_map[key] = len(alleles_seen) > 1 if alleles_seen else None

# Map variance back onto df_filtered
df_filtered = df_filtered.copy()
df_filtered["is_variable"] = [
    variance_map.get((str(c), str(p)), None)
    for c, p in zip(df_filtered["chrom"], df_filtered["pos"])
]

n_variable   = (df_filtered["is_variable"] == True).sum()
n_invariant  = (df_filtered["is_variable"] == False).sum()
n_unresolved = df_filtered["is_variable"].isna().sum()

print("\nGenotype Variance Summary (filtered VNTRs > 6 bp):")
print("-" * 40)
print(f"  Variable   (>1 unique allele) : {n_variable:>10,}")
print(f"  Invariant  (1 unique allele)  : {n_invariant:>10,}")
print(f"  Unresolved (all missing GTs)  : {n_unresolved:>10,}")

# Per-chromosome breakdown
print(f"\n{'Chromosome':<12} {'Variable':>10} {'Invariant':>10} {'Unresolved':>12}")
print("-" * 46)
for chrom in present_chroms:
    sub = df_filtered[df_filtered["chrom"] == chrom]
    v = (sub["is_variable"] == True).sum()
    i = (sub["is_variable"] == False).sum()
    u = sub["is_variable"].isna().sum()
    print(f"{chrom:<12} {v:>10,} {i:>10,} {u:>12,}")

# ---------------------------------------------------------------------------
# 6. Per-chromosome histograms of consensus sequence length (filtered)
# ---------------------------------------------------------------------------
n_chroms = len(present_chroms)
ncols = 4
nrows = -(-n_chroms // ncols)  # ceiling division

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
axes = axes.flatten()

for i, chrom in enumerate(present_chroms):
    ax = axes[i]
    data = df_filtered[df_filtered["chrom"] == chrom]["cons_len"]
    ax.hist(data, bins=50, range=(7, 100), color="steelblue", edgecolor="white", linewidth=0.4)
    ax.set_xlim(7, 100)
    ax.set_title(chrom, fontsize=11)
    ax.set_xlabel("Consensus Length (bp)", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.tick_params(labelsize=8)

# Hide any unused subplot panels
for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle(
    "VNTR Consensus Sequence Length Distribution per Chromosome\n"
    "(vamosExpanded v3.0 GRCh38 oriMotifs — > 6 bp)",
    fontsize=13, y=1.01,
)
plt.tight_layout()

plot_path = os.path.join(out_dir, "cons_len_histograms.png")
plt.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to: {plot_path}")

# ---------------------------------------------------------------------------
# 7. Further filter: > 6 bp AND variable (>1 unique allele)
# ---------------------------------------------------------------------------
df_variable = df_filtered[df_filtered["is_variable"] == True].copy()

print(f"\nFurther filtered (> 6 bp + variable): {len(df_variable):>10,} VNTRs")

# Per-chromosome counts
variable_counts = df_variable["chrom"].value_counts()
present_variable_chroms = [c for c in chrom_order if c in df_variable["chrom"].values]

print(f"\n{'Chromosome':<12} {'Variable VNTRs':>16}")
print("-" * 30)
variable_total = 0
for chrom in present_variable_chroms:
    count = variable_counts.get(chrom, 0)
    print(f"{chrom:<12} {count:>16,}")
    variable_total += count
print("-" * 30)
print(f"{'Total':<12} {variable_total:>16,}")

# ---------------------------------------------------------------------------
# 8. Per-chromosome histograms — variable VNTRs only
# ---------------------------------------------------------------------------
n_chroms_v = len(present_variable_chroms)
ncols = 4
nrows = -(-n_chroms_v // ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
axes = axes.flatten()

for i, chrom in enumerate(present_variable_chroms):
    ax = axes[i]
    data = df_variable[df_variable["chrom"] == chrom]["cons_len"]
    ax.hist(data, bins=50, range=(7, 100), color="steelblue", edgecolor="white", linewidth=0.4)
    ax.set_xlim(7, 100)
    ax.set_title(chrom, fontsize=11)
    ax.set_xlabel("Consensus Length (bp)", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.tick_params(labelsize=8)

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle(
    "VNTR Consensus Sequence Length Distribution per Chromosome\n"
    "(vamosExpanded v3.0 GRCh38 oriMotifs — > 6 bp, variable only)",
    fontsize=13, y=1.01,
)
plt.tight_layout()

plot_path_v = os.path.join(out_dir, "cons_len_histograms_variable.png")
plt.savefig(plot_path_v, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to: {plot_path_v}")

# ---------------------------------------------------------------------------
# 9. Export variable VNTRs to CSV with per-sample haplotype copy numbers
# ---------------------------------------------------------------------------
variable_export_keys = set(
    zip(df_variable["chrom"].astype(str), df_variable["pos"].astype(str))
)

tsv_col_names = [
    "chrom", "start", "end", "motifs", "cluster_type", "vntr_type",
    "cons_len", "cons_seq", "col_8", "col_9", "col_10", "col_11",
    "genomic_context", "bio_context_1", "bio_context_2", "bio_context_3",
]

# Read TSV rows for variable VNTRs into memory
print("\nReading TSV metadata for variable VNTRs...")
tsv_data = {}
with open(tsv_path, "r") as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        key = (fields[0], fields[1])
        if key in variable_export_keys:
            tsv_data[key] = fields

# Read sample names from VCF header
sample_names = []
with open(vcf_path, "r") as f:
    for line in f:
        if line.startswith("#CHROM"):
            sample_names = line.rstrip("\n").split("\t")[9:]
            break

# Build CSV header: TSV cols | sample_hap1, sample_hap2 (interleaved) | summary stats
sample_cols = [col for s in sample_names for col in (f"{s}_hap1", f"{s}_hap2")]
csv_header = tsv_col_names + sample_cols + [
    "mean_hap1", "mean_hap2", "var_hap1", "var_hap2",
    "mean_combined", "var_combined", "h_exp", "h_obs",
]

csv_path = os.path.join(out_dir, "variable_vntrs.csv")
print(f"Streaming VCF and writing CSV (second pass — may take a few minutes)...")

n_written = 0
summary_rows = []  # collect TSV metadata + stats only (for Excel export)

with open(csv_path, "w", newline="") as out_f:
    writer = csv.writer(out_f)
    writer.writerow(csv_header)

    with open(vcf_path, "r") as vcf_f:
        for line in vcf_f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            key = (fields[0], fields[1])
            if key not in variable_export_keys or key not in tsv_data:
                continue

            # Parse ALTANNO: allele index (0-based) → copy number (element count)
            allele_copy_num = {}
            for part in fields[7].split(";"):
                if part.startswith("ALTANNO="):
                    for idx, entry in enumerate(part[8:].split(",")):
                        allele_copy_num[idx] = len(entry.split("-"))
                    break

            # Per-sample hap1/hap2 copy numbers
            hap1_cols, hap2_cols = [], []
            hap1_vals, hap2_vals, combined_vals = [], [], []
            allele_counts = Counter()  # copy-number allele frequencies (both haplotypes)
            n_complete = 0             # individuals with both haplotypes called
            n_het = 0                  # individuals where cn1 != cn2

            for gt in fields[9:]:
                a1_str, a2_str = gt.split("/")
                cn1 = allele_copy_num.get(int(a1_str)) if a1_str != "." else None
                cn2 = allele_copy_num.get(int(a2_str)) if a2_str != "." else None

                if cn1 is not None and cn2 is not None:
                    hap1, hap2 = max(cn1, cn2), min(cn1, cn2)
                    combined_vals.append(cn1 + cn2)
                    n_complete += 1
                    if cn1 != cn2:
                        n_het += 1
                elif cn1 is not None:
                    hap1, hap2 = cn1, ""
                elif cn2 is not None:
                    hap1, hap2 = cn2, ""
                else:
                    hap1, hap2 = "", ""

                if cn1 is not None:
                    allele_counts[cn1] += 1
                if cn2 is not None:
                    allele_counts[cn2] += 1

                hap1_cols.append(hap1)
                hap2_cols.append(hap2)
                if hap1 != "":
                    hap1_vals.append(hap1)
                if hap2 != "":
                    hap2_vals.append(hap2)

            mean_hap1     = round(sum(hap1_vals) / len(hap1_vals), 4) if hap1_vals else ""
            mean_hap2     = round(sum(hap2_vals) / len(hap2_vals), 4) if hap2_vals else ""
            var_hap1      = round(statistics.variance(hap1_vals), 4) if len(hap1_vals) > 1 else ""
            var_hap2      = round(statistics.variance(hap2_vals), 4) if len(hap2_vals) > 1 else ""
            mean_combined = round(sum(combined_vals) / len(combined_vals), 4) if combined_vals else ""
            var_combined  = round(statistics.variance(combined_vals), 4) if len(combined_vals) > 1 else ""

            total_alleles = sum(allele_counts.values())
            h_exp = round(1 - sum((c / total_alleles) ** 2 for c in allele_counts.values()), 4) \
                    if total_alleles > 0 else ""
            h_obs = round(n_het / n_complete, 4) if n_complete > 0 else ""

            sample_vals = [v for pair in zip(hap1_cols, hap2_cols) for v in pair]
            writer.writerow(tsv_data[key] + sample_vals + [mean_hap1, mean_hap2, var_hap1, var_hap2, mean_combined, var_combined, h_exp, h_obs])
            summary_rows.append(tsv_data[key] + [mean_hap1, mean_hap2, var_hap1, var_hap2, mean_combined, var_combined, h_exp, h_obs])
            n_written += 1

print(f"CSV saved to: {csv_path}  ({n_written:,} rows)")

# ---------------------------------------------------------------------------
# 10. Export summary Excel (TSV metadata + stats, no per-sample columns)
# ---------------------------------------------------------------------------
summary_col_names = tsv_col_names + ["mean_hap1", "mean_hap2", "var_hap1", "var_hap2", "mean_combined", "var_combined", "h_exp", "h_obs"]
df_summary = pd.DataFrame(summary_rows, columns=summary_col_names)

xlsx_path = os.path.join(out_dir, "variable_vntrs_summary.xlsx")
df_summary.drop(columns=["motifs"]).to_excel(xlsx_path, index=False)
print(f"Summary Excel saved to: {xlsx_path}  ({len(df_summary):,} rows)")

# ---------------------------------------------------------------------------
# 11. Export copy-number-polymorphic VNTRs (var_hap1 > 0 or var_hap2 > 0)
# ---------------------------------------------------------------------------
cn_var_mask = (
    (pd.to_numeric(df_summary["var_hap1"], errors="coerce") > 0) |
    (pd.to_numeric(df_summary["var_hap2"], errors="coerce") > 0)
)
df_cn_variable = df_summary[cn_var_mask].copy()

xlsx_cn_path = os.path.join(out_dir, "copy_number_variable_vntrs_summary.xlsx")
df_cn_variable.drop(columns=["motifs"]).to_excel(xlsx_cn_path, index=False)
print(f"Copy-number-variable Excel saved to: {xlsx_cn_path}  ({len(df_cn_variable):,} rows)")

# ---------------------------------------------------------------------------
# 12. Per-chromosome histograms — copy-number-variable VNTRs only
# ---------------------------------------------------------------------------
df_cn_variable["cons_len"] = pd.to_numeric(df_cn_variable["cons_len"], errors="coerce")
present_cn_chroms = [c for c in chrom_order if c in df_cn_variable["chrom"].values]

n_chroms_cn = len(present_cn_chroms)
ncols = 4
nrows = -(-n_chroms_cn // ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 4))
axes = axes.flatten()

for i, chrom in enumerate(present_cn_chroms):
    ax = axes[i]
    data = df_cn_variable[df_cn_variable["chrom"] == chrom]["cons_len"].dropna()
    ax.hist(data, bins=50, range=(7, 100), color="steelblue", edgecolor="white", linewidth=0.4)
    ax.set_xlim(7, 100)
    ax.set_title(chrom, fontsize=11)
    ax.set_xlabel("Consensus Length (bp)", fontsize=9)
    ax.set_ylabel("Count", fontsize=9)
    ax.tick_params(labelsize=8)

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle(
    "VNTR Consensus Sequence Length Distribution per Chromosome\n"
    "(vamosExpanded v3.0 GRCh38 oriMotifs — copy-number-variable only)",
    fontsize=13, y=1.01,
)
plt.tight_layout()

plot_path_cn = os.path.join(out_dir, "cons_len_histograms_cn_variable.png")
plt.savefig(plot_path_cn, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to: {plot_path_cn}")
