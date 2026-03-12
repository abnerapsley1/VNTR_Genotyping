# VNTR Genotyping

A genome-wide tool for quantifying VNTR (Variable Number Tandem Repeat) read coverage from short-read alignment files (BAM/CRAM). Coverage at each VNTR is expressed as a density ratio relative to the host gene, providing a sample-comparable metric that accounts for differences in sequencing depth and regional GC bias.

---

## How it works

For each VNTR region and each input alignment file, the tool:

1. Counts the unique reads overlapping the VNTR (primary chromosome + any alt scaffold coordinates).
2. Counts the unique reads overlapping the host gene region, excluding reads already counted in any VNTR within that gene.
3. Computes the **density ratio**: `(VNTR reads / VNTR length) / (gene reads / gene length)`

A ratio of 1.0 indicates average coverage relative to the gene. Values above 1.0 indicate enrichment (more copies); values below 1.0 indicate depletion (fewer copies).

When no GTF annotation is provided, raw read counts are returned instead of density ratios.

---

## Installation

### Prerequisites

- Python ≥ 3.9
- `pysam` and `pandas` (installed automatically)
- A reference FASTA is required for CRAM input without an embedded reference

### Install the Python package

```bash
pip install src/VNTR_GenotypingTool/Python
```

This also installs the `vntr-count` command-line tool.

### Install the R package

The R package wraps the Python package via [reticulate](https://rstudio.github.io/reticulate/). Install the Python package first, then:

```r
devtools::install("src/VNTR_GenotypingTool/R/vntrGenotyping")
```

---

## Usage

### Command line

```bash
# All built-in VNTRs, with density ratio normalization
vntr-count --default \
    -i sample.cram \
    -g gencode.v38.annotation.gtf.gz \
    -T hg38.fa \
    -o results.csv

# Subset to VNTRs in specific genes
vntr-count --gene SLC6A3 MAOA \
    -i sample1.cram sample2.cram \
    -g gencode.v38.annotation.gtf.gz \
    -T hg38.fa \
    -o results.csv

# Specific VNTR IDs
vntr-count --vntr TR30 TR22 TR17 \
    -i sample.cram \
    -g gencode.v38.annotation.gtf.gz \
    -T hg38.fa \
    -o results.csv

# Custom BED file (no GTF — returns raw read counts)
vntr-count -r my_vntrs.bed \
    -i sample.bam \
    -o counts.csv

# Custom BED combined with built-in VNTRs
vntr-count -r my_vntrs.bed --gene SLC6A3 \
    -i sample.cram \
    -g annotation.gtf.gz \
    -T hg38.fa \
    -o results.csv
```

**Key options:**

| Option | Description |
|--------|-------------|
| `--default` | Use all VNTRs in the built-in annotation |
| `--gene GENE [...]` | Subset built-in VNTRs by gene name |
| `--vntr ID [...]` | Subset built-in VNTRs by VNTR ID |
| `-r / --regions BED` | Custom BED file (combinable with above) |
| `-i FILE [...]` | Input BAM/CRAM/SAM file(s) |
| `-g / --gtf GTF` | Gene annotation GTF/GTF.gz (enables density ratio) |
| `-T / --reference FASTA` | Reference FASTA (required for CRAM) |
| `-o / --output CSV` | Output file (default: `read_counts.csv`) |
| `--no-alt-contigs` | Disable alt-scaffold normalization |
| `--norm-gene-types TYPE [...]` | Biotypes eligible for nearest-gene assignment (default: `protein_coding`) |
| `--psl PSL` | Custom UCSC altSeqLiftOverPsl file (default: bundled GRCh38) |

---

### Python

```python
import vntr_genotyping as vg

# All built-in VNTRs — returns a pandas DataFrame
df = vg.count_vntrs(
    "sample.cram",
    default=True,
    gtf="gencode.v38.annotation.gtf.gz",
    reference="hg38.fa",
)

# Subset to specific genes, multiple samples, save CSV
df = vg.count_vntrs(
    ["sample1.cram", "sample2.cram"],
    gene=["SLC6A3", "MAOA"],
    gtf="gencode.v38.annotation.gtf.gz",
    reference="hg38.fa",
    output_csv="results.csv",
)

# Custom BED, raw counts (no GTF)
df = vg.count_vntrs(
    "sample.bam",
    regions="my_vntrs.bed",
)

print(df)
#    sample         metric    TR30    TR22    TR17
# 0  sample  density_ratio  1.2341  0.9812  1.0543
```

**`count_vntrs()` parameters** mirror the CLI options:

| Parameter | Type | Description |
|-----------|------|-------------|
| `input_files` | `str` or `list[str]` | BAM/CRAM/SAM path(s) |
| `default` | `bool` | Use all built-in VNTRs |
| `gene` | `list[str]` | Filter by gene name |
| `vntr` | `list[str]` | Filter by VNTR ID |
| `regions` | `str` | Custom BED file path |
| `gtf` | `str` | GTF annotation path |
| `reference` | `str` | Reference FASTA path |
| `psl` | `str` | altSeqLiftOverPsl path |
| `no_alt_contigs` | `bool` | Disable alt-scaffold normalization |
| `norm_gene_types` | `list[str]` | Biotypes for nearest-gene assignment |
| `output_csv` | `str` | Optional path to save CSV |

Returns a `pandas.DataFrame` with columns `sample`, `metric`, and one column per VNTR.

---

### R

```r
library(vntrGenotyping)

# All built-in VNTRs
df <- count_vntrs(
  "sample.cram",
  default   = TRUE,
  gtf       = "gencode.v38.annotation.gtf.gz",
  reference = "hg38.fa"
)

# Subset to specific genes, multiple samples, save CSV
df <- count_vntrs(
  c("sample1.cram", "sample2.cram"),
  gene       = c("SLC6A3", "MAOA"),
  gtf        = "gencode.v38.annotation.gtf.gz",
  reference  = "hg38.fa",
  output_csv = "results.csv"
)

# Custom BED, raw counts (no GTF)
df <- count_vntrs("sample.bam", regions = "my_vntrs.bed")

head(df)
#     sample        metric    TR30    TR22    TR17
# 1   sample density_ratio  1.2341  0.9812  1.0543
```

The R function accepts the same parameters as the Python API. The returned object is a standard R `data.frame`.

**Configuring the Python environment:**

By default, reticulate uses whichever Python has `vntr_genotyping` installed. To specify a particular environment:

```r
reticulate::use_condaenv("myenv", required = TRUE)
# or
reticulate::use_virtualenv("~/.venvs/vntr", required = TRUE)
library(vntrGenotyping)
```

---

## Output format

Results are returned as a wide-format table with one row per sample:

| sample | metric | TR30 | TR22 | TR17 | ... |
|--------|--------|------|------|------|-----|
| sample1 | density_ratio | 1.234 | 0.981 | 1.054 | ... |
| sample2 | density_ratio | 1.102 | 1.203 | 0.876 | ... |

- `metric` is `density_ratio` when a GTF is provided, or `read_count` otherwise.
- Cells where normalization could not be computed (e.g., zero gene reads) contain `NA`.

---

## Custom BED file format

Custom BED files follow standard BED format with two optional extensions:

```
chrom  start  end  name  gene  [alt_chrom  alt_start  alt_end  ...]
```

- **Columns 1–3** (required): chromosome, 0-based start, end.
- **Column 4** (recommended): VNTR identifier.
- **Column 5** (optional): gene name for normalization, or a coordinate string
  (e.g. `chr5:1318298-1509772`) to use an explicit normalization window.
  Use `.` to leave unspecified; the nearest protein-coding gene will be assigned automatically.
- **Columns 6–8, 9–11, ...** (optional): alt-contig triplets (`chrom start end`)
  for including reads aligned to alternate scaffolds.

---

## Built-in annotation

The built-in VNTR annotation (`data/VNTRAnnotations/v1_VNTR_Annotations_03_07_2026.bed`)
covers GRCh38 (hg38) coordinates and contains ~60,500 genome-wide VNTR regions.
Alternate scaffold coordinates are pre-computed from the UCSC `altSeqLiftOverPsl.txt.gz`
alignment for the GRCh38 reference assembly.
