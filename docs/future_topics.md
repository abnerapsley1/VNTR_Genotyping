# Future Topics

A running list of known issues, design gaps, and planned improvements to address in future development.

---

## ~~Gene-level normalization: alternate contig support~~ — Resolved

Gene background read counts now include alt-contig reads via the same PSL-based
block-level mapping used for VNTRs. The `--psl` flag in `count_reads.py` defaults
to the stored `altSeqLiftOverPsl.txt.gz`; use `--no-alt-contigs` to disable for
primary-only references.

---

## Tool naming

**Context:** The tool currently has working names for all three interfaces (`vntr_genotyping` Python package, `vntrGenotyping` R package, `vntr-count` CLI entry point) chosen for clarity during development. A proper name should be decided before any public release.

**Things to consider:**
- A single, memorable name that works across Python, R, and the CLI (e.g. the Python/R package names and CLI command should ideally share a common root).
- Whether the name should reflect the method (density ratio / read counting), the target (VNTRs / tandem repeats), or be more brand-like.
- Namespace conflicts — check PyPI and CRAN for existing packages before committing to a name.
- Once chosen, the name will need to be updated in `pyproject.toml`, the R `DESCRIPTION`, `__init__.py`, `cli.py`, and the README.

---

## Multi-reference support: GRCh38 and CHM13

**Context:** The tool currently assumes GRCh38 (hg38) coordinates in both the built-in BED file (`dat1_regions.bed`) and any user-supplied BED files. CHM13 (T2T-CHM13) is an increasingly used reference that has different chromosome coordinates and contig naming conventions.

**What needs to happen:**
- The built-in BED file (and any future bundled BED files) would need a CHM13 coordinate set alongside the GRCh38 one, or maintained as separate versioned BED files.
- A `--reference-build` (or similar) flag should be added to `count_reads.py` so the tool selects the appropriate built-in BED and, if normalization is enabled, expects a matching GTF.
- Chromosome naming differences (e.g., CHM13 uses `chr5` but contig names and coordinates differ) must be accounted for — alternate contig IDs will also differ between builds.
- User documentation should clearly state which reference build a custom BED file's coordinates correspond to.

---

## ~~`assign_nearest_genes`: ambiguous multi-gene overlap handling~~ — Resolved

`assign_nearest_genes` now uses a three-level tie-breaker (max bp overlap →
protein_coding preferred → alphabetical) and emits a warning on any surviving
ties. `parse_gtf` stores `gene_type` for all genes.

Additional improvements also implemented:
- **`--norm-gene-types`** flag (default: `protein_coding`) restricts which biotypes
  are eligible as the nearest-gene fallback. Use `any` to restore the old
  type-agnostic behaviour, or list multiple types (e.g. `protein_coding lncRNA`).
- **Coordinate-format normalization** — BED column 5 now accepts a coordinate string
  (`chrN:start-end`) as a direct normalization window instead of a gene name.
  PSL alt-contig expansion is applied automatically unless `--no-alt-contigs` is set.

---
