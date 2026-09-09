# Ground-truth spec

**Source:** `data/1000Genomes/vamosExpanded_v3.0_GRCh38_oriMotifs.vcf` (10 GB,
416 diploid samples). This is the same call set the old tool used
(`src/1000Genomes_Analysis/1000Genomes_Analysis.py:10`). Not committed (size);
present locally. The new pipeline (`pipeline/ground_truth/parse_vamos.py`) parses
it into `work/ground_truth/`.

## Per-allele copy number

= number of `-`-separated tokens in that allele's `ALTANNO` entry (INFO field).
`ALTANNO` is a comma-separated list, one entry per allele.

## Per-sample genotype

`GT` is `a1/a2`. `cn_i = token_count(ALTANNO[a_i - 1])` (see 1-based note below);
`.` = uncalled. Convention: `hap_hi = max(cn1,cn2)`, `hap_lo = min`,
`combined = cn1 + cn2`. **Combined is defined only when BOTH haplotypes are
called** (`n_haps_called == 2`); the model table filters to these.

## Motif bp lengths

`RU` (INFO) is a comma-separated list of motif sequences. Allele bp length =
Σ len(RU[token]) over that allele's ALTANNO tokens. `parse_vamos.py` emits both
motif-unit and bp lengths.

## ⚠ vamos GT is 1-based — the OLD pipeline has an off-by-one bug

Verified over 50,000 records on 2026-08-31: `max GT index == len(ALTANNO)` in
every record, never exceeds; GT allele `0` never occurs. So **GT value `k` maps
to `ALTANNO[k-1]`**.

The old `1000Genomes_Analysis.py` uses `allele_copy_num[int(gt)]` with a dict
keyed `0..N-1`, so it:
(a) reads the **wrong allele's** token count for every call, and
(b) maps `GT == N` (the highest allele, present in EVERY record) to `None`,
silently dropping that sample as uncalled.

**Consequence:** the old `variable_vntrs.csv`, `cn_variable_genotypes_*.csv`, and
all `vntr_correlations*.csv` ground truth are corrupted. Do **not** reuse them.
The new `parse_vamos.py` uses `ALTANNO[gt-1]` correctly.

## BED ↔ vamos join

The bundled BED stores `start == vamos POS` and `end == vamos INFO.END` with
**zero offset** — despite the BED header claiming "0-based half-open", `start` is
effectively 1-based here. Join is exact on `(chrom, start, end)`; `parse_vamos.py`
also has an overlap fallback and reports the match rate (full run: 60,504/60,504
exact, 0 fallback).

For 0-based read fetching, `loci_targets.parquet` carries
`array_start0 = start - 1`, `array_end0 = end`. Some BED rows are degenerate
(`start == end`, e.g. `chr1:2682949-2682949`) — extraction returns
`note="degenerate_interval"`.

## Sample-column naming

Columns are `APR-phase1__aprNNN_hap1/APR-phase1__aprNNN_hap2`. Clean sample id =
strip trailing `_hap[12]` and everything after `/`. The manifest's
`vamos_column` field must be the **exact** header string; get the list via
`python -m pipeline.cli list-vamos`.

## Old-pipeline filter (for reference)

The old pipeline restricted to motif length > 6 bp, and the 60,519-locus bundled
BED is the `var_hap1 > 0 | var_hap2 > 0` (copy-number-variable) subset of the
full vamos catalog. The new pipeline inherits that BED as its locus set.
