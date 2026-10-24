# vntrlen (draft)

Absolute **combined** VNTR length (both alleles added together) from **short-read**
whole-genome sequencing, with a 90% interval and a confidence tier for every call.

* `vntrlen genotype` — genotype samples at repeats that have a model
* `vntrlen build` — build models for your own repeats from samples with long-read truth
* `vntrlen models` — list every reference-repeat model available, or install the shipped models

Status: **draft** (`0.1.0.dev0`), branch `vntr-package`.

## Install

```bash
conda env create -f environment.yml      # includes bwa + samtools for FASTQ input
conda activate vntrlen
vntrlen models install --bundle vntrlen-shipped-GRCh38-v1.tar.gz   # or with no --bundle once a URL is published
```

Data live in `$VNTRLEN_HOME` (default `~/.vntrlen`): models, downloaded references and a feature cache.

## Shipped models

* Reference: GRCh38 as used by the 1000 Genomes NYGC 30x data
  (`GRCh38_full_analysis_set_plus_decoy_hla`; downloaded on first use if `--reference` is omitted).
* 60,504 VNTRs (the default BED). Every repeat has its own calibration curve: the 11 short-read
  signals are combined per repeat by ridge regression, and a monotone (isotonic) curve with
  Theil-Sen extensions maps the combined score to length.
* Test people (52, never used in training): mean absolute error 0.83 repeat units over all repeats.
* Trained on 212 people; 90% intervals calibrated on 26; per-repeat performance from 52 held-out
  test people (`vntrlen models list --tsv models.tsv`, columns `perf_*`).
* Ground truth: vamos long-read genotypes (TRCompDB v3.0, Gu & Chaisson).

## Genotyping

```bash
vntrlen genotype --reads sample1.cram sample2.bam --out results/
vntrlen genotype --samples sheet.tsv --bed my_repeats.bed --reference ref.fa --out results/
```

Inputs are checked before anything runs: BAM/CRAM headers must match the reference (names and
lengths), reads must be short (median <= 300 bp; long-read platforms are rejected), and paired-
vs single-end is detected. FASTQ input is aligned with BWA-MEM. Single-end samples are supported,
but their calls are flagged `single_end_unvalidated` (the models were trained on paired-end data).

`results/genotypes.tsv` — one row per sample and repeat:

| column | meaning |
|---|---|
| `sample_id`, `locus_id` | sample; repeat as `chrom:start-end` (1-based) |
| `chrom`, `start`, `end`, `period` | repeat coordinates (BED, 0-based start) and motif length |
| `est_units`, `est_bp` | point estimate of the combined length, repeat units and bp |
| `ci90_lo/hi_units`, `ci90_lo/hi_bp` | 90% interval (whole repeat units; bp = units x bp per unit) |
| `tier` | `PASS`, `LOW_CONF` (interval wider than the estimate, or signal outside the curve's training range), `NA` |
| `flags` | `single_end_unvalidated`, `signal_outside_training_range` |
| `model_set`, `model_version` | which model produced the call |

Also `genotypes_matrix.tsv` (repeats x samples), `skipped_repeats.tsv` (BED repeats with no model),
`run_info.json`.

## Building models for your repeats

```bash
vntrlen build --samples sheet.tsv --bed my_repeats.bed --truth truth.tsv --name my_set --out build_report/
```

* `truth.tsv`: `sample_id`, `locus_id` (or `chrom`/`start`/`end`), `length`, `unit` (`units` or `bp`),
  optional `allele1`/`allele2` (enable the allele-count and heterozygosity features).
* `--annotation`: optional per-repeat values (`segdup_identity`, `te_frac_locus`, ...) overriding the
  bundled TRCompDB lookup.
* Each repeat needs ground truth in >= 20 people (50+ recommended).
* Each repeat gets its own calibration curve. The shipped cross-repeat correction models are reused
  and the 90% intervals recalibrated on your data; with >= 50 people and >= 1,000 repeats the
  correction models are retrained (`--retrain auto|never|always`).
* All reported performance is 5-fold cross-validated over people (families kept together via the
  sheet's `family_id`): `build_report/performance_per_repeat.tsv` and `performance_summary.html`.
* The new model set is saved to the store and used automatically by `vntrlen genotype` (the newest
  user set wins for a repeat; pick one with `--model-set`).

## File formats

**BED** (tab-separated, no header): `chrom  start  end  motif_length  motifs[  name]` — `motifs` is a
comma-separated list of motif sequences; the first is treated as the consensus.

**Sample sheet** (tab-separated, header): `sample_id` plus `alignment` (BAM/CRAM), or `fastq1`
(+ `fastq2` for paired-end); optional `family_id`. Alternatively pass `--reads` with BAM/CRAM paths
(sample IDs from the `@RG SM` tag, else the file name).
