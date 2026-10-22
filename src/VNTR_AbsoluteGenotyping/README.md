# VNTR Absolute Genotyping

A pipeline + tool to estimate **absolute VNTR copy number / array length** from
short-read alignments (BAM/CRAM), for a user-supplied coordinate list, with
calibrated uncertainty.

This directory is self-contained and does **not** modify the existing
`VNTR_GenotypingTool/` (the density-ratio tool). It supersedes it methodologically.

---

## Why this exists

The density-ratio tool assumes normalized read count scales linearly with
combined allele length. Validation against long-read ground truth showed this
overestimates absolute copy number for ~97% of loci (median 4.4x) and correlates
per-locus for <1%. See `../../` project memory (`tool-performance-diagnosis`).

Root causes: mapping saturation (short reads stop tracking dosage as arrays
grow), paralog read recruitment, and a hard-coded `2 * ref_copies` scaling with
no per-locus calibration.

---

## Approach (two-tier, panel-anchored)

1. **Rich signal extraction** — one pass per locus emitting depth (MAPQ-tiered),
   in-repeat-read counts, spanning read-pair insert-size shifts, motif-composition
   profiles, soft-clip pileups. Not just the depth ratio.
2. **Per-locus empirical calibration** — for every locus in the bundled
   annotation, fit `combined_length = f(signals)` against the vamos long-read
   panel (~70 samples with short-read BAMs). Ship the curves.
3. **Calibration transfer** — a model trained across the ~60k bundled loci that
   predicts the calibration from locus covariates (mappability, GC, motif
   entropy, segdup, reference copy number, ...), for **novel user loci** not in
   the panel. With an out-of-distribution check that downgrades confidence.
4. **Global residual model** — gradient-boosted regression over all
   (locus, sample) rows, using the Phase-2 prediction as a feature, with quantile
   heads for prediction intervals.
5. **Per-allele module (opt-in)** — spanning-pair insert-size mixture +
   motif-composition read clustering; emitted only when confidence-gated.

Primary output: **combined array length** (motif units + bp) with a CI and a
`PASS / LOW_CONF / NA` tier. Per-allele is a gated secondary output.

---

## Ground truth

`data/1000Genomes/vamosExpanded_v3.0_GRCh38_oriMotifs.vcf` (same call set the
old validation used). Per-allele copy number = token count of the allele's
`ALTANNO` entry. **vamos GT is 1-based**: GT value `k` -> `ALTANNO[k-1]`
(the old pipeline got this wrong — do not reuse its GT CSVs).

---

## Data flow

```
config/samples.tsv  (you fill: sample -> URL -> vamos column)
config/references.tsv (you fill: reference FASTA per alignment set)
          |
          v
[ground_truth] parse_vamos.py ....... vamos VCF  -> gt_long.parquet (locus x sample lengths)
[covariates]   build_locus_covariates .. BED + mutRate + anova + reference -> locus_covariates.parquet
          |
          v
[run] run_all.py  (per sample, one at a time):
        download the CRAM (parallel HTTP ranges; next sample prefetched)
        or remote-slice per locus (runner.transfer_mode)
        -> extract signals for ~60k loci + random background windows
        -> features/<sample>.parquet ; delete the local CRAM
          |
          v
[model] splits.py -> per_locus_calibration.py -> calibration_transfer.py
        -> global_residual.py -> evaluate.py  (frozen locus x sample test set)
          |
          v
      model artifacts + evaluation report
```

---

## Compute target

32 CPU cores (the 70-sample run used a box with 2x RTX 6000 Ada). Feature extraction is CPU-parallel over loci
(32 workers). Global residual model + any NN components use GPU.

## Setup

```bash
mamba env update -f src/VNTR_AbsoluteGenotyping/environment.yml   # adds samtools, xgboost, pyarrow, ...
```

## Quickstart

```bash
# 0. environment (adds pyarrow, scikit-learn, joblib, libcurl-htslib, ...)
mamba env update -f src/VNTR_AbsoluteGenotyping/environment.yml
conda activate vntr_analysis
cd src/VNTR_AbsoluteGenotyping

# 1. ground truth + covariates  (local files only; ~2 min total)
python -m pipeline.cli ground-truth
python -m pipeline.cli covariates

# 2. fill in the manifests (copy the templates first)
cp config/samples.template.tsv    config/samples.tsv
cp config/references.template.tsv config/references.tsv
python -m pipeline.cli list-vamos > /tmp/vamos_columns.txt   # 416 column headers
$EDITOR config/samples.tsv        # one row per short-read sample
$EDITOR config/references.tsv     # FASTA per alignment set (CRAM needs this)

# 3. feature extraction — one sample at a time (download mode: ~10 min/sample)
python -m pipeline.cli extract                 # all pending samples
#   or: python -m pipeline.cli extract-one --sample HG00733

# 4. build + evaluate the model chain
python -m pipeline.cli all-model
cat work/model/evaluation.json
```

`work/` holds every artifact; feature extraction is resumable (skips samples with
an existing `features/<id>.parquet`).

## Status

Run end-to-end on real data (2026-10-08): 70 short-read samples, 60.5k loci.
On held-out individuals at bundled loci: mean abs error 1.02 motif units (vs
1.76 for the panel mean), within-locus r 0.76; on held-out loci: 2.77 units
(vs 9.40 naive), r 0.59; 90% intervals cover 91%.  Per-allele output is
effectively combined ÷ 2 — short reads can't separate most alleles.  See
`HANDOFF.md` §2 and `PROGRESS.md`.
