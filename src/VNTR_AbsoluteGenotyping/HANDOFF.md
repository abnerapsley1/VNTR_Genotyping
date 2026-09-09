# HANDOFF — read this first

Written 2026-09-09 to hand this project to a fresh session with no chat history.
Everything needed to continue is in this directory + `docs/context/`.

---

## 1. What we're building and why

**Goal:** a tool that estimates **absolute VNTR array length** from short-read
alignments (BAM/CRAM) for a **user-supplied coordinate list** (BED, like the
existing `data/VNTRAnnotations/v1_VNTR_Annotations_03_07_2026.bed`), with
**calibrated uncertainty**.

- Primary output: **combined array length** (both alleles summed) in motif units
  and bp, with a 90% interval and a `PASS / LOW_CONF / NA` confidence tier.
- Secondary, opt-in, confidence-gated: **per-allele** lengths.
- Method-agnostic by intent (naive / regression / ML / simulation all on the table).

**Why a new tool:** the existing density-ratio tool (`src/VNTR_GenotypingTool/`)
assumes normalized read count scales linearly with combined length. Validation
against long-read ground truth showed it **overestimates absolute copy number for
97% of loci (median 4.4x)** and produces a usable per-locus rank correlation for
only ~0.6% of loci. Full diagnosis: `docs/context/01-old-tool-diagnosis.md`.
The new pipeline is self-contained and does not touch the old tool.

---

## 2. Current state (2026-09-09)

**All pipeline stages are BUILT.** 25 modules, one CLI:
`python -m pipeline.cli <command>`.

| Stage | Status |
|---|---|
| `ground-truth` (parse vamos VCF) | ✅ **ran on the full 10 GB VCF** — 60,504/60,504 loci matched exact; 25.2 M-row `work/ground_truth/gt_long.parquet` |
| `covariates` (per-locus covariate table) | ✅ **ran** — `work/covariates/locus_covariates.parquet`, 60,504 × 78 cols |
| feature extraction (`signals.py`, `nuisance.py`) | ✅ built; validated on the toy BAM `data/Samples/test_sample.bam` only — **not yet run on a real 30x BAM** |
| runner (`run/*.py`) | ✅ built; **not yet exercised against a live remote CRAM** |
| model chain (`model/*.py`) | ✅ built; **full 7-step chain validated end-to-end on synthetic data** (mechanics only — synthetic signal was near-noise-free) |

`PROGRESS.md` has the detailed checklist.

### The one thing blocking real results

The `config/samples.tsv` manifest is not filled in. ~70 samples have short-read
BAMs available (HPRC + HGSVC2). Someone needs to fill in, per sample: the
download URL, the index URL, and which of the 416 vamos VCF sample-columns it
corresponds to.

---

## 3. Exact next steps

```bash
cd src/VNTR_AbsoluteGenotyping

# 0. environment  (adds scikit-learn, joblib, libcurl-htslib, optional xgboost/lightgbm)
mamba env update -f environment.yml
conda activate vntr_analysis

# 1. regenerate ground truth + covariates if work/ is absent (local files only, ~2 min)
python -m pipeline.cli ground-truth
python -m pipeline.cli covariates

# 2. fill the manifests
cp config/samples.template.tsv    config/samples.tsv
cp config/references.template.tsv config/references.tsv
python -m pipeline.cli list-vamos > /tmp/vamos_columns.txt   # 416 exact column headers
#   -> edit config/samples.tsv    (one row per short-read sample; ~70 rows)
#   -> edit config/references.tsv (FASTA per alignment set; required for CRAM)

# 3. feature extraction — one sample at a time, remote-sliced, only the index on disk
python -m pipeline.cli extract               # all pending; resumable; log at work/features/_run_log.tsv
#   run this on the NVIDIA box (2x A100, 32 cores). ~60k remote range-reads per sample.

# 4. build + score the model chain
python -m pipeline.cli all-model             # split -> assemble -> calibrate -> transfer -> global -> per-allele -> evaluate
cat work/model/evaluation.json
```

### After the first real evaluation, the tuning knobs are:

- `pipeline/model/per_locus_calibration.py::_SIGNAL_SHORTLIST` — which 1-D signal
  each per-locus calibration uses.
- `pipeline/model/global_residual.py` — the `interval_width_rel > 1.0` threshold
  that splits `PASS` vs `LOW_CONF`.
- `pipeline/model/per_allele.py` — the per-allele confidence threshold.
- `config/config.yaml` — flank windows, MAPQ tiers, `min_panel_samples_per_locus`,
  holdout fractions, target space (`log1p_units` default).
- Consider `export VNTR_GBM=xgboost` to use the A100s for the GBMs.

### Phase-5 ideas not yet built (see `docs/context/02-approach-and-alternatives.md`)

Graph-realignment front-end (biggest principled fix for the mismapping that sank
the old tool); simulation-based per-locus calibration curves for novel loci with
no panel data.

---

## 4. Architecture in one paragraph

**Two-tier, panel-anchored calibration.** (1) Rich signal extraction — one
BAM pass per locus emitting depth (MAPQ-tiered) × flank windows, in-repeat-read
counts, spanning read-pair insert-size mixture, soft-clip pileups, motif
composition. (2) **Per-locus empirical calibration** (`per_locus_calibration.py`)
— for each bundled locus, fit `combined_length = f(signal)` (isotonic + Theil-Sen)
against the ~70-sample vamos panel; generalizes across *samples*. (3)
**Calibration transfer** (`calibration_transfer.py`) — a gradient-boosted model
over locus covariates + signals with **no locus identity**, trained across the
60k bundled loci, for *novel* user-BED loci; carries an out-of-distribution score.
(4) **Global residual model** (`global_residual.py`) — stacked quantile GBM over
everything incl. the tier-2/3 predictions, emitting the point estimate, 90%
interval, and confidence tier. (5) **Per-allele** (`per_allele.py`) — spanning-pair
insert-size 2-component split, gated hard on confidence.

Evaluation (`evaluate.py`) reports per-quadrant metrics — the frozen test split
holds out whole **loci** and whole **samples** independently, so
`test_novel_locus` measures generalization to an arbitrary user BED — always
alongside the naive `density_ratio × 2 × ref_units` baseline.

---

## 5. Gotchas a fresh session must know

- **vamos GT is 1-based.** GT value `k` -> `ALTANNO[k-1]`. The OLD
  `src/1000Genomes_Analysis/1000Genomes_Analysis.py` has an off-by-one here that
  corrupts its ground truth CSVs — **do not reuse `vntr_correlations*.csv`,
  `cn_variable_genotypes_*.csv`, `variable_vntrs.csv`.** The new
  `parse_vamos.py` does it correctly. Full spec: `docs/context/03-ground-truth-spec.md`.
- **Bundled BED coordinate convention:** `start == vamos POS`, `end == vamos END`,
  with **zero offset** — `start` is effectively 1-based despite the file header.
  `parse_vamos.py` writes `loci_targets.parquet` with `array_start0 = start - 1`
  for 0-based read fetching. Some BED rows are degenerate (`start == end`);
  extraction returns `note="degenerate_interval"` for those.
- **Ground truth file:** `data/1000Genomes/vamosExpanded_v3.0_GRCh38_oriMotifs.vcf`
  (10 GB, 416 diploid samples, columns `APR-phase1__aprNNN_hap1/..._hap2`). Same
  file the old tool used. Not committed (too large); present locally.
- **Runner needs no `samtools` CLI** — it uses pysam/htslib HTTP-range reads
  directly. But pysam must be built against libcurl-enabled htslib (the
  `environment.yml` pins `htslib` from bioconda to ensure this).
- `work/` and the filled `config/samples.tsv` / `config/references.tsv` are
  git-ignored. `work/ground_truth/` and `work/covariates/` regenerate in ~2 min
  from local files; `work/features/` requires the real extraction run.

---

## 6. Data assets (local, not in git)

| path | what | size |
|---|---|---|
| `data/1000Genomes/vamosExpanded_v3.0_GRCh38_oriMotifs.vcf` | ground-truth call set | 10 GB |
| `data/1000Genomes/vamosExpanded_v3.0_GRCh38_oriMotifs.tsv` | per-locus annotation (motifs, context) | 584 MB |
| `data/1000Genomes/mutRate_GRCh38_masterFile.tsv` | per-locus covariates (segdup, GC, entropy, mutation rate, ...) | 1.3 GB |
| `data/1000Genomes/anova.{annoLen,ed2Major}.tsv` | population-structure stats per locus | ~200-250 MB |
| `data/VNTRAnnotations/v1_VNTR_Annotations_03_07_2026.bed` | the 60,519-locus bundled annotation | small, in git |
| `data/Samples/test_sample.bam` | toy BAM (chr5 DAT1, ~0.3x) for smoke tests | small, in git |

---

## 7. The April validation work (`src/1000Genomes_Analysis/`, untracked)

Prior exploratory validation of the OLD tool (Apr 2026): HPRC + HGSVC2 whole-genome
runs, per-VNTR correlation analysis, ~15 diagnostic figures. Its conclusions are
summarized in `docs/context/01-old-tool-diagnosis.md`. **Its ground-truth CSVs are
corrupted by the 1-based-GT bug** (see §5) so the absolute-scale numbers may be
worse than reality — but the qualitative findings (rank works for a small
minority, scale is uncalibrated per-locus, normalization window doesn't matter)
still hold. Kept for reference; the scripts are reusable, the derived CSVs are not.
