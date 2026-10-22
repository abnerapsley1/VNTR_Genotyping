# HANDOFF — read this first

Written 2026-09-09; **updated 2026-10-08 after the first real 70-sample run.**
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

## 2. Current state (2026-10-08)

**The pipeline has run end-to-end on real data.** Ground truth, covariates,
feature extraction for all 70 short-read samples (34 HGSVC2 + 36 HPRC, 0
failures), and the full model chain.  Code changes from that run are
**uncommitted** on branch `absolute-genotyping-pipeline` (see §3).

### Headline results (2026-10-08 run: 92% of loci train / 3% val / 5% test; 48 train / 7 val / 15 test people; panel stats exclude the 22 held-out people)

| question | metric | **model** | baseline |
|---|---|---|---|
| known locus, new individual (`test_novel_sample`, 55k loci) | mean abs err | **1.02 u** | panel mean 1.80 u |
| | within-locus r | **0.79** | panel mean 0 |
| | 90% interval coverage | **0.91** (0.91–0.92 in every core-spread bin) | |
| novel locus (`test_novel_locus`, 3k loci) | mean abs err | **2.46 u** | naive depth ratio, much worse |
| | within-locus r | **0.69** | ~0 |
| | 90% interval coverage | **0.91** overall, but 0.70–0.80 at high-spread loci | |
| rare-allele carriers (person far outside the panel core) | right direction | **70%** (known) / 79% (novel) | |
| | 90% coverage | **0.17** (known) — intervals miss extreme alleles | |
| per-allele, both alleles within ±2 u | | 0.88–0.92 | (old split: 0.05) |

Locus variability is described by **core spread** (p95 − p5 of the panel's
combined length; a few extreme people can't move it) plus an **outlier count**
(people > 2 × max(core, 1) beyond the core) — not SD, which outliers dominate.
`evaluate.py` reports performance by core-spread bin and by site type
(stable / broad spread / stable core + outliers / broad spread + outliers).

`u` = motif units.  **Read the metrics correctly:** most VNTRs barely vary
between people, so pooled correlation (~0.99) and median error (~0) look
excellent for *any* method incl. "predict the panel mean". The informative
metrics are **within-locus r** (pooled correlation of locus-demeaned prediction
vs truth — does it rank individuals?) and **mean abs error vs the panel-mean
baseline**.  `evaluate.py` reports both.

### What changed in the 2026-10-07/08 session (and why)

1. **Download mode** (`runner.transfer_mode: download`, now default).  For the
   dense 60k-locus BED, per-locus remote slicing re-reads overlapping CRAM
   slices: ~300 KB/locus ≈ 18 GB/sample, more than the 15–30 GB CRAM.  The
   runner now pulls the CRAM with 16 parallel HTTP ranges (`remote.fetch_parallel`;
   single EBI connections throttle to ~1 MB/s), extracts locally, deletes it, and
   prefetches the next sample.  ~10 min/sample, peak scratch ~60 GB.
2. **Leakage fixes** — the first real run's novel-locus score (r 0.76) was
   inflated: per-locus calibration was fit on held-out loci too, and the global
   model saw panel-derived covariates.  Now: per-locus calibration only for
   train-fold loci; novel path uses `calibration_transfer.NOVEL_SAFE_COVARIATES`
   (BED/reference-derivable only); stacked inputs (`plc_pred`, `transfer_pred`)
   are **5-fold cross-fitted over samples** on train rows (the stacker had
   learned to copy in-sample calibration curves).
3. **Global model redesign** (`global_residual.py`): two paths — *known* (bundled
   loci, all features, anchor = per-locus calibration or panel mean) and *novel*
   (novel-safe features, anchor = 2 × reference units).  Each predicts the
   residual `log1p(y) − log1p(anchor)` with a **squared-error** GBM (the old
   quantile-median point estimate was the main reason the shipped model lost to
   its own components).  Intervals: 5%/95% quantile GBMs, reported in **whole
   motif units** and conformal-calibrated on the path's validation rows
   (zero-width float intervals had broken coverage).
4. **Per-allele** (`per_allele.py`): the extraction-time split had a sign error
   (TLEN is in reference coordinates → shift ≈ −Δallele) and always split
   homozygous loci.  More fundamentally ~85% of genotypes are homozygous and the
   median het allele difference is ~22 bp vs ~100 bp insert-size SD, so
   **short reads cannot resolve alleles at most loci**.  The estimator now calls
   homozygous (combined ÷ 2) unless the spanning-pair gap clearly exceeds the
   single-allele expectation; on this data the tuned threshold never calls het.
   Treat per-allele output as a homozygous-assumption estimate.
5. **Per-locus signal shortlist** extended with `spanning_tlen_shift`,
   `depth_ratio_d0_w2000`, `anchored_irr_count_norm` (out-of-fold within-locus
   r 0.685 → 0.714).
6. **Evaluation**: constant predictions score r = 0 (were silently dropped);
   added within-locus r, mean abs err, panel-mean baseline; NaN-safe; core-spread
   / site-type breakdown and a rare-allele-carrier check.
7. **Locus holdout cut to 5% test + 3% val** (was 20% + 10%): users query the
   bundled loci, and a held-out locus loses its per-locus calibration.
8. **Panel statistics exclude held-out people** (`model/panel_stats.py`,
   applied in `assemble`): panel min/max otherwise contained a held-out rare-
   allele carrier's own answer.  Adds `panel_core90_units`, `panel_n_outliers`.
9. **Per-group interval calibration**: known path by core-spread bin, novel path
   by predicted-width quintile.

---

## 3. Exact next steps

```bash
cd src/VNTR_AbsoluteGenotyping && conda activate vntr_analysis
python -m pipeline.cli all-model        # ~45 min on 32 cores; reruns everything from features
cat work/model/evaluation.json
```

1. **Review + commit** the uncommitted changes above.
1b. **Ship refit**: once evaluation is settled, refit on ALL loci and all 70
   people (keep a calibration set for the intervals) and report the held-out
   numbers above as the model's performance.
1c. **Rare-allele carriers**: direction is right ~70% of the time but size is
   underestimated and intervals cover only 17%.  Consider a "possible rare
   expansion/contraction" flag when y_pred leaves the panel core, and check
   whether the extreme panel calls (e.g. chr20:45320191, one person at ~1,140
   copies vs median 8) are real or long-read artifacts.
1d. **Novel-path interval groups**: width quintiles don't track true locus
   variability (coverage 0.70–0.80 at high-spread novel loci); find a
   reference-derivable grouping (e.g. reference length, motif length).
2. **Tier calibration** — `PASS` is ~91% (known) / ~83% (novel) of calls, set
   by `model.pass_max_interval_width_rel` (1.0) + OOD p90.  PASS rows have ~30%
   lower error, but the threshold was not tuned; pick it on validation rows for
   a target error rate.
3. **`signals.py` per-allele fields** still carry the old sign in
   `pa_units_a1/a2` (kept unchanged so all 70 samples stay consistent; GBMs are
   sign-agnostic).  Fix at the next full re-extraction.
4. **Novel-locus motif dependence** — in-repeat-read signals use vamos motif sets
   (population-derived).  A truly novel user locus needs user-supplied motifs;
   decide the BED contract.
5. **Phase 5** (`docs/context/02-approach-and-alternatives.md`,
   `04-improvement-ideas-and-roadmap.md`): apply the GC-depth curve (computed,
   unused), mappability-masked flanks, graph realignment, simulation-based
   calibration for novel loci, benchmarking vs adVNTR / ExpansionHunter.

Knobs: `per_locus_calibration._SIGNAL_SHORTLIST`; `config.yaml` →
`model.stack_cv_folds`, `model.n_jobs`, `model.pass_max_interval_width_rel`,
`runner.download_streams`; `export VNTR_GBM=xgboost` for GPU GBMs.

---

## 4. Architecture in one paragraph

(1) **Signal extraction** — one BAM pass per locus: MAPQ-tiered depth ×
flank windows, in-repeat reads, spanning-pair insert sizes, soft clips, motif
composition, plus per-sample nuisance stats.  (2) **Per-locus calibration**
(`per_locus_calibration.py`) — for each train-fold bundled locus, isotonic +
Theil-Sen `combined_length = f(best shortlist signal)`; cross-fitted on train
rows.  (3) **Calibration transfer** (`calibration_transfer.py`) — GBM over
signals + novel-safe covariates, no locus identity; cross-fitted; OOD score.
(4) **Global model** (`global_residual.py`) — known/novel paths, residual-vs-
anchor squared-error GBM + quantile GBMs, whole-unit conformal 90% interval,
`PASS/LOW_CONF/NA` tier.  (5) **Per-allele** (`per_allele.py`) — homozygous
default around the combined estimate.

Evaluation (`evaluate.py`) holds out whole **loci** and whole **samples**
independently: `test_novel_sample` = known locus/new person, `test_novel_locus`
= arbitrary user BED.  Always alongside naive and panel-mean baselines.

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
  from local files; `work/features/` (1.4 GB, 70 samples) took ~12 h to extract.
- **Manifest provenance:** `config/samples.tsv` was built from
  `data/1000Genomes/{hgsvc2,hprc}_sr_download_globus.txt` (ENA Globus paths
  `/vol1/run/...` → `https://ftp.sra.ebi.ac.uk/vol1/run/...`); sample ID parsed
  from the CRAM filename and matched to exactly one vamos column for all 70.
  Several HGSVC2 short-read samples have their long-read calls under
  `HGSVC-phase3__…_h1`, `HGSVC-phase2__…_CCS_h1` or `HPRC-phase2__…_mat` columns —
  same individuals (CRAM `@RG SM` verified).  All CRAMs are aligned to
  `GRCh38_full_analysis_set_plus_decoy_hla.fa` (cached under
  `work/scratch/_refcache/`).
- `gt_long.vamos_sample` stores the LEFT half of the vamos column with only a
  `_hap1` suffix stripped (e.g. `HPRC-phase2__HG00733_mat`); join via
  `assemble._clean_vamos_col`, not by bare sample ID.
- Comparisons ACROSS quadrants are unreliable: per-locus statistics depend on
  how many samples a quadrant has per locus (48 vs 15).  Compare methods within
  a quadrant.

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

## 7. The April validation work (`src/1000Genomes_Analysis/`)

Prior exploratory validation of the OLD tool (Apr 2026): HPRC + HGSVC2 whole-genome
runs, per-VNTR correlation analysis, ~15 diagnostic figures. Its conclusions are
summarized in `docs/context/01-old-tool-diagnosis.md`. **Its ground-truth CSVs are
corrupted by the 1-based-GT bug** (see §5) so the absolute-scale numbers may be
worse than reality — but the qualitative findings (rank works for a small
minority, scale is uncalibrated per-locus, normalization window doesn't matter)
still hold. Kept for reference; the scripts are reusable, the derived CSVs are not.
