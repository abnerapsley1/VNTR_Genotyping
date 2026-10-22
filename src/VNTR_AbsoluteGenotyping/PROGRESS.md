# Build progress

Legend: [x] done & tested  ·  [~] built, needs real-data validation  ·  [ ] not started

## Foundation
- [x] Directory + README + config templates (`config/config.yaml`, `samples.tsv`, `references.tsv`)
- [x] `pipeline/paths.py` — repo-root + config loader
- [x] `pipeline/cli.py` — single entrypoint (`python -m pipeline.cli <cmd>`)
- [x] `environment.yml` (additive: pyarrow, scikit-learn, joblib, htslib, optional xgboost/lightgbm)

## Ground truth  — DONE, ran over full 10 GB VCF
- [x] `pipeline/ground_truth/parse_vamos.py`
      1-based GT `ALTANNO[gt-1]`; coord join offset (0,0); motif-unit + bp lengths
- [x] full run: 60,504/60,504 loci matched (exact), 25.2M-row `gt_long.parquet`,
      `gt_locus_summary.parquet`, `loci_targets.parquet` (0-based array coords + motifs)

## Covariates  — DONE
- [x] `pipeline/covariates/build_locus_covariates.py`
      gt summary + vamos TSV (exact) + mutRate masterFile (stream/overlap) + anova (overlap)
      -> `work/covariates/locus_covariates.parquet`  (60,504 loci x 78 cols)
      match rates: vamos TSV 100%, mutRate 86%, anova 72–99%

## Feature extraction  — DONE, 70/70 real samples (2026-10-07/08)
- [x] `pipeline/features/schema.py` — 96-col contract (depth MAPQ tiers x flank windows,
      IRR, spanning-pair insert mixture, soft-clip, motif composition, per-allele)
- [x] `pipeline/features/signals.py` — one-fetch/one-pass extractor incl. spanning-pair
      2-component per-allele estimate; degenerate-interval + no-reads guards
- [x] `pipeline/features/nuisance.py` — genome-wide depth / insert dist / read len / GC curve
- [x] ran on all 70 HGSVC2 + HPRC 30x CRAMs: 0 failures, ~60.4k/60.5k loci ok per sample

## Runner  — DONE
- [x] `pipeline/run/manifest.py` — read samples.tsv / references.tsv
- [x] `pipeline/run/remote.py` — pysam/htslib remote range reads; reference + index caching
- [x] `pipeline/run/process_sample.py` — 1 sample: cache index -> nuisance -> 32-worker
      per-locus extraction -> `features/<sid>.parquet` (+ nuisance.json); deletes index
- [x] `pipeline/run/run_all.py` — resumable manifest loop, one sample at a time, `_run_log.tsv`
- [x] remote CRAM slicing works (pysam/libcurl), but moves more bytes than the file for
      the dense bundled BED -> added `download` mode (parallel ranges, prefetch, delete)

## Model  — DONE on real data (70 samples); see HANDOFF §2 for results
- [x] `pipeline/model/splits.py` — frozen locus x sample folds -> 4 test quadrants
- [x] `pipeline/model/assemble.py` — features + gt + covariates + folds -> `model_table.parquet`
- [x] `pipeline/model/per_locus_calibration.py` — Phase 1 isotonic+TheilSen per locus
- [x] `pipeline/model/calibration_transfer.py` — Phase 2 GBM (no locus id) + OOD score
- [x] `pipeline/model/global_residual.py` — Phase 3 known/novel residual GBM + conformal interval + PASS/LOW_CONF/NA tier
- [x] `pipeline/model/evaluate.py` — per-quadrant metrics vs naive baseline + per-locus r
- [x] `pipeline/model/per_allele.py` — homozygous-default per-allele estimate around the combined call
- [x] `pipeline/model/_gbm.py` — sklearn HGB default; `VNTR_GBM=xgboost|lightgbm` to swap
- [x] synthetic run: all 7 steps pass; naive baseline reproduces slope~0.2 / perLocus_r~0.12,
      model recovers slope~1 (mechanics only — synthetic signal is near-noise-free)

- [x] leakage fixes: train-fold-only per-locus calibration, novel-safe covariate
      allowlist, 5-fold sample cross-fitting of stacked inputs
- [x] global model: known/novel paths, residual-vs-anchor squared-error GBM,
      whole-unit conformal 90% intervals (coverage 0.91 on test)
- [x] per-allele: homozygous-default estimator (both alleles ±2u: 0.05 -> ~0.9)
- [x] evaluate: within-locus r, mean abs err, panel-mean baseline, r=0 for constant preds

## Remaining
- [ ] commit the 2026-10-07/08 changes
- [ ] tune PASS tier threshold on validation rows
- [ ] fix pa_units sign in signals.py at the next full re-extraction
- [ ] decide motif contract for truly novel user loci
- [ ] Phase-5: GC-depth correction, mappability mask, graph realignment,
      simulation calibration, benchmark vs adVNTR / ExpansionHunter
