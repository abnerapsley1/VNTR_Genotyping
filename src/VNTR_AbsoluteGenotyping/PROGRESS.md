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

## Feature extraction  — built; validated on local toy BAM
- [x] `pipeline/features/schema.py` — 96-col contract (depth MAPQ tiers x flank windows,
      IRR, spanning-pair insert mixture, soft-clip, motif composition, per-allele)
- [x] `pipeline/features/signals.py` — one-fetch/one-pass extractor incl. spanning-pair
      2-component per-allele estimate; degenerate-interval + no-reads guards
- [x] `pipeline/features/nuisance.py` — genome-wide depth / insert dist / read len / GC curve
- [~] not yet run on a real 30x BAM (needs manifest + network)

## Runner  — built; not exercised against a live remote CRAM
- [x] `pipeline/run/manifest.py` — read samples.tsv / references.tsv
- [x] `pipeline/run/remote.py` — pysam/htslib remote range reads; reference + index caching
- [x] `pipeline/run/process_sample.py` — 1 sample: cache index -> nuisance -> 32-worker
      per-locus extraction -> `features/<sid>.parquet` (+ nuisance.json); deletes index
- [x] `pipeline/run/run_all.py` — resumable manifest loop, one sample at a time, `_run_log.tsv`
- [~] needs a real manifest row to confirm htslib remote CRAM slicing on the target host

## Model  — built; full chain validated end-to-end on synthetic data
- [x] `pipeline/model/splits.py` — frozen locus x sample folds -> 4 test quadrants
- [x] `pipeline/model/assemble.py` — features + gt + covariates + folds -> `model_table.parquet`
- [x] `pipeline/model/per_locus_calibration.py` — Phase 1 isotonic+TheilSen per locus
- [x] `pipeline/model/calibration_transfer.py` — Phase 2 GBM (no locus id) + OOD score
- [x] `pipeline/model/global_residual.py` — Phase 3 stacked quantile GBM + PASS/LOW_CONF/NA tier
- [x] `pipeline/model/evaluate.py` — per-quadrant metrics vs naive baseline + per-locus r
- [x] `pipeline/model/per_allele.py` — evaluate/gate the spanning-pair per-allele estimate
- [x] `pipeline/model/_gbm.py` — sklearn HGB default; `VNTR_GBM=xgboost|lightgbm` to swap
- [x] synthetic run: all 7 steps pass; naive baseline reproduces slope~0.2 / perLocus_r~0.12,
      model recovers slope~1 (mechanics only — synthetic signal is near-noise-free)

## Remaining before real results
- [ ] user fills `config/samples.tsv` (70 rows) + `config/references.tsv`
- [ ] `python -m pipeline.cli list-vamos` to get exact vamos column headers
- [ ] run `pipeline.cli extract` on the NVIDIA box (needs htslib+libcurl; `mamba env update`)
- [ ] `pipeline.cli all-model` -> read `work/model/evaluation.json`
- [ ] tune: signal shortlist for per-locus calibration, PASS-tier interval-width threshold,
      per-allele confidence threshold — all against real held-out metrics
- [ ] consider Phase-5 items (graph realignment front-end, simulation augmentation)
