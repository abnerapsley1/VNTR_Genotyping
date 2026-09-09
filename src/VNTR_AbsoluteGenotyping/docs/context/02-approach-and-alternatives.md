# Approach and alternatives considered

## What's physically recoverable from short reads

Reads (100–150 bp) are shorter than most VNTR arrays, so this is **dosage
inference**, not measurement.

| Target | Feasibility |
|---|---|
| Combined array length (L1+L2) | Realistic for most loci with unique flanks + OK mappability. **The primary output.** |
| Per-allele (L1, L2) | Only a subset: heterozygous loci separable by motif composition, flank phasing SNPs, or a bimodal spanning-pair insert-size signal. **Opt-in, confidence-gated secondary.** |
| Exact motif-string alleles | Only arrays shorter than ~1 fragment. Small fraction of VNTRs. |

Recurring identifiability fact: **depth alone never splits two alleles** — it
gives the sum. Per-allele needs the spanning-pair insert mixture, motif/SNP read
phasing, or assembly bubble-ordering.

## The signal inventory (extracted by `pipeline/features/signals.py`)

The old tool used one signal (normalized depth ratio) and discarded the rest.
The new layer extracts, in one pass per locus:

1. Depth, MAPQ-stratified {0,1,20,60}, array vs multiple flank windows vs genome-wide
2. In-repeat-read (IRR) count — reads that are ~all motif k-mers
3. Motif-composition profile — frequency vector of motif variants among locus reads
4. Spanning read-pair insert-size distribution — mean shift = length change;
   bimodality = het alleles (this is the main per-allele signal)
5. Anchored-pair counts (one mate in unique flank, one in repeat)
6. Soft-clip / junction pileups at both repeat boundaries
7. Locus covariates — GC, mappability, segdup, reference copy number, motif
   entropy, mutation rate (from `mutRate_GRCh38_masterFile.tsv` etc.)

## Method families

**A. Mechanistic forward model** (ExpansionHunter/GangSTR style). Write
`P(signals | L1, L2, nuisances)` and invert. Pro: no training data, extrapolates
to any BED. Con: the per-locus mapping-efficiency / paralog-contamination
parameters are exactly what sank the old tool.

**B. Supervised ML regression.** Featurize (locus×sample) → gradient-boosted
trees → long-read truth, quantile heads for CIs. Pro: learns the corrections.
Con: ~70 samples; generalization to *novel loci* is the risk unless the model
leans on transferable covariates rather than locus identity.

**C. Simulation-based calibration.** Synthesize alleles at known L, simulate reads
with the sample's empirical profile, run identical extraction → per-locus
`signal→L` curve. Pro: captures mapping saturation / paralog recruitment with no
ground truth; works on arbitrary BED. Con: simulator-realism risk; expensive.

**D. Graph-realignment front-end.** Realign locus + MAPQ-0 reads to a small graph
(flanks + motif cycle + optional paralog nodes) before extraction. Directly fixes
the mismapping root cause. Not an estimator — a preprocessing step. Deferred to
Phase 5 but the extraction layer is designed to accept it.

**E. Local assembly.** de-Bruijn over locus reads → motif truth + coverage-based
CN + partial phasing. A per-allele / motif helper, not a length estimator.

## Chosen approach: two-tier, panel-anchored (B-centric, with A/C/D as roadmap)

Rationale:

- With ~70 samples you get ~70 (length, signal) points **per bundled locus** —
  enough to fit a per-locus calibration for all 60k. This directly removes the
  4.4× multiplicative bias and the per-locus slope distortion. Nothing else gives
  that much locus-specific correction for that little modeling risk.
- The 60k bundled loci then become **60k training units** for a model that
  predicts the calibration from locus covariates — a well-posed problem even
  with few *samples* — which is what makes an arbitrary user BED work.
- Incremental: Tier-1 alone is a shippable, measurable win; each later tier is
  additive, not a rewrite.

Implemented as:

1. `per_locus_calibration.py` — isotonic + Theil-Sen per locus, panel-fit,
   generalizes across samples. Tier for "bundled locus in panel".
2. `calibration_transfer.py` — GBM over covariates + signals, **no locus id**,
   trained across bundled loci, + an out-of-distribution score. Tier for "novel
   user locus".
3. `global_residual.py` — stacked quantile GBM over everything (incl. tier-1/2
   predictions) → point estimate + 90% interval + `PASS/LOW_CONF/NA`.
4. `per_allele.py` — spanning-pair insert-size 2-component split, gated.

Not worth pursuing (from the old-tool diagnosis): further normalization-window
tuning; any global (non-per-locus) rescale.

## Evaluation design

Frozen split holds out whole **loci** and whole **samples** independently →
quadrants:

| quadrant | loci | samples | tests |
|---|---|---|---|
| `train` | train | train | fitting |
| `test_novel_sample` | train | test | new individual, known locus |
| `test_novel_locus` | test | train | **new locus — the arbitrary-user-BED case** |
| `test_novel_both` | test | test | both |

`evaluate.py` reports, per quadrant: pooled + median-per-locus Pearson/Spearman,
MAE (units), median relative error, calibration slope, |err|≤1/≤2 units, 90%
interval coverage, fraction of loci with per-locus r>0.5/0.7 — always next to the
naive `density_ratio × 2 × ref_units` baseline.

## Honest expectations

- Tier-1 on bundled CN-variable loci: per-locus Pearson likely 0.5–0.8, slope ≈ 1
  by construction, abs error ~1–3 motif units for typical loci, wider for
  long / low-mappability arrays.
- 40–70% of CN-variable loci get a useful combined length; the rest return wide
  CIs or `NA`.
- Per-allele fires on maybe 20–40% of het loci. Frame it as a minority feature.

## Open questions still worth deciding

- Real usable panel size (samples with BOTH good-depth short-read BAM AND
  assembly/vamos calls) — determines how far tier-2/3 ML can carry.
- Target sequencing profile (30x PCR-free 150bp vs 15x / 100bp / PCR+).
- CHM13 in scope? (currently GRCh38 only; vamos CHM13 files exist locally).
- Acceptable per-sample runtime for on-the-fly simulation on novel BEDs.
- Assume user loci are genuinely motif-structured VNTRs, or arbitrary CNV regions?
- Appetite for external deps (vg / GraphAligner / ExpansionHunter) as components.
