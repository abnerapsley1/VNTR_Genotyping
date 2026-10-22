# Improvement ideas / roadmap

Original ranked list from the 2026-08-31 diagnosis, annotated with what the new
pipeline already implements.

1. **Per-VNTR empirical calibration.** ✅ implemented as
   `pipeline/model/per_locus_calibration.py` (isotonic + Theil-Sen per locus,
   fit on the panel, generalizes across samples) and
   `pipeline/model/calibration_transfer.py` (covariate model for novel loci).
   Emits `NA` / `LOW_CONF` when the fit is weak or the locus is out-of-distribution.

2. **MAPQ handling / paralog subtraction.** ~partially: `signals.py` records
   depth at MAPQ tiers {0,1,20,60} and `mapq0_frac` as separate features so the
   model can use them. A principled paralog-expectation subtraction is not done —
   see Phase-5 graph realignment below.

3. **Explicit reference copy number.** ~partially: covariates include a
   `ref_units_est` from `array_span_bp / period`; a directly-measured `ref_copies`
   from the vamos annotation would be better. `cov_mr_size_bp` is also available.

4. **Model mapping saturation.** ✅ in spirit: the per-locus calibration is
   isotonic (monotone, non-linear), so it fits a saturating `signal→length`
   curve rather than assuming linearity. Simulation-derived curves (below) would
   extend this into CN ranges the panel doesn't cover.

5. **GC-bias correction.** ~partially: `nuisance.py` estimates a per-sample
   GC-depth curve; it is stored but not yet applied as an explicit correction to
   the depth signals. TODO.

6. **Tighten background to unique-mappable bases.** not done — the flank-depth
   windows in `signals.py` are raw genomic windows. A mappability track mask would
   improve them.

7. **Flank-anchored read-pair counting.** ✅ `signals.py` emits `n_anchored_pairs`
   and the spanning read-pair insert-size mixture.

## Phase 5 (not built) — highest-value next additions

- **Graph-realignment front-end (method D).** Realign locus-neighborhood reads +
  MAPQ-0 reads to a small graph (unique flanks + motif cycle + optional paralog
  nodes) before signal extraction. Directly attacks the mismapping that caused
  the old tool's 4.4× overestimation. Design `signals.py` to accept a
  pre-realigned read stream. Candidate tools: `vg giraffe`, `GraphAligner`, or
  ExpansionHunter's graph aligner.

- **Simulation-based per-locus calibration (method C).** For loci with thin or no
  panel data (novel user BED), synthesize alleles at a grid of known lengths
  using the locus's motif set, simulate reads with the sample's empirical
  profile, run identical extraction, fit `signal→length`. Cross-validate against
  the panel where they overlap. Precompute for bundled loci; run on-the-fly
  (with a `--fast` covariate-model fallback) for novel BEDs.

- **Apply the GC-depth curve** as an explicit correction (item 5).

- **Mappability-masked flank windows** (item 6) — needs a precomputed
  GRCh38 mappability track added to `data/`.

## Benchmarking

Once real evaluation exists, compare on the frozen test set against:
- the naive `density_ratio × 2 × ref_units` baseline (already in `evaluate.py`),
- adVNTR-NN,
- ExpansionHunter (VNTR mode).
Either as competitors to beat or components to borrow.

## Deferred: mate rescue for in-repeat reads (noted 2026-10-22)

In-repeat read signals (`irr_count_norm`, `partial_repeat_frac`, `anchored_irr_count_norm`)
only count reads the aligner placed on the repeat's reference coordinates. Pure-repeat reads
that were placed at a similar repeat elsewhere (MAPQ 0) or left unmapped are missed, and
misplaced reads from similar repeats can be counted. Per-repeat calibration absorbs a steady
misplacement rate, but not one that varies between people (segmental duplications, common motifs).

Proposed: mate rescue, as in ExpansionHunter/GangSTR — when one read of a pair is firmly placed
in this repeat's unique flank and its partner is pure repeat, count the partner for this repeat
wherever it was aligned (including unmapped reads). This changes features, so test on a subset of
repeats (extract, refit, compare test error) before a full re-extraction and retrain.
