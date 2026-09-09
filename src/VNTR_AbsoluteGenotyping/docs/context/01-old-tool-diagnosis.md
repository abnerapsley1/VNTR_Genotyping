# Old density-ratio tool — validation diagnosis

Source: `src/1000Genomes_Analysis/` April-2026 work + a re-analysis on 2026-08-31.
Compares the old tool's `predicted_copies` (`= density_ratio × 2 × vntr_len/period`)
to long-read ground-truth combined copy number (HPRC + HGSVC2, ≤70 sample pairs
per VNTR, pairs kept only where both haplotypes are called).

> Caveat: the April ground-truth CSVs are corrupted by a 1-based-GT off-by-one
> (see `03-ground-truth-spec.md`). Absolute-scale magnitudes may be somewhat
> worse than reality; the qualitative conclusions below still hold.

## Rank / relative performance — works for a small minority

- Genome-wide median per-VNTR Pearson r ≈ **−0.017** (essentially zero for the
  typical VNTR).
- Only ~391 / 60,500 VNTRs (**0.6%**) reach r > 0.5; ~29 reach r > 0.8.
- The r > 0.5 subset is enriched for: high expected heterozygosity, high CN
  variance, long total repeat length, segDup context (28% vs 8% background),
  chromosomes 15 and 21.

## Absolute `predicted_copies` — not usable as-is

- Overestimates true CN for **97.3%** of VNTRs; median ratio tool/truth =
  **4.4×** (IQR 2.2–9.2×). Median tool output ≈ 57 copies vs median true ≈ 9.5.
- Root cause is the **density ratio running ~6× hot** (inflated numerator), not
  the `2 × vntr_len/period` term (that alone would give ≈ 9). Drivers:
  - no MAPQ filter in `get_read_names` — paralog / MAPQ-0 reads counted;
  - short reference repeat span — reads from long alleles pile onto it;
  - gene-background depth biased low by introns/UTRs.
- **Dynamic-range distortion:** OLS slope (truth→tool) among well-correlating
  high-variance VNTRs centers ~0.53 (≈2× compression = short-read mapping
  saturation), heavy-tailed, uncalibrated per locus.
- **Correlation and absolute scale are decoupled** — e.g. chr7:78715845 has
  r = 0.94 but slope = 43, MAE = 48 copies.
- segDup is **not** the dominant absolute-bias driver — the overestimation is
  near-universal.

## Normalization method is irrelevant

Gene vs 5 kb vs 10 kb local-window backgrounds give **identical** r and MAE. The
bottleneck is VNTR read recruitment, not the background choice. Do not spend more
time tuning the normalization window.

## Implications for the new tool

1. Per-locus calibration is mandatory — a single global rescale cannot fix a
   per-locus slope that ranges from <0.5 to >4.
2. The signal layer must go beyond the depth ratio (add IRR counts, spanning-pair
   insert shifts, motif composition).
3. Confidence gating is core — the old tool's worst property is emitting
   confident nonsense. `NA` must beat a confident wrong number.
4. Only a minority of VNTRs will ever yield a useful absolute length from short
   reads (realistically 40–70% of CN-variable loci). Say so honestly.
