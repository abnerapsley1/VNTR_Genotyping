"""
schema.py — the per-(locus, sample) feature contract.

`FEATURE_COLUMNS` is the ordered, canonical column list written to
features/<sample>.parquet.  Signal extraction (`signals.py`) returns a dict;
`row_from_signals()` projects it onto this schema (missing -> None) so every
sample's parquet has identical columns regardless of extraction hiccups.

Naming
------
  d{T}          MAPQ tier, e.g. d0, d1, d20, d60  (reads with mapq >= T)
  w{W}          flank window in bp, e.g. w500, w1000, w2000, w5000
  *_norm        divided by the reference flank depth (w2000, d20) for that sample
"""

from __future__ import annotations

MAPQ_TIERS = (0, 1, 20, 60)
FLANK_WINDOWS = (500, 1000, 2000, 5000)
REF_FLANK_W = 2000            # window used for *_norm normalisation
REF_FLANK_MAPQ = 20

_ID_COLS = [
    "sample_id", "locus_id", "chrom", "array_start0", "array_end0", "period",
    "n_reads_fetched", "extraction_ok", "note",
]

_DEPTH_COLS = (
    [f"array_depth_d{t}" for t in MAPQ_TIERS]
    + [f"array_reads_d{t}" for t in MAPQ_TIERS]
    + [f"flank_depth_d{t}_w{w}" for t in MAPQ_TIERS for w in FLANK_WINDOWS]
    + [f"depth_ratio_d{t}_w{w}" for t in MAPQ_TIERS for w in FLANK_WINDOWS]
    + [
        "mapq0_frac", "array_depth_cv", "array_depth_end_mid_ratio",
        "array_gc", "flank_gc_w2000",
        "array_depth_d20_norm",          # array_depth_d20 / ref flank depth
    ]
)

_IRR_COLS = [
    "irr_count", "irr_frac", "irr_count_norm",
    "partial_repeat_count", "partial_repeat_frac",
    "anchored_irr_count", "anchored_irr_count_norm",
    "mate_unmapped_count", "mate_unmapped_count_norm",
]

_SPAN_COLS = [
    "n_spanning", "n_spanning_norm",
    "spanning_tlen_median", "spanning_tlen_mad",
    "spanning_tlen_shift", "spanning_tlen_p10", "spanning_tlen_p90",
    "spanning_tlen_bimodality", "spanning_len_gap",
    "n_anchored_pairs", "n_anchored_pairs_norm",
]

_CLIP_COLS = [
    "clip_left_count", "clip_right_count", "clip_count_norm",
    "clip_left_bp_mean", "clip_right_bp_mean",
    "discordant_count", "discordant_count_norm",
]

_MOTIF_COLS = [
    "mc_n_reads_used", "mc_n_motifs_observed", "mc_entropy",
    "mc_frac_top_motif", "mc_frac_ref_consensus", "mc_divergence_from_ref",
    "mc_mean_unit_mismatch",
]

# opt-in per-allele estimate from the spanning read-pair insert-size mixture
# (only informative when arrays are shorter than ~1 fragment and the two
# alleles differ).  Downstream gating decides whether to surface it.
_PERALLELE_COLS = [
    "pa_spanning_n", "pa_bp_shift_a1", "pa_bp_shift_a2",
    "pa_units_a1", "pa_units_a2", "pa_mix_separation", "pa_confidence",
]

FEATURE_COLUMNS: list[str] = (
    _ID_COLS + _DEPTH_COLS + _IRR_COLS + _SPAN_COLS + _CLIP_COLS + _MOTIF_COLS
    + _PERALLELE_COLS
)

# columns that are categorical / string (everything else numeric)
STRING_COLUMNS = {"sample_id", "locus_id", "chrom", "note"}


def row_from_signals(sig: dict) -> dict:
    """Project a signals dict onto the fixed schema (missing keys -> None)."""
    return {c: sig.get(c) for c in FEATURE_COLUMNS}


def empty_row(sample_id: str, locus_id: str, note: str) -> dict:
    r = {c: None for c in FEATURE_COLUMNS}
    r.update(sample_id=sample_id, locus_id=locus_id,
             extraction_ok=False, note=note, n_reads_fetched=0)
    return r
