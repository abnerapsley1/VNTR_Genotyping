"""Inventory of the inputs to (1) per-locus calibration, (2) the
calibration-transfer model and (3) the extra inputs of the final known-repeat
model, with descriptions, where each value comes from, and how much each was used.

data_origin:
  short reads                     measured from the person's BAM/CRAM
  reference / repeat BED          GRCh38 sequence or the repeat coordinates
  long-read annotation            fixed per-repeat descriptors from long-read resources
                                  (vamos motif catalog; mutation-rate table; ANOVA tables)
  long-read panel genotypes       statistics of the vamos panel's repeat lengths,
                                  recomputed with validation/test people removed
  model estimate                  output of an earlier model stage (short-read input,
                                  fitted against training people's long-read truth)

    python scripts/feature_inventory.py   ->  work/model/feature_inventory.csv
"""
import collections
import re
import sys
from pathlib import Path

import joblib
import polars as pl

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from pipeline.model.per_locus_calibration import _SIGNAL_SHORTLIST  # noqa: E402

M = HERE / "work" / "model"
MQ = {"0": "all mapped reads (MAPQ >= 0)", "1": "MAPQ >= 1", "20": "MAPQ >= 20", "60": "MAPQ >= 60"}
NORM = " Divided by the locus's reference flank depth (MAPQ >= 20, 2 kb flanks) to remove sequencing-depth differences."

FIXED = {
    "period": ("locus", "annotation BED", "Motif period of the locus from the bundled annotation", "bp"),
    "mapq0_frac": ("read depth", "per locus x person", "Fraction of reads overlapping the array with mapping quality 0 (ambiguous placement)", "fraction"),
    "array_depth_cv": ("read depth", "per locus x person", "Coefficient of variation of MAPQ >= 20 depth along the array", "ratio"),
    "array_depth_end_mid_ratio": ("read depth", "per locus x person", "Mean depth over the array's outer fifths divided by depth over its middle (a dip in the middle suggests reads cannot span a long array)", "ratio"),
    "array_gc": ("GC content", "reference", "GC fraction of the reference array sequence", "fraction"),
    "flank_gc_w2000": ("GC content", "reference", "GC fraction of the 2 kb upstream reference flank", "fraction"),
    "array_depth_d20_norm": ("read depth", "per locus x person", "Mean MAPQ >= 20 depth over the array, normalised by the reference flank depth (identical to depth_ratio_d20_w2000)", "ratio"),
    "irr_count": ("in-repeat reads", "per locus x person", "Number of reads overlapping the array whose 5-mers are >= 80% from the locus's motif set (reads lying wholly in the repeat)", "reads"),
    "irr_frac": ("in-repeat reads", "per locus x person", "In-repeat reads as a fraction of all reads overlapping the array", "fraction"),
    "irr_count_norm": ("in-repeat reads", "per locus x person", "In-repeat read count." + NORM, "reads / depth"),
    "partial_repeat_count": ("in-repeat reads", "per locus x person", "Reads overlapping the array whose motif 5-mer fraction is 30-80% (reads straddling a repeat edge)", "reads"),
    "partial_repeat_frac": ("in-repeat reads", "per locus x person", "Partial-repeat reads as a fraction of reads overlapping the array", "fraction"),
    "anchored_irr_count": ("in-repeat reads", "per locus x person", "In-repeat reads whose mate is unmapped or on another chromosome", "reads"),
    "anchored_irr_count_norm": ("in-repeat reads", "per locus x person", "Anchored in-repeat read count." + NORM, "reads / depth"),
    "mate_unmapped_count": ("read pairs", "per locus x person", "Reads overlapping the array whose mate is unmapped", "reads"),
    "mate_unmapped_count_norm": ("read pairs", "per locus x person", "Mate-unmapped read count." + NORM, "reads / depth"),
    "n_spanning": ("read pairs", "per locus x person", "Read pairs whose outer ends extend >= 30 bp beyond both array edges (fragment spans the whole array)", "pairs"),
    "n_spanning_norm": ("read pairs", "per locus x person", "Spanning pair count." + NORM, "pairs / depth"),
    "spanning_tlen_median": ("read pairs", "per locus x person", "Median template length of spanning pairs (reference coordinates)", "bp"),
    "spanning_tlen_mad": ("read pairs", "per locus x person", "Median absolute deviation of spanning-pair template length", "bp"),
    "spanning_tlen_shift": ("read pairs", "per locus x person", "Median spanning template length minus the person's genome-wide median insert; an allele D bp longer than the reference shortens it by ~D", "bp"),
    "spanning_tlen_p10": ("read pairs", "per locus x person", "10th percentile of spanning-pair template length", "bp"),
    "spanning_tlen_p90": ("read pairs", "per locus x person", "90th percentile of spanning-pair template length", "bp"),
    "spanning_tlen_bimodality": ("read pairs", "per locus x person", "Sarle's bimodality coefficient of spanning template lengths (two distinct alleles push it up)", "coefficient"),
    "spanning_len_gap": ("read pairs", "per locus x person", "Largest gap between sorted spanning template lengths within the 10th-90th percentile range", "bp"),
    "n_anchored_pairs": ("read pairs", "per locus x person", "Pairs with exactly one read overlapping the array (one end anchored in flank)", "pairs"),
    "n_anchored_pairs_norm": ("read pairs", "per locus x person", "Anchored pair count." + NORM, "pairs / depth"),
    "clip_left_count": ("soft clipping", "per locus x person", "Reads soft-clipped by >= 10 bp at the array's left boundary", "reads"),
    "clip_right_count": ("soft clipping", "per locus x person", "Reads soft-clipped by >= 10 bp at the array's right boundary", "reads"),
    "clip_count_norm": ("soft clipping", "per locus x person", "Left + right clipped reads." + NORM, "reads / depth"),
    "clip_left_bp_mean": ("soft clipping", "per locus x person", "Mean clipped length of left-boundary clipped reads", "bp"),
    "clip_right_bp_mean": ("soft clipping", "per locus x person", "Mean clipped length of right-boundary clipped reads", "bp"),
    "discordant_count": ("read pairs", "per locus x person", "Pairs overlapping the array with mates on another chromosome or with templates > 10x the median insert", "pairs"),
    "discordant_count_norm": ("read pairs", "per locus x person", "Discordant pair count." + NORM, "pairs / depth"),
    "mc_n_reads_used": ("motif composition", "per locus x person", "Repeat-rich reads (motif 5-mer fraction >= 30%) tiled into motif-length windows", "reads"),
    "mc_n_motifs_observed": ("motif composition", "per locus x person", "Number of distinct vamos motifs assigned to read windows", "count"),
    "mc_entropy": ("motif composition", "per locus x person", "Shannon entropy of motif usage across read windows", "bits"),
    "mc_frac_top_motif": ("motif composition", "per locus x person", "Fraction of read windows assigned to the most common motif", "fraction"),
    "mc_frac_ref_consensus": ("motif composition", "per locus x person", "Fraction of read windows assigned to the consensus motif", "fraction"),
    "mc_divergence_from_ref": ("motif composition", "per locus x person", "Total-variation distance between observed motif usage and a uniform distribution over the locus's motifs", "0-1"),
    "mc_mean_unit_mismatch": ("motif composition", "per locus x person", "Mean per-base mismatch rate of read windows against their assigned motif", "fraction"),
    "pa_spanning_n": ("per-allele split", "per locus x person", "Number of spanning pairs available for the two-allele split", "pairs"),
    "pa_bp_shift_a1": ("per-allele split", "per locus x person", "Mean (template length - median insert) of the lower group of a 2-means split of spanning pairs", "bp"),
    "pa_bp_shift_a2": ("per-allele split", "per locus x person", "Mean (template length - median insert) of the upper group of the 2-means split", "bp"),
    "pa_units_a1": ("per-allele split", "per locus x person", "Reference span plus pa_bp_shift_a1, divided by period (extraction-time allele estimate; sign convention inverted, see Methods)", "repeat units"),
    "pa_units_a2": ("per-allele split", "per locus x person", "Reference span plus pa_bp_shift_a2, divided by period (as above)", "repeat units"),
    "pa_mix_separation": ("per-allele split", "per locus x person", "Gap between the two split means divided by the pooled within-group SD", "SD units"),
    "pa_confidence": ("per-allele split", "per locus x person", "Heuristic 0-1 confidence from pair count, separation and group balance", "0-1"),
    "nuis_depth_mean": ("per-person sequencing", "per person", "Mean MAPQ >= 20 depth over 3,000 random 2 kb autosomal windows", "x"),
    "nuis_depth_median": ("per-person sequencing", "per person", "Median MAPQ >= 20 depth over the random windows", "x"),
    "nuis_depth_sd": ("per-person sequencing", "per person", "Standard deviation of window depth", "x"),
    "nuis_insert_median": ("per-person sequencing", "per person", "Median insert size of proper pairs", "bp"),
    "nuis_insert_mad": ("per-person sequencing", "per person", "Median absolute deviation of insert size", "bp"),
    "nuis_insert_p10": ("per-person sequencing", "per person", "10th percentile of insert size", "bp"),
    "nuis_insert_p90": ("per-person sequencing", "per person", "90th percentile of insert size", "bp"),
    "nuis_read_length_mode": ("per-person sequencing", "per person", "Most common read length (constant 150 bp in this cohort)", "bp"),
    "nuis_mapq0_frac_global": ("per-person sequencing", "per person", "Genome-wide fraction of reads with MAPQ 0 in the random windows", "fraction"),
    "nuis_n_windows_used": ("per-person sequencing", "per person", "Number of random windows successfully sampled (constant 3,000)", "count"),
    "nuis_n_reads_sampled": ("per-person sequencing", "per person", "Number of reads in the random windows", "reads"),
    "cov_period": ("locus covariate", "annotation BED", "Motif period from the annotation", "bp"),
    "cov_ref_units_est": ("locus covariate", "annotation BED", "Reference repeat count: array span / period", "repeat units"),
    "cov_array_span_bp": ("locus covariate", "annotation BED", "Reference array length", "bp"),
    "cov_log_array_span_bp": ("locus covariate", "annotation BED", "Natural log of the reference array length", "log bp"),
    "cov_consensus_gc": ("locus covariate", "vamos annotation", "GC fraction of the consensus motif", "fraction"),
    "cov_consensus_base_entropy": ("locus covariate", "vamos annotation", "Shannon entropy of base composition of the consensus motif", "bits"),
    "cov_consensus_len": ("locus covariate", "vamos annotation", "Length of the consensus motif", "bp"),
    "cov_is_segdup": ("locus covariate", "vamos annotation", "1 if the locus lies in a segmental duplication", "0/1"),
    "cov_mr_segdup_identity": ("locus covariate", "mutation-rate table", "Sequence identity of the overlapping segmental duplication (-1 if none)", "fraction"),
    "cov_mr_segdup_cnv": ("locus covariate", "mutation-rate table", "Copy-number variability flag of the overlapping segmental duplication", "code"),
    "cov_mr_cons_gc": ("locus covariate", "mutation-rate table", "GC of the consensus motif (mutation-rate table)", "fraction"),
    "cov_mr_all_gc": ("locus covariate", "mutation-rate table", "GC of the whole locus sequence", "fraction"),
    "cov_mr_te_frac_locus": ("locus covariate", "mutation-rate table", "Fraction of the locus overlapped by transposable elements", "fraction"),
    "cov_mr_te_frac_element": ("locus covariate", "mutation-rate table", "Fraction of the overlapping transposable element covered by the locus", "fraction"),
    "cov_mr_dist_to_recomb": ("locus covariate", "mutation-rate table", "Distance to the nearest recombination hotspot", "bp"),
    "cov_mr_size_bp": ("locus covariate", "mutation-rate table", "Locus size in the mutation-rate table", "bp"),
    "cov_mr_period": ("locus covariate", "mutation-rate table", "Period in the mutation-rate table", "bp"),
}


ANOVA = ("ANOVA table across continental groups of the long-read panel (external, precomputed)")
for tag, what in (("anno_len", "annotated repeat length"), ("ed2_major", "edit distance to the major allele")):
    FIXED.update({
        f"cov_{tag}_grand_mean": ("locus covariate", "ANOVA table", f"Mean {what} across the long-read panel", "units"),
        f"cov_{tag}_between_var": ("locus covariate", "ANOVA table", f"Between-population variance of {what}", "variance"),
        f"cov_{tag}_within_var": ("locus covariate", "ANOVA table", f"Within-population variance of {what}", "variance"),
        f"cov_{tag}_f": ("locus covariate", "ANOVA table", f"ANOVA F statistic for population differences in {what}", "F"),
        f"cov_{tag}_p": ("locus covariate", "ANOVA table", f"ANOVA p-value for population differences in {what}", "p"),
        f"cov_{tag}_effect_size": ("locus covariate", "ANOVA table", f"Effect size of population differences in {what}", "effect size"),
    })
PANEL = "long-read panel"
FIXED.update({
    "cov_n_called": ("panel statistic", PANEL, "Panel people with both alleles called at this repeat", "people"),
    "cov_n_alleles": ("panel statistic", PANEL, "Distinct alleles carried by panel people", "count"),
    "cov_h_exp": ("panel statistic", PANEL, "Expected heterozygosity from panel allele-length frequencies", "0-1"),
    "cov_h_obs": ("panel statistic", PANEL, "Share of panel people whose two alleles differ", "0-1"),
    "cov_mean_combined_units": ("panel statistic", PANEL, "Mean combined length (both alleles) in the panel", "repeat units"),
    "cov_var_combined_units": ("panel statistic", PANEL, "Variance of combined length in the panel", "units^2"),
    "cov_min_combined_units": ("panel statistic", PANEL, "Smallest combined length in the panel", "repeat units"),
    "cov_max_combined_units": ("panel statistic", PANEL, "Largest combined length in the panel", "repeat units"),
    "cov_mean_combined_bp": ("panel statistic", PANEL, "Mean combined length in the panel", "bp"),
    "cov_var_combined_bp": ("panel statistic", PANEL, "Variance of combined length in the panel", "bp^2"),
    "cov_n_distinct_hap_units": ("panel statistic", PANEL, "Distinct haplotype lengths in the panel", "count"),
    "cov_panel_cn_range_units": ("panel statistic", PANEL, "Max minus min combined length in the panel", "repeat units"),
    "cov_panel_median_units": ("panel statistic", PANEL, "Median combined length in the panel", "repeat units"),
    "cov_panel_p5_units": ("panel statistic", PANEL, "5th percentile of combined length in the panel", "repeat units"),
    "cov_panel_p95_units": ("panel statistic", PANEL, "95th percentile of combined length in the panel", "repeat units"),
    "cov_panel_core90_units": ("panel statistic", PANEL, "Core spread: p95 - p5 of combined length in the panel", "repeat units"),
    "cov_panel_n_outliers": ("panel statistic", PANEL, "Panel people beyond the core by more than 2 x max(core spread, 1)", "people"),
    "cov_panel_outlier_frac": ("panel statistic", PANEL, "Outlier people as a share of called panel people", "fraction"),
    "cov_n_motifs": ("locus covariate", "vamos annotation", "Number of motifs in the repeat's vamos motif set", "count"),
    "cov_motif_len_mean": ("locus covariate", "vamos annotation", "Mean length of the repeat's vamos motifs", "bp"),
    "cov_motif_len_min": ("locus covariate", "vamos annotation", "Shortest vamos motif", "bp"),
    "cov_motif_len_max": ("locus covariate", "vamos annotation", "Longest vamos motif", "bp"),
    "cov_motif_len_cv": ("locus covariate", "vamos annotation", "Coefficient of variation of vamos motif lengths", "ratio"),
    "cov_mr_motif_entropy": ("locus covariate", "mutation-rate table", "Entropy of motif usage across alleles", "bits"),
    "cov_mr_snp_density": ("locus covariate", "mutation-rate table", "SNP density in the surrounding window", "per bp"),
    "cov_mr_std_len_units": ("locus covariate", "mutation-rate table", "Standard deviation of allele length across the study's assemblies", "repeat units"),
    "cov_mr_homo_len": ("locus covariate", "mutation-rate table", "Share of people with two equal-length alleles", "fraction"),
    "cov_mr_hete_len": ("locus covariate", "mutation-rate table", "Share of people with two different-length alleles", "fraction"),
    "cov_mr_homo_comp": ("locus covariate", "mutation-rate table", "Share of people with two alleles of identical motif composition", "fraction"),
    "cov_mr_hete_comp": ("locus covariate", "mutation-rate table", "Share of people with two alleles of different motif composition", "fraction"),
    "cov_mr_num_pattern": ("locus covariate", "mutation-rate table", "Number of distinct motif-composition patterns across alleles", "count"),
    "cov_mr_ave_block": ("locus covariate", "mutation-rate table", "Average length of uninterrupted motif blocks", "units"),
    "cov_mr_mut_rate": ("locus covariate", "mutation-rate table", "Estimated mutation rate of the repeat", "per generation"),
    "cov_mr_mut_rate_normlen": ("locus covariate", "mutation-rate table", "Mutation rate normalised for repeat length", "per generation"),
    "cov_reciprocal_overlap": ("locus covariate", "annotation join", "Reciprocal overlap between our repeat coordinates and the matched mutation-rate-table record (join quality)", "0-1"),
    "plc_pred": ("model estimate", "per-locus calibration", "Estimate from the repeat's own calibration curve (short-read signal fitted to training people's long-read truth; cross-fitted for training people)", "repeat units"),
    "plc_available": ("model estimate", "per-locus calibration", "1 if the repeat has a calibration curve", "0/1"),
    "plc_resid_sd": ("model estimate", "per-locus calibration", "Typical error of the repeat's calibration curve on training people", "repeat units"),
    "transfer_pred": ("model estimate", "calibration transfer", "Estimate from the cross-repeat transfer model (cross-fitted for training people)", "repeat units"),
    "transfer_ood": ("model estimate", "calibration transfer", "Novelty score: distance from the repeat's annotation to the nearest training repeat", "distance"),
})

# where each source's values come from, and caveats
ORIGIN = {
    "per locus x person": "short reads", "per person": "short reads",
    "annotation BED": "reference / repeat BED", "reference": "reference / repeat BED", "locus": "reference / repeat BED",
    "annotation join": "reference / repeat BED",
    "vamos annotation": "long-read annotation", "mutation-rate table": "long-read annotation",
    "ANOVA table": "long-read annotation",
    PANEL: "long-read panel genotypes",
    "per-locus calibration": "model estimate", "calibration transfer": "model estimate",
}
POPULATION_MR = {"cov_mr_all_gc", "cov_mr_motif_entropy", "cov_mr_std_len_units", "cov_mr_homo_len", "cov_mr_hete_len",
                 "cov_mr_homo_comp", "cov_mr_hete_comp", "cov_mr_num_pattern", "cov_mr_ave_block", "cov_mr_mut_rate",
                 "cov_mr_mut_rate_normlen", "cov_mr_snp_density"}


def origin_note(f: str, cat: str, src: str) -> str:
    if cat in ("in-repeat reads", "motif composition") or f.startswith("partial_repeat"):
        return "read counts from short reads; which reads count as repeat reads uses the vamos motif set"
    if src == PANEL:
        return "recomputed from vamos genotypes with validation and test people removed; known repeats only"
    if src == "ANOVA table" or f in POPULATION_MR:
        return "population summary precomputed by the source study over its long-read samples, which may include our validation/test people"
    if src in ("vamos annotation", "mutation-rate table"):
        return "per-repeat descriptor from a long-read resource; same for every person"
    if src == "per-locus calibration":
        return "input is short reads; curve fitted against training people's long-read lengths"
    if src == "calibration transfer":
        return "input is short reads + annotation; model trained against training people's long-read lengths"
    return ""


def describe(f: str):
    if f in FIXED:
        return FIXED[f]
    m = re.fullmatch(r"array_depth_d(\d+)", f)
    if m:
        return ("read depth", "per locus x person", f"Mean per-base depth over the array, {MQ[m[1]]}", "x")
    m = re.fullmatch(r"array_reads_d(\d+)", f)
    if m:
        return ("read depth", "per locus x person", f"Number of reads overlapping the array, {MQ[m[1]]}", "reads")
    m = re.fullmatch(r"flank_depth_d(\d+)_w(\d+)", f)
    if m:
        return ("read depth", "per locus x person", f"Mean per-base depth over {m[2]} bp on each side of the array, {MQ[m[1]]}", "x")
    m = re.fullmatch(r"depth_ratio_d(\d+)_w(\d+)", f)
    if m:
        return ("read depth", "per locus x person", f"Array depth / flank depth ({m[2]} bp flanks), {MQ[m[1]]}; ~1 when the person's array matches the reference length", "ratio")
    raise KeyError(f)


def split_gain(model, feats):
    """Share of total split gain per feature (sklearn HistGradientBoosting or the xgboost wrapper)."""
    gain = collections.Counter()
    if hasattr(model, "model") and hasattr(model.model, "get_booster"):          # _XGBEarlyStop
        for k, g in model.model.get_booster().get_score(importance_type="total_gain").items():
            gain[feats[int(k[1:])] if k.startswith("f") and k[1:].isdigit() else k] += float(g)
    else:
        for it in model._predictors:
            for tree in it:
                n = tree.nodes
                for fi, g, leaf in zip(n["feature_idx"], n["gain"], n["is_leaf"]):
                    if not leaf:
                        gain[feats[fi]] += float(g)
    tot = sum(gain.values()) or 1.0
    return {f: gain[f] / tot for f in feats}


def main() -> None:
    rows = []
    sel = dict(pl.read_parquet(M / "per_locus_calibration.parquet")["signal"].value_counts().drop_nulls().iter_rows())
    n_fit = sum(sel.values())
    for f in _SIGNAL_SHORTLIST:
        cat, src, desc, unit = describe(f)
        rows.append(dict(model="1_per_locus_calibration", feature=f, category=cat, source=src,
                         data_origin=ORIGIN[src], origin_note=origin_note(f, cat, src), units=unit,
                         description=desc, usage_metric="share of calibrated loci choosing this signal",
                         usage_value=round(sel.get(f, 0) / n_fit, 4), loci_choosing=sel.get(f, 0)))
    tm = joblib.load(M / "transfer_model.joblib")
    feats = tm["features"]
    gain = split_gain(tm["model"], feats)
    for f in feats:
        cat, src, desc, unit = describe(f)
        rows.append(dict(model="2_calibration_transfer", feature=f, category=cat, source=src,
                         data_origin=ORIGIN[src], origin_note=origin_note(f, cat, src), units=unit,
                         description=desc, usage_metric="share of total split gain",
                         usage_value=round(gain[f], 5), loci_choosing=None))
    # final known-repeat model: inputs beyond the transfer model's (the new-repeat model uses exactly the 116 above)
    gm = joblib.load(M / "global_model.joblib")["known"]
    kfeats = gm["features"]
    kgain = split_gain(gm["point"], kfeats)
    for f in [f for f in kfeats if f not in set(feats)]:
        cat, src, desc, unit = describe(f)
        rows.append(dict(model="3_final_known_repeats_extra", feature=f, category=cat, source=src,
                         data_origin=ORIGIN[src], origin_note=origin_note(f, cat, src), units=unit,
                         description=desc, usage_metric="share of total split gain in the final known-repeat point model",
                         usage_value=round(kgain[f], 5), loci_choosing=None))
    out = pl.DataFrame(rows, infer_schema_length=None)
    out.write_csv(M / "feature_inventory.csv")
    print(out.group_by("model").len())
    print(out.group_by("model", "data_origin").len().sort("model", "data_origin"))
    print(out.filter(pl.col("model") == "2_calibration_transfer").sort("usage_value", descending=True)
          .select("feature", "category", "usage_value").head(10))
    print(f"wrote {M / 'feature_inventory.csv'}")


if __name__ == "__main__":
    main()
