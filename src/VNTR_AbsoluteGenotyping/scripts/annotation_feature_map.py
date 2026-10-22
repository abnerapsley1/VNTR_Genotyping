"""Map every field of the per-repeat annotation files to the model features derived from it,
and record whether each feature is used by the current default model.

    python scripts/annotation_feature_map.py  ->  docs/annotation_feature_map.csv
"""
import sys
from pathlib import Path

import joblib
import polars as pl

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from pipeline.model.assemble import EXTERNAL_POPULATION_FEATURES  # noqa: E402
from pipeline.model.panel_stats import PANEL_COLUMNS  # noqa: E402

M = HERE / "work" / "model"
BED = "repeat BED (v1_VNTR_Annotations_03_07_2026.bed)"
VCF = "vamos VCF (vamosExpanded_v3.0_GRCh38_oriMotifs.vcf)"
TSV = "vamos TSV (vamosExpanded_v3.0_GRCh38_oriMotifs.tsv)"
MR = "mutation-rate table (mutRate_GRCh38_masterFile.tsv)"
A1 = "ANOVA, allele length (anova.annoLen.tsv)"
A2 = "ANOVA, edit distance to major allele (anova.ed2Major.tsv)"

# (source file, source field, model feature, how it is derived)
ROWS = [
    (BED, "chrom, start, end", "(locus definition)", "defines the 60,504 repeats and their coordinates; all files are joined to these"),
    (BED, "column 6 (period)", "period, cov_period", "motif length as given"),
    (BED, "start, end", "cov_array_span_bp", "end - start"),
    (BED, "start, end", "cov_log_array_span_bp", "log(end - start)"),
    (BED, "start, end, period", "cov_ref_units_est", "(end - start) / period"),
    (BED, "start, end", "array_gc, flank_gc_w2000 (with GRCh38 sequence)", "GC of reference array / 2 kb flank"),
    (BED, "gene", "(not a feature)", "kept for reporting only"),
    (VCF, "GT + ALTANNO + RU, per person", "TRAINING TARGET combined_units / combined_bp", "sum of both alleles' motif tokens (units) or motif lengths (bp)"),
    (VCF, "RU (motif list)", "irr_*, partial_repeat_*, anchored_irr_*, mc_* (short-read signals)", "which reads count as repeat reads: 5-mers of the RU motifs; mc_* assign read windows to RU motifs"),
    (VCF, "RU[0]", "mc_frac_ref_consensus", "first RU motif treated as the consensus"),
    (VCF, "POS, INFO END", "cov_reciprocal_overlap", "overlap between the VCF record and our BED repeat (join quality)"),
    (VCF, "ALTANNO", "cov_n_alleles", "distinct alleles carried by panel people (recomputed without validation/test people)"),
]
for c in PANEL_COLUMNS:
    if c != "n_alleles":
        ROWS.append((VCF, "GT + ALTANNO + RU, panel people", f"cov_{c}",
                     "panel summary recomputed without validation/test people (panel_stats.py)"))
ROWS += [
    (TSV, "motifs (column 4)", "cov_n_motifs", "number of motifs"),
    (TSV, "motifs", "cov_motif_len_mean / _min / _max", "motif length statistics"),
    (TSV, "motifs", "cov_motif_len_cv", "(max - min) / mean motif length"),
    (TSV, "consensus (column 8)", "cov_consensus_gc", "GC of the consensus motif"),
    (TSV, "consensus", "cov_consensus_base_entropy", "Shannon entropy of the consensus motif's bases"),
    (TSV, "consensus", "cov_consensus_len", "consensus motif length"),
    (TSV, "genomic_context", "cov_is_segdup", "1 if the context mentions segDup"),
    (TSV, "cluster_type, vntr_type, genomic_context, bio_context_1-3", "(not features)", "text categories; not used by the model"),
    (TSV, "period, sc1-sc4", "(not read)", ""),
]
MRMAP = {"segDupIdentity": "cov_mr_segdup_identity", "segDupCNV": "cov_mr_segdup_cnv", "consGC": "cov_mr_cons_gc",
         "allGC": "cov_mr_all_gc", "op1": "cov_mr_te_frac_locus", "op2": "cov_mr_te_frac_element",
         "entropy": "cov_mr_motif_entropy", "snpDensity": "cov_mr_snp_density", "stdLenAll": "cov_mr_std_len_units",
         "homoLenAll": "cov_mr_homo_len", "heteLenAll": "cov_mr_hete_len", "homoComAll": "cov_mr_homo_comp",
         "heteComAll": "cov_mr_hete_comp", "numPattern": "cov_mr_num_pattern", "aveBlock": "cov_mr_ave_block",
         "dist2Recom": "cov_mr_dist_to_recomb", "rate": "cov_mr_mut_rate", "rateNormAfter": "cov_mr_mut_rate_normlen",
         "size": "cov_mr_size_bp", "period": "cov_mr_period"}
for src, dst in MRMAP.items():
    ROWS.append((MR, src, dst, "value of the best-overlapping record (reciprocal overlap >= 0.5)"))
ROWS.append((MR, "coding, exon, transcribed, telomere, hotspot, type", "(not features)", "text categories; read but not used by the model"))
ROWS.append((MR, "all other columns (tree, rate components, purity, alleles, recombination hotspot coordinates, ...)", "(not read)", ""))
for path, tag in ((A1, "anno_len"), (A2, "ed2_major")):
    for src, suf in (("grandMean", "grand_mean"), ("betweenGroupVariance_mean", "between_var"),
                     ("withinGroupVariance_mean", "within_var"), ("f", "f"), ("p", "p"), ("effectSize", "effect_size")):
        ROWS.append((path, src, f"cov_{tag}_{suf}", "value of the best-overlapping record (reciprocal overlap >= 0.5)"))
    ROWS.append((path, "per-population means (EUR AMR AFR SAS EAS MDE), withinGroupVariance, betweenGroupVariance, df",
                 "(not read)", ""))


def main():
    g = joblib.load(M / "global_model.joblib")
    known, new = set(g["known"]["features"]), set(g["novel"]["features"])
    out = []
    for src, field, feat, how in ROWS:
        names = [f.strip() for f in feat.replace("/", ",").split(",")]
        names = [n if n.startswith(("cov_", "array_", "flank_", "period", "irr", "partial", "anchored", "mc_"))
                 else n for n in names]
        if feat.startswith("("):
            status = "not a model feature" if "TARGET" not in feat else "training target"
        elif feat.startswith("TRAINING"):
            status = "training target"
        elif "short-read signals" in feat or feat.startswith(("array_gc", "mc_frac")):
            status = "used (known and new repeats)"
        else:
            n = feat.split()[0].rstrip(",")
            if n.startswith("cov_motif_len_mean"):
                n = "cov_motif_len_mean"
            if n in EXTERNAL_POPULATION_FEATURES:
                status = "excluded (external population summary; may include held-out people)"
            elif n in new and n in known:
                status = "used (known and new repeats)"
            elif n in known:
                status = "used (known repeats only)"
            else:
                status = "not used"
        out.append(dict(source_file=src, source_field=field, model_feature=feat, derivation=how, status=status))
    df = pl.DataFrame(out)
    df.write_csv(HERE / "docs" / "annotation_feature_map.csv")
    print(df.group_by("source_file", "status").len().sort("source_file", "status"))


if __name__ == "__main__":
    main()
