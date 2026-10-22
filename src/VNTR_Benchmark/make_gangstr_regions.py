"""GangSTR regions file for the bundled VNTR loci.

GangSTR reference rows: chrom, start, end (1-based, inclusive), period, motif.
The bundled annotation already stores 1-based starts equal to vamos POS, so
coordinates are used as is.  The motif is vamos's consensus motif for the
locus and the period is its length (GangSTR documents motifs up to 20 bp; all
loci are written, and loci it cannot call are reported as no-calls).

Writes src/VNTR_Benchmark/work/gangstr/regions.bed (sorted) and
regions_index.tsv mapping each row back to locus_id.
"""
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "src" / "VNTR_Benchmark" / "work" / "gangstr"

cols = ["chrom", "start", "end", "motifs", "cluster_type", "vntr_type", "period", "consensus",
        "sc1", "sc2", "sc3", "sc4", "genomic_context", "bio1", "bio2", "bio3"]
tsv = pl.read_csv(ROOT / "data" / "1000Genomes" / "vamosExpanded_v3.0_GRCh38_oriMotifs.tsv",
                  separator="\t", has_header=False, infer_schema_length=0,
                  truncate_ragged_lines=True, new_columns=cols).select("chrom", "start", "end", "consensus")
tsv = tsv.with_columns(locus_id=pl.col("chrom") + ":" + pl.col("start") + "-" + pl.col("end"))
loci = pl.read_parquet(ROOT / "src" / "VNTR_AbsoluteGenotyping" / "work" / "ground_truth" / "loci_targets.parquet")
d = (loci.select("locus_id", "chrom", "array_start0", "array_end0")
     .join(tsv.select("locus_id", "consensus"), on="locus_id", how="left")
     .with_columns(start=pl.col("array_start0") + 1, end=pl.col("array_end0"),
                   motif=pl.col("consensus").str.to_uppercase())
     .with_columns(period=pl.col("motif").str.len_chars())
     .filter(pl.col("motif").is_not_null() & (pl.col("period") > 0)))
chrom_order = {f"chr{c}": i for i, c in enumerate(list(range(1, 23)) + ["X", "Y"])}
d = d.with_columns(o=pl.col("chrom").replace_strict(chrom_order, default=99)).sort("o", "start")
OUT.mkdir(parents=True, exist_ok=True)
d.select("chrom", "start", "end", "period", "motif").write_csv(OUT / "regions.bed", separator="\t", include_header=False)
d.select("locus_id", "chrom", "start", "end", "period", "motif").write_csv(OUT / "regions_index.tsv", separator="\t")
print(f"{d.height:,} of {loci.height:,} loci written; period <= 20 bp: {(d['period'] <= 20).mean():.1%}")
