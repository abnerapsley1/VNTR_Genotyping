"""
build_locus_covariates.py — one covariate row per bundled-BED locus.

Sources (each optional except the first; a *_joined flag column records which hit):
  1. work/ground_truth/gt_locus_summary.parquet   (from parse_vamos; REQUIRED)
       period, ref_units_est, panel h_exp/h_obs, panel combined-length mean/var/range
  2. inputs.vamos_tsv     exact (chrom,start,end) join
       cluster_type, vntr_type, consensus motif, genomic_context, bio_context_*,
       motif list -> n_motifs / motif-length stats / motif GC / motif entropy
  3. inputs.mutrate_tsv   reciprocal-overlap join (streamed; ~1.3 GB)
       segDupIdentity, segDupCNV, consGC, allGC, TE overlap (op1/op2), coding/exon/
       transcribed, motif entropy, snpDensity, length het/homozygosity, numPattern,
       aveBlock, dist2Recom, telomere, mutation rate, hotspot
  4. inputs.anova_annolen_tsv / inputs.anova_ed2major_tsv  reciprocal-overlap join
       population-structure effect size / between-group variance / F / p

Output: work/covariates/locus_covariates.parquet

Usage
-----
    python -m pipeline.covariates.build_locus_covariates
"""

from __future__ import annotations

import argparse
import math
from bisect import bisect_left, bisect_right
from pathlib import Path

import polars as pl

from ..paths import load_config

_AUX_MIN_OVERLAP = 0.5   # reciprocal-overlap threshold for streamed aux sources


# ---------------------------------------------------------------------------
# small sequence helpers
# ---------------------------------------------------------------------------

def _gc(seq: str) -> float | None:
    if not seq:
        return None
    seq = seq.upper()
    n = sum(c in "ACGT" for c in seq)
    return (seq.count("G") + seq.count("C")) / n if n else None


def _shannon_bases(seq: str) -> float | None:
    if not seq:
        return None
    seq = seq.upper()
    counts = {b: seq.count(b) for b in "ACGT"}
    tot = sum(counts.values())
    if tot == 0:
        return None
    h = 0.0
    for c in counts.values():
        if c:
            p = c / tot
            h -= p * math.log2(p)
    return h  # 0..2 bits


def _num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# interval index for overlap joins
# ---------------------------------------------------------------------------

class IntervalIndex:
    """Per-chrom sorted-start index; query() returns best-overlap locus row-idx."""

    def __init__(self, loci: pl.DataFrame):
        self.by_chrom: dict[str, tuple[list[int], list[int], list[int]]] = {}
        rows = loci.select(["chrom", "start", "end"]).with_row_index("i")
        for chrom, grp in rows.group_by("chrom"):
            g = grp.sort("start")
            self.by_chrom[chrom[0]] = (
                g["start"].to_list(), g["end"].to_list(), g["i"].to_list()
            )

    def best(self, chrom: str, s: int, e: int, min_ro: float) -> int | None:
        hit = self.by_chrom.get(chrom)
        if hit is None:
            return None
        starts, ends, idxs = hit
        lo = bisect_left(starts, s - 200_000)
        hi = bisect_right(starts, e + 200_000)
        best_i, best_ro = None, 0.0
        span_q = max(1, e - s)
        for j in range(lo, hi):
            inter = min(e, ends[j]) - max(s, starts[j])
            if inter <= 0:
                continue
            ro = min(inter / span_q, inter / max(1, ends[j] - starts[j]))
            if ro >= min_ro and ro > best_ro:
                best_i, best_ro = idxs[j], ro
        return best_i


# ---------------------------------------------------------------------------
# source 2: vamos annotation TSV (exact join)
# ---------------------------------------------------------------------------

_VAMOS_TSV_COLS = [
    "chrom", "start", "end", "motifs", "cluster_type", "vntr_type",
    "period", "consensus", "sc1", "sc2", "sc3", "sc4",
    "genomic_context", "bio_context_1", "bio_context_2", "bio_context_3",
]


def load_vamos_tsv(path: Path) -> pl.DataFrame:
    df = pl.read_csv(
        path, separator="\t", has_header=False, new_columns=_VAMOS_TSV_COLS,
        infer_schema_length=0, truncate_ragged_lines=True,
    )
    motif_feats = []
    for motifs, cons in zip(df["motifs"].to_list(), df["consensus"].to_list()):
        mlist = [m for m in (motifs or "").split(",") if m]
        mlens = [len(m) for m in mlist]
        motif_feats.append(dict(
            n_motifs=len(mlist),
            motif_len_mean=(sum(mlens) / len(mlens)) if mlens else None,
            motif_len_min=min(mlens) if mlens else None,
            motif_len_max=max(mlens) if mlens else None,
            consensus_gc=_gc(cons or ""),
            consensus_base_entropy=_shannon_bases(cons or ""),
            consensus_len=len(cons or "") or None,
        ))
    mf = pl.DataFrame(motif_feats)
    df = df.with_columns(
        pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64)
    ).with_columns(
        locus_key=pl.col("chrom") + ":" + pl.col("start").cast(pl.Utf8)
        + "-" + pl.col("end").cast(pl.Utf8)
    )
    keep = ["locus_key", "cluster_type", "vntr_type", "genomic_context",
            "bio_context_1", "bio_context_2", "bio_context_3"]
    return pl.concat([df.select(keep), mf], how="horizontal")


# ---------------------------------------------------------------------------
# source 3: mutRate masterFile  (streamed overlap join)
# ---------------------------------------------------------------------------

_MUTRATE_KEEP = {
    "segDupIdentity": "mr_segdup_identity",
    "segDupCNV": "mr_segdup_cnv",
    "consGC": "mr_cons_gc",
    "allGC": "mr_all_gc",
    "op1": "mr_te_frac_locus",
    "op2": "mr_te_frac_element",
    "coding": "mr_coding",
    "exon": "mr_exon",
    "transcribed": "mr_transcribed",
    "entropy": "mr_motif_entropy",
    "snpDensity": "mr_snp_density",
    "stdLenAll": "mr_std_len_units",
    "homoLenAll": "mr_homo_len",
    "heteLenAll": "mr_hete_len",
    "homoComAll": "mr_homo_comp",
    "heteComAll": "mr_hete_comp",
    "numPattern": "mr_num_pattern",
    "aveBlock": "mr_ave_block",
    "dist2Recom": "mr_dist_to_recomb",
    "telomere": "mr_telomere",
    "rate": "mr_mut_rate",
    "rateNormAfter": "mr_mut_rate_normlen",
    "hotspot": "mr_hotspot",
    "size": "mr_size_bp",
    "period": "mr_period",
    "type": "mr_type",
}
_CATEGORICAL_MR = {"mr_coding", "mr_exon", "mr_transcribed", "mr_telomere",
                   "mr_hotspot", "mr_type"}


def join_streamed_tsv(path: Path, idx: IntervalIndex, n_loci: int,
                      col_map: dict[str, str], prefix_tag: str,
                      categorical: set[str] | None = None) -> pl.DataFrame:
    """Stream a big TSV with a header; assign each row to its best-overlap locus."""
    categorical = categorical or set()
    out: list[dict | None] = [None] * n_loci
    n_hit = 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_i = {c: k for k, c in enumerate(header)}
        need = {src: col_i[src] for src in col_map if src in col_i}
        ci_chr = col_i.get("chr", col_i.get("chrom", 0))
        ci_s = col_i.get("start", 1)
        ci_e = col_i.get("end", 2)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= ci_e:
                continue
            try:
                s, e = int(f[ci_s]), int(f[ci_e])
            except ValueError:
                continue
            li = idx.best(f[ci_chr], s, e, _AUX_MIN_OVERLAP)
            if li is None or out[li] is not None:
                continue
            row = {}
            for src, dst in col_map.items():
                if src not in need:
                    continue
                raw = f[need[src]]
                row[dst] = raw if dst in categorical else _num(raw)
            out[li] = row
            n_hit += 1
    print(f"  {prefix_tag}: matched {n_hit:,}/{n_loci:,} loci")
    filled = [r if r is not None else {} for r in out]
    return pl.DataFrame(filled, infer_schema_length=n_loci)


# ---------------------------------------------------------------------------
# source 4: anova files (overlap join, small enough to read whole)
# ---------------------------------------------------------------------------

def join_anova(path: Path, idx: IntervalIndex, n_loci: int, tag: str) -> pl.DataFrame:
    df = pl.read_csv(path, separator="\t", infer_schema_length=0)
    out: list[dict] = [{} for _ in range(n_loci)]
    seen = set()
    for r in df.iter_rows(named=True):
        try:
            s, e = int(r["start"]), int(r["end"])
        except (ValueError, TypeError, KeyError):
            continue
        li = idx.best(r["chr"], s, e, _AUX_MIN_OVERLAP)
        if li is None or li in seen:
            continue
        seen.add(li)
        out[li] = {
            f"{tag}_grand_mean": _num(r.get("grandMean")),
            f"{tag}_between_var": _num(r.get("betweenGroupVariance_mean")),
            f"{tag}_within_var": _num(r.get("withinGroupVariance_mean")),
            f"{tag}_f": _num(r.get("f")),
            f"{tag}_p": _num(r.get("p")),
            f"{tag}_effect_size": _num(r.get("effectSize")),
        }
    print(f"  anova/{tag}: matched {len(seen):,}/{n_loci:,} loci")
    return pl.DataFrame(out, infer_schema_length=n_loci)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def build(cfg) -> Path:
    gt_dir = cfg.path("outputs.ground_truth_dir")
    summ_path = gt_dir / "gt_locus_summary.parquet"
    if not summ_path.exists():
        raise FileNotFoundError(
            f"{summ_path} missing — run `python -m pipeline.ground_truth.parse_vamos --run` first"
        )
    base = pl.read_parquet(summ_path)
    n = base.height
    print(f"base loci: {n:,}")
    base = base.with_columns(
        locus_key=pl.col("chrom") + ":" + pl.col("start").cast(pl.Utf8)
        + "-" + pl.col("end").cast(pl.Utf8),
        array_span_bp=(pl.col("end") - pl.col("start")),
        panel_cn_range_units=(pl.col("max_combined_units") - pl.col("min_combined_units")),
    )
    idx = IntervalIndex(base)

    parts = [base]

    tsv_path = cfg.path("inputs.vamos_tsv")
    if tsv_path.exists():
        print("source: vamos annotation TSV (exact join)")
        vt = load_vamos_tsv(tsv_path)
        base_keys = base.select(["locus_key"]).with_row_index("row")
        vt = base_keys.join(vt, on="locus_key", how="left").drop("row", "locus_key")
        print(f"  vamos TSV: matched "
              f"{vt['cluster_type'].is_not_null().sum():,}/{n:,} loci")
        parts.append(vt)
    else:
        print("source: vamos annotation TSV — SKIPPED (not found)")

    mr_path = cfg.path("inputs.mutrate_tsv")
    if mr_path.exists():
        print("source: mutRate masterFile (streamed overlap join)")
        parts.append(join_streamed_tsv(mr_path, idx, n, _MUTRATE_KEEP,
                                       "mutRate", _CATEGORICAL_MR))
    else:
        print("source: mutRate masterFile — SKIPPED (not found)")

    for key, tag in [("inputs.anova_annolen_tsv", "anno_len"),
                     ("inputs.anova_ed2major_tsv", "ed2_major")]:
        p = cfg.path(key)
        if p.exists():
            print(f"source: {p.name} (overlap join)")
            parts.append(join_anova(p, idx, n, tag))
        else:
            print(f"source: {key} — SKIPPED (not found)")

    cov = pl.concat(parts, how="horizontal")
    # derived
    cov = cov.with_columns(
        motif_len_cv=(pl.col("motif_len_max") - pl.col("motif_len_min"))
        / pl.col("motif_len_mean"),
        log_array_span_bp=pl.col("array_span_bp").clip(1).log(),
        is_segdup=pl.when(pl.col("genomic_context").str.contains("segDup"))
        .then(1).when(pl.col("genomic_context").is_not_null()).then(0).otherwise(None),
    )

    out_dir = cfg.ensure_dir("outputs.covariates_dir")
    out_path = out_dir / "locus_covariates.parquet"
    cov.write_parquet(out_path)
    print(f"\nwrote {out_path}  ({cov.height:,} loci x {cov.width} cols)")
    nullpct = (cov.null_count().transpose(include_header=True, header_name="col",
                                          column_names=["nulls"])
               .with_columns(pct=(pl.col("nulls") / cov.height * 100).round(1))
               .sort("pct", descending=True))
    print("top columns by % missing:")
    print(nullpct.head(12))
    return out_path


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    args = ap.parse_args(argv)
    build(load_config(args.config))


if __name__ == "__main__":
    main()
