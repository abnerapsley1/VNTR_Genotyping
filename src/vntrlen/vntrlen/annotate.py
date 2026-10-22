"""Per-repeat annotation features (the cov_* model inputs that are not panel statistics).

Computed from the BED and motifs:
  cov_period, cov_ref_units_est, cov_array_span_bp, cov_log_array_span_bp,
  cov_n_motifs, cov_motif_len_mean/min/max/cv,
  cov_consensus_gc, cov_consensus_base_entropy, cov_consensus_len (first motif),
  cov_reciprocal_overlap (1: the user's coordinates define the repeat)
Looked up in the bundled TRCompDB annotation table (best overlap >= 0.5):
  cov_mr_segdup_identity, cov_mr_segdup_cnv, cov_mr_cons_gc, cov_mr_te_frac_locus,
  cov_mr_te_frac_element, cov_mr_dist_to_recomb, cov_mr_size_bp, cov_mr_period,
  cov_is_segdup (segdup identity >= 0)
A user annotation table (--annotation) overrides any of these per repeat; values
neither supplied nor found stay missing (the models accept missing values).
"""
from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

import numpy as np
import polars as pl

TRCOMPDB_COLS = {"segDupIdentity": "cov_mr_segdup_identity", "segDupCNV": "cov_mr_segdup_cnv",
                 "consGC": "cov_mr_cons_gc", "op1": "cov_mr_te_frac_locus", "op2": "cov_mr_te_frac_element",
                 "dist2Recom": "cov_mr_dist_to_recomb", "size": "cov_mr_size_bp", "period": "cov_mr_period"}
MIN_OVERLAP = 0.5


def _gc(s: str) -> float | None:
    return (s.count("G") + s.count("C")) / len(s) if s else None


def _entropy(s: str) -> float | None:
    if not s:
        return None
    c = Counter(s)
    return -sum(v / len(s) * math.log2(v / len(s)) for v in c.values())


def bed_covariates(bed: pl.DataFrame) -> pl.DataFrame:
    rows = []
    for r in bed.iter_rows(named=True):
        motifs = r["motifs"].split(",")
        lens = [len(m) for m in motifs]
        span = r["end"] - r["start"]
        mean = sum(lens) / len(lens)
        rows.append(dict(
            locus_id=r["locus_id"], cov_period=float(r["period"]), cov_array_span_bp=float(span),
            cov_log_array_span_bp=math.log(max(span, 1)), cov_ref_units_est=span / r["period"],
            cov_n_motifs=float(len(motifs)), cov_motif_len_mean=mean, cov_motif_len_min=float(min(lens)),
            cov_motif_len_max=float(max(lens)), cov_motif_len_cv=(max(lens) - min(lens)) / mean,
            cov_consensus_gc=_gc(motifs[0]), cov_consensus_base_entropy=_entropy(motifs[0]),
            cov_consensus_len=float(len(motifs[0])), cov_reciprocal_overlap=1.0))
    return pl.DataFrame(rows)


def trcompdb_lookup(bed: pl.DataFrame, table: Path | None) -> pl.DataFrame:
    """Best-overlapping TRCompDB record per repeat (reciprocal overlap >= 0.5)."""
    out = bed.select("locus_id")
    if table is None or not Path(table).exists():
        return out.with_columns([pl.lit(None, dtype=pl.Float64).alias(c) for c in TRCOMPDB_COLS.values()])
    t = pl.read_parquet(table)
    hits = {}
    for (chrom,), g in bed.group_by(["chrom"]):
        tc = t.filter(pl.col("chrom") == chrom).sort("start")
        if tc.height == 0:
            continue
        ts, te = tc["start"].to_numpy(), tc["end"].to_numpy()
        for r in g.iter_rows(named=True):
            s, e = r["start"], r["end"]
            i0 = max(0, int(np.searchsorted(ts, s - 100_000)))
            i1 = int(np.searchsorted(ts, e))
            best, best_ro = None, MIN_OVERLAP
            for j in range(i0, i1):
                inter = min(e, te[j]) - max(s, ts[j])
                if inter <= 0:
                    continue
                ro = min(inter / (e - s), inter / max(1, te[j] - ts[j]))
                if ro >= best_ro:
                    best, best_ro = j, ro
            if best is not None:
                hits[r["locus_id"]] = {v: tc[k][best] for k, v in TRCOMPDB_COLS.items()}
    rows = [dict(locus_id=l, **hits.get(l, {})) for l in out["locus_id"]]
    df = pl.DataFrame(rows, schema={"locus_id": pl.Utf8, **{c: pl.Float64 for c in TRCOMPDB_COLS.values()}})
    return df


def annotate(bed: pl.DataFrame, trcompdb: Path | None = None, user_table: Path | None = None) -> pl.DataFrame:
    cov = bed_covariates(bed).join(trcompdb_lookup(bed, trcompdb), on="locus_id", how="left")
    cov = cov.with_columns(cov_is_segdup=pl.when(pl.col("cov_mr_segdup_identity").is_null()).then(None)
                           .when(pl.col("cov_mr_segdup_identity") >= 0).then(1.0).otherwise(0.0))
    if user_table:
        u = pl.read_csv(user_table, separator="\t", infer_schema_length=0)
        if "locus_id" not in u.columns:
            u = u.with_columns(locus_id=(pl.col("chrom") + ":" + (pl.col("start").cast(pl.Int64) + 1)
                                         .cast(pl.Utf8) + "-" + pl.col("end")))
        u = u.rename({c: (c if c.startswith("cov_") or c == "locus_id" else f"cov_{c}") for c in u.columns})
        keep = [c for c in u.columns if c in cov.columns and c != "locus_id"]
        u = u.select(["locus_id"] + keep).with_columns([pl.col(c).cast(pl.Float64, strict=False) for c in keep])
        cov = cov.join(u, on="locus_id", how="left", suffix="_user").with_columns(
            [pl.coalesce(f"{c}_user", c).alias(c) for c in keep]).drop([f"{c}_user" for c in keep])
    return cov
