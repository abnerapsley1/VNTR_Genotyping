"""Repeat BED files.

Columns (tab-separated, no header; lines starting with '#' are ignored):
  1 chrom    2 start (0-based)    3 end (exclusive)    4 motif length (bp)
  5 motifs   comma-separated motif sequence(s); the first is treated as the consensus
  6 name     optional

Each repeat gets locus_id = "chrom:start+1-end" (1-based, inclusive), the ID used
throughout vntrlen's outputs and model store.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

SCHEMA = {"chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64, "period": pl.Int64,
          "motifs": pl.Utf8, "name": pl.Utf8}


def locus_id(chrom: str, start0: int, end: int) -> str:
    return f"{chrom}:{start0 + 1}-{end}"


def read_bed(path: str | Path) -> pl.DataFrame:
    rows = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                raise ValueError(f"{path}:{n}: need at least 5 columns "
                                 f"(chrom, start, end, motif length, motifs); got {len(f)}")
            chrom, s, e = f[0], int(f[1]), int(f[2])
            if e <= s:
                raise ValueError(f"{path}:{n}: end <= start")
            motifs = [m.strip().upper() for m in f[4].split(",") if m.strip()]
            if not motifs or not all(set(m) <= set("ACGTN") for m in motifs):
                raise ValueError(f"{path}:{n}: column 5 must list motif sequences (A/C/G/T)")
            period = int(float(f[3])) if f[3] not in ("", ".") else len(motifs[0])
            rows.append((chrom, s, e, period, ",".join(motifs), f[5] if len(f) > 5 else None))
    df = pl.DataFrame(rows, schema=SCHEMA, orient="row")
    df = df.with_columns(locus_id=(pl.col("chrom") + ":" + (pl.col("start") + 1).cast(pl.Utf8)
                                   + "-" + pl.col("end").cast(pl.Utf8)))
    dup = df.filter(pl.col("locus_id").is_duplicated())
    if dup.height:
        raise ValueError(f"{path}: duplicate repeats, e.g. {dup['locus_id'][0]}")
    return df


def to_targets(bed: pl.DataFrame) -> list[dict]:
    """Rows in the form the signal extractor expects."""
    return [dict(locus_id=r["locus_id"], chrom=r["chrom"], array_start0=r["start"],
                 array_end0=r["end"], period=r["period"], motifs=r["motifs"].split(","),
                 consensus=r["motifs"].split(",")[0])
            for r in bed.iter_rows(named=True)]
