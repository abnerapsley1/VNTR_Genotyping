"""Per-sample feature extraction (same signals the models were trained on).

For each sample: genome-wide sequencing profile (depth, fragment size, read length),
then the per-repeat read signals in parallel workers.  Single-end samples get their
read-pair signals set to missing (0 spanning pairs would wrongly mean "long array").
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
import pysam

from .features.nuisance import estimate_sample_nuisance
from .features.schema import FEATURE_COLUMNS, STRING_COLUMNS, empty_row, row_from_signals
from .features.signals import extract_locus_signals

PAIR_FEATURES = [
    "mate_unmapped_count", "mate_unmapped_count_norm", "anchored_irr_count", "anchored_irr_count_norm",
    "n_spanning", "n_spanning_norm", "spanning_tlen_median", "spanning_tlen_mad", "spanning_tlen_shift",
    "spanning_tlen_p10", "spanning_tlen_p90", "spanning_tlen_bimodality", "spanning_len_gap",
    "n_anchored_pairs", "n_anchored_pairs_norm", "discordant_count", "discordant_count_norm",
    "pa_spanning_n", "pa_bp_shift_a1", "pa_bp_shift_a2", "pa_units_a1", "pa_units_a2",
    "pa_mix_separation", "pa_confidence",
]
PAIR_NUISANCE = ["insert_median", "insert_mad", "insert_p10", "insert_p90"]

_ALN = _FASTA = None
_NUIS: dict = {}


def _open(path: str, fasta: str):
    return pysam.AlignmentFile(path, "rc" if path.endswith(".cram") else "rb",
                               reference_filename=str(fasta))


def _resolve(aln, chrom: str) -> str | None:
    refs = set(aln.references)
    for c in (chrom, chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"):
        if c in refs:
            return c
    return None


def _winit(path, fasta, nuis):
    global _ALN, _FASTA, _NUIS
    _ALN = _open(path, fasta)
    _FASTA = pysam.FastaFile(str(fasta))
    _NUIS = nuis


def _work(chunk):
    rows = []
    for lc in chunk:
        c = _resolve(_ALN, lc["chrom"])
        if c is None:
            rows.append(empty_row(_NUIS["sample_id"], lc["locus_id"], f"contig_absent:{lc['chrom']}"))
            continue
        try:
            sig = extract_locus_signals(_ALN, dict(lc, chrom=c), _NUIS, fasta=_FASTA)
        except Exception as exc:  # noqa: BLE001 - one bad locus must not stop the sample
            sig = empty_row(_NUIS["sample_id"], lc["locus_id"], f"error:{type(exc).__name__}"[:120])
        rows.append(row_from_signals(sig))
    return rows


def cache_key(alignment: str, targets: list[dict]) -> str:
    h = hashlib.sha1(str(Path(alignment).resolve()).encode())
    h.update(str(Path(alignment).stat().st_mtime_ns if Path(alignment).exists() else 0).encode())
    for t in targets:
        h.update(f"{t['locus_id']}|{','.join(t['motifs'])}".encode())
    return h.hexdigest()[:16]


def extract_sample(sample_id: str, alignment: str, fasta: str, targets: list[dict],
                   paired: bool, workdir: Path, workers: int = 8) -> pl.DataFrame:
    """Return one row per target repeat with the model's read-level features,
    plus nuis_* sample-profile columns.  Cached in workdir by input identity."""
    workdir.mkdir(parents=True, exist_ok=True)
    key = cache_key(alignment, targets)
    out = workdir / f"{sample_id}.{key}.features.parquet"
    if out.exists():
        return pl.read_parquet(out)
    t0 = time.time()
    with _open(alignment, fasta) as aln:
        nuis = estimate_sample_nuisance(aln, sample_id, fasta=pysam.FastaFile(str(fasta)))
    nuis["sample_id"] = sample_id
    if not paired:
        for k in PAIR_NUISANCE:
            nuis[k] = None
    k = max(1, math.ceil(len(targets) / (workers * 8)))
    chunks = [targets[i:i + k] for i in range(0, len(targets), k)]
    rows = []
    with ProcessPoolExecutor(workers, initializer=_winit, initargs=(alignment, fasta, nuis)) as ex:
        for part in ex.map(_work, chunks):
            rows.extend(part)
    df = pl.DataFrame(rows, schema={c: (pl.Utf8 if c in STRING_COLUMNS else pl.Float64)
                                    for c in FEATURE_COLUMNS}, strict=False)
    if not paired:
        df = df.with_columns([pl.lit(None, dtype=pl.Float64).alias(c) for c in PAIR_FEATURES])
    for name, v in nuis.items():
        if name not in ("sample_id", "gc_depth_curve"):
            df = df.with_columns(pl.lit(v, dtype=pl.Float64).alias(f"nuis_{name}"))
    df = df.with_columns(sample_id=pl.lit(sample_id), single_end=pl.lit(not paired))
    df.write_parquet(out)
    (workdir / f"{sample_id}.{key}.nuisance.json").write_text(json.dumps(nuis, indent=1, default=str))
    ok = int((df["extraction_ok"] == 1.0).sum())
    print(f"[{sample_id}] {df.height:,} repeats extracted ({ok:,} with reads) in {time.time()-t0:.0f}s",
          flush=True)
    return df
