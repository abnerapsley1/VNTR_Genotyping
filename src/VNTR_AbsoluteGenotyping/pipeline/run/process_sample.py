"""
process_sample.py — extract the full feature table for ONE sample.

    python -m pipeline.run.process_sample --sample HG00733
    python -m pipeline.run.process_sample --row 3          # 0-based manifest row

Steps: cache reference + index -> open remote alignment -> genome-wide nuisance
-> parallel per-locus signal extraction (32 workers, each its own htslib handle)
-> write features/<sample_id>.parquet and <sample_id>.nuisance.json.
Nothing large is written to disk; only the small index is cached.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
import pysam

from ..features.nuisance import estimate_sample_nuisance
from ..features.schema import FEATURE_COLUMNS, STRING_COLUMNS, empty_row, row_from_signals
from ..features.signals import extract_locus_signals
from ..paths import load_config
from .manifest import load_references, load_samples
from .remote import ensure_index, ensure_reference, open_alignment, resolve_contig

# per-worker globals
_ALN: pysam.AlignmentFile | None = None
_FASTA: pysam.FastaFile | None = None
_NUIS: dict = {}
_PARAMS: dict = {}


def _winit(data_url, index_path, ref_path, nuis, params):
    global _ALN, _FASTA, _NUIS, _PARAMS
    _ALN = open_alignment(data_url, Path(index_path) if index_path else None,
                          Path(ref_path) if ref_path else None)
    _FASTA = pysam.FastaFile(ref_path) if ref_path and Path(ref_path).exists() else None
    _NUIS = nuis
    _PARAMS = params


def _work(chunk: list[dict]) -> list[dict]:
    rows = []
    for lc in chunk:
        rc = lc.get("_contig")
        if rc is None:
            rows.append(empty_row(_NUIS["sample_id"], lc["locus_id"],
                                  f"contig_absent:{lc['chrom']}"))
            continue
        loc = dict(lc, chrom=rc)
        try:
            sig = extract_locus_signals(
                _ALN, loc, _NUIS, fasta=_FASTA,
                irr_kmer=_PARAMS.get("irr_kmer", 5),
                irr_motif_frac=_PARAMS.get("irr_motif_frac", 0.80),
            )
        except Exception as exc:                       # noqa: BLE001 - keep the run alive
            sig = empty_row(_NUIS["sample_id"], lc["locus_id"],
                            f"error:{type(exc).__name__}:{exc}"[:180])
        row = row_from_signals(sig)
        rows.append(row)
    return rows


def _chunks(seq, n):
    k = math.ceil(len(seq) / n)
    for i in range(0, len(seq), k):
        yield seq[i:i + k]


def process(cfg, sample: dict) -> Path:
    sid = sample["sample_id"]
    feat_dir = cfg.ensure_dir("outputs.features_dir")
    out_path = feat_dir / f"{sid}.parquet"
    if cfg.get("runner.resume", True) and out_path.exists():
        print(f"[{sid}] features exist -> skip ({out_path})")
        return out_path

    scratch = cfg.ensure_dir("outputs.scratch_dir") / sid
    scratch.mkdir(parents=True, exist_ok=True)
    refs = load_references(cfg)
    ref_key = sample.get("align_reference")
    ref_path = None
    if ref_key and ref_key in refs:
        r = refs[ref_key]
        print(f"[{sid}] caching reference {ref_key} ...")
        ref_path = str(ensure_reference(r["fasta_url"], r.get("fai_url"), scratch.parent / "_refcache"))
    elif sample["data_url"].lower().endswith(".cram"):
        raise SystemExit(f"[{sid}] CRAM needs align_reference; none in references.tsv for '{ref_key}'")

    print(f"[{sid}] caching index ...")
    index_path = str(ensure_index(sample["data_url"], sample.get("index_url"), scratch))

    t0 = time.time()
    aln = open_alignment(sample["data_url"], Path(index_path),
                         Path(ref_path) if ref_path else None)
    print(f"[{sid}] estimating genome-wide nuisance ...")
    nuis = estimate_sample_nuisance(
        aln, sid,
        n_windows=int(cfg.get("features.nuisance_n_windows", 3000)),
        window_bp=int(cfg.get("features.nuisance_window_bp", 2000)),
        seed=int(cfg.get("model.random_seed", 20260831)),
        fasta=pysam.FastaFile(ref_path) if ref_path and Path(ref_path).exists() else None,
    )
    print(f"[{sid}]   depth_mean={nuis['depth_mean']:.2f} "
          f"insert_median={nuis['insert_median']} rl={nuis['read_length_mode']} "
          f"mapq0={nuis['mapq0_frac_global']}")

    targets = pl.read_parquet(cfg.path("outputs.ground_truth_dir") / "loci_targets.parquet")
    loci = targets.to_dicts()
    # resolve contigs once
    cmap: dict[str, str | None] = {}
    for lc in loci:
        c = lc["chrom"]
        if c not in cmap:
            cmap[c] = resolve_contig(aln, c)
        lc["_contig"] = cmap[c]
    aln.close()
    n_missing = sum(v is None for v in cmap.values())
    if n_missing:
        print(f"[{sid}] WARNING: {n_missing} contig(s) not found in alignment: "
              f"{[k for k,v in cmap.items() if v is None][:10]}")

    params = dict(irr_kmer=int(cfg.get("features.irr_kmer", 5)),
                  irr_motif_frac=float(cfg.get("features.irr_motif_frac", 0.80)))
    workers = int(cfg.get("runner.workers", os.cpu_count() or 8))
    print(f"[{sid}] extracting {len(loci):,} loci on {workers} workers ...")

    rows: list[dict] = []
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_winit,
        initargs=(sample["data_url"], index_path, ref_path, nuis, params),
    ) as ex:
        done = 0
        for part in ex.map(_work, list(_chunks(loci, workers * 8))):
            rows.extend(part)
            done += len(part)
            if done % 8000 < len(part):
                print(f"[{sid}]   {done:,}/{len(loci):,} "
                      f"({time.time()-t0:.0f}s)")

    df = pl.DataFrame(rows, schema={c: (pl.Utf8 if c in STRING_COLUMNS else pl.Float64)
                                    for c in FEATURE_COLUMNS},
                      strict=False)
    df = df.with_columns(pl.lit(sid).alias("sample_id"))
    df.write_parquet(out_path)
    (feat_dir / f"{sid}.nuisance.json").write_text(json.dumps(nuis, indent=2))
    ok = df.filter(pl.col("extraction_ok") == 1.0).height
    print(f"[{sid}] wrote {out_path}  ({df.height:,} loci; {ok:,} extraction_ok; "
          f"{time.time()-t0:.0f}s)")

    # cleanup transient index (keep reference cache for other samples)
    try:
        Path(index_path).unlink(missing_ok=True)
    except OSError:
        pass
    return out_path


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--sample", help="sample_id from the manifest")
    g.add_argument("--row", type=int, help="0-based manifest row index")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    samples = load_samples(cfg)
    if args.sample:
        hit = [s for s in samples if s["sample_id"] == args.sample]
        if not hit:
            raise SystemExit(f"sample '{args.sample}' not in manifest")
        sample = hit[0]
    else:
        sample = samples[args.row]
    process(cfg, sample)


if __name__ == "__main__":
    main()
