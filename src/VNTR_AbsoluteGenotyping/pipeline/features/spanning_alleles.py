"""
spanning_alleles.py — direct per-read allele lengths at short repeats.

When a repeat allele plus anchors on both sides fits inside one read, a read
aligned across the whole array measures that haplotype's array length
exactly: reference span + inserted - deleted bases (from the CIGAR) within the
anchored window.  A person's reads then form one (same-length alleles) or two
(different alleles) clusters of lengths.

Loci: known loci whose panel allele lengths (p95 over non-held-out people) are
<= MAX_ALLELE_BP, so both alleles usually fit in a 150 bp read with anchors.
Reads: primary, not duplicate/QC-fail, MAPQ >= 1, aligned from <= start - ANCHOR
to >= end + ANCHOR.  Reads of long alleles may fail to align end-to-end
(allelic dropout); the scoring model treats read counts as evidence, not proof.

CRAMs are read from the kept local copy when present; otherwise the CRAM is
downloaded (MD5-checked) to a scratch directory on the system disk, read, and
deleted.  (Reading in place from the AWS mirror dropped connections under load
and silently lost regions.)  A person whose share of loci with any spanning
read is below MIN_COVERED is rejected and retried.

    python -m pipeline.features.spanning_alleles --fold val test [--jobs 24]
    python -m pipeline.features.spanning_alleles --fold train --local-only   # per-locus offsets
-> work/features/spanning/<sample>.parquet  (locus_id, n_reads, lengths[list])
"""

from __future__ import annotations

import argparse
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

from ..paths import load_config
from ..run.manifest import load_references, load_samples
from ..run.process_sample import _md5
from ..run.remote import ensure_index, ensure_reference, fetch_parallel, open_alignment, resolve_contig

MAX_ALLELE_BP = 120
ANCHOR = 10
MIN_COVERED = 0.95
STAGE = Path.home() / "vntr_scratch" / "span_cram"

_ALN = None


def _init(url, idx, ref):
    global _ALN
    _ALN = open_alignment(url, Path(idx), Path(ref))


def _lengths(read, s: int, e: int) -> int | None:
    """Array length on this read's haplotype, or None if it does not span [s-ANCHOR, e+ANCHOR)."""
    if read.reference_start > s - ANCHOR or read.reference_end is None or read.reference_end < e + ANCHOR:
        return None
    ins = dele = 0
    pos = read.reference_start
    for op, n in read.cigartuples:
        if op in (0, 7, 8):                      # M = X consume reference
            pos += n
        elif op == 2:                            # deletion
            lo, hi = max(pos, s - ANCHOR), min(pos + n, e + ANCHOR)
            dele += max(0, hi - lo)
            pos += n
        elif op == 3:
            return None
        elif op == 1 and s - ANCHOR <= pos <= e + ANCHOR:
            ins += n
    return (e - s) + ins - dele


def _work(chunk):
    out = []
    for lid, chrom, s, e in chunk:
        c = resolve_contig(_ALN, chrom)
        if c is None:
            continue
        L = []
        try:
            for r in _ALN.fetch(c, s - ANCHOR, e + ANCHOR):
                if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality < 1:
                    continue
                x = _lengths(r, s, e)
                if x is not None:
                    L.append(x)
        except (OSError, ValueError):
            continue                      # region unreadable: no row (counted as uncovered)
        out.append((lid, len(L), L))
    return out


def target_loci(cfg) -> pl.DataFrame:
    mdir = cfg.path("outputs.model_dir")
    from ..model.assemble import _clean_vamos_col
    smap = {s["sample_id"]: _clean_vamos_col(s["vamos_column"]) for s in load_samples(cfg) if s.get("vamos_column")}
    sp = pl.read_parquet(mdir / "sample_split.parquet")
    held = [smap[s] for s in sp.filter(pl.col("fold") != "train")["sample_id"] if s in smap]
    gt = pl.read_parquet(cfg.path("outputs.ground_truth_dir") / "gt_long.parquet",
                         columns=["locus_id", "vamos_sample", "a1_bp", "a2_bp"]).filter(~pl.col("vamos_sample").is_in(held))
    bp = pl.concat([gt.select("locus_id", b=pl.col("a1_bp")), gt.select("locus_id", b=pl.col("a2_bp"))]).drop_nulls()
    p95 = bp.group_by("locus_id").agg(p95=pl.col("b").quantile(0.95, "nearest"))
    known = pl.read_parquet(mdir / "locus_split.parquet").filter(pl.col("fold") == "train")["locus_id"]
    idx = pl.read_csv(Path(cfg.repo_root) / "src" / "VNTR_Benchmark" / "work" / "danbing" / "tr_index.tsv", separator="\t")
    return (idx.join(p95, on="locus_id").filter((pl.col("p95") <= MAX_ALLELE_BP) & pl.col("locus_id").is_in(known))
            .select("locus_id", "chrom", "array_start0", "array_end0"))


def run_sample(cfg, sample: dict, loci: pl.DataFrame, jobs: int, out: Path) -> str:
    sid = sample["sample_id"]
    t0 = time.time()
    refs = load_references(cfg)
    r = refs[sample["align_reference"]]
    cache = cfg.ensure_dir("outputs.scratch_dir")
    ref = ensure_reference(r["fasta_url"], r.get("fai_url"), cache / "_refcache")
    kept = cfg.path("runner.kept_cram_dir") / f"{sid}.cram"
    staged = None
    if kept.exists() and Path(f"{kept}.crai").exists():
        url, idx = str(kept), Path(f"{kept}.crai")
    else:
        STAGE.mkdir(parents=True, exist_ok=True)
        staged = STAGE / f"{sid}.cram"
        fetch_parallel(sample["data_url"] + ".crai", Path(f"{staged}.crai"), n_streams=1)
        fetch_parallel(sample["data_url"], staged, n_streams=8)
        if sample.get("md5") and _md5(staged) != sample["md5"]:
            staged.unlink()
            raise OSError("MD5 mismatch")
        url, idx = str(staged), Path(f"{staged}.crai")
    rows = list(loci.iter_rows())
    k = max(1, len(rows) // (jobs * 16))
    chunks = [rows[i:i + k] for i in range(0, len(rows), k)]
    res = []
    with ProcessPoolExecutor(jobs, initializer=_init, initargs=(url, str(idx), str(ref))) as ex:
        for part in ex.map(_work, chunks):
            res.extend(part)
    if staged is not None:
        staged.unlink(missing_ok=True)
        Path(f"{staged}.crai").unlink(missing_ok=True)
    df = pl.DataFrame(res, schema={"locus_id": pl.Utf8, "n_reads": pl.Int32, "lengths": pl.List(pl.Int32)}, orient="row")
    covered = float((df["n_reads"] > 0).sum()) / max(loci.height, 1)
    if covered < MIN_COVERED:
        raise OSError(f"only {covered:.0%} of loci have spanning reads; not saved")
    df.write_parquet(out)
    return f"{sid}: {df.height:,} loci, median {int(df['n_reads'].median() or 0)} spanning reads, " \
           f"{covered:.0%} covered, {'kept' if staged is None else 'downloaded'}, {time.time()-t0:.0f}s"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fold", nargs="+", default=["val", "test"])
    ap.add_argument("--jobs", type=int, default=24)
    ap.add_argument("--local-only", action="store_true", help="only people whose CRAM is kept locally")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    loci = target_loci(cfg)
    sp = pl.read_parquet(cfg.path("outputs.model_dir") / "sample_split.parquet")
    want = set(sp.filter(pl.col("fold").is_in(args.fold))["sample_id"])
    outdir = cfg.ensure_dir("outputs.features_dir") / "spanning"
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"{loci.height:,} short loci; {len(want)} people", flush=True)
    for s in load_samples(cfg):
        if s["sample_id"] not in want:
            continue
        out = outdir / f"{s['sample_id']}.parquet"
        if out.exists():
            continue
        kept = cfg.path("runner.kept_cram_dir") / f"{s['sample_id']}.cram"
        if args.local_only and not kept.exists():
            continue
        for attempt in range(3):
            try:
                print(run_sample(cfg, s, loci, args.jobs, out), flush=True)
                break
            except Exception as exc:  # noqa: BLE001 - keep going, retry remote hiccups
                print(f"{s['sample_id']}: attempt {attempt + 1} failed: {exc}", flush=True)
                time.sleep(30)


if __name__ == "__main__":
    main()
