"""
nuisance.py — genome-wide per-sample sequencing profile.

    nuis = estimate_sample_nuisance(aln, sample_id, n_windows=3000, window_bp=2000)

Samples random autosomal windows (fixed seed) and accumulates:
  depth_mean / depth_median / depth_sd        (per-base, MAPQ>=20)
  insert_median / insert_mad / insert_p10/p90 (proper pairs)
  read_length_mode
  mapq0_frac_global
  gc_depth_curve   : dict {gc_bin_center: relative_depth}  (for GC-bias correction)

These feed signal normalisation and become sample-level model features.
"""

from __future__ import annotations

import random

import numpy as np

_AUTOSOMES = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}


def estimate_sample_nuisance(aln, sample_id: str, n_windows: int = 3000,
                             window_bp: int = 2000, seed: int = 20260831,
                             fasta=None) -> dict:
    rng = random.Random(seed)
    contigs = [
        (name, aln.get_reference_length(name))
        for name in aln.references
        if name in _AUTOSOMES and (aln.get_reference_length(name) or 0) > 5 * window_bp
    ]
    if not contigs:
        contigs = [(n, aln.get_reference_length(n)) for n in aln.references[:5]]
    total = sum(L for _, L in contigs)

    depths: list[float] = []
    inserts: list[int] = []
    readlens: list[int] = []
    n_reads = 0
    n_mapq0 = 0
    gc_pairs: list[tuple[float, float]] = []   # (gc, depth) per window, if fasta

    for _ in range(n_windows):
        pick = rng.uniform(0, total)
        acc = 0
        chrom, L = contigs[-1]
        for name, length in contigs:
            if acc + length >= pick:
                chrom, L = name, length
                break
            acc += length
        start = rng.randint(0, max(1, L - window_bp - 1))
        end = start + window_bp

        cov = np.zeros(window_bp + 1, dtype=np.int32)
        try:
            it = aln.fetch(chrom, start, end)
        except (ValueError, OSError):
            continue
        for r in it:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            n_reads += 1
            if (r.mapping_quality or 0) == 0:
                n_mapq0 += 1
            if (r.mapping_quality or 0) >= 20 and r.reference_end is not None:
                for bs, be in r.get_blocks():
                    lo = max(bs, start) - start
                    hi = min(be, end) - start
                    if hi > lo:
                        cov[lo] += 1
                        cov[hi] -= 1
            rl = r.query_length or r.infer_read_length()
            if rl:
                readlens.append(rl)
            if (r.is_paired and r.is_proper_pair and not r.is_read2
                    and r.template_length):
                tl = abs(r.template_length)
                if 0 < tl < 2000:
                    inserts.append(tl)
        w = np.cumsum(cov[:window_bp])
        wd = float(w.mean())
        depths.append(wd)
        if fasta is not None:
            try:
                seq = fasta.fetch(chrom, start, end).upper()
                nb = sum(c in "ACGT" for c in seq)
                if nb:
                    gc_pairs.append(((seq.count("G") + seq.count("C")) / nb, wd))
            except (KeyError, ValueError):
                pass

    dep = np.asarray(depths, dtype=float)
    ins = np.asarray(inserts, dtype=float)
    rl = np.asarray(readlens, dtype=float)
    med_ins = float(np.median(ins)) if ins.size else None

    nuis = dict(
        sample_id=sample_id,
        depth_mean=float(dep.mean()) if dep.size else None,
        depth_median=float(np.median(dep)) if dep.size else None,
        depth_sd=float(dep.std()) if dep.size else None,
        insert_median=med_ins,
        insert_mad=float(np.median(np.abs(ins - med_ins))) if ins.size else None,
        insert_p10=float(np.percentile(ins, 10)) if ins.size else None,
        insert_p90=float(np.percentile(ins, 90)) if ins.size else None,
        read_length_mode=int(np.bincount(rl.astype(int)).argmax()) if rl.size else None,
        mapq0_frac_global=(n_mapq0 / n_reads) if n_reads else None,
        n_windows_used=int(dep.size),
        n_reads_sampled=n_reads,
    )
    nuis["gc_depth_curve"] = _gc_curve(gc_pairs, nuis["depth_mean"])
    return nuis


def _gc_curve(gc_pairs: list[tuple[float, float]], depth_mean) -> dict | None:
    if len(gc_pairs) < 100 or not depth_mean:
        return None
    arr = np.asarray(gc_pairs)
    bins = np.linspace(0.25, 0.65, 9)
    out: dict[str, float] = {}
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (arr[:, 0] >= lo) & (arr[:, 0] < hi)
        if m.sum() >= 10:
            out[f"{(lo + hi) / 2:.3f}"] = float(arr[m, 1].mean() / depth_mean)
    return out or None
