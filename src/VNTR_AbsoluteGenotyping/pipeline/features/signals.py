"""
signals.py — per-locus signal extraction from a (possibly remote) BAM/CRAM.

    sig = extract_locus_signals(aln, locus, sample_nuisance, fasta=None)

`aln`     : open pysam.AlignmentFile (local path or remote URL)
`locus`   : dict(chrom, array_start0, array_end0, period, motifs=[...], consensus)
`sample_nuisance` : dict from features.nuisance.estimate_sample_nuisance()
`fasta`   : optional pysam.FastaFile for reference GC (else consensus-motif GC)

Returns a dict keyed by (a superset of) schema.FEATURE_COLUMNS.

One fetch + one pass over reads.  All depth work is done with numpy
difference arrays over the fetched span, stratified by MAPQ tier.
"""

from __future__ import annotations

import numpy as np

from .schema import (FLANK_WINDOWS, MAPQ_TIERS, REF_FLANK_MAPQ, REF_FLANK_W,
                     empty_row)

_PAD = max(FLANK_WINDOWS) + 6000     # fetched on each side of the array
_CLIP_MIN_BP = 10
_SPAN_MARGIN = 30                    # mate must clear the array by this many bp
_EPS = 1e-6


# ---------------------------------------------------------------------------
# motif k-mer machinery
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def build_motif_kmers(motifs: list[str], k: int) -> frozenset[str]:
    kms: set[str] = set()
    for m in motifs:
        if not m:
            continue
        for seq in (m + m, _revcomp(m + m)):          # wrap-around
            for i in range(len(seq) - k + 1):
                kms.add(seq[i:i + k])
    return frozenset(kms)


def _kmer_frac(seq: str, kmers: frozenset[str], k: int) -> float:
    if not seq or len(seq) < k:
        return 0.0
    n = len(seq) - k + 1
    hit = sum(1 for i in range(n) if seq[i:i + k] in kmers)
    return hit / n


def _best_motif(window: str, motifs: list[str]) -> tuple[int, int]:
    """(index of nearest motif, hamming distance).  window length == period."""
    best_i, best_d = -1, len(window) + 1
    for i, m in enumerate(motifs):
        if len(m) != len(window):
            continue
        d = sum(a != b for a, b in zip(window, m))
        if d < best_d:
            best_i, best_d = i, d
            if d == 0:
                break
    return best_i, best_d


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def extract_locus_signals(aln, locus: dict, sample_nuisance: dict,
                          fasta=None, irr_kmer: int = 5,
                          irr_motif_frac: float = 0.80) -> dict:
    sid = sample_nuisance.get("sample_id", "?")
    lid = locus["locus_id"]
    chrom = locus["chrom"]
    a0 = int(locus["array_start0"])
    a1 = int(locus["array_end0"])
    period = int(locus["period"]) if locus.get("period") else 0
    motifs = [m for m in (locus.get("motifs") or []) if m]
    consensus = locus.get("consensus") or (motifs[0] if motifs else "")

    if a1 <= a0:
        return empty_row(sid, lid, "degenerate_interval")

    try:
        contig_len = aln.get_reference_length(chrom)
    except (KeyError, ValueError):
        contig_len = None
    if contig_len is None:
        return empty_row(sid, lid, f"contig_absent:{chrom}")

    fs = max(0, a0 - _PAD)
    fe = min(contig_len, a1 + _PAD)
    span = fe - fs
    arr_lo, arr_hi = a0 - fs, a1 - fs        # array slice within the span frame

    # difference arrays: one per MAPQ tier
    diff = {t: np.zeros(span + 1, dtype=np.int32) for t in MAPQ_TIERS}

    ins_med = float(sample_nuisance.get("insert_median") or 350.0)
    ins_max_ok = 10.0 * ins_med

    n_fetched = 0
    array_read_mapq: list[int] = []
    irr_count = partial_count = 0
    anchored_irr = mate_unmapped = 0
    span_tlens: list[int] = []
    n_anchored_pairs = 0
    clip_left = clip_right = 0
    clip_left_bp: list[int] = []
    clip_right_bp: list[int] = []
    discordant = 0
    mc_motif_idx: list[int] = []
    mc_mismatch: list[float] = []
    mc_reads_used = 0

    kmers = build_motif_kmers(motifs or [consensus], irr_kmer) if (motifs or consensus) else frozenset()

    try:
        it = aln.fetch(chrom, fs, fe)
    except (ValueError, OSError) as exc:
        return empty_row(sid, lid, f"fetch_error:{type(exc).__name__}")

    for r in it:
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        n_fetched += 1
        rs, re = r.reference_start, r.reference_end
        if re is None:
            continue
        mq = r.mapping_quality or 0

        # ---- coverage (reference-aligned blocks) ----
        for bs, be in r.get_blocks():
            lo = max(bs, fs) - fs
            hi = min(be, fe) - fs
            if hi <= lo:
                continue
            for t in MAPQ_TIERS:
                if mq >= t:
                    diff[t][lo] += 1
                    diff[t][hi] -= 1

        overlaps_array = rs < a1 and re > a0

        # ---- soft-clip at array boundaries ----
        cig = r.cigartuples or []
        if cig:
            left_clip = cig[0][1] if cig[0][0] == 4 else 0
            right_clip = cig[-1][1] if cig[-1][0] == 4 else 0
            if left_clip >= _CLIP_MIN_BP and abs(rs - a1) <= 20:
                clip_right += 1                # clip just right of array, pointing in
                clip_right_bp.append(left_clip)
            if left_clip >= _CLIP_MIN_BP and abs(rs - a0) <= period + 20 and rs >= a0 - 5:
                clip_left += 1
                clip_left_bp.append(left_clip)
            if right_clip >= _CLIP_MIN_BP and abs(re - a0) <= 20:
                clip_left += 1
                clip_left_bp.append(right_clip)
            if right_clip >= _CLIP_MIN_BP and abs(re - a1) <= period + 20 and re <= a1 + 5:
                clip_right += 1
                clip_right_bp.append(right_clip)

        # ---- pair-level signals (count each pair once, on read1) ----
        if r.is_paired:
            mate_unm = r.mate_is_unmapped
            if overlaps_array and mate_unm:
                mate_unmapped += 1
            if not r.is_read2:                 # read1 (or unflagged) drives pair stats
                if mate_unm:
                    pass
                elif r.next_reference_id != r.reference_id:
                    if overlaps_array:
                        discordant += 1
                else:
                    m0 = r.next_reference_start
                    tlen = abs(r.template_length or 0)
                    m1 = m0 + (r.infer_read_length() or r.query_length or 150)
                    pair_lo, pair_hi = min(rs, m0), max(re, m1)
                    spans = (pair_lo <= a0 - _SPAN_MARGIN) and (pair_hi >= a1 + _SPAN_MARGIN)
                    if spans and 0 < tlen < ins_max_ok:
                        span_tlens.append(tlen)
                    elif spans and tlen >= ins_max_ok and overlaps_array:
                        discordant += 1
                    # anchored: this read in array, mate fully outside in flank
                    this_in = rs < a1 and re > a0
                    mate_in = m0 < a1 and m1 > a0
                    if this_in ^ mate_in:
                        n_anchored_pairs += 1

        if not overlaps_array:
            continue

        array_read_mapq.append(mq)
        seq = r.query_sequence or ""
        # restrict to the portion of the read inside the array for motif work
        if kmers:
            frac = _kmer_frac(seq, kmers, irr_kmer)
            if frac >= irr_motif_frac:
                irr_count += 1
                if r.is_paired and (r.mate_is_unmapped or
                                    r.next_reference_id != r.reference_id):
                    anchored_irr += 1
            elif frac >= 0.30:
                partial_count += 1

            # ---- motif composition on repeat-rich reads ----
            if frac >= 0.30 and period and len(seq) >= period and motifs:
                mc_reads_used += 1
                step = period
                for i in range(0, len(seq) - period + 1, step):
                    mi, dist = _best_motif(seq[i:i + period], motifs)
                    if mi >= 0:
                        mc_motif_idx.append(mi)
                        mc_mismatch.append(dist / period)

    if n_fetched == 0:
        return empty_row(sid, lid, "no_reads")

    # ---- assemble coverage vectors ----
    cov = {t: np.cumsum(diff[t][:span]) for t in MAPQ_TIERS}
    ref_flank_depth = _flank_depth(cov[REF_FLANK_MAPQ], arr_lo, arr_hi, REF_FLANK_W, span)
    norm = ref_flank_depth if ref_flank_depth > _EPS else np.nan

    out: dict = dict(
        sample_id=sid, locus_id=lid, chrom=chrom,
        array_start0=a0, array_end0=a1, period=period or None,
        n_reads_fetched=n_fetched, extraction_ok=True, note=None,
    )

    arr_slice = {t: cov[t][arr_lo:arr_hi] for t in MAPQ_TIERS}
    for t in MAPQ_TIERS:
        ad = float(arr_slice[t].mean()) if arr_slice[t].size else 0.0
        out[f"array_depth_d{t}"] = ad
        out[f"array_reads_d{t}"] = int(sum(1 for q in array_read_mapq if q >= t))
        for w in FLANK_WINDOWS:
            fd = _flank_depth(cov[t], arr_lo, arr_hi, w, span)
            out[f"flank_depth_d{t}_w{w}"] = fd
            out[f"depth_ratio_d{t}_w{w}"] = (ad / fd) if fd > _EPS else None

    d20 = arr_slice[20]
    out["array_depth_d20_norm"] = (float(d20.mean()) / norm) if d20.size else None
    out["mapq0_frac"] = (
        float(np.mean([q == 0 for q in array_read_mapq])) if array_read_mapq else None
    )
    out["array_depth_cv"] = (
        float(d20.std() / d20.mean()) if d20.size and d20.mean() > _EPS else None
    )
    if d20.size >= 6:
        ends = np.concatenate([d20[: max(1, d20.size // 5)], d20[-max(1, d20.size // 5):]])
        mid = d20[d20.size // 5: -d20.size // 5] if d20.size >= 10 else d20
        out["array_depth_end_mid_ratio"] = (
            float(ends.mean() / mid.mean()) if mid.mean() > _EPS else None
        )
    out["array_gc"] = _array_gc(fasta, chrom, a0, a1, consensus)
    out["flank_gc_w2000"] = _array_gc(fasta, chrom, max(0, a0 - REF_FLANK_W), a0, consensus)

    # ---- IRR ----
    n_arr = max(1, len(array_read_mapq))
    out["irr_count"] = irr_count
    out["irr_frac"] = irr_count / n_arr
    out["irr_count_norm"] = irr_count / norm if np.isfinite(norm) else None
    out["partial_repeat_count"] = partial_count
    out["partial_repeat_frac"] = partial_count / n_arr
    out["anchored_irr_count"] = anchored_irr
    out["anchored_irr_count_norm"] = anchored_irr / norm if np.isfinite(norm) else None
    out["mate_unmapped_count"] = mate_unmapped
    out["mate_unmapped_count_norm"] = mate_unmapped / norm if np.isfinite(norm) else None

    # ---- spanning pairs ----
    st = np.asarray(span_tlens, dtype=float)
    out["n_spanning"] = int(st.size)
    out["n_spanning_norm"] = (st.size / norm) if np.isfinite(norm) else None
    if st.size:
        med = float(np.median(st))
        out["spanning_tlen_median"] = med
        out["spanning_tlen_mad"] = float(np.median(np.abs(st - med)))
        out["spanning_tlen_shift"] = med - ins_med
        out["spanning_tlen_p10"] = float(np.percentile(st, 10))
        out["spanning_tlen_p90"] = float(np.percentile(st, 90))
        out["spanning_tlen_bimodality"] = _bimodality(st)
        out["spanning_len_gap"] = _largest_central_gap(st)
    out["n_anchored_pairs"] = n_anchored_pairs
    out["n_anchored_pairs_norm"] = (n_anchored_pairs / norm) if np.isfinite(norm) else None

    # ---- per-allele from spanning-pair insert-size mixture ----
    ref_span_bp = a1 - a0
    _perallele_from_spanning(out, st, ins_med, ref_span_bp, period)

    # ---- clips / discordant ----
    out["clip_left_count"] = clip_left
    out["clip_right_count"] = clip_right
    out["clip_count_norm"] = ((clip_left + clip_right) / norm) if np.isfinite(norm) else None
    out["clip_left_bp_mean"] = float(np.mean(clip_left_bp)) if clip_left_bp else None
    out["clip_right_bp_mean"] = float(np.mean(clip_right_bp)) if clip_right_bp else None
    out["discordant_count"] = discordant
    out["discordant_count_norm"] = (discordant / norm) if np.isfinite(norm) else None

    # ---- motif composition ----
    if mc_motif_idx:
        idxs = np.asarray(mc_motif_idx)
        counts = np.bincount(idxs, minlength=len(motifs)).astype(float)
        p = counts / counts.sum()
        nz = p[p > 0]
        out["mc_n_reads_used"] = mc_reads_used
        out["mc_n_motifs_observed"] = int((counts > 0).sum())
        out["mc_entropy"] = float(-(nz * np.log2(nz)).sum())
        out["mc_frac_top_motif"] = float(p.max())
        cons_i = motifs.index(consensus) if consensus in motifs else 0
        out["mc_frac_ref_consensus"] = float(p[cons_i]) if cons_i < len(p) else None
        unif = np.full_like(p, 1.0 / len(p))
        out["mc_divergence_from_ref"] = float(0.5 * np.abs(p - unif).sum())  # TV distance
        out["mc_mean_unit_mismatch"] = float(np.mean(mc_mismatch)) if mc_mismatch else None

    return out


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _flank_depth(cov: np.ndarray, arr_lo: int, arr_hi: int, w: int, span: int) -> float:
    left_lo, left_hi = max(0, arr_lo - w), max(0, arr_lo)
    right_lo, right_hi = min(span, arr_hi), min(span, arr_hi + w)
    segs = []
    if left_hi > left_lo:
        segs.append(cov[left_lo:left_hi])
    if right_hi > right_lo:
        segs.append(cov[right_lo:right_hi])
    if not segs:
        return 0.0
    cat = np.concatenate(segs)
    return float(cat.mean()) if cat.size else 0.0


def _bimodality(x: np.ndarray) -> float | None:
    """Sarle's bimodality coefficient; > ~0.555 suggests bimodal."""
    n = x.size
    if n < 8:
        return None
    m = x.mean()
    s = x.std()
    if s < _EPS:
        return 0.0
    z = (x - m) / s
    g1 = float((z ** 3).mean())
    g2 = float((z ** 4).mean()) - 3.0
    denom = g2 + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return (g1 ** 2 + 1.0) / denom if abs(denom) > _EPS else None


def _largest_central_gap(x: np.ndarray) -> float | None:
    if x.size < 6:
        return None
    xs = np.sort(x)
    lo = int(0.1 * xs.size)
    hi = int(0.9 * xs.size)
    if hi - lo < 2:
        return None
    return float(np.max(np.diff(xs[lo:hi])))


def _perallele_from_spanning(out: dict, st: np.ndarray, ins_med: float,
                             ref_span_bp: int, period: int) -> None:
    """1-D 2-component split of spanning insert sizes -> per-allele bp/units.

    Each spanning pair's insert ~= genomic_insert + (allele_bp - ref_span_bp).
    So (tlen - ins_med) estimates (allele_bp - ref_span_bp) for that haplotype.
    """
    out["pa_spanning_n"] = int(st.size)
    if st.size < 10:
        return
    shift = np.sort(st - ins_med)                    # ~ allele_bp - ref_span_bp
    # 2-means (1-D, exact via best split point)
    best_j, best_ss = None, np.inf
    csum = np.cumsum(shift)
    csq = np.cumsum(shift ** 2)
    n = shift.size
    for j in range(2, n - 2):
        m1 = csum[j - 1] / j
        m2 = (csum[-1] - csum[j - 1]) / (n - j)
        ss = (csq[j - 1] - j * m1 * m1) + ((csq[-1] - csq[j - 1]) - (n - j) * m2 * m2)
        if ss < best_ss:
            best_ss, best_j = ss, j
    if best_j is None:
        return
    g1, g2 = shift[:best_j], shift[best_j:]
    m1, m2 = float(g1.mean()), float(g2.mean())
    pooled_sd = float(np.sqrt(best_ss / n)) + 1e-6
    sep = abs(m2 - m1) / pooled_sd
    a1_bp = max(0.0, ref_span_bp + m1)
    a2_bp = max(0.0, ref_span_bp + m2)
    out["pa_bp_shift_a1"] = m1
    out["pa_bp_shift_a2"] = m2
    out["pa_units_a1"] = (a1_bp / period) if period else None
    out["pa_units_a2"] = (a2_bp / period) if period else None
    out["pa_mix_separation"] = sep
    # confidence: needs enough pairs, clear separation, both clusters populated
    bal = min(len(g1), len(g2)) / max(len(g1), len(g2))
    out["pa_confidence"] = float(
        min(1.0, st.size / 40) * min(1.0, sep / 2.0) * bal
    )


def _array_gc(fasta, chrom: str, s: int, e: int, consensus: str) -> float | None:
    if fasta is not None and e > s:
        try:
            seq = fasta.fetch(chrom, s, e).upper()
            n = sum(c in "ACGT" for c in seq)
            if n:
                return (seq.count("G") + seq.count("C")) / n
        except (KeyError, ValueError):
            pass
    if consensus:
        cu = consensus.upper()
        n = sum(c in "ACGT" for c in cu)
        if n:
            return (cu.count("G") + cu.count("C")) / n
    return None
