"""
parse_vamos.py — build the copy-number ground truth from the vamos VCF.

For every locus in the bundled BED that can be matched to a vamos record,
emit per-sample haplotype lengths in BOTH motif units and bp:

    combined_units = tokens(ALTANNO[a1-1]) + tokens(ALTANNO[a2-1])
    combined_bp    = sum(len(RU[t]) for t in ALTANNO[a1-1]) + (same for a2)

vamos GT is 1-based: GT value k indexes ALTANNO[k-1] and RU is the motif list
(INFO field).  GT '.' = uncalled haplotype.  (The legacy 1000Genomes_Analysis.py
used a 0-based index here and is wrong — do not reuse its outputs.)

Outputs (parquet, under outputs.ground_truth_dir):
    gt_long.parquet          one row per (locus, sample)
    gt_locus_summary.parquet one row per locus (panel-level stats + match info)
    gt_unmatched_loci.tsv    BED loci with no vamos record (for inspection)

Usage
-----
    # Diagnose the BED<->VCF coordinate join (fast-ish; streams the VCF once):
    python -m pipeline.ground_truth.parse_vamos --report

    # Full extraction:
    python -m pipeline.ground_truth.parse_vamos --run
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from bisect import bisect_left, bisect_right
from collections import Counter
from pathlib import Path

import polars as pl

from ..paths import load_config

_SAMPLE_CLEAN_RE = re.compile(r"_hap[12]\b")


# ---------------------------------------------------------------------------
# BED loading
# ---------------------------------------------------------------------------

def load_bed(bed_path: Path) -> list[dict]:
    """Read the VNTR BED. Returns list of dicts: chrom,start,end,name,gene,period."""
    loci: list[dict] = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, start, end = f[0], int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 and f[3] not in (".", "") else f"{chrom}:{start}-{end}"
            gene = f[4] if len(f) > 4 and f[4] not in (".", "") else None
            period = None
            if len(f) > 5 and f[5] not in (".", ""):
                try:
                    period = int(float(f[5]))
                except ValueError:
                    period = None
            loci.append(
                dict(chrom=chrom, start=start, end=end,
                     locus_id=f"{chrom}:{start}-{end}", name=name,
                     gene=gene, period=period)
            )
    return loci


def _clean_sample(col: str) -> str:
    """'APR-phase1__apr001_hap1/APR-phase1__apr001_hap2' -> 'APR-phase1__apr001'."""
    left = col.split("/")[0]
    return _SAMPLE_CLEAN_RE.sub("", left).rstrip("_")


# ---------------------------------------------------------------------------
# Coordinate join helpers
# ---------------------------------------------------------------------------

# Offset hypotheses tried by --report.  Convention:  bed.start == vcf_pos + dp
# and  bed.end == vcf_end + de.  The standard 1-based-POS / 0-based-BED-start
# relationship is  bed.start == vcf_pos - 1  ->  dp = -1, de = 0.
_OFFSET_HYPOTHESES = {
    "bed.start = POS-1 , bed.end = END   (standard)": (-1, 0),
    "bed.start = POS   , bed.end = END":              (0, 0),
    "bed.start = POS-1 , bed.end = END-1":            (-1, -1),
    "bed.start = POS+1 , bed.end = END":              (1, 0),
}


def _build_chrom_interval_index(loci: list[dict]) -> dict[str, tuple[list[int], list[int], list[int]]]:
    """Per chrom: (sorted starts, matching ends, matching loci-array-indices)."""
    by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for i, lc in enumerate(loci):
        by_chrom.setdefault(lc["chrom"], []).append((lc["start"], lc["end"], i))
    out: dict[str, tuple[list[int], list[int], list[int]]] = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        out[chrom] = ([r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows])
    return out


def _reciprocal_overlap(a0: int, a1: int, b0: int, b1: int) -> float:
    inter = max(0, min(a1, b1) - max(a0, b0))
    if inter == 0:
        return 0.0
    return min(inter / (a1 - a0), inter / (b1 - b0))


def _iter_vcf_records(vcf_path: Path):
    """Yield (chrom, pos_1based, end, ru_list, altanno_list, gt_fields, raw_gt_list).

    Streams a plain-text (optionally .gz) VCF.  Header '#CHROM' line sets samples.
    """
    opener = _open_maybe_gzip(vcf_path)
    samples: list[str] = []
    with opener as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                yield ("__SAMPLES__", samples)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos = f[0], int(f[1])
            info = f[7]
            end = None
            ru: list[str] = []
            altanno: list[str] = []
            for part in info.split(";"):
                if part.startswith("END="):
                    end = int(part[4:])
                elif part.startswith("RU="):
                    ru = part[3:].split(",")
                elif part.startswith("ALTANNO="):
                    altanno = part[8:].split(",")
            if end is None:
                continue
            yield (chrom, pos, end, ru, altanno, f[9:])


def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path)


# ---------------------------------------------------------------------------
# Allele length maths
# ---------------------------------------------------------------------------

def _allele_lengths(altanno: list[str], ru: list[str]) -> tuple[list[int], list[int]]:
    """Per-allele (units, bp) for each ALTANNO entry.  1-based GT indexes these."""
    ru_len = [len(m) for m in ru]
    units: list[int] = []
    bp: list[int] = []
    for entry in altanno:
        if entry == "":
            units.append(0)
            bp.append(0)
            continue
        toks = entry.split("-")
        units.append(len(toks))
        total = 0
        for t in toks:
            try:
                ti = int(t)
            except ValueError:
                continue
            if 0 <= ti < len(ru_len):
                total += ru_len[ti]
        bp.append(total)
    return units, bp


def _parse_gt(gt_field: str) -> tuple[str | None, str | None]:
    """'3/7' -> ('3','7'); './.' -> (None,None); '3/.' -> ('3',None)."""
    core = gt_field.split(":")[0].replace("|", "/")
    if "/" not in core:
        return (None, None)
    a, b = core.split("/", 1)
    return (a if a != "." else None, b if b != "." else None)


# ---------------------------------------------------------------------------
# Report mode: test coordinate-join hypotheses
# ---------------------------------------------------------------------------

def run_report(cfg) -> None:
    bed_path = cfg.path("inputs.regions_bed", must_exist=True)
    vcf_path = cfg.path("inputs.vamos_vcf", must_exist=True)
    loci = load_bed(bed_path)
    print(f"BED loci: {len(loci):,}   VCF: {vcf_path.name}")

    # Build a set per hypothesis: (chrom, derived_start, derived_end) from VCF records.
    hyp_hits = {h: 0 for h in _OFFSET_HYPOTHESES}
    bed_keys = {(lc["chrom"], lc["start"], lc["end"]) for lc in loci}
    idx = _build_chrom_interval_index(loci)

    n_rec = 0
    overlap_only = 0
    t0 = time.time()
    for rec in _iter_vcf_records(vcf_path):
        if rec[0] == "__SAMPLES__":
            print(f"VCF samples: {len(rec[1]):,}  e.g. {rec[1][0]}")
            continue
        chrom, pos, end = rec[0], rec[1], rec[2]
        n_rec += 1
        for h, (dp, de) in _OFFSET_HYPOTHESES.items():
            if (chrom, pos + dp, end + de) in bed_keys:
                hyp_hits[h] += 1
        # overlap check (against best hypothesis-independent): does this record
        # reciprocally overlap some BED locus >= threshold without an exact hit?
        if chrom in idx and not any(
            (chrom, pos + dp, end + de) in bed_keys
            for dp, de in _OFFSET_HYPOTHESES.values()
        ):
            starts, ends, _ = idx[chrom]
            lo = bisect_left(starts, pos - 100_000)
            hi = bisect_right(starts, end + 100_000)
            for j in range(lo, hi):
                if _reciprocal_overlap(pos, end, starts[j], ends[j]) >= cfg.get(
                    "ground_truth.min_reciprocal_overlap", 0.9
                ):
                    overlap_only += 1
                    break
        if n_rec % 500_000 == 0:
            print(f"  ...{n_rec:,} records  ({time.time()-t0:.0f}s)")

    print(f"\nVCF records: {n_rec:,}")
    print("Exact-join hits by offset hypothesis (higher = correct convention):")
    for h, c in sorted(hyp_hits.items(), key=lambda kv: -kv[1]):
        print(f"  {c:>8,}  {h}")
    print(f"\nRecords that ONLY match by reciprocal overlap "
          f">= {cfg.get('ground_truth.min_reciprocal_overlap', 0.9)}: {overlap_only:,}")
    best = max(hyp_hits, key=hyp_hits.get)
    print(f"\nBest hypothesis: {best}  ->  {hyp_hits[best]:,}/{len(loci):,} BED loci "
          f"({100*hyp_hits[best]/len(loci):.1f}%)")


# ---------------------------------------------------------------------------
# Run mode: full extraction
# ---------------------------------------------------------------------------

def run_extract(cfg, offset: tuple[int, int] | None, limit: int | None) -> None:
    bed_path = cfg.path("inputs.regions_bed", must_exist=True)
    vcf_path = cfg.path("inputs.vamos_vcf", must_exist=True)
    out_dir = cfg.ensure_dir("outputs.ground_truth_dir")
    require_both = bool(cfg.get("ground_truth.require_both_haps", True))
    min_ro = float(cfg.get("ground_truth.min_reciprocal_overlap", 0.9))

    loci = load_bed(bed_path)
    if limit:
        loci = loci[:limit]
    loci_by_i = loci
    exact_lookup: dict[tuple, int] = {}
    dp, de = offset if offset is not None else _resolve_offset(cfg)
    # Convention: bed.start == vcf_pos + dp  ->  vcf_pos == bed.start - dp.
    # exact_lookup is keyed by the (chrom, vcf_pos, vcf_end) that maps to this locus.
    for i, lc in enumerate(loci):
        exact_lookup[(lc["chrom"], lc["start"] - dp, lc["end"] - de)] = i
    print(f"join offset: bed.start = POS{dp:+d}, bed.end = END{de:+d}")
    idx = _build_chrom_interval_index(loci)

    matched: dict[int, dict] = {}          # locus array-index -> match meta
    long_rows: list[tuple] = []
    samples: list[str] = []
    clean_samples: list[str] = []

    n_rec = 0
    t0 = time.time()
    for rec in _iter_vcf_records(vcf_path):
        if rec[0] == "__SAMPLES__":
            samples = rec[1]
            clean_samples = [_clean_sample(s) for s in samples]
            continue
        chrom, pos, end, ru, altanno, gts = rec
        n_rec += 1
        if n_rec % 500_000 == 0:
            print(f"  ...{n_rec:,} VCF records, {len(matched):,} loci matched "
                  f"({time.time()-t0:.0f}s)")

        li = exact_lookup.get((chrom, pos, end))
        match_type, ro = "exact", 1.0
        if li is None:
            li, ro = _best_overlap(chrom, pos, end, idx, min_ro)
            match_type = "overlap"
            if li is None:
                continue
        if li in matched and matched[li]["reciprocal_overlap"] >= ro:
            continue  # keep the better (or first exact) match

        units, bp = _allele_lengths(altanno, ru)

        def _len(idx_str: str | None) -> tuple[int | None, int | None]:
            if idx_str is None:
                return (None, None)
            k = int(idx_str) - 1  # 1-based GT -> 0-based ALTANNO
            if 0 <= k < len(units):
                return (units[k], bp[k])
            return (None, None)

        rows_this_locus = []
        allele_unit_counts: Counter = Counter()
        n_complete = n_het = 0
        for cs, gt in zip(clean_samples, gts):
            a1, a2 = _parse_gt(gt)
            u1, b1 = _len(a1)
            u2, b2 = _len(a2)
            n_called = int(u1 is not None) + int(u2 is not None)
            if u1 is not None:
                allele_unit_counts[u1] += 1
            if u2 is not None:
                allele_unit_counts[u2] += 1
            comb_u = comb_b = None
            if n_called == 2:
                comb_u, comb_b = u1 + u2, b1 + b2
                n_complete += 1
                if u1 != u2:
                    n_het += 1
            elif n_called == 1 and not require_both:
                comb_u = u1 if u1 is not None else u2
                comb_b = b1 if b1 is not None else b2
            hi_u = max(x for x in (u1, u2) if x is not None) if n_called else None
            lo_u = min(x for x in (u1, u2) if x is not None) if n_called == 2 else None
            rows_this_locus.append(
                (loci_by_i[li]["locus_id"], cs, gt.split(":")[0],
                 u1, u2, hi_u, lo_u, comb_u, b1, b2, comb_b, n_called)
            )

        total_alleles = sum(allele_unit_counts.values())
        h_exp = (1 - sum((c / total_alleles) ** 2 for c in allele_unit_counts.values())
                 ) if total_alleles else None
        h_obs = (n_het / n_complete) if n_complete else None

        # drop any earlier rows for this locus if we're replacing an overlap match
        if li in matched:
            long_rows[:] = [r for r in long_rows if r[0] != loci_by_i[li]["locus_id"]]
        long_rows.extend(rows_this_locus)
        matched[li] = dict(
            match_type=match_type, reciprocal_overlap=ro,
            vcf_pos=pos, vcf_end=end, n_alleles=len(altanno),
            h_exp=h_exp, h_obs=h_obs, n_complete=n_complete,
            ru=list(ru),
        )

    print(f"\nVCF records scanned: {n_rec:,}")
    print(f"BED loci matched:    {len(matched):,}/{len(loci):,} "
          f"({100*len(matched)/len(loci):.1f}%)")
    mt = Counter(m["match_type"] for m in matched.values())
    print(f"  by type: {dict(mt)}")

    _write_outputs(cfg, out_dir, loci, matched, long_rows, clean_samples)


def _resolve_offset(cfg) -> tuple[int, int]:
    """Coordinate join offset  (bed.start == vcf_pos + dp, bed.end == vcf_end + de).

    Empirically (parse_vamos --report) the bundled BED stores start == vamos POS
    and end == vamos END with NO offset -- the BED 'start' column is effectively
    1-based (it was generated from the vamos TSV), despite the file header.
    So the default is (0, 0).  Override with --offset if --report disagrees.
    """
    raw = cfg.get("ground_truth.join_offset")
    if raw:
        dp, de = (int(x) for x in str(raw).split(","))
        return (dp, de)
    return (0, 0)


def _best_overlap(chrom, pos, end, idx, min_ro):
    if chrom not in idx:
        return None, 0.0
    starts, ends, arr_i = idx[chrom]
    lo = bisect_left(starts, pos - 100_000)
    hi = bisect_right(starts, end + 100_000)
    best_i, best_ro = None, 0.0
    for j in range(lo, hi):
        ro = _reciprocal_overlap(pos, end, starts[j], ends[j])
        if ro >= min_ro and ro > best_ro:
            best_i, best_ro = arr_i[j], ro
    return best_i, best_ro


def _write_outputs(cfg, out_dir: Path, loci, matched, long_rows, clean_samples) -> None:
    long_schema = ["locus_id", "vamos_sample", "gt_raw",
                   "a1_units", "a2_units", "hap_hi_units", "hap_lo_units",
                   "combined_units", "a1_bp", "a2_bp", "combined_bp", "n_haps_called"]
    long_df = pl.DataFrame(long_rows, schema=long_schema, orient="row")
    long_path = out_dir / "gt_long.parquet"
    long_df.write_parquet(long_path)
    print(f"\nwrote {long_path}  ({long_df.height:,} rows)")

    # Per-locus summary (join panel stats onto every BED locus, matched or not).
    summ_rows = []
    per_locus = long_df.group_by("locus_id").agg(
        n_called=pl.col("combined_units").is_not_null().sum(),
        mean_combined_units=pl.col("combined_units").mean(),
        var_combined_units=pl.col("combined_units").var(),
        min_combined_units=pl.col("combined_units").min(),
        max_combined_units=pl.col("combined_units").max(),
        mean_combined_bp=pl.col("combined_bp").mean(),
        var_combined_bp=pl.col("combined_bp").var(),
        n_distinct_hap_units=pl.col("hap_hi_units").n_unique(),
    )
    stats = {r["locus_id"]: r for r in per_locus.to_dicts()}
    li_by_id = {lc["locus_id"]: i for i, lc in enumerate(loci)}
    for lc in loci:
        i = li_by_id[lc["locus_id"]]
        m = matched.get(i)
        s = stats.get(lc["locus_id"], {})
        ref_units = (lc["end"] - lc["start"]) / lc["period"] if lc["period"] else None
        summ_rows.append(dict(
            locus_id=lc["locus_id"], chrom=lc["chrom"], start=lc["start"],
            end=lc["end"], period=lc["period"], gene=lc["gene"],
            ref_units_est=ref_units,
            matched=m is not None,
            match_type=(m["match_type"] if m else None),
            reciprocal_overlap=(m["reciprocal_overlap"] if m else None),
            n_alleles=(m["n_alleles"] if m else None),
            h_exp=(m["h_exp"] if m else None),
            h_obs=(m["h_obs"] if m else None),
            n_called=s.get("n_called"),
            mean_combined_units=s.get("mean_combined_units"),
            var_combined_units=s.get("var_combined_units"),
            min_combined_units=s.get("min_combined_units"),
            max_combined_units=s.get("max_combined_units"),
            mean_combined_bp=s.get("mean_combined_bp"),
            var_combined_bp=s.get("var_combined_bp"),
            n_distinct_hap_units=s.get("n_distinct_hap_units"),
        ))
    summ_df = pl.DataFrame(summ_rows)
    summ_path = out_dir / "gt_locus_summary.parquet"
    summ_df.write_parquet(summ_path)
    print(f"wrote {summ_path}  ({summ_df.height:,} loci; "
          f"{summ_df['matched'].sum():,} matched)")

    unmatched = summ_df.filter(~pl.col("matched")).select(
        ["locus_id", "chrom", "start", "end", "period"]
    )
    um_path = out_dir / "gt_unmatched_loci.tsv"
    unmatched.write_csv(um_path, separator="\t")
    print(f"wrote {um_path}  ({unmatched.height:,} unmatched loci)")

    # Feature-extraction targets: 0-based array coords + motif list per locus.
    # BED 'start' == vamos POS (1-based); 0-based half-open array is [start-1, end).
    tgt_rows = []
    for lc in loci:
        i = li_by_id[lc["locus_id"]]
        m = matched.get(i)
        if not m:
            continue
        ru = m.get("ru") or []
        tgt_rows.append(dict(
            locus_id=lc["locus_id"], chrom=lc["chrom"],
            array_start0=lc["start"] - 1, array_end0=lc["end"],
            period=lc["period"], motifs=ru,
            consensus=(ru[0] if ru else None),
        ))
    tgt_df = pl.DataFrame(tgt_rows)
    tgt_path = out_dir / "loci_targets.parquet"
    tgt_df.write_parquet(tgt_path)
    print(f"wrote {tgt_path}  ({tgt_df.height:,} feature-extraction targets)")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None, help="path to config.yaml")
    ap.add_argument("--report", action="store_true",
                    help="diagnose the BED<->VCF coordinate join and exit")
    ap.add_argument("--run", action="store_true", help="run full extraction")
    ap.add_argument("--limit", type=int, default=None,
                    help="process only the first N BED loci (dev)")
    ap.add_argument("--offset", default=None,
                    help="override join offset as 'dp,de' (vcf_pos+dp==bed.start)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if args.report:
        run_report(cfg)
        return
    if args.run:
        off = None
        if args.offset:
            dp, de = (int(x) for x in args.offset.split(","))
            off = (dp, de)
        run_extract(cfg, off, args.limit)
        return
    ap.error("pass --report or --run")


if __name__ == "__main__":
    main()
