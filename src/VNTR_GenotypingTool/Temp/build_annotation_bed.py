#!/usr/bin/env python3
"""
build_annotation_bed.py - Build a finalized VNTR annotation BED file.

Pipeline:
  1. Load VNTR coordinates from an input file (XLSX, CSV, TSV, or BED).
  2. Annotate each VNTR with a gene name using a GENCODE/Ensembl GTF.
     Priority: 1) most bp overlap  2) protein_coding  3) alphabetical name.
  3. Optionally add alternate-contig coordinate columns using a UCSC
     altSeqLiftOverPsl file. Skip this step with --no-alt-contigs when
     the reference genome contains only primary chromosome sequences.

Usage:
    python build_annotation_bed.py \\
        --vntr-list  <vntrs.xlsx|vntrs.csv|vntrs.tsv|vntrs.bed> \\
        --gtf        <gencode.gtf.gz> \\
        --reference  GRCh38 \\
        --out        <output.bed> \\
        [--psl       <altSeqLiftOverPsl.txt.gz>] \\
        [--no-alt-contigs]

Input file formats:
    XLSX / CSV / TSV  : must have header columns  chrom (or chr), start, end
                        optional column cons_len (or period) carries the repeat unit length
    BED (headerless)  : columns 1-3 used as chrom, start, end; period set to '.'

Output BED columns:
    chrom  start  end  name  gene_name  period  [alt_chrom alt_start alt_end ...]

    period (column 6): repeat unit length in bp ('.' if unknown).

Coordinates are 0-based half-open [start, end), consistent with BED format.
"""

import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict

try:
    import pandas as pd
    _HAVE_PANDAS = True
except ImportError:
    _HAVE_PANDAS = False

_PRIMARY_CHROM_RE = re.compile(r"^chr([0-9]+|[XYM])$")


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------

def _parse_period(value):
    """Return period as a positive int, or None if missing/invalid."""
    try:
        p = int(float(str(value)))
        return p if p > 0 else None
    except (ValueError, TypeError):
        return None


def load_vntr_list(path):
    """
    Load VNTR coordinates from XLSX, CSV, TSV, or headerless BED.
    Returns a list of (chrom, start, end, period) where period is an int
    or None if the repeat unit length is not available in the source file.
    """
    low = path.lower()

    # XLSX
    if low.endswith(".xlsx") or low.endswith(".xls"):
        if not _HAVE_PANDAS:
            sys.exit("pandas is required to read XLSX files: pip install pandas openpyxl")
        df = pd.read_excel(path, dtype={"start": int, "end": int})
        cols = {c.lower(): c for c in df.columns}
        chrom_col = cols.get("chrom") or cols.get("chr")
        if not chrom_col or "start" not in cols or "end" not in cols:
            sys.exit(f"XLSX must have columns: chrom (or chr), start, end. Found: {list(df.columns)}")
        period_col = cols.get("cons_len") or cols.get("period")
        records = []
        for _, r in df.iterrows():
            period = _parse_period(r[period_col]) if period_col else None
            records.append((str(r[chrom_col]), int(r[cols["start"]]), int(r[cols["end"]]), period))
        return records

    # CSV
    if low.endswith(".csv"):
        if not _HAVE_PANDAS:
            sys.exit("pandas is required to read CSV files: pip install pandas")
        df = pd.read_csv(path)
        cols = {c.lower(): c for c in df.columns}
        chrom_col = cols.get("chrom") or cols.get("chr")
        if not chrom_col or "start" not in cols or "end" not in cols:
            sys.exit(f"CSV must have columns: chrom (or chr), start, end. Found: {list(df.columns)}")
        period_col = cols.get("cons_len") or cols.get("period")
        records = []
        for _, r in df.iterrows():
            period = _parse_period(r[period_col]) if period_col else None
            records.append((str(r[chrom_col]), int(r[cols["start"]]), int(r[cols["end"]]), period))
        return records

    # TSV or BED — peek at first non-comment line to detect header
    opener = gzip.open if low.endswith(".gz") else open
    with opener(path, "rt") as fh:
        lines = [l for l in fh if not l.startswith("#")]

    if not lines:
        sys.exit(f"No data rows found in {path}")

    first = lines[0].rstrip("\n").split("\t")
    # Headerless if column 2 (start) looks numeric
    has_header = not (len(first) > 1 and first[1].lstrip("-").isdigit())

    if has_header:
        header = [c.lower() for c in first]
        try:
            ci = header.index("chrom") if "chrom" in header else header.index("chr")
            si = header.index("start")
            ei = header.index("end")
        except ValueError:
            sys.exit(f"TSV header must contain chrom/chr, start, end. Found: {first}")
        pi = next((i for i, h in enumerate(header) if h in ("cons_len", "period")), None)
        records = []
        for line in lines[1:]:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ci, si, ei):
                continue
            period = _parse_period(parts[pi]) if pi is not None and pi < len(parts) else None
            records.append((parts[ci], int(parts[si]), int(parts[ei]), period))
        return records
    else:
        # Headerless BED — no period information available
        records = []
        for line in lines:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                records.append((parts[0], int(parts[1]), int(parts[2]), None))
            except ValueError:
                continue
        return records


# ---------------------------------------------------------------------------
# GTF parsing and gene annotation
# ---------------------------------------------------------------------------

def _gtf_attr(attrs_str, key):
    for token in attrs_str.split(";"):
        token = token.strip()
        if token.startswith(key):
            parts = token.split('"')
            if len(parts) >= 2:
                return parts[1]
    return None


def parse_gtf(gtf_file):
    """
    Parse GTF; return genes grouped by primary chromosome.
    Each entry: (start, end, gene_name, is_protein_coding).
    Sorted by start within each chromosome.
    """
    genes_by_chrom = defaultdict(list)
    opener = gzip.open if gtf_file.endswith(".gz") else open

    print(f"Parsing GTF: {gtf_file} ...")
    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            if not _PRIMARY_CHROM_RE.match(chrom):
                continue
            start     = int(parts[3]) - 1   # GTF 1-based → 0-based
            end       = int(parts[4])
            gene_name = _gtf_attr(parts[8], "gene_name") or _gtf_attr(parts[8], "gene_id")
            gene_type = _gtf_attr(parts[8], "gene_type") or _gtf_attr(parts[8], "gene_biotype") or ""
            is_coding = (gene_type == "protein_coding")
            if gene_name:
                genes_by_chrom[chrom].append((start, end, gene_name, is_coding))

    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda x: x[0])

    total = sum(len(v) for v in genes_by_chrom.values())
    print(f"  Loaded {total:,} gene records on primary chromosomes.")
    return genes_by_chrom


def find_overlapping_genes(chrom, vntr_start, vntr_end, genes_by_chrom, starts_by_chrom):
    genes  = genes_by_chrom.get(chrom)
    starts = starts_by_chrom.get(chrom)
    if not genes:
        return []
    idx = bisect.bisect_left(starts, vntr_end)
    overlapping = []
    for i in range(idx - 1, -1, -1):
        g_start, g_end, g_name, is_coding = genes[i]
        if g_end <= vntr_start:
            continue
        overlap_bp = min(vntr_end, g_end) - max(vntr_start, g_start)
        overlapping.append((overlap_bp, is_coding, g_name))
    return overlapping


def assign_gene(overlapping):
    if not overlapping:
        return None
    best = sorted(overlapping, key=lambda x: (-x[0], -x[1], x[2]))[0]
    return best[2]


# ---------------------------------------------------------------------------
# PSL parsing and alt-contig coordinate mapping
# ---------------------------------------------------------------------------

def parse_psl(psl_file):
    """
    Parse UCSC altSeqLiftOverPsl (22-column BIN-prepended PSL format).

    Column layout (0-indexed, BIN prepended at position 0):
      9=strand  10=qName  11=qSize  12=qStart  13=qEnd
      14=tName  15=tSize  16=tStart 17=tEnd
      18=blockCount  19=blockSizes  20=qStarts  21=tStarts

    Only records where tName matches a primary chromosome are kept.
    Returns dict { primary_chrom: [record, ...] } sorted by t_start.
    """
    records_by_chrom = defaultdict(list)
    opener = gzip.open if psl_file.endswith(".gz") else open

    with opener(psl_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 22:
                continue
            strand      = cols[9]
            alt_name    = cols[10]
            q_size      = int(cols[11])
            t_name      = cols[14]
            t_start     = int(cols[16])
            t_end       = int(cols[17])
            block_sizes = [int(x) for x in cols[19].rstrip(",").split(",") if x]
            q_starts    = [int(x) for x in cols[20].rstrip(",").split(",") if x]
            t_starts    = [int(x) for x in cols[21].rstrip(",").split(",") if x]

            if not _PRIMARY_CHROM_RE.match(t_name):
                continue

            blocks = list(zip(t_starts, q_starts, block_sizes))
            records_by_chrom[t_name].append({
                "alt_name": alt_name,
                "q_size":   q_size,
                "strand":   strand,
                "t_start":  t_start,
                "t_end":    t_end,
                "blocks":   blocks,
            })

    for chrom in records_by_chrom:
        records_by_chrom[chrom].sort(key=lambda r: r["t_start"])

    return dict(records_by_chrom)


def _map_vntr_to_alt(vntr_start, vntr_end, record):
    """
    Map VNTR [vntr_start, vntr_end) to alt scaffold coordinates via PSL blocks.

    Case 1 — Aligned blocks:
        For each PSL block that overlaps the VNTR, map the covered portion.
        + strand: alt_pos = q_blk + (primary_pos - t_blk)
        - strand: alt_pos = q_size - q_blk - (primary_pos - t_blk)
                  (returns forward-strand coordinates of the alt scaffold)

    Case 2 — Alt-scaffold insertions at VNTR boundaries:
        When consecutive blocks are adjacent on the primary (no primary gap)
        but have a gap in the alt scaffold, and that junction falls within
        [vntr_start, vntr_end], the gap contains extra sequence on the alt
        (e.g. additional VNTR repeat copies). Both boundaries of the gap are
        included so the returned interval spans the full alt VNTR extent.

    Returns (alt_start, alt_end) as 0-based half-open, or None if no blocks
    overlap the VNTR. The interval is the minimal span covering all blocks
    and any insertion gaps touching the VNTR.
    """
    q_size = record["q_size"]
    strand = record["strand"]
    blocks = record["blocks"]

    alt_starts = []
    alt_ends   = []

    for i, (t_blk, q_blk, size) in enumerate(blocks):
        blk_t_end = t_blk + size

        # Case 1: block overlaps the VNTR
        if blk_t_end > vntr_start and t_blk < vntr_end:
            dp_start = max(vntr_start, t_blk)   - t_blk
            dp_end   = min(vntr_end,   blk_t_end) - t_blk
            if strand == "+":
                alt_starts.append(q_blk + dp_start)
                alt_ends.append(q_blk   + dp_end)
            else:
                alt_starts.append(q_size - q_blk - dp_end)
                alt_ends.append(q_size   - q_blk - dp_start)

        # Case 2: alt insertion at a VNTR boundary
        if i + 1 < len(blocks):
            next_t_blk, next_q_blk, _ = blocks[i + 1]
            if (vntr_start <= blk_t_end <= vntr_end   # junction within VNTR
                    and next_t_blk == blk_t_end        # adjacent on primary
                    and next_q_blk > q_blk + size):    # insertion in alt
                if strand == "+":
                    alt_starts.append(q_blk + size)
                    alt_ends.append(next_q_blk)
                else:
                    alt_starts.append(q_size - next_q_blk)
                    alt_ends.append(q_size - (q_blk + size))

    if not alt_starts:
        return None
    return min(alt_starts), max(alt_ends)


def get_alt_coords(vntr_chrom, vntr_start, vntr_end, records_by_chrom):
    """
    Find all alt scaffolds with PSL alignment overlapping [vntr_start, vntr_end)
    on vntr_chrom. Returns sorted list of (alt_name, alt_start, alt_end).
    """
    records = records_by_chrom.get(vntr_chrom, [])
    if not records:
        return []

    results = {}
    for rec in records:
        if rec["t_end"] <= vntr_start or rec["t_start"] >= vntr_end:
            continue
        coords = _map_vntr_to_alt(vntr_start, vntr_end, rec)
        if coords is None:
            continue
        alt_name = rec["alt_name"]
        if alt_name in results:
            existing = results[alt_name]
            results[alt_name] = (min(existing[0], coords[0]),
                                 max(existing[1], coords[1]))
        else:
            results[alt_name] = coords

    return sorted(results.items())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vntr-list", required=True, metavar="FILE",
                        help="VNTR coordinate list (XLSX, CSV, TSV, or BED)")
    parser.add_argument("--gtf", required=True, metavar="GTF",
                        help="GENCODE/Ensembl GTF or GTF.gz for gene annotation")
    parser.add_argument("--reference", required=True, metavar="LABEL",
                        help="Reference genome label used in output header (e.g. GRCh38, CHM13)")
    parser.add_argument("--out", required=True, metavar="BED",
                        help="Output BED file path")
    parser.add_argument("--psl", metavar="PSL",
                        help="UCSC altSeqLiftOverPsl.txt.gz (required unless --no-alt-contigs)")
    parser.add_argument("--no-alt-contigs", action="store_true",
                        help="Skip alternate-contig mapping; use when the reference "
                             "contains only primary chromosome sequences")
    args = parser.parse_args()

    if not args.no_alt_contigs and not args.psl:
        parser.error("--psl is required unless --no-alt-contigs is set")

    # ------------------------------------------------------------------
    # Load VNTRs
    # ------------------------------------------------------------------
    print(f"\nLoading VNTR list: {args.vntr_list} ...")
    vntr_coords = load_vntr_list(args.vntr_list)
    print(f"  {len(vntr_coords):,} VNTRs loaded.")

    # ------------------------------------------------------------------
    # Parse GTF
    # ------------------------------------------------------------------
    print()
    genes_by_chrom = parse_gtf(args.gtf)
    starts_by_chrom = {
        chrom: [g[0] for g in genes]
        for chrom, genes in genes_by_chrom.items()
    }

    # ------------------------------------------------------------------
    # Parse PSL (if alternate-contig mapping is enabled)
    # ------------------------------------------------------------------
    records_by_chrom = {}
    if not args.no_alt_contigs:
        print(f"\nParsing PSL: {args.psl} ...")
        records_by_chrom = parse_psl(args.psl)
        n_psl = sum(len(v) for v in records_by_chrom.values())
        print(f"  {n_psl} alignment records across {len(records_by_chrom)} primary chromosomes.")

    # ------------------------------------------------------------------
    # Annotate each VNTR
    # ------------------------------------------------------------------
    print("\nAnnotating VNTRs ...")
    n_no_overlap     = 0
    n_single_overlap = 0
    n_multi_overlap  = 0
    n_with_alt       = 0
    n_no_alt         = 0
    alt_count_dist   = defaultdict(int)

    bed_rows = []
    for chrom, start, end, period in vntr_coords:
        name = f"{chrom}:{start}-{end}"

        # Gene annotation
        overlapping = find_overlapping_genes(chrom, start, end,
                                             genes_by_chrom, starts_by_chrom)
        n_hits = len(overlapping)
        if n_hits == 0:
            n_no_overlap += 1
            gene = "."
        elif n_hits == 1:
            n_single_overlap += 1
            gene = overlapping[0][2]
        else:
            n_multi_overlap += 1
            gene = assign_gene(overlapping)

        # Alternate-contig coordinates
        alt_cols = []
        if not args.no_alt_contigs:
            alts = get_alt_coords(chrom, start, end, records_by_chrom)
            if alts:
                n_with_alt += 1
                alt_count_dist[len(alts)] += 1
                for alt_name, (alt_start, alt_end) in alts:
                    alt_cols.extend([alt_name, str(alt_start), str(alt_end)])
            else:
                n_no_alt += 1
                alt_count_dist[0] += 1

        bed_rows.append((chrom, start, end, name, gene, period, alt_cols))

    # ------------------------------------------------------------------
    # Write output BED
    # ------------------------------------------------------------------
    n_with_period = sum(1 for *_, p, _ in bed_rows if p is not None)

    if args.no_alt_contigs:
        alt_header = "# Alternate contig columns not included (--no-alt-contigs)."
    else:
        alt_header = (
            "# Alternate contig columns (7-9, 10-12, ...) list equivalent regions on\n"
            "# alternate reference sequences. Reads from all contigs are unioned before counting."
        )

    header_lines = [
        f"# Copy-number-variable VNTR regions - {args.reference}",
        f"# Source: {args.vntr_list}",
        "# Format: chrom\tstart\tend\tname\tgene_name\tperiod\t[alt_chrom\talt_start\talt_end ...]",
        "# BED coordinates are 0-based half-open [start, end)",
        "#",
        "# gene_name (column 5): gene with greatest bp overlap with the VNTR.",
        "#   Priority: 1) most bp overlap  2) protein_coding  3) alphabetical name",
        "#   '.' = no overlapping gene found on primary chromosome.",
        "#",
        "# period (column 6): repeat unit length in bp ('.' if unknown).",
        "#   Used to convert density ratios to predicted copy numbers:",
        "#   predicted_copies = density_ratio * 2 * (VNTR_length / period)",
        "#",
        alt_header,
    ]

    with open(args.out, "w") as fh:
        for line in header_lines:
            fh.write(line + "\n")
        for chrom, start, end, name, gene, period, alt_cols in bed_rows:
            period_str = str(period) if period is not None else "."
            parts = [chrom, str(start), str(end), name, gene, period_str] + alt_cols
            fh.write("\t".join(parts) + "\n")

    # ------------------------------------------------------------------
    # Summary report
    # ------------------------------------------------------------------
    total = len(bed_rows)
    print(f"\n{'='*55}")
    print(f"  Reference              : {args.reference}")
    print(f"  Total VNTRs written    : {total:,}")
    print(f"  With period (cons_len) : {n_with_period:,}  ({100*n_with_period/total:.1f}%)")
    print(f"  No gene overlap        : {n_no_overlap:,}  ({100*n_no_overlap/total:.1f}%)")
    print(f"  Single gene overlap    : {n_single_overlap:,}  ({100*n_single_overlap/total:.1f}%)")
    print(f"  Multi-gene overlap     : {n_multi_overlap:,}  ({100*n_multi_overlap/total:.1f}%)")
    if not args.no_alt_contigs:
        print(f"  With alt contig coords : {n_with_alt:,}  ({100*n_with_alt/total:.1f}%)")
        print(f"  Without alt coords     : {n_no_alt:,}  ({100*n_no_alt/total:.1f}%)")
        print(f"\n  Alt scaffold count per VNTR:")
        for n in sorted(alt_count_dist):
            if n == 0:
                continue
            print(f"    {n} alt scaffold(s) : {alt_count_dist[n]:,} VNTRs")
    print(f"{'='*55}")
    print(f"\nBED written to: {args.out}")


if __name__ == "__main__":
    main()
