"""
_mapping.py — PSL-based alt-contig coordinate mapping for vntr_genotyping.
"""

import gzip
import re
from collections import defaultdict

_PRIMARY_CHROM_RE = re.compile(r"^chr([0-9]+|[XYM])$")


# ---------------------------------------------------------------------------
# PSL parsing
# ---------------------------------------------------------------------------

def parse_psl(psl_file):
    """
    Parse a UCSC altSeqLiftOverPsl file (22-column BIN-prepended PSL format).

    Column layout (0-indexed, BIN column at position 0):
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


# ---------------------------------------------------------------------------
# Block-level coordinate mapping
# ---------------------------------------------------------------------------

def _map_region_to_alt(region_start, region_end, record):
    """
    Map [region_start, region_end) to alt scaffold coordinates via PSL blocks.

    Case 1 — Aligned blocks: map the covered portion of each overlapping block.
    Case 2 — Alt insertions: when consecutive blocks are adjacent on the primary
              but have a gap in the alt scaffold, and that junction falls within
              [region_start, region_end], include the gap boundaries so that
              extra copies on the alt scaffold are fully captured.

    Returns (alt_start, alt_end) as 0-based half-open, or None if no overlap.
    """
    q_size = record["q_size"]
    strand = record["strand"]
    blocks = record["blocks"]

    alt_starts = []
    alt_ends   = []

    for i, (t_blk, q_blk, size) in enumerate(blocks):
        blk_t_end = t_blk + size

        # Case 1: block overlaps the region
        if blk_t_end > region_start and t_blk < region_end:
            dp_start = max(region_start, t_blk)    - t_blk
            dp_end   = min(region_end,   blk_t_end) - t_blk
            if strand == "+":
                alt_starts.append(q_blk + dp_start)
                alt_ends.append(q_blk   + dp_end)
            else:
                alt_starts.append(q_size - q_blk - dp_end)
                alt_ends.append(q_size   - q_blk - dp_start)

        # Case 2: alt insertion at a region boundary
        if i + 1 < len(blocks):
            next_t_blk, next_q_blk, _ = blocks[i + 1]
            if (region_start <= blk_t_end <= region_end   # junction within region
                    and next_t_blk == blk_t_end            # adjacent on primary
                    and next_q_blk > q_blk + size):        # insertion in alt
                if strand == "+":
                    alt_starts.append(q_blk + size)
                    alt_ends.append(next_q_blk)
                else:
                    alt_starts.append(q_size - next_q_blk)
                    alt_ends.append(q_size - (q_blk + size))

    if not alt_starts:
        return None
    return min(alt_starts), max(alt_ends)


def get_alt_coords_for_region(chrom, start, end, records_by_chrom):
    """
    Return a list of alt-contig region dicts for [start, end) on chrom,
    suitable for passing directly to get_read_names as alt_regions.
    """
    records = records_by_chrom.get(chrom, [])
    results = {}

    for rec in records:
        if rec["t_end"] <= start or rec["t_start"] >= end:
            continue
        coords = _map_region_to_alt(start, end, rec)
        if coords is None:
            continue
        alt_name = rec["alt_name"]
        if alt_name in results:
            existing = results[alt_name]
            results[alt_name] = (min(existing[0], coords[0]),
                                 max(existing[1], coords[1]))
        else:
            results[alt_name] = coords

    return [
        {"chrom": name, "start": s, "end": e}
        for name, (s, e) in sorted(results.items())
    ]
