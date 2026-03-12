"""
_io.py — BED, GTF, and alignment file I/O for vntr_genotyping.
"""

import gzip
import os
import re
import sys
from collections import defaultdict
from importlib.resources import files

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Matches standard primary chromosomes: chr1-22, chrX, chrY, chrM
_PRIMARY_CHROM_RE = re.compile(r"^chr([0-9]+|[XYM])$")

_DATA_DIR = files("vntr_genotyping") / "data"

DEFAULT_BED = str(_DATA_DIR / "v1_VNTR_Annotations_03_07_2026.bed")
DEFAULT_PSL = str(_DATA_DIR / "altSeqLiftOverPsl.txt.gz")
DEFAULT_GTF = str(_DATA_DIR / "gencode.v38.genes.gtf.gz")


# ---------------------------------------------------------------------------
# BED parsing
# ---------------------------------------------------------------------------

def parse_regions(bed_file):
    """
    Parse a BED file. Returns a list of dicts:
        {chrom, start, end, name, gene, period, alt_regions}
    BED coordinates are 0-based half-open [start, end).
    Column 5 (gene) is optional; set to None if absent or '.'.
    Column 6 (period) is optional; a positive integer repeat unit length in bp,
        or None if absent or '.'. Auto-detected: if column 6 is a positive integer
        it is read as the period and alt-contig triplets start at column 7;
        otherwise alt-contig triplets start at column 6 (legacy format).
    Columns 7-9, 10-12, ... (or 6-8, 9-11, ... for legacy) are alt-contig triplets.
    """
    regions = []
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                sys.exit(f"ERROR: Invalid BED line (need >= 3 columns): {line!r}")
            gene = parts[4].strip() if len(parts) > 4 else None
            if gene in (".", "", "NA"):
                gene = None

            # Auto-detect period column: if parts[5] is a positive integer, treat
            # it as the repeat unit length; otherwise fall through to alt-contigs.
            period = None
            idx = 5
            if len(parts) > 5:
                try:
                    p = int(parts[5])
                    if p > 0:
                        period = p
                        idx = 6   # alt-contigs start at col 7 (0-indexed: 6)
                except ValueError:
                    pass          # not an integer — old format, alt-contigs at col 6

            alt_regions = []
            while idx + 2 < len(parts):
                try:
                    alt_regions.append({
                        "chrom": parts[idx].strip(),
                        "start": int(parts[idx + 1]),
                        "end":   int(parts[idx + 2]),
                    })
                except ValueError:
                    break
                idx += 3

            regions.append({
                "chrom":       parts[0],
                "start":       int(parts[1]),
                "end":         int(parts[2]),
                "name":        parts[3] if len(parts) > 3 else f"{parts[0]}:{parts[1]}-{parts[2]}",
                "gene":        gene,
                "period":      period,
                "alt_regions": alt_regions,
            })
    if not regions:
        sys.exit(f"ERROR: No regions found in {bed_file}")
    return regions


# ---------------------------------------------------------------------------
# Region selection and merging
# ---------------------------------------------------------------------------

def load_default_regions(gene_filter=None, vntr_filter=None, verbose=True):
    """
    Load regions from the bundled default BED file.

    gene_filter: set of gene names to include (None = no filter by gene)
    vntr_filter: set of VNTR names to include (None = no filter by name)

    When both filters are provided, their union is returned (OR logic).
    When neither is provided, all regions are returned.
    """
    if not os.path.exists(DEFAULT_BED):
        sys.exit(
            f"ERROR: Built-in BED file not found at {DEFAULT_BED}. "
            "The package may be corrupted — try reinstalling with: pip install --force-reinstall vntr-genotyping"
        )

    all_regions = parse_regions(DEFAULT_BED)

    if gene_filter is None and vntr_filter is None:
        return all_regions

    selected = []
    unmatched_genes = set(gene_filter) if gene_filter else set()
    unmatched_vntrs = set(vntr_filter) if vntr_filter else set()

    for r in all_regions:
        by_gene = gene_filter is not None and r["gene"] in gene_filter
        by_vntr = vntr_filter is not None and r["name"] in vntr_filter
        if by_gene or by_vntr:
            selected.append(r)
            unmatched_genes.discard(r["gene"])
            unmatched_vntrs.discard(r["name"])

    if verbose:
        for g in sorted(unmatched_genes):
            print(f"  WARNING: --gene '{g}' did not match any region in the built-in BED.")
        for v in sorted(unmatched_vntrs):
            print(f"  WARNING: --vntr '{v}' did not match any region in the built-in BED.")

    return selected


def merge_regions(default_regions, custom_regions, verbose=True):
    """
    Merge default and custom region lists.

    If any custom region shares a name with a default region, it is renamed
    '<name>_UserSupplied' and the user is warned.
    """
    default_names = {r["name"] for r in default_regions}
    collisions = []

    renamed_custom = []
    for r in custom_regions:
        if r["name"] in default_names:
            collisions.append(r["name"])
            r = dict(r)
            r["name"] = r["name"] + "_UserSupplied"
        renamed_custom.append(r)

    if collisions and verbose:
        print(
            f"\nWARNING: {len(collisions)} VNTR name(s) appear in both the built-in BED "
            "and your custom BED. Custom entries have been renamed with '_UserSupplied':"
        )
        for name in collisions:
            print(f"  '{name}' -> '{name}_UserSupplied'")
        print()

    return default_regions + renamed_custom


def build_regions(default=False, gene_filter=None, vntr_filter=None, regions_file=None,
                  verbose=True):
    """
    Resolve the final region list from user selection parameters.

    default:      use all built-in VNTRs
    gene_filter:  set/list of gene names (subset of built-in BED)
    vntr_filter:  set/list of VNTR names (subset of built-in BED)
    regions_file: path to a custom BED file (combinable with any of the above)
    verbose:      if False, suppress progress messages
    """
    if not default and gene_filter is None and vntr_filter is None and not regions_file:
        sys.exit(
            "ERROR: No regions specified. Provide at least one of: "
            "default=True, gene_filter, vntr_filter, or regions_file."
        )

    if default and (gene_filter is not None or vntr_filter is not None):
        sys.exit(
            "ERROR: default=True cannot be combined with gene_filter or vntr_filter. "
            "Use default=True to select all built-in regions, or use "
            "gene_filter/vntr_filter to select a subset."
        )

    gf = set(gene_filter) if gene_filter else None
    vf = set(vntr_filter) if vntr_filter else None

    default_regions = []
    if default:
        default_regions = load_default_regions(verbose=verbose)
        if verbose:
            print(f"Loaded {len(default_regions)} region(s) from built-in BED (all).")
    elif gf is not None or vf is not None:
        default_regions = load_default_regions(gf, vf, verbose=verbose)
        if verbose:
            print(f"Selected {len(default_regions)} region(s) from built-in BED.")

    custom_regions = []
    if regions_file:
        custom_regions = parse_regions(regions_file)
        if verbose:
            print(f"Loaded {len(custom_regions)} region(s) from custom BED: {regions_file}")

    if not default_regions and not custom_regions:
        sys.exit(
            "ERROR: No regions were selected. "
            "Check that your gene_filter/vntr_filter values match entries in the built-in BED."
        )

    if default_regions and custom_regions:
        regions = merge_regions(default_regions, custom_regions, verbose=verbose)
        if verbose:
            print(f"Combined total: {len(regions)} region(s).")
        return regions

    return default_regions or custom_regions


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf(gtf_file):
    """
    Parse a GTF (or GTF.gz) file and extract gene-level features.
    Returns dict: {gene_name: {'chrom', 'start', 'end', 'gene_type'}}
    Coordinates are 0-based half-open (BED convention).
    Primary chromosome entries are preferred over alt/patch contigs.
    """
    genes = {}
    opener = gzip.open if gtf_file.endswith(".gz") else open

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1
            end   = int(parts[4])

            gene_name = (_gtf_attr(parts[8], "gene_name")
                         or _gtf_attr(parts[8], "gene_id"))
            gene_type = (_gtf_attr(parts[8], "gene_type")
                         or _gtf_attr(parts[8], "gene_biotype")
                         or "unknown")

            if gene_name:
                existing  = genes.get(gene_name)
                is_primary = bool(_PRIMARY_CHROM_RE.match(chrom))
                if existing is None or (
                    is_primary and not _PRIMARY_CHROM_RE.match(existing["chrom"])
                ):
                    genes[gene_name] = {
                        "chrom":     chrom,
                        "start":     start,
                        "end":       end,
                        "gene_type": gene_type,
                    }

    if not genes:
        sys.exit(
            f"ERROR: No gene-level features found in {gtf_file}. "
            "Ensure the file uses GTF format with a 'gene' feature type."
        )
    return genes


def _gtf_attr(attrs_str, key):
    """Extract the value of a single attribute from a GTF attribute string."""
    for token in attrs_str.split(";"):
        token = token.strip()
        if token.startswith(key):
            parts = token.split('"')
            if len(parts) >= 2:
                return parts[1]
    return None


# ---------------------------------------------------------------------------
# Alignment file helpers
# ---------------------------------------------------------------------------

def get_sample_name(filepath):
    """Derive a sample name from a BAM/CRAM/SAM file path."""
    bn = os.path.basename(filepath)
    for ext in (".cram", ".bam", ".sam"):
        if bn.lower().endswith(ext):
            return bn[: -len(ext)]
    return bn


def open_alignment_file(filepath, reference=None):
    """Open a SAM/BAM/CRAM file with pysam. Exits on error."""
    import pysam
    kwargs = {}
    if reference:
        kwargs["reference_filename"] = reference
    mode = "r" if filepath.lower().endswith(".sam") else "rb"
    try:
        return pysam.AlignmentFile(filepath, mode, **kwargs)
    except Exception as exc:
        sys.exit(f"ERROR: Could not open {filepath}: {exc}")
