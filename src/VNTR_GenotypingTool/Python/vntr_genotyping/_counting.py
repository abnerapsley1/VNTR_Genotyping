"""
_counting.py — Read counting, normalization, and the primary count_vntrs() API.
"""

import os
import re
import sys
from collections import defaultdict

import pandas as pd

from ._io import (
    DEFAULT_GTF,
    DEFAULT_PSL,
    build_regions,
    get_sample_name,
    open_alignment_file,
    parse_gtf,
)
from ._mapping import get_alt_coords_for_region, parse_psl


# ---------------------------------------------------------------------------
# Coordinate-format normalization window helper
# ---------------------------------------------------------------------------

# Matches "chrN:start-end" used as an explicit normalization window in BED col 5
_COORD_RE = re.compile(r"^(chr[^\s:]+):([0-9]+)-([0-9]+)$")


def _parse_coord_norm(s):
    """
    If s matches 'chrN:start-end', return (chrom, start, end); else return None.
    """
    if s is None:
        return None
    m = _COORD_RE.match(s)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


# ---------------------------------------------------------------------------
# Nearest-gene assignment
# ---------------------------------------------------------------------------

def assign_nearest_genes(regions, gene_dict, allowed_types=None, verbose=True):
    """
    For each region with gene=None, find the nearest eligible gene on the same
    chromosome and assign it in-place.

    allowed_types : set of gene_type strings eligible for normalization,
                    e.g. {"protein_coding"} (default) or {"protein_coding","lncRNA"}.
                    Pass {"any"} to allow all gene types.

    Tie-breaking (distance == 0, multiple overlapping genes):
      1. Most base-pair overlap
      2. protein_coding preferred
      3. Alphabetical gene name (deterministic final tiebreaker)
    """
    if allowed_types is None:
        allowed_types = {"protein_coding"}
    use_any = "any" in allowed_types

    for r in regions:
        if r["gene"] is not None:
            continue

        chrom = r["chrom"]
        candidates = []   # (dist, -overlap_bp, not_protein_coding, name)

        for name, g in gene_dict.items():
            if g["chrom"] != chrom:
                continue
            if not use_any and g.get("gene_type") not in allowed_types:
                continue

            if r["start"] < g["end"] and r["end"] > g["start"]:
                overlap_bp = min(r["end"], g["end"]) - max(r["start"], g["start"])
                dist = 0
            else:
                overlap_bp = 0
                dist = min(abs(r["start"] - g["end"]), abs(r["end"] - g["start"]))

            is_coding = (g.get("gene_type") == "protein_coding")
            candidates.append((dist, -overlap_bp, not is_coding, name))

        if not candidates:
            if verbose:
                print(f"  WARNING: {r['name']}: no eligible gene found on {chrom}. "
                      "Normalization will be skipped (output: NA).")
            continue

        candidates.sort()
        best = candidates[0]

        ties = [c for c in candidates
                if c[0] == best[0] and c[1] == best[1] and c[2] == best[2]
                and c[3] != best[3]]
        if ties and verbose:
            tied_names = [best[3]] + [c[3] for c in ties]
            print(f"  WARNING: {r['name']}: normalization gene tied between "
                  f"{tied_names}. Using '{best[3]}' (alphabetical first).")

        dist_bp   = best[0]
        gene_name = best[3]
        r["gene"] = gene_name
        if verbose:
            label = "overlapping" if dist_bp == 0 else f"nearest eligible (distance {dist_bp:,} bp)"
            print(f"  {r['name']}: no gene in column 5 — assigned {label} gene '{gene_name}'")


# ---------------------------------------------------------------------------
# Read counting
# ---------------------------------------------------------------------------

def get_read_names(aln, chrom, start, end, alt_regions=None):
    """
    Return the set of unique, mapped QNAMEs overlapping [start, end) on chrom,
    unioned with reads from any alternate contig regions in alt_regions.

    alt_regions: list of dicts with keys 'chrom', 'start', 'end', or None.
    Returns an empty set if all regions are absent from the index.
    """
    names = set()
    for region in [{"chrom": chrom, "start": start, "end": end}] + (alt_regions or []):
        try:
            for read in aln.fetch(region["chrom"], region["start"], region["end"]):
                if not read.is_unmapped:
                    names.add(read.query_name)
        except ValueError:
            pass
    return names


# ---------------------------------------------------------------------------
# Density ratio metric
# ---------------------------------------------------------------------------

def density_ratio(vntr_reads, vntr_len, gene_reads, gene_len):
    """
    Relative read density of a VNTR vs. its host gene (excluding VNTRs).

        ratio = (vntr_reads / vntr_len) / (gene_reads / gene_len)

    Returns None if any denominator is zero.
    1.0 = average coverage; >1 = enriched; <1 = depleted.
    """
    if gene_reads == 0 or vntr_len == 0 or gene_len == 0:
        return None
    return (vntr_reads * gene_len) / (gene_reads * vntr_len)


# ---------------------------------------------------------------------------
# Primary public API
# ---------------------------------------------------------------------------

def count_vntrs(
    input_files,
    *,
    default=False,
    gene=None,
    vntr=None,
    regions=None,
    gtf=DEFAULT_GTF,
    reference=None,
    psl=DEFAULT_PSL,
    no_alt_contigs=False,
    norm_gene_types=("protein_coding",),
    norm_method="gene",
    norm_window=5000,
    output_csv=None,
    verbose=None,
):
    """
    Count reads overlapping VNTR regions from one or more BAM/CRAM files.

    Parameters
    ----------
    input_files : str or list of str
        Path(s) to SAM/BAM/CRAM alignment files.
    default : bool
        Use all VNTRs from the built-in BED file.
    gene : list of str, optional
        Subset built-in VNTRs to those assigned to these gene(s).
    vntr : list of str, optional
        Subset built-in VNTRs to these VNTR name(s).
    regions : str, optional
        Path to a custom BED file (combinable with default/gene/vntr).
    gtf : str or None, optional
        Path to a GTF or GTF.gz annotation file. Defaults to the bundled
        GENCODE v38 gene annotation. Used only when norm_method='gene'.
        Pass gtf=None to disable normalization and return raw read counts.
        Ignored when norm_method='local'.
    reference : str, optional
        Path to a reference FASTA. Required for CRAM files without an
        embedded reference sequence.
    psl : str, optional
        Path to UCSC altSeqLiftOverPsl.txt(.gz). Used to include alt-contig
        reads in gene background normalization. Defaults to the bundled
        GRCh38 file. Only used when norm_method='gene'.
    no_alt_contigs : bool
        Disable alt-contig normalization. When norm_method='gene', disables
        PSL-based alt-contig expansion of gene regions. When norm_method='local',
        disables extending each VNTR's alt-contig coordinates by norm_window.
    norm_gene_types : list of str
        Gene biotype(s) eligible as the normalization gene when a VNTR has
        no overlapping gene. Default: ["protein_coding"]. Pass ["any"] to
        allow all biotypes. Only used when norm_method='gene'.
    norm_method : str
        Normalization strategy. One of:
          'gene'  (default) — normalize against the host gene region, excluding
                  VNTR reads from the gene background. Requires a GTF.
          'local' — normalize against a fixed-size flanking window around each
                  VNTR (norm_window bp on each side). No GTF required.
    norm_window : int
        Flank size in bp used for local normalization (default: 5000).
        The background window spans [VNTR_start - norm_window, VNTR_end + norm_window]
        with VNTR reads excluded. Effective background length = 2 * norm_window.
        Only used when norm_method='local'.
    output_csv : str, optional
        If provided, write results to this CSV file path in addition to
        returning the DataFrame.
    verbose : bool or None, optional
        If True, print progress messages to stdout. If False, suppress all
        output. If None (default), automatically set to True when output_csv
        is None (interactive use) and False when output_csv is provided.

    Returns
    -------
    pandas.DataFrame
        Wide-format DataFrame with one row per sample.
        Columns: sample, metric, <vntr_1>, <vntr_2>, ...
        metric = 'predicted_copies' when normalization is enabled and at least
            one VNTR has a period (repeat unit length) in the BED file.
            predicted_copies = density_ratio * 2 * (VNTR_length / period).
            VNTRs without a period produce NA in this mode.
        metric = 'density_ratio' when normalization is enabled but no VNTRs
            have a period.
        metric = 'read_count' when normalization is disabled.
    """
    if isinstance(input_files, str):
        input_files = [input_files]

    # When verbose is not explicitly set, suppress output if writing to a file.
    if verbose is None:
        verbose = (output_csv is None)

    # ------------------------------------------------------------------
    # Build region list
    # ------------------------------------------------------------------
    region_list = build_regions(
        default=default,
        gene_filter=gene,
        vntr_filter=vntr,
        regions_file=regions,
        verbose=verbose,
    )

    # ------------------------------------------------------------------
    # Validate norm_method
    # ------------------------------------------------------------------
    if norm_method not in ("gene", "local"):
        sys.exit(f"ERROR: norm_method must be 'gene' or 'local', got '{norm_method}'.")

    # ------------------------------------------------------------------
    # Normalization setup
    # ------------------------------------------------------------------
    normalize        = (norm_method == "local") or (gtf is not None)
    gene_dict        = {}
    gene_to_vntrs    = defaultdict(list)
    gene_eff_len     = {}
    gene_alt_regions = {}
    gene_excl_vntrs  = {}   # all VNTRs overlapping each norm region (by coords)

    if norm_method == "gene" and normalize:
        if verbose:
            if gtf == DEFAULT_GTF:
                print(f"\nUsing bundled GENCODE v38 gene annotation for normalization.")
            print(f"\nParsing GTF: {gtf} ...")
        gene_dict = parse_gtf(gtf)
        if verbose:
            print(f"  Loaded {len(gene_dict):,} genes.")

        allowed_types = set(norm_gene_types)
        assign_nearest_genes(region_list, gene_dict, allowed_types, verbose=verbose)

        norm_regions = {}
        for r in region_list:
            key = r["gene"]
            if key is None or key in norm_regions:
                continue
            coord = _parse_coord_norm(key)
            if coord:
                norm_regions[key] = {
                    "chrom": coord[0], "start": coord[1], "end": coord[2],
                    "gene_type": None,
                }
            elif key in gene_dict:
                norm_regions[key] = gene_dict[key]
            elif verbose:
                print(f"  WARNING: Normalization target '{key}' not found in GTF "
                      "and is not a valid coordinate string. "
                      "VNTRs assigned to this target will output NA.")

        for r in region_list:
            if r["gene"] and r["gene"] in norm_regions:
                gene_to_vntrs[r["gene"]].append(r)

        for norm_key in gene_to_vntrs:
            g = norm_regions[norm_key]
            # All VNTRs (built-in or custom) whose primary coords overlap this region
            overlapping = [
                r for r in region_list
                if r["chrom"] == g["chrom"]
                and r["start"] < g["end"]
                and r["end"]   > g["start"]
            ]
            gene_excl_vntrs[norm_key] = overlapping
            excl_bp = sum(
                min(r["end"], g["end"]) - max(r["start"], g["start"])
                for r in overlapping
            )
            gene_eff_len[norm_key] = max((g["end"] - g["start"]) - excl_bp, 1)

        use_alt  = not no_alt_contigs
        psl_path = psl if use_alt else None

        if use_alt and psl_path and os.path.exists(psl_path):
            if verbose:
                print(f"\nParsing PSL for alt-contig gene normalization: {psl_path} ...")
            psl_records = parse_psl(psl_path)
            if verbose:
                n_psl = sum(len(v) for v in psl_records.values())
                print(f"  {n_psl} alignment records across {len(psl_records)} primary chromosomes.")
            for norm_key in gene_to_vntrs:
                g    = norm_regions[norm_key]
                alts = get_alt_coords_for_region(
                    g["chrom"], g["start"], g["end"], psl_records
                )
                if alts:
                    gene_alt_regions[norm_key] = alts
            if verbose:
                print(f"  {len(gene_alt_regions)} of {len(gene_to_vntrs)} normalization "
                      f"region(s) have alt-contig coordinates.")
        elif use_alt and psl_path and not os.path.exists(psl_path):
            if verbose:
                print(
                    f"\nWARNING: PSL file not found at {psl_path}. "
                    "Alt-contig gene normalization will be skipped. "
                    "Pass no_alt_contigs=True to suppress this warning, or "
                    "supply the correct path with psl=<path>."
                )

    elif norm_method == "local":
        use_alt = not no_alt_contigs
        if verbose:
            print(f"\nUsing local {norm_window:,} bp window normalization "
                  f"({'primary + alt contigs' if use_alt else 'primary chromosome only'}).")

    # ------------------------------------------------------------------
    # Per-sample processing
    # ------------------------------------------------------------------
    vntr_names = [r["name"] for r in region_list]

    # Determine output metric. When a GTF is provided and at least one VNTR
    # has a period, convert density ratios to predicted copy numbers.
    has_period = normalize and any(r.get("period") for r in region_list)
    if has_period:
        metric_label = "predicted_copies"
        if verbose:
            n_with_period = sum(1 for r in region_list if r.get("period"))
            n_without = len(region_list) - n_with_period
            print(f"\nPeriod (repeat unit length) available for {n_with_period} of "
                  f"{len(region_list)} VNTRs — outputting predicted copy numbers.")
            if n_without:
                print(f"  {n_without} VNTR(s) without period will output NA.")
    elif normalize:
        metric_label = "density_ratio"
    else:
        metric_label = "read_count"

    rows = []

    for input_file in input_files:
        sample = get_sample_name(input_file)
        if verbose:
            print(f"\nProcessing {sample} ...")

        with open_alignment_file(input_file, reference) as aln:

            if normalize:
                vntr_read_cache = {
                    r["name"]: get_read_names(
                        aln, r["chrom"], r["start"], r["end"], r["alt_regions"]
                    )
                    for r in region_list
                }

                if norm_method == "gene":
                    gene_excl_cache = {}
                    for norm_key in gene_to_vntrs:
                        g          = norm_regions[norm_key]
                        gene_reads = get_read_names(
                            aln, g["chrom"], g["start"], g["end"],
                            gene_alt_regions.get(norm_key),
                        )
                        # Exclude ALL VNTRs overlapping this region by coordinates,
                        # not just those whose gene assignment matches norm_key.
                        # This ensures custom-BED VNTRs never inflate gene backgrounds.
                        excl_reads = set()
                        for v in gene_excl_vntrs.get(norm_key, []):
                            excl_reads |= vntr_read_cache[v["name"]]
                        gene_excl_cache[norm_key] = len(gene_reads - excl_reads)

                values = []
                for r in region_list:
                    vntr_reads = vntr_read_cache[r["name"]]
                    vntr_count = len(vntr_reads)
                    vntr_len   = r["end"] - r["start"]

                    if norm_method == "gene":
                        gene_name = r["gene"]
                        if gene_name and gene_name in gene_excl_cache:
                            ratio = density_ratio(
                                vntr_count, vntr_len,
                                gene_excl_cache[gene_name],
                                gene_eff_len[gene_name],
                            )
                        else:
                            ratio = None

                    else:  # local
                        win_start  = max(0, r["start"] - norm_window)
                        win_end    = r["end"] + norm_window

                        # All other VNTRs (built-in or custom) whose primary coords
                        # overlap this window — exclude their reads from the background.
                        win_neighbors = [
                            other for other in region_list
                            if other["name"] != r["name"]
                            and other["chrom"] == r["chrom"]
                            and other["start"] < win_end
                            and other["end"]   > win_start
                        ]
                        excl_reads = set(vntr_reads)
                        for other in win_neighbors:
                            excl_reads |= vntr_read_cache[other["name"]]
                        # Subtract neighbor bp from effective window length
                        neighbor_bp = sum(
                            min(other["end"], win_end) - max(other["start"], win_start)
                            for other in win_neighbors
                        )

                        win_alts   = None
                        if use_alt and r["alt_regions"]:
                            win_alts = [
                                {
                                    "chrom": a["chrom"],
                                    "start": max(0, a["start"] - norm_window),
                                    "end":   a["end"] + norm_window,
                                }
                                for a in r["alt_regions"]
                            ]
                        window_reads    = get_read_names(aln, r["chrom"], win_start, win_end, win_alts)
                        window_excl     = len(window_reads - excl_reads)
                        window_eff_len  = max(2 * norm_window - neighbor_bp, 1)
                        ratio = density_ratio(vntr_count, vntr_len, window_excl, window_eff_len)

                    if has_period:
                        period = r.get("period")
                        if ratio is not None and period:
                            n_ref = vntr_len / period
                            val   = ratio * 2.0 * n_ref
                        else:
                            val = None
                        display = f"{val:.2f} copies" if val is not None else "NA"
                    else:
                        val = ratio
                        display = f"{val:.6f}" if val is not None else "NA"

                    if verbose:
                        print(f"  {r['name']}: {vntr_count} VNTR reads | {metric_label} = {display}")
                    values.append(val)

            else:
                values = []
                for r in region_list:
                    count = len(get_read_names(
                        aln, r["chrom"], r["start"], r["end"], r["alt_regions"]
                    ))
                    if verbose:
                        print(f"  {r['name']}: {count} reads")
                    values.append(count)

        rows.append({"sample": sample, "metric": metric_label, **dict(zip(vntr_names, values))})

    df = pd.DataFrame(rows, columns=["sample", "metric"] + vntr_names)

    if output_csv is not None:
        df.to_csv(output_csv, index=False)
        if verbose:
            print(f"\nResults written to {output_csv}")

    return df
