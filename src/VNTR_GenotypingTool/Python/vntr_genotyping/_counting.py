"""
_counting.py — Read counting, normalization, and the primary count_vntrs() API.
"""

import bisect
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

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

    # ------------------------------------------------------------------
    # Pre-index eligible genes by chromosome.
    # Two sorted structures per chromosome replace the original O(N_regions
    # × N_genes) nested loop with O(N_regions × K) where K is the number
    # of overlapping genes per region (typically 0–5).
    #
    #   chrom_genes  sorted by (start, ...)  → overlap detection + right neighbor
    #   chrom_ends   sorted by (end,   ...)  → left neighbor via bisect_right
    #
    # The backward scan for overlapping genes terminates early once
    # gene_start < r_start - max_gene_len_on_chrom, because no gene starting
    # that far left can have an end > r_start (by definition of max_gene_len).
    # ------------------------------------------------------------------
    chrom_genes = defaultdict(list)  # chrom → [(start, end, name, gene_type), ...]
    chrom_ends  = defaultdict(list)  # chrom → [(end,  name, gene_type), ...]

    for name, g in gene_dict.items():
        if not use_any and g.get("gene_type") not in allowed_types:
            continue
        g_type = g.get("gene_type", "unknown")
        chrom_genes[g["chrom"]].append((g["start"], g["end"], name, g_type))
        chrom_ends[g["chrom"]].append((g["end"], name, g_type))

    for entries in chrom_genes.values():
        entries.sort()
    for entries in chrom_ends.values():
        entries.sort()

    chrom_starts   = {c: [e[0] for e in v] for c, v in chrom_genes.items()}
    chrom_end_vals = {c: [e[0] for e in v] for c, v in chrom_ends.items()}
    chrom_max_len  = {c: max(e[1] - e[0] for e in v) for c, v in chrom_genes.items()}

    for r in regions:
        if r["gene"] is not None:
            continue

        chrom  = r["chrom"]
        genes  = chrom_genes.get(chrom)
        starts = chrom_starts.get(chrom)

        if not genes:
            if verbose:
                print(f"  WARNING: {r['name']}: no eligible gene found on {chrom}. "
                      "Normalization will be skipped (output: NA).")
            continue

        r_start    = r["start"]
        r_end      = r["end"]
        max_len    = chrom_max_len[chrom]
        candidates = []   # (dist, -overlap_bp, not_protein_coding, name)

        # ---- Overlapping genes ------------------------------------------------
        # idx = first gene with start >= r_end; genes[0:idx] all have start < r_end
        # and are the only ones that can overlap. Scan backward; break once
        # g_start < r_start - max_len (g_end ≤ g_start + max_len < r_start).
        idx = bisect.bisect_left(starts, r_end)
        for i in range(idx - 1, -1, -1):
            g_start, g_end, name, g_type = genes[i]
            if g_start < r_start - max_len:
                break
            if g_end > r_start:
                overlap_bp = min(r_end, g_end) - max(r_start, g_start)
                is_coding  = (g_type == "protein_coding")
                candidates.append((0, -overlap_bp, not is_coding, name))

        # ---- Right neighbor ---------------------------------------------------
        # genes[idx] is the first gene starting at or after r_end.
        if idx < len(genes):
            g_start, _, name, g_type = genes[idx]
            dist_right = g_start - r_end
            candidates.append((dist_right, 0, g_type != "protein_coding", name))
            # Collect any ties at exactly the same distance (same start position).
            for j in range(idx + 1, len(genes)):
                if genes[j][0] - r_end != dist_right:
                    break
                candidates.append((dist_right, 0, genes[j][3] != "protein_coding", genes[j][2]))

        # ---- Left neighbor ----------------------------------------------------
        # bisect_right on end values finds the insertion point after all ends
        # ≤ r_start; the entry just before it is the gene whose end is closest
        # to r_start from the left (i.e. minimum distance).
        end_vals = chrom_end_vals[chrom]
        ends     = chrom_ends[chrom]
        j = bisect.bisect_right(end_vals, r_start) - 1
        if j >= 0:
            left_end, left_name, left_gtype = ends[j]
            dist_left = r_start - left_end
            candidates.append((dist_left, 0, left_gtype != "protein_coding", left_name))
            # Collect ties (same end value → same distance).
            for k in range(j - 1, -1, -1):
                if ends[k][0] != left_end:
                    break
                candidates.append((dist_left, 0, ends[k][2] != "protein_coding", ends[k][1]))

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
# Parallel I/O helpers
# ---------------------------------------------------------------------------

# Module-level shared state populated once per worker process by _worker_init.
# Only non-empty inside ProcessPoolExecutor worker processes.
_SHARED_DATA: dict = {}


def _worker_init(shared: dict) -> None:
    """Populate per-process shared state at worker startup (called once per worker)."""
    global _SHARED_DATA
    _SHARED_DATA = shared


def _fetch_chunk(args):
    """
    Fetch VNTR read names for a batch of regions, opening a private file handle.
    Each thread/process gets its own handle to avoid pysam thread-safety issues.

    args = (input_file, reference, chunk)
    Returns {name: read_set} for every region in chunk.
    """
    input_file, reference, chunk = args
    from ._io import open_alignment_file as _open_aln
    result = {}
    with _open_aln(input_file, reference) as aln:
        for r in chunk:
            result[r["name"]] = get_read_names(
                aln, r["chrom"], r["start"], r["end"], r["alt_regions"]
            )
    return result


def _fetch_window_chunk(args):
    """
    Fetch local-normalization window read names for a batch of VNTRs, opening a
    private file handle.

    args = (input_file, reference, chunk, norm_window, use_alt)
    Returns {name: read_set} for every region in chunk.
    """
    input_file, reference, chunk, norm_window, use_alt = args
    from ._io import open_alignment_file as _open_aln
    result = {}
    with _open_aln(input_file, reference) as aln:
        for r in chunk:
            win_start = max(0, r["start"] - norm_window)
            win_end   = r["end"] + norm_window
            win_alts  = None
            if use_alt and r["alt_regions"]:
                win_alts = [
                    {
                        "chrom": a["chrom"],
                        "start": max(0, a["start"] - norm_window),
                        "end":   a["end"] + norm_window,
                    }
                    for a in r["alt_regions"]
                ]
            result[r["name"]] = get_read_names(aln, r["chrom"], win_start, win_end, win_alts)
    return result


def _parallel_fetch(fetch_fn, input_file, reference, region_list, extra_args, workers):
    """
    Split region_list into up to `workers` equal chunks and dispatch via
    ThreadPoolExecutor, with each thread running fetch_fn on its chunk.

    fetch_fn signature: (input_file, reference, chunk, *extra_args) -> {name: read_set}
    Returns a merged {name: read_set} dict.
    """
    n = len(region_list)
    chunk_size = max(1, (n + workers - 1) // workers)
    chunks = [region_list[i:i + chunk_size] for i in range(0, n, chunk_size)]
    args_list = [(input_file, reference, chunk) + extra_args for chunk in chunks]
    cache = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for partial in pool.map(fetch_fn, args_list):
            cache.update(partial)
    return cache


# ---------------------------------------------------------------------------
# Per-sample worker (top-level for ProcessPoolExecutor pickling)
# ---------------------------------------------------------------------------

def _process_sample_task(task: dict):
    """
    Process a single alignment file and return (row_dict, log_lines).

    When running inside a ProcessPoolExecutor worker, shared read-only data is
    taken from the module-level _SHARED_DATA (populated once by _worker_init).
    When called directly on the main process (serial path), all data must be
    present in *task* itself.
    """
    shared = _SHARED_DATA if _SHARED_DATA else task

    input_file   = task["input_file"]
    vntr_threads = task["vntr_threads"]
    verbose      = task["verbose"]

    region_list      = shared["region_list"]
    reference        = shared["reference"]
    normalize        = shared["normalize"]
    norm_method      = shared["norm_method"]
    norm_window      = shared["norm_window"]
    use_alt          = shared["use_alt"]
    gene_to_vntrs    = shared["gene_to_vntrs"]
    norm_regions     = shared["norm_regions"]
    gene_excl_vntrs  = shared["gene_excl_vntrs"]
    gene_eff_len     = shared["gene_eff_len"]
    gene_alt_regions = shared["gene_alt_regions"]
    has_period       = shared["has_period"]
    metric_label     = shared["metric_label"]
    vntr_names       = shared["vntr_names"]

    log = []
    sample = get_sample_name(input_file)

    if verbose:
        log.append(f"\nProcessing {sample} ...")

    if normalize:
        # ---- Build VNTR read cache ----
        if vntr_threads > 1:
            vntr_read_cache = _parallel_fetch(
                _fetch_chunk, input_file, reference, region_list, (), vntr_threads
            )
        else:
            with open_alignment_file(input_file, reference) as aln:
                vntr_read_cache = {
                    r["name"]: get_read_names(
                        aln, r["chrom"], r["start"], r["end"], r["alt_regions"]
                    )
                    for r in region_list
                }

        if norm_method == "gene":
            # Gene background reads — one fetch per normalization key; fewer calls
            # than VNTR fetches so kept serial on the main process thread.
            gene_excl_cache = {}
            with open_alignment_file(input_file, reference) as aln:
                for norm_key in gene_to_vntrs:
                    g          = norm_regions[norm_key]
                    gene_reads = get_read_names(
                        aln, g["chrom"], g["start"], g["end"],
                        gene_alt_regions.get(norm_key),
                    )
                    # Exclude ALL VNTRs overlapping this region by coordinates,
                    # not just those whose gene assignment matches norm_key.
                    excl_reads = set()
                    for v in gene_excl_vntrs.get(norm_key, []):
                        excl_reads |= vntr_read_cache[v["name"]]
                    gene_excl_cache[norm_key] = len(gene_reads - excl_reads)

        else:  # local
            # ---- Build window read cache ----
            if vntr_threads > 1:
                window_read_cache = _parallel_fetch(
                    _fetch_window_chunk, input_file, reference, region_list,
                    (norm_window, use_alt), vntr_threads
                )
            else:
                window_read_cache = {}
                with open_alignment_file(input_file, reference) as aln:
                    for r in region_list:
                        win_start = max(0, r["start"] - norm_window)
                        win_end   = r["end"] + norm_window
                        win_alts  = None
                        if use_alt and r["alt_regions"]:
                            win_alts = [
                                {
                                    "chrom": a["chrom"],
                                    "start": max(0, a["start"] - norm_window),
                                    "end":   a["end"] + norm_window,
                                }
                                for a in r["alt_regions"]
                            ]
                        window_read_cache[r["name"]] = get_read_names(
                            aln, r["chrom"], win_start, win_end, win_alts
                        )

        # ---- Compute per-VNTR metrics (all from cache — no more I/O) ----
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

                # All other VNTRs whose primary coords overlap this window —
                # exclude their reads from the background.
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

                window_reads   = window_read_cache[r["name"]]
                window_excl    = len(window_reads - excl_reads)
                window_eff_len = max(2 * norm_window - neighbor_bp, 1)
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
                log.append(f"  {r['name']}: {vntr_count} VNTR reads | {metric_label} = {display}")
            values.append(val)

    else:
        # Raw read counts — no normalization
        values = []
        with open_alignment_file(input_file, reference) as aln:
            for r in region_list:
                count = len(get_read_names(
                    aln, r["chrom"], r["start"], r["end"], r["alt_regions"]
                ))
                if verbose:
                    log.append(f"  {r['name']}: {count} reads")
                values.append(count)

    row = {"sample": sample, "metric": metric_label, **dict(zip(vntr_names, values))}
    return row, log


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
    chrom=None,
    gtf=DEFAULT_GTF,
    reference=None,
    psl=DEFAULT_PSL,
    no_alt_contigs=False,
    norm_gene_types=("protein_coding",),
    norm_method="gene",
    norm_window=5000,
    output_csv=None,
    verbose=None,
    workers=1,
    vntr_workers=None,
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
    chrom : str or list of str, optional
        Restrict genotyping to VNTRs on these chromosome(s). Accepts a single
        chromosome name (e.g. 'chr5') or a list (e.g. ['chr5', 'chr6']).
        Applied after all region selection; VNTRs on other chromosomes are
        silently dropped. Exits with an error if no VNTRs remain after filtering.
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
    workers : int or None, optional
        Number of parallel workers (default: 1 — fully serial).
        Pass None to use all available CPU cores (os.cpu_count()).
        With multiple input files, workers are assigned one process per sample
        via ProcessPoolExecutor, giving near-linear speedup with sample count.
        With a single input file, workers are assigned as threads for parallel
        VNTR read fetching within that sample.
    vntr_workers : int or None, optional
        Number of threads used for parallel VNTR read fetching within each
        sample (default: None — auto).
        When None, the value is derived from *workers*: 1 thread per process
        pool worker (no over-subscription), or *workers* threads when running
        a single sample serially.
        Set explicitly to use both parallelism strategies simultaneously —
        e.g. workers=4, vntr_workers=4 spawns 4 sample processes each using
        4 VNTR-fetch threads (16 total concurrent file handles per run).
        Pass 0 to use all available CPU cores.

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
    # Optional chromosome filter
    # ------------------------------------------------------------------
    if chrom is not None:
        chrom_set = {chrom} if isinstance(chrom, str) else set(chrom)
        region_list = [r for r in region_list if r["chrom"] in chrom_set]
        if not region_list:
            sys.exit(
                f"ERROR: No VNTRs found on chromosome(s) {sorted(chrom_set)} "
                "after filtering. Check that chromosome names match the BED file "
                "(e.g. 'chr5' not '5')."
            )
        if verbose:
            found = sorted({r["chrom"] for r in region_list})
            missing = sorted(chrom_set - set(found))
            print(f"Chromosome filter: {len(region_list)} VNTR(s) retained on {found}.")
            if missing:
                for c in missing:
                    print(f"  WARNING: no VNTRs found on '{c}' — skipped.")

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
    norm_regions     = {}   # populated in gene mode; empty dict otherwise
    use_alt          = False

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

    # Shared read-only data passed to every worker.
    shared_data = {
        "reference":       reference,
        "region_list":     region_list,
        "normalize":       normalize,
        "norm_method":     norm_method,
        "norm_window":     norm_window,
        "use_alt":         use_alt,
        "gene_to_vntrs":   dict(gene_to_vntrs),
        "norm_regions":    norm_regions,
        "gene_excl_vntrs": gene_excl_vntrs,
        "gene_eff_len":    gene_eff_len,
        "gene_alt_regions": gene_alt_regions,
        "has_period":      has_period,
        "metric_label":    metric_label,
        "vntr_names":      vntr_names,
    }

    n_sample_workers = os.cpu_count() if workers is None else max(1, int(workers))
    n_samples = len(input_files)

    # Resolve VNTR thread count.
    # vntr_workers=None (auto): use n_sample_workers for single-sample serial runs;
    # default to 1 inside each process pool worker to avoid over-subscription.
    # vntr_workers set explicitly: honor it in all paths so both strategies can
    # run simultaneously (e.g. workers=4, vntr_workers=4 → 16 concurrent handles).
    _vntr_explicit = vntr_workers is not None
    if _vntr_explicit:
        n_vntr_workers = os.cpu_count() if vntr_workers == 0 else max(1, int(vntr_workers))
    else:
        n_vntr_workers = n_sample_workers  # used only in serial/single-sample path

    rows = []

    if n_sample_workers > 1 and n_samples > 1:
        # ------------------------------------------------------------------
        # Process-level parallelism: one process per sample.
        # Shared data is sent once per worker via the initializer, not once
        # per task, to avoid redundant pickling overhead.
        # vntr_threads defaults to 1 per worker unless vntr_workers is
        # explicitly set, in which case both strategies run simultaneously.
        # ------------------------------------------------------------------
        vntr_threads_in_pool = n_vntr_workers if _vntr_explicit else 1
        if verbose:
            print(f"\nDispatching {n_samples} samples across "
                  f"{min(n_sample_workers, n_samples)} worker process(es)"
                  + (f", {vntr_threads_in_pool} VNTR thread(s) each ..."
                     if vntr_threads_in_pool > 1 else " ..."))

        tasks = [
            {"input_file": f, "vntr_threads": vntr_threads_in_pool, "verbose": verbose}
            for f in input_files
        ]
        with ProcessPoolExecutor(
            max_workers=min(n_sample_workers, n_samples),
            initializer=_worker_init,
            initargs=(shared_data,),
        ) as pool:
            results = list(pool.map(_process_sample_task, tasks))

        for row, log_lines in results:
            if verbose:
                for line in log_lines:
                    print(line)
            rows.append(row)

    else:
        # ------------------------------------------------------------------
        # Serial sample loop.
        # Use threads for VNTR reads when: (a) vntr_workers is explicitly set,
        # or (b) auto mode with a single sample (inherit workers count).
        # ------------------------------------------------------------------
        if _vntr_explicit:
            vntr_threads = n_vntr_workers
        else:
            vntr_threads = n_sample_workers if n_samples == 1 else 1
        for input_file in input_files:
            task = {
                "input_file": input_file,
                "vntr_threads": vntr_threads,
                "verbose": verbose,
                **shared_data,
            }
            row, log_lines = _process_sample_task(task)
            if verbose:
                for line in log_lines:
                    print(line)
            rows.append(row)

    df = pd.DataFrame(rows, columns=["sample", "metric"] + vntr_names)

    if output_csv is not None:
        df.to_csv(output_csv, index=False)
        if verbose:
            print(f"\nResults written to {output_csv}")

    return df
