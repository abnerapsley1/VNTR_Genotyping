"""Input checks: reference identity, short reads only, single- vs paired-end."""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import numpy as np
import pysam

LONG_READ_PLATFORMS = ("PACBIO", "ONT", "NANOPORE", "OXFORD")
MAX_SHORT_READ_LEN = 300


class InputError(Exception):
    """A user input that vntrlen cannot use (message explains why)."""


def reference_contigs(fasta: str | Path) -> dict[str, int]:
    fai = Path(f"{fasta}.fai")
    if not fai.exists():
        pysam.faidx(str(fasta))
    out = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def reference_fingerprint(fasta: str | Path) -> str:
    """Short ID of a reference: hash of its primary contig names and lengths
    (chr1-22, X, Y, M with or without 'chr'), so decoy/ALT additions do not change it
    but a different build does."""
    contigs = reference_contigs(fasta)
    prim = sorted((k, v) for k, v in contigs.items()
                  if k.removeprefix("chr") in {*map(str, range(1, 23)), "X", "Y", "M", "MT"})
    if not prim:
        prim = sorted(contigs.items())
    return hashlib.sha1(repr(prim).encode()).hexdigest()[:12]


def check_alignment_reference(path: str, fasta: str | Path) -> None:
    """Stop if the BAM/CRAM was aligned to a different reference than `fasta`:
    every @SQ in the alignment must exist in the reference with the same length,
    and the same MD5 where both declare one."""
    ref = reference_contigs(fasta)
    with pysam.AlignmentFile(path, "rc" if path.endswith(".cram") else "rb",
                             reference_filename=str(fasta), check_sq=False) as a:
        sq = a.header.to_dict().get("SQ", [])
    if not sq:
        raise InputError(f"{path}: no @SQ lines; cannot verify the reference (is it aligned?)")
    missing = [s["SN"] for s in sq if s["SN"] not in ref]
    wrong = [s["SN"] for s in sq if s["SN"] in ref and int(s["LN"]) != ref[s["SN"]]]
    if wrong or len(missing) > 0.5 * len(sq):
        ex = (wrong or missing)[:3]
        raise InputError(f"{path} was aligned to a different reference than {fasta} "
                         f"(contig names or lengths differ, e.g. {ex}). Re-align, or pass the "
                         f"reference it was aligned to.")
    prim = [s for s in sq if s["SN"].removeprefix("chr") in {"1", "2", "X"}]
    if not prim:
        raise InputError(f"{path}: none of chr1/chr2/chrX in the header; vntrlen needs a "
                         "whole-genome human alignment")


def _read_lengths_bam(path: str, fasta, n: int = 20000) -> tuple[np.ndarray, float, str]:
    lens, paired = [], 0
    with pysam.AlignmentFile(path, "rc" if path.endswith(".cram") else "rb",
                             reference_filename=str(fasta) if fasta else None) as a:
        plat = " ".join(rg.get("PL", "") for rg in a.header.to_dict().get("RG", [])).upper()
        for r in a.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            lens.append(r.infer_read_length() or r.query_length or 0)
            paired += r.is_paired
            if len(lens) >= n:
                break
    return np.asarray(lens), paired / max(1, len(lens)), plat


def _read_lengths_fastq(path: str, n: int = 20000) -> np.ndarray:
    op = gzip.open if str(path).endswith(".gz") else open
    lens = []
    with op(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lens.append(len(line.strip()))
                if len(lens) >= n:
                    break
    return np.asarray(lens)


def check_short_reads(sample, fasta=None) -> bool:
    """Raise InputError for long-read data; return True if the sample is paired-end."""
    if sample.alignment:
        lens, frac_paired, plat = _read_lengths_bam(sample.alignment, fasta)
        if any(p in plat for p in LONG_READ_PLATFORMS):
            raise InputError(f"{sample.sample_id}: read group platform '{plat.strip()}' is long-read; "
                             "vntrlen supports short-read (Illumina-type) data only")
        paired = frac_paired > 0.5
    else:
        lens = _read_lengths_fastq(sample.fastq1)
        paired = sample.fastq2 is not None
    if lens.size == 0:
        raise InputError(f"{sample.sample_id}: no reads found")
    med = float(np.median(lens))
    if med > MAX_SHORT_READ_LEN:
        raise InputError(f"{sample.sample_id}: median read length {med:.0f} bp looks like long-read "
                         f"data; vntrlen supports short reads (<= {MAX_SHORT_READ_LEN} bp) only")
    return paired
