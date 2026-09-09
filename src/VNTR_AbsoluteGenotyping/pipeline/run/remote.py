"""
remote.py — open BAM/CRAM by remote region-slicing, plus reference/index caching.

htslib (via pysam) does HTTP-range / S3 range reads on a remote BAM/CRAM as
long as an index is available and, for CRAM, a local reference FASTA is given.
Nothing is fully downloaded except the (small) index and the reference (cached).
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import urllib.request
from pathlib import Path

import pysam


# ---------------------------------------------------------------------------
# reference + index caching
# ---------------------------------------------------------------------------

def _cache_name(url: str) -> str:
    return hashlib.sha1(url.encode()).hexdigest()[:16]


def fetch_to_cache(url: str, cache_dir: Path, suffix: str = "") -> Path:
    """Download `url` into cache_dir (skip if present). Local paths pass through."""
    if "://" not in url or url.startswith("file://"):
        p = Path(url.replace("file://", ""))
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{_cache_name(url)}{suffix or Path(url).suffix}"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    if url.startswith("s3://"):
        subprocess.run(["aws", "s3", "cp", "--no-progress", url, str(tmp)], check=True)
    else:
        urllib.request.urlretrieve(url, tmp)  # noqa: S310 (trusted manifest URLs)
    tmp.rename(dest)
    return dest


def ensure_reference(fasta_url: str, fai_url: str | None, cache_dir: Path) -> Path:
    fa = fetch_to_cache(fasta_url, cache_dir, suffix=".fa")
    fai = Path(str(fa) + ".fai")
    if not fai.exists():
        if fai_url and fai_url not in (".", ""):
            got = fetch_to_cache(fai_url, cache_dir, suffix=".fa.fai")
            if got != fai:
                fai.write_bytes(got.read_bytes())
        else:
            pysam.faidx(str(fa))
    return fa


def ensure_index(data_url: str, index_url: str | None, cache_dir: Path) -> Path:
    if not index_url or index_url in (".", ""):
        low = data_url.lower()
        index_url = data_url + (".crai" if low.endswith(".cram") else ".bai")
    suffix = ".crai" if index_url.lower().endswith(".crai") else ".bai"
    return fetch_to_cache(index_url, cache_dir, suffix=suffix)


# ---------------------------------------------------------------------------
# alignment handle
# ---------------------------------------------------------------------------

def open_alignment(data_url: str, index_path: Path | None,
                   reference_fasta: Path | None) -> pysam.AlignmentFile:
    is_cram = data_url.lower().endswith(".cram")
    mode = "rc" if is_cram else "rb"
    kw: dict = {}
    if index_path is not None:
        kw["index_filename"] = str(index_path)
    if reference_fasta is not None:
        kw["reference_filename"] = str(reference_fasta)
    # htslib env knobs that make remote reads more robust
    os.environ.setdefault("HTS_S3_ADDRESSING_STYLE", "path")
    return pysam.AlignmentFile(data_url, mode, **kw)


def resolve_contig(aln: pysam.AlignmentFile, chrom: str) -> str | None:
    """Map a 'chr'-style locus contig onto whatever the alignment uses."""
    refs = set(aln.references)
    if chrom in refs:
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    if alt in refs:
        return alt
    if chrom == "chrM" and "MT" in refs:
        return "MT"
    if chrom == "MT" and "chrM" in refs:
        return "chrM"
    return None
