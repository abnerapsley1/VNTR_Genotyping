"""
remote.py — open BAM/CRAM by remote region-slicing, plus reference/index caching.

htslib (via pysam) does HTTP-range / S3 range reads on a remote BAM/CRAM as
long as an index is available and, for CRAM, a local reference FASTA is given.
In `runner.transfer_mode: remote` nothing is fully downloaded except the (small)
index and the reference (cached).  `download` mode instead pulls the whole
alignment with `fetch_parallel` — cheaper for the dense bundled BED, where
per-locus slicing re-reads overlapping CRAM slices and moves more bytes than
the file itself.
"""

from __future__ import annotations

import hashlib
import http.client
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


def fetch_parallel(url: str, dest: Path, n_streams: int = 16,
                   retries: int = 30, chunk_bytes: int = 64 << 20) -> Path:
    """Download `url` to `dest` with parallel HTTP range requests (skip if present).

    Single EBI connections get throttled to ~1 MB/s; many concurrent ranges
    recover most of the link.  Each chunk is retried independently and its
    byte count verified; the file is renamed into place only when complete.
    Finished chunks are logged beside the partial file, so an interrupted
    download resumes instead of starting over.
    """
    from concurrent.futures import ThreadPoolExecutor
    import threading
    import time

    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):     # the size request can be dropped too
        try:
            with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"),  # noqa: S310
                                        timeout=60) as r:
                size = int(r.headers["Content-Length"])
            break
        except (OSError, http.client.HTTPException):
            if attempt == retries - 1:
                raise
            time.sleep(min(120, 10 * (attempt + 1)))
    tmp = dest.with_suffix(dest.suffix + ".part")
    done_log = dest.with_suffix(dest.suffix + ".part.done")
    done: set[int] = set()
    if tmp.exists() and tmp.stat().st_size == size and done_log.exists():
        done = {int(x) for x in done_log.read_text().split()}
    else:
        with open(tmp, "wb") as fh:
            fh.truncate(size)
        done_log.write_text("")
    lock = threading.Lock()

    def _get(rng: tuple[int, int]) -> None:
        s, e = rng
        want = e - s + 1
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e}"})
                with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "r+b") as fh:  # noqa: S310
                    # Only a 206 for exactly this range may be written: a server that
                    # ignores Range answers 200 with the whole file, which written at
                    # offset s would silently overwrite neighbouring, finished chunks.
                    cr = r.headers.get("Content-Range", "")
                    if r.status != 206 or not cr.startswith(f"bytes {s}-{e}/"):
                        raise OSError(f"bad range response {r.status} {cr!r}")
                    fh.seek(s)
                    got = 0
                    while got < want and (buf := r.read(min(1 << 20, want - got))):
                        fh.write(buf)
                        got += len(buf)
                if got == want:
                    with lock, open(done_log, "a") as fh:
                        fh.write(f"{s}\n")
                    return
            except (OSError, http.client.HTTPException):   # resets, timeouts, short reads
                pass
            time.sleep(min(120, 10 * (attempt + 1)))
        raise OSError(f"range {s}-{e} of {url} failed after {retries} attempts")

    ranges = [(s, min(s + chunk_bytes, size) - 1) for s in range(0, size, chunk_bytes)
              if s not in done]
    with ThreadPoolExecutor(n_streams) as ex:
        list(ex.map(_get, ranges))
    tmp.rename(dest)
    done_log.unlink(missing_ok=True)
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
