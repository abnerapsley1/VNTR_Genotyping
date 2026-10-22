"""Downloads: default reference and the shipped model bundle (fetched on first use).

Parallel HTTP range requests with per-chunk retries and resume, then an MD5/SHA-256
check where the source publishes one.
"""
from __future__ import annotations

import hashlib
import http.client
import json
import tarfile
import urllib.request
from importlib import resources
from pathlib import Path

from . import paths

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



def file_hash(path: Path, algo: str = "sha256") -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as fh:
        while b := fh.read(16 << 20):
            h.update(b)
    return h.hexdigest()


def registry() -> dict:
    """Known references and model bundles (vntrlen/data/registry.json)."""
    return json.loads(resources.files("vntrlen").joinpath("data/registry.json").read_text())


def resolve_reference(ref: str | None) -> Path:
    """A FASTA path, or a registry name (default: the 1000 Genomes GRCh38 reference),
    downloaded with its index on first use."""
    reg = registry()["references"]
    name = ref or registry()["default_reference"]
    if name not in reg:
        p = Path(name).expanduser()
        if not p.exists():
            raise SystemExit(f"reference {name} not found (not a file, not one of {sorted(reg)})")
        return p
    info = reg[name]
    dest = paths.references_dir() / Path(info["fasta_url"]).name
    if not dest.exists():
        print(f"downloading reference {name} ({info.get('size_gb', '?')} GB, one-time) ...", flush=True)
        fetch_parallel(info["fasta_url"], dest)
        if info.get("fasta_md5") and file_hash(dest, "md5") != info["fasta_md5"]:
            dest.unlink()
            raise SystemExit(f"{name}: downloaded FASTA failed its MD5 check; please retry")
    fai = Path(f"{dest}.fai")
    if not fai.exists():
        if info.get("fai_url"):
            fetch_parallel(info["fai_url"], fai, n_streams=1)
        else:
            import pysam
            pysam.faidx(str(dest))
    return dest


def install_bundle(source: str | None = None) -> Path:
    """Install the shipped model bundle from a local .tar.gz or the registry URL."""
    info = registry()["model_bundle"]
    if source and Path(source).exists():
        tgz = Path(source)
    else:
        url = source or info.get("url")
        if not url:
            raise SystemExit("no download URL for the shipped models is configured yet; install a bundle "
                             "file with: vntrlen models install --bundle PATH.tar.gz")
        tgz = paths.cache_dir() / Path(url).name
        fetch_parallel(url, tgz)
        if info.get("sha256") and file_hash(tgz) != info["sha256"]:
            tgz.unlink()
            raise SystemExit("model bundle failed its checksum; please retry")
    with tarfile.open(tgz) as tf:
        tf.extractall(paths.models_dir(), filter="data")
    print(f"installed model bundle {tgz.name} into {paths.models_dir()}")
    return paths.models_dir()
