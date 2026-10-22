"""danbing-tk genotyping for the 290-sample benchmark, unattended.

  1. waits for the RPGG built from the training individuals' assemblies
     (build_rpgg_290.sh logs "RPGG_290 built")
  2. genotypes the test individuals first (CRAMs kept in data/Benchmark/cram),
     then the training individuals (needed for the per-locus length model),
     whose CRAMs are fetched from the AWS mirror, MD5-checked, and deleted
     after alignment
  3. runs 2 samples at a time while our extraction / GangSTR services still
     use the CPU, 4 at a time once they have finished

Reads reach danbing-tk as name-paired FASTQ: `samtools collate -f` (fast mode,
primary reads only -- the same reads `samtools fastq` emits by default) keeps
mates adjacent without writing a full temporary copy; its small spill files go
to /dev/shm.  Staged CRAMs go to the system disk to keep /data free.
Options follow the danbing-tk README (-gc 85 3 -kf 4 1 -cth 45 -k 21).

Output: work/danbing/align_290/<sample>.tr.kmers (+ .done with seconds taken)

    python -u src/VNTR_Benchmark/run_danbing_290.py
"""
import subprocess
import sys
import threading
import time
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402
from pipeline.run.process_sample import _md5  # noqa: E402
from pipeline.run.remote import fetch_parallel  # noqa: E402

DB = ROOT / "src" / "VNTR_Benchmark" / "work" / "danbing"
RPGG = DB / "rpgg_290" / "pan"
OUT = DB / "align_290"
RPGG_LOG = ROOT / "src" / "VNTR_Benchmark" / "work" / "logs" / "build_rpgg_290.log"
KEPT = ROOT / "data" / "Benchmark" / "cram"
STAGE = Path.home() / "vntr_scratch" / "danbing_cram"
TMP = Path("/dev/shm/danbing_tmp")
REF = ROOT / "data" / "Benchmark" / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
SAMTOOLS = "/data/miniconda3/envs/vntr_benchmark/bin/samtools"
DANBING = ROOT / "tools" / "danbing-tk" / "bin" / "danbing-tk"
MAX_WAITING = 2         # downloaded training CRAMs not yet being aligned (~17 GB each)
BUSY_SERVICES = ["vntr-extract.service", "vntr-gangstr.service"]


def log(msg: str) -> None:
    print(f"{time.strftime('%F %T')} {msg}", flush=True)


def slots() -> tuple[int, int]:
    """(concurrent samples, danbing threads each) given what else is running."""
    busy = subprocess.run(["systemctl", "--user", "is-active", *BUSY_SERVICES],
                          capture_output=True, text=True).stdout.split()
    return (2, 16) if "active" in busy else (4, 8)


def done(s: str) -> bool:
    return (OUT / f"{s}.done").exists()


def align(s: str, cram: Path, threads: int) -> None:
    out = OUT / s
    t0 = time.time()
    cmd = (f"{SAMTOOLS} collate -f -r 2000000 -@2 -u -O --reference {REF} {cram} {TMP / s} | "
           f"{SAMTOOLS} fastq -n - 2>/dev/null | "
           f"{DANBING} -gc 85 3 -kf 4 1 -cth 45 -k 21 -qs {RPGG} -fq /dev/stdin -p {threads} -o {out} "
           f"> /dev/null 2> {out}.log")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    (OUT / f"{s}.done").write_text(f"{time.time() - t0:.0f}\n")


def main() -> None:
    log("waiting for RPGG_290")
    while not (RPGG_LOG.exists() and "RPGG_290 built" in RPGG_LOG.read_text(errors="ignore")
               and Path(f"{RPGG}.kmerDBi.vv").exists()):
        time.sleep(300)
    log("RPGG_290 ready")

    cfg = load_config(None)
    split = pl.read_parquet(cfg.path("outputs.model_dir") / "sample_split.parquet")
    test = sorted(split.filter(pl.col("fold") == "test")["sample_id"])
    train = (DB / "genomes_290.txt").read_text().split()
    manifest = {r["sample_id"]: r for r in load_samples(cfg)}
    order = [s for s in dict.fromkeys(test + train) if not done(s)]
    for d in (OUT, STAGE, TMP):
        d.mkdir(parents=True, exist_ok=True)
    log(f"{len(order)} samples to genotype ({len(test)} test, {len(train)} training)")

    ready: dict[str, Path] = {}          # sample -> CRAM ready to align
    staged: set[str] = set()             # downloaded copies to delete afterwards
    failed: dict[str, int] = {}
    running: set[str] = set()            # samples being aligned
    lock = threading.Condition()

    def downloader():
        for s in order:
            if (KEPT / f"{s}.cram").exists() and (KEPT / f"{s}.cram.crai").exists():
                with lock:
                    ready[s] = KEPT / f"{s}.cram"
                    lock.notify_all()
                continue
            with lock:
                lock.wait_for(lambda: len(staged - running) < MAX_WAITING)
            row = manifest[s]
            dest = STAGE / f"{s}.cram"
            try:
                fetch_parallel(row["data_url"] + ".crai", Path(f"{dest}.crai"), n_streams=1)
                for attempt in range(3):
                    fetch_parallel(row["data_url"], dest, n_streams=8)
                    if _md5(dest) == row["md5"]:
                        break
                    dest.unlink(missing_ok=True)
                    log(f"[{s}] MD5 mismatch on attempt {attempt + 1}; re-downloading")
                else:
                    raise OSError("MD5 mismatch three times")
            except Exception as exc:  # noqa: BLE001 - skip the sample, keep going
                log(f"[{s}] download FAILED: {exc}")
                with lock:
                    failed[s] = 99
                    lock.notify_all()
                continue
            with lock:
                staged.add(s)
                ready[s] = dest
                lock.notify_all()
            log(f"[{s}] staged")

    def worker(s: str, cram: Path, threads: int):
        try:
            align(s, cram, threads)
            log(f"[{s}] done in {(OUT / f'{s}.done').read_text().strip()} s ({threads} threads)")
        except subprocess.CalledProcessError as exc:
            log(f"[{s}] FAILED ({exc.returncode}); see {OUT / s}.log")
            with lock:
                failed[s] = failed.get(s, 0) + 1
        finally:
            for f in TMP.glob(f"{s}*"):
                f.unlink(missing_ok=True)
            with lock:
                if s in staged and (done(s) or failed.get(s, 0) >= 2):
                    cram.unlink(missing_ok=True)
                    Path(f"{cram}.crai").unlink(missing_ok=True)
                    staged.discard(s)
                running.discard(s)
                lock.notify_all()

    threading.Thread(target=downloader, daemon=True).start()
    while True:
        with lock:
            left = [s for s in order if not done(s) and failed.get(s, 0) < 2]
            if not left and not running:
                break
            n_slots, threads = slots()
            nxt = next((s for s in left if s in ready and s not in running), None)
            if nxt is None or len(running) >= n_slots:
                lock.wait(timeout=60)
                continue
            running.add(nxt)
        log(f"[{nxt}] aligning ({threads} threads, {len(running)} running)")
        threading.Thread(target=worker, args=(nxt, ready[nxt], threads)).start()
        time.sleep(5)
    bad = [s for s in order if not done(s)]
    log(f"all danbing-tk alignments finished; {len(bad)} failed: {bad}")


if __name__ == "__main__":
    main()
