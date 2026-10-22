"""Download (and keep) CRAM + CRAI files for the GangSTR / danbing-tk comparison,
using the pipeline's parallel range fetcher.  Default: the 15 test individuals;
danbing-tk's length model also needs the 48 training individuals.

    python src/VNTR_Benchmark/fetch_crams.py                 # test
    python src/VNTR_Benchmark/fetch_crams.py --folds train   # training
"""
import argparse
import hashlib
import time
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))

import polars as pl  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402
from pipeline.run.remote import fetch_parallel  # noqa: E402

OUT = ROOT / "data" / "Benchmark" / "cram"


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while buf := fh.read(16 << 20):
            h.update(buf)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="+", default=["test"], choices=["test", "val", "train"])
    ap.add_argument("--streams", type=int, default=8)
    ap.add_argument("--only", nargs="+", default=None, help="restrict to these sample IDs")
    args = ap.parse_args()
    cfg = load_config(None)
    split = pl.read_parquet(cfg.path("outputs.model_dir") / "sample_split.parquet")
    test = set(split.filter(pl.col("fold").is_in(args.folds))["sample_id"])
    if args.only:
        test &= set(args.only)
    md5_tsv = ROOT / "src" / "VNTR_Benchmark" / "work" / "cram_md5_expected.tsv"
    expected = {}
    if md5_tsv.exists():
        for line in md5_tsv.read_text().splitlines()[1:]:
            f = line.split("\t")
            expected[f[0]] = f[2]
    check = ROOT / "src" / "VNTR_Benchmark" / "work" / "cram_md5_check.tsv"
    verified = ({l.split("\t")[0] for l in check.read_text().splitlines() if l.endswith("\tOK")}
                if check.exists() else set())
    failed = []
    for s in load_samples(cfg):
        if s["sample_id"] not in test:
            continue
        sid, url = s["sample_id"], s["data_url"]
        idx = s.get("index_url") or url + ".crai"
        for attempt in range(3):          # each attempt resumes the partial file
            try:
                print(f"[{sid}] index ...", flush=True)
                fetch_parallel(idx, OUT / f"{sid}.cram.crai", n_streams=1)
                print(f"[{sid}] cram ...", flush=True)
                fetch_parallel(url, OUT / f"{sid}.cram", n_streams=args.streams)
                if sid in expected and sid not in verified and md5(OUT / f"{sid}.cram") != expected[sid]:
                    (OUT / f"{sid}.cram").unlink()
                    raise OSError("MD5 mismatch against ENA; deleted, re-downloading")
                print(f"[{sid}] done {(OUT / f'{sid}.cram').stat().st_size / 1e9:.1f} GB", flush=True)
                break
            except OSError as exc:
                print(f"[{sid}] attempt {attempt + 1} interrupted: {exc}", flush=True)
                time.sleep(300)
        else:
            failed.append(sid)
    print(f"finished; failed: {failed or 'none'}", flush=True)


if __name__ == "__main__":
    main()
