"""
run_all.py — iterate the sample manifest, extracting features ONE SAMPLE AT A TIME.

    python -m pipeline.run.run_all                 # all pending samples
    python -m pipeline.run.run_all --only HG00733 HG00512
    python -m pipeline.run.run_all --limit 3       # first 3 pending

Resumable: samples whose features/<sample_id>.parquet already exists are skipped
(unless runner.resume is false).  Per-sample failures are logged to
work/features/_run_log.tsv and do not stop the loop.
"""

from __future__ import annotations

import argparse
import time
import traceback
from datetime import datetime, timezone

from ..paths import load_config
from .manifest import load_samples
from .process_sample import process


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--only", nargs="+", default=None, help="restrict to these sample_ids")
    ap.add_argument("--limit", type=int, default=None, help="process at most N pending")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    samples = load_samples(cfg)
    if args.only:
        want = set(args.only)
        samples = [s for s in samples if s["sample_id"] in want]

    feat_dir = cfg.ensure_dir("outputs.features_dir")
    log_path = feat_dir / "_run_log.tsv"
    if not log_path.exists():
        log_path.write_text("utc\tsample_id\tstatus\tseconds\tnote\n")

    pending = [s for s in samples
               if not (cfg.get("runner.resume", True)
                       and (feat_dir / f"{s['sample_id']}.parquet").exists())]
    if args.limit:
        pending = pending[: args.limit]
    print(f"{len(samples)} in manifest | {len(pending)} to process")

    for i, s in enumerate(pending, 1):
        sid = s["sample_id"]
        print(f"\n===== [{i}/{len(pending)}] {sid} =====")
        t0 = time.time()
        try:
            process(cfg, s)
            status, note = "ok", ""
        except KeyboardInterrupt:
            raise
        except Exception as exc:                       # noqa: BLE001
            status, note = "FAIL", f"{type(exc).__name__}: {exc}"
            traceback.print_exc()
        dt = time.time() - t0
        with open(log_path, "a") as fh:
            fh.write(f"{datetime.now(timezone.utc).isoformat()}\t{sid}\t{status}"
                     f"\t{dt:.0f}\t{note}\n")
        print(f"----- {sid}: {status} ({dt:.0f}s) -----")

    print(f"\ndone. log: {log_path}")


if __name__ == "__main__":
    main()
