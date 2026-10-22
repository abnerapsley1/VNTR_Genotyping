"""Sample inputs.

Either pass alignment files directly (`--reads a.cram b.bam ...`; sample IDs come
from the @RG SM tag, else the file name), or a sample sheet (`--samples sheet.tsv`),
a tab-separated file with a header:

  sample_id    required
  alignment    BAM/CRAM path            } give either an alignment,
  fastq1       FASTQ (read 1 / single)  } or fastq1 (+ fastq2 for paired-end)
  fastq2       FASTQ read 2 (optional)
  family_id    optional; relatives are kept together during cross-validation
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl
import pysam


@dataclass
class Sample:
    sample_id: str
    alignment: str | None = None
    fastq1: str | None = None
    fastq2: str | None = None
    family_id: str | None = None

    @property
    def needs_alignment(self) -> bool:
        return self.alignment is None


def _sm_tag(path: str) -> str | None:
    try:
        with pysam.AlignmentFile(path, "rc" if path.endswith(".cram") else "rb", check_sq=False) as a:
            sms = {rg.get("SM") for rg in a.header.to_dict().get("RG", []) if rg.get("SM")}
        return sms.pop() if len(sms) == 1 else None
    except (ValueError, OSError):
        return None


def from_paths(paths: list[str]) -> list[Sample]:
    out = []
    for p in paths:
        sid = _sm_tag(p) or Path(p).name.split(".")[0]
        out.append(Sample(sample_id=sid, alignment=p))
    return _check(out)


def from_sheet(path: str | Path) -> list[Sample]:
    df = pl.read_csv(path, separator="\t", infer_schema_length=0, comment_prefix="#")
    if "sample_id" not in df.columns:
        raise ValueError(f"{path}: needs a sample_id column")
    base = Path(path).resolve().parent

    def resolve(v):
        if v is None or v in ("", "."):
            return None
        p = Path(v)
        return str(p if p.is_absolute() or "://" in v else base / p)

    out = []
    for r in df.iter_rows(named=True):
        s = Sample(sample_id=r["sample_id"], alignment=resolve(r.get("alignment")),
                   fastq1=resolve(r.get("fastq1")), fastq2=resolve(r.get("fastq2")),
                   family_id=r.get("family_id") or None)
        if not s.alignment and not s.fastq1:
            raise ValueError(f"{path}: sample {s.sample_id} has neither alignment nor fastq1")
        if s.alignment and s.fastq1:
            raise ValueError(f"{path}: sample {s.sample_id} has both alignment and FASTQ; give one")
        out.append(s)
    return _check(out)


def _check(samples: list[Sample]) -> list[Sample]:
    ids = [s.sample_id for s in samples]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise ValueError(f"duplicate sample IDs: {sorted(dup)[:5]}")
    for s in samples:
        for p in (s.alignment, s.fastq1, s.fastq2):
            if p and "://" not in p and not Path(p).exists():
                raise FileNotFoundError(f"{s.sample_id}: {p} not found")
    return samples
