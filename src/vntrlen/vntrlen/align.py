"""Align FASTQ input to the reference (BWA-MEM, as the 1000 Genomes NYGC data were).

Requires `bwa` and `samtools` on PATH (both in the conda environment).  The BWA index
(~1 h, ~5 GB for GRCh38) is built next to the reference the first time it is needed.
Output: <workdir>/<sample>.cram (+ .crai), coordinate-sorted, duplicates marked.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from .checks import InputError


def _need(tool: str) -> str:
    p = shutil.which(tool)
    if not p:
        raise InputError(f"'{tool}' not found on PATH; it is needed to align FASTQ input "
                         "(install with: conda install -c bioconda bwa samtools)")
    return p


def ensure_bwa_index(fasta: str | Path) -> None:
    if all(Path(f"{fasta}.{ext}").exists() for ext in ("bwt", "pac", "sa", "ann", "amb")):
        return
    print(f"building BWA index for {fasta} (one-time, about an hour for a human genome) ...", flush=True)
    subprocess.run([_need("bwa"), "index", str(fasta)], check=True)


def align_sample(sample, fasta: str | Path, workdir: Path, threads: int = 8) -> str:
    bwa, samtools = _need("bwa"), _need("samtools")
    ensure_bwa_index(fasta)
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / f"{sample.sample_id}.cram"
    if out.exists() and Path(f"{out}.crai").exists():
        return str(out)
    rg = f"@RG\\tID:{sample.sample_id}\\tSM:{sample.sample_id}\\tPL:ILLUMINA"
    reads = f"{sample.fastq1}" + (f" {sample.fastq2}" if sample.fastq2 else "")
    tmp = workdir / f"{sample.sample_id}.tmp"
    cmd = (f"{bwa} mem -t {threads} -K 100000000 -Y -R '{rg}' {fasta} {reads} | "
           f"{samtools} fixmate -m -u -@ {threads} - - | "
           f"{samtools} sort -u -@ {threads} -T {tmp} - | "
           f"{samtools} markdup -@ {threads} --reference {fasta} -O cram - {out}.part")
    print(f"[{sample.sample_id}] aligning {'paired' if sample.fastq2 else 'single'}-end FASTQ ...", flush=True)
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    Path(f"{out}.part").rename(out)
    subprocess.run([samtools, "index", str(out)], check=True)
    return str(out)
