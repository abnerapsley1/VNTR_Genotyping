"""Genotype benchmark CRAMs with danbing-tk align against the RPGG built from the
48 training individuals' assemblies.

danbing-tk pairs reads on the fly, so mates must be adjacent: the CRAM is
collated by read name and streamed as FASTQ.  FASTQ (not FASTA) is required:
v1.3.2.5's quality-aware baiting indexes base qualities and segfaults on FASTA.
Options follow the danbing-tk README (-gc 85 3 -kf 4 1 -cth 45 -k 21).

Output: work/danbing/align/<sample>.tr.kmers (per-locus k-mer counts)

    python src/VNTR_Benchmark/run_danbing_align.py HG00513 [...] [--threads 24]
"""
import argparse
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CRAM = ROOT / "data" / "Benchmark" / "cram"
REF = ROOT / "data" / "Benchmark" / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
RPGG = ROOT / "src" / "VNTR_Benchmark" / "work" / "danbing" / "rpgg" / "pan"
OUT = ROOT / "src" / "VNTR_Benchmark" / "work" / "danbing" / "align"
TMP = ROOT / "src" / "VNTR_Benchmark" / "work" / "danbing" / "tmp"
SAMTOOLS = "/data/miniconda3/envs/vntr_benchmark/bin/samtools"
DANBING = ROOT / "tools" / "danbing-tk" / "bin" / "danbing-tk"


def align(sample: str, threads: int) -> None:
    out = OUT / sample
    if Path(f"{out}.tr.kmers").exists() and Path(f"{out}.done").exists():
        print(f"[{sample}] already done", flush=True)
        return
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    cmd = (f"{SAMTOOLS} collate -@4 -u -O --reference {REF} {CRAM / (sample + '.cram')} {TMP / sample} | "
           f"{SAMTOOLS} fastq -n - 2>/dev/null | "
           f"{DANBING} -gc 85 3 -kf 4 1 -cth 45 -k 21 -qs {RPGG} -fq /dev/stdin -p {threads} -o {out} "
           f"> /dev/null 2> {out}.log")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    Path(f"{out}.done").write_text(f"{time.time() - t0:.0f}\n")
    print(f"[{sample}] {time.time() - t0:.0f}s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--threads", type=int, default=24)
    args = ap.parse_args()
    for s in args.samples:
        align(s, args.threads)


if __name__ == "__main__":
    main()
