"""Run GangSTR 2.5 on one or more benchmark CRAMs over all bundled VNTR loci.

Per sample: infer sex from chrY coverage (needed for chrX/chrY calls), split
the regions file into chunks, run GangSTR on the chunks in parallel (GangSTR
is single-threaded), concatenate, and parse the VCF into
work/gangstr/calls/<sample>.tsv with one row per locus:
    locus_id, a1_copies, a2_copies, ci, q, read classes, filter

    python src/VNTR_Benchmark/run_gangstr.py HG00513 [HG00621 ...] [--jobs 32]
"""
import argparse
import gzip
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "src" / "VNTR_Benchmark" / "work" / "gangstr"
CRAM = ROOT / "data" / "Benchmark" / "cram"
REF = ROOT / "data" / "Benchmark" / "ref" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
BIN = "/data/miniconda3/envs/vntr_benchmark/bin/"


def infer_sex(cram: Path) -> str:
    """Male if chrY's male-specific region carries reads at a meaningful depth."""
    def count(region):
        out = subprocess.run([BIN + "samtools", "view", "-c", "-q", "30", "-T", str(REF),
                              str(cram), region], capture_output=True, text=True, check=True)
        return int(out.stdout.strip())
    y = count("chrY:2781480-10000000")      # male-specific Y, ~7.2 Mb
    a = count("chr7:30000000-37218520")     # autosomal control, same width
    return "M" if a and y / a > 0.10 else "F"


def split_regions(n: int) -> list[Path]:
    rows = (WORK / "regions.bed").read_text().splitlines()
    k = -(-len(rows) // n)
    paths = []
    for i in range(n):
        p = WORK / "chunks" / f"regions.{i:03d}.bed"
        p.parent.mkdir(parents=True, exist_ok=True)
        if rows[i * k:(i + 1) * k]:
            p.write_text("\n".join(rows[i * k:(i + 1) * k]) + "\n")
            paths.append(p)
    return paths


def run_chunk(sample: str, sex: str, chunk: Path) -> Path:
    out = WORK / "raw" / sample / chunk.stem
    out.parent.mkdir(parents=True, exist_ok=True)
    if Path(str(out) + ".vcf").exists():
        return Path(str(out) + ".vcf")
    with open(str(out) + ".log", "w") as log:
        subprocess.run([BIN + "GangSTR", "--bam", str(CRAM / f"{sample}.cram"), "--ref", str(REF),
                        "--regions", str(chunk), "--out", str(out), "--bam-samps", sample,
                        "--samp-sex", sex], stdout=log, stderr=subprocess.STDOUT, check=False)
    return Path(str(out) + ".vcf")


def parse(sample: str, vcfs: list[Path]) -> int:
    idx = {}
    for line in (WORK / "regions_index.tsv").read_text().splitlines()[1:]:
        lid, chrom, start, *_ = line.split("\t")
        idx[(chrom, start)] = lid
    out = WORK / "calls" / f"{sample}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(out, "w") as fo:
        fo.write("locus_id\tfilter\tperiod\tref_copies\ta1_copies\ta2_copies\tci\tq\trc\n")
        for v in vcfs:
            if not v.exists():
                continue
            for line in open(v):
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
                cn = fmt.get("REPCN", ".").split(",")
                a1, a2 = (cn + ["."])[:2] if cn != ["."] else (".", ".")
                fo.write("\t".join([idx.get((f[0], f[1]), f"{f[0]}:{f[1]}"), f[6], info.get("PERIOD", "."),
                                    info.get("REF", "."), a1, a2, fmt.get("REPCI", "."),
                                    fmt.get("Q", "."), fmt.get("RC", ".")]) + "\n")
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("samples", nargs="+")
    ap.add_argument("--jobs", type=int, default=32)
    args = ap.parse_args()
    chunks = split_regions(args.jobs * 4)
    for s in args.samples:
        t0 = time.time()
        sex = infer_sex(CRAM / f"{s}.cram")
        with ThreadPoolExecutor(args.jobs) as ex:
            vcfs = list(ex.map(lambda c: run_chunk(s, sex, c), chunks))
        n = parse(s, vcfs)
        print(f"[{s}] sex={sex} {n:,} records, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
