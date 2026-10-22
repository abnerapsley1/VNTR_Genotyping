"""Phased SNP haplotypes around our loci, from the 1kGP 2022 phased panel
(Byrska-Bishop et al. 2022; 3,202 people; SNV/INDEL/SV, SHAPEIT-phased with pedigrees).

Per chromosome: parallel-range download from EBI (single connections are
throttled) with MD5 verification against the EBI manifest and piece-wise
repair of corrupted transfers, keep our 290 people and biallelic SNPs with MAF >= 1% among them
within +/-100 kb of our loci, then delete the full file.

    python src/VNTR_Benchmark/fetch_phased_snps.py [--jobs 4]
-> data/1000Genomes/phased_snps/<chr>.bcf (+ .csi)
"""
import argparse
import hashlib
import subprocess
import urllib.request
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))
from pipeline.run.remote import fetch_parallel  # noqa: E402

E = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
     "20220422_3202_phased_SNV_INDEL_SV/")
OUT = ROOT / "data" / "1000Genomes" / "phased_snps"
BCF = "/data/miniconda3/envs/vntr_benchmark/bin/bcftools"
CHUNK = 64 << 20


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while b := fh.read(16 << 20):
            h.update(b)
    return h.hexdigest()


def get_range(url: str, s: int, e: int) -> bytes:
    for attempt in range(10):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e - 1}"})
            with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
                b = r.read()
            if len(b) == e - s:
                return b
        except OSError:
            pass
    raise OSError(f"range {s}-{e} failed")


def verified_download(url: str, dest: Path, want: str) -> None:
    """Download; on MD5 mismatch fetch a second copy, re-fetch the 64 MB pieces where the
    copies differ, keep each piece's majority version, until the MD5 matches."""
    fetch_parallel(url, dest, n_streams=12)
    if md5(dest) == want:
        return
    alt = dest.with_name(dest.name + ".b")
    fetch_parallel(url, alt, n_streams=12)
    size = dest.stat().st_size
    for rnd in range(4):
        with open(dest, "r+b") as a, open(alt, "r+b") as b:
            for s in range(0, size, CHUNK):
                e = min(size, s + CHUNK)
                a.seek(s); x = a.read(e - s)
                b.seek(s); y = b.read(e - s)
                if x == y and rnd == 0:
                    continue
                if x != y or rnd > 0:
                    z = get_range(url, s, e)
                    good = x if x in (y, z) else (y if y == z else z)
                    a.seek(s); a.write(good)
                    b.seek(s); b.write(z)
        if md5(dest) == want:
            alt.unlink()
            print(f"  repaired {dest.name} (round {rnd + 1})", flush=True)
            return
    alt.unlink(missing_ok=True)
    raise OSError(f"{dest.name}: MD5 still wrong after repair")


def one(c: str) -> str:
    out = OUT / f"{c}.bcf"
    if Path(f"{out}.csi").exists():
        return f"{c} already done"
    name = f"1kGP_high_coverage_Illumina.{c}.filtered.SNV_INDEL_SV_phased_panel{'.v2' if c == 'chrX' else ''}.vcf.gz"
    full = OUT / "full" / name
    for attempt in range(3):
        if attempt < 2 and c not in USE_2020:
            url, want = E + name, MANIFEST[name]
        else:
            name = f"CCDG_14151_B01_GRM_WGS_2020-08-05_{c}.filtered.shapeit2-duohmm-phased.vcf.gz"
            full = OUT / "full" / name
            url, want = E2020 + name, MANIFEST_2020[name]
            print(f"{c}: falling back to the 2020 phased release", flush=True)
        try:
            verified_download(url, full, want)
        except OSError as exc:
            print(f"{c} attempt {attempt + 1} download failed: {exc}", flush=True)
            full.unlink(missing_ok=True)
            continue
        cmd = (f"{BCF} view -S {OUT/'samples.txt'} --force-samples -v snps -m2 -M2 -T {OUT/f'regions.{c}.tsv'} -Ou {full} | "
               f"{BCF} view -i 'MAF>=0.01' -Ob -o {out}.tmp - && mv {out}.tmp {out} && {BCF} index -f {out}")
        r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True)
        if r.returncode == 0:
            full.unlink()
            n = subprocess.run([BCF, "index", "-n", str(out)], capture_output=True, text=True).stdout.strip()
            return f"{c} ok: {n} SNPs"
        print(f"{c} attempt {attempt + 1} failed: {r.stderr[-300:]}", flush=True)
        full.unlink(missing_ok=True)
    return f"{c} FAILED"


MANIFEST: dict = {}
# fallback: the 2020 release of the same panel (SHAPEIT2-duohmm, same 3,202 people), used for a
# chromosome whose 2022 file fails its EBI MD5 even after piece-wise repair (seen for chr6, chr9)
E2020 = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/"
         "20201028_3202_phased/")
MANIFEST_2020: dict = {}
USE_2020: set = set()


def main():
    MANIFEST.update(dict(line.split() for line in urllib.request.urlopen(  # noqa: S310
        E + "20220804_manifest.txt").read().decode().splitlines() if line.strip()))
    MANIFEST_2020.update({ln.split("\t")[0]: ln.split("\t")[2] for ln in urllib.request.urlopen(  # noqa: S310
        E2020 + "phased-manifest_July2021.tsv").read().decode().splitlines() if ln.count("\t") >= 2})
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--use-2020", nargs="*", default=[], help="chromosomes to take from the 2020 release directly")
    args = ap.parse_args()
    USE_2020.update(args.use_2020)
    (OUT / "full").mkdir(parents=True, exist_ok=True)
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    chroms = [c for c in chroms if (OUT / f"regions.{c}.tsv").exists()]
    with ThreadPoolExecutor(args.jobs) as ex:
        for msg in ex.map(one, chroms):
            print(msg, flush=True)
    print("all chromosomes done", flush=True)


if __name__ == "__main__":
    main()
