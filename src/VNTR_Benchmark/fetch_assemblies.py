"""Download the haplotype-resolved assemblies of the 48 TRAINING individuals for
danbing-tk's RPGG (test/validation individuals are deliberately excluded, so
danbing-tk never sees a held-out person's alleles).

Writes data/Benchmark/assemblies/{sample}.{0,1}.fa (+ .fai), the layout
danbing-tk's GoodPanGenomeGraph pipeline expects, and assemblies.tsv listing
the source of every haplotype.

Sources, chosen to match the vamos call set each individual came from:
  HPRC phase 2  -> HPRC release 2 hifiasm assemblies (trio: pat = 0, mat = 1;
                   Hi-C phased: hap1 = 0, hap2 = 1)
  HPRC phase 1  -> HPRC year-1 assemblies (paternal = 0, maternal = 1)
  HGSVC phase 3 -> HGSVC3 Verkko assemblies (hap1 = 0, hap2 = 1)
  HGSVC phase 2 -> HGSVC2 phased CCS assemblies (h1 = 0, h2 = 1)

    python src/VNTR_Benchmark/fetch_assemblies.py
"""
import csv
import re
import gzip
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "VNTR_AbsoluteGenotyping"))

import polars as pl  # noqa: E402
from pipeline.paths import load_config  # noqa: E402
from pipeline.run.manifest import load_samples  # noqa: E402
from pipeline.run.remote import fetch_parallel  # noqa: E402

OUT = ROOT / "data" / "Benchmark" / "assemblies"
S3 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/"
EBI = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
HPRC_R2_INDEX = ("https://raw.githubusercontent.com/human-pangenomics/"
                 "hprc_intermediate_assembly/main/data_tables/assemblies_release2_v1.0.index.csv")
HGSVC3_BATCHES = ["20230818_verkko_batch1", "20230927_verkko_batch2", "20240201_verkko_batch3"]
HGSVC2_DIR = {"NA12878": "20200628_HHU_assembly-results_CCS_v12", "HG03125": "20200717_HHU_assembly-results_CCS_v12"}
HGSVC2_TAG = {"NA12878": "giab", "HG03125": "hgsvc"}


def hgsvc3_batches() -> dict[str, str]:
    """sample -> Verkko batch folder, read from the HGSVC3 FTP listings."""
    out = {}
    for b in HGSVC3_BATCHES:
        with urllib.request.urlopen(f"{EBI}HGSVC3/working/{b}/assemblies/") as r:  # noqa: S310
            for sid in re.findall(r'href="([A-Z]{2}[0-9]+)/"', r.read().decode()):
                out[sid] = b
    return out


def sources(sample: str, cohort: str, r2: dict, hgsvc3_batch: dict) -> list[str]:
    if cohort == "HPRC-phase2":
        return [r2[(sample, "pat")], r2[(sample, "mat")]]
    if cohort == "HPRC-phase1":
        d = f"{S3}working/HPRC_PLUS/{sample}/assemblies/year1_f1_assembly_v2_genbank/"
        return [f"{d}{sample}.{h}.f1_assembly_v2_genbank.fa.gz" for h in ("paternal", "maternal")]
    if cohort == "HGSVC-phase3":
        if sample not in hgsvc3_batch:
            raise ValueError(f"{sample} not in any HGSVC3 Verkko batch")
        d = f"{EBI}HGSVC3/working/{hgsvc3_batch[sample]}/assemblies/{sample}/"
        return [f"{d}{sample}.vrk-ps-sseq.asm-hap{h}.fasta.gz" for h in (1, 2)]
    if cohort == "HGSVC-phase2":
        d = f"{EBI}HGSVC2/release/v1.0/assemblies/{HGSVC2_DIR[sample]}/assemblies/phased/"
        return [f"{d}v12_{sample}_{HGSVC2_TAG[sample]}_pbsq2-ccs_1000-pereg.h{h}-un.racon-p2.fasta"
                for h in (1, 2)]
    raise ValueError(f"no assembly source for {sample} ({cohort})")


def main() -> None:
    cfg = load_config(None)
    split = pl.read_parquet(cfg.path("outputs.model_dir") / "sample_split.parquet")
    train = set(split.filter(pl.col("fold") == "train")["sample_id"])
    with urllib.request.urlopen(HPRC_R2_INDEX) as r:  # noqa: S310
        rows = list(csv.DictReader(r.read().decode().splitlines()))
    r2 = {}
    for row in rows:
        # trio-phased assemblies are named _pat_/_mat_, Hi-C-phased ones _hap1_/_hap2_
        name = row["assembly_name"]
        hap = ("pat" if "_pat_" in name or "_hap1_" in name else
               "mat" if "_mat_" in name or "_hap2_" in name else None)
        if hap:
            r2[(row["sample_id"], hap)] = row["assembly"].replace("s3://human-pangenomics/", S3)

    OUT.mkdir(parents=True, exist_ok=True)
    batch = hgsvc3_batches()
    plan, skipped = [], []
    for s in load_samples(cfg):
        if s["sample_id"] in train:
            cohort = s["vamos_column"].split("__")[0]
            try:
                urls = sources(s["sample_id"], cohort, r2, batch)
            except (ValueError, KeyError):
                skipped.append((s["sample_id"], cohort))
                continue
            plan += [(s["sample_id"], h, url) for h, url in enumerate(urls)]
    if skipped:
        print(f"no public assembly source, excluded from the RPGG: {skipped}", flush=True)
    with open(OUT / "assemblies.tsv", "w") as fh:
        fh.write("sample\thap\turl\n")
        fh.writelines(f"{a}\t{b}\t{c}\n" for a, b, c in plan)
    print(f"{len(plan)} haplotypes for {len(plan) // 2} training individuals", flush=True)

    for sid, hap, url in plan:
        fa = OUT / f"{sid}.{hap}.fa"
        if fa.exists() and (OUT / f"{sid}.{hap}.fa.fai").exists():
            continue
        raw = OUT / f"{sid}.{hap}.download{'.gz' if url.endswith('.gz') else ''}"
        fetch_parallel(url, raw, n_streams=8)
        tmp = fa.with_suffix(".fa.part")
        if url.endswith(".gz"):
            with gzip.open(raw, "rb") as src, open(tmp, "wb") as dst:
                shutil.copyfileobj(src, dst, 16 << 20)
            raw.unlink()
        else:
            raw.rename(tmp)
        tmp.rename(fa)
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
        print(f"[{sid}.{hap}] {fa.stat().st_size / 1e9:.2f} GB  <- {url}", flush=True)


if __name__ == "__main__":
    main()
