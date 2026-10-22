"""vntrlen command line.

  vntrlen genotype --reads a.cram b.cram --out results/          (shipped GRCh38 models, default BED)
  vntrlen genotype --samples sheet.tsv --bed my.bed --reference ref.fa --out results/
  vntrlen build    --samples sheet.tsv --bed my.bed --truth truth.tsv --name my_set --out build_report/
  vntrlen models   list [--reference ref.fa] [--region chr1:1-2000000] [--tsv out.tsv]
  vntrlen models   install [--bundle vntrlen-shipped-GRCh38.tar.gz]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .checks import InputError


def _samples(args):
    from . import samples
    if bool(args.reads) == bool(args.samples):
        raise SystemExit("give either --reads FILE ... or --samples SHEET.tsv")
    return samples.from_paths(args.reads) if args.reads else samples.from_sheet(args.samples)


def _common(p, bed_required=False):
    p.add_argument("--reference", help="reference FASTA, or a registry name (default: GRCh38_1kGP, "
                                       "downloaded on first use)")
    p.add_argument("--bed", required=bed_required, help="repeat BED: chrom, start, end, motif length, motifs[, name]"
                   + ("" if bed_required else " (default: every repeat in the shipped models)"))
    g = p.add_argument_group("samples (short reads only)")
    g.add_argument("--reads", nargs="+", help="BAM/CRAM files")
    g.add_argument("--samples", help="sample sheet (sample_id + alignment, or fastq1[/fastq2]; optional family_id)")
    p.add_argument("--out", required=True, type=Path, help="output directory")
    p.add_argument("--workers", type=int, default=8, help="parallel processes for extraction (default 8)")
    p.add_argument("--threads", type=int, default=8, help="threads for alignment of FASTQ input (default 8)")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="vntrlen", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", action="version", version=f"vntrlen {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("genotype", help="combined repeat length per sample and repeat")
    _common(g)
    g.add_argument("--model-set", help="use only this model set (default: newest user set, else shipped)")

    b = sub.add_parser("build", help="build models for your own repeats from samples with ground truth")
    _common(b, bed_required=True)
    b.add_argument("--truth", required=True, help="ground truth: sample_id, locus_id (or chrom/start/end), "
                                                  "length, unit (units|bp)[, allele1, allele2]")
    b.add_argument("--annotation", help="optional per-repeat annotation overriding the TRCompDB lookup")
    b.add_argument("--name", required=True, help="name for the new model set")
    b.add_argument("--retrain", choices=["auto", "never", "always"], default="auto",
                   help="retrain the correction models (auto: with >= 50 people and >= 1,000 repeats)")

    m = sub.add_parser("models", help="list or install models")
    msub = m.add_subparsers(dest="mcmd", required=True)
    ml = msub.add_parser("list", help="every reference-repeat model pair available")
    ml.add_argument("--reference", help="only models for this reference (FASTA or registry name)")
    ml.add_argument("--region", help="chrom or chrom:start-end")
    ml.add_argument("--tsv", type=Path, help="write the full table here")
    mi = msub.add_parser("install", help="install the shipped GRCh38 model bundle")
    mi.add_argument("--bundle", help="local bundle .tar.gz (default: download)")

    args = ap.parse_args(argv)
    try:
        if args.cmd == "genotype":
            from . import genotype
            genotype.run(args.reference, args.bed, _samples(args), args.out, args.model_set,
                         args.workers, args.threads)
        elif args.cmd == "build":
            from . import build
            build.run(args.reference, args.bed, _samples(args), args.truth, args.out, args.name,
                      args.annotation, args.retrain, args.workers, args.threads)
        elif args.mcmd == "install":
            from . import download
            download.install_bundle(args.bundle)
        else:
            from . import checks, download, store
            fp = checks.reference_fingerprint(download.resolve_reference(args.reference)) if args.reference else None
            df = store.query(fp, args.region)
            if args.tsv:
                df.write_csv(args.tsv, separator="\t", float_precision=4)
                print(f"wrote {df.height:,} rows to {args.tsv}")
            summ = df.group_by("reference_fingerprint", "model_set", "source").len().sort("model_set") \
                if df.height else df
            print(summ)
    except InputError as exc:
        print(f"vntrlen: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
