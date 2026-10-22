"""
cli.py — one entrypoint for the whole pipeline.

    python -m pipeline.cli <command> [--config path]

Commands (usual order):
    ground-truth     parse the vamos VCF -> gt_long / gt_locus_summary / loci_targets
    covariates       build per-locus covariate table
    list-vamos       print the 416 vamos sample-column headers (for samples.tsv)
    extract          run feature extraction for all pending manifest samples
    extract-one      run feature extraction for a single sample  (--sample ID)
    split            build frozen locus / sample train-val-test folds
    assemble         join features + ground truth + covariates -> model_table
    calibrate        Phase 1: per-locus empirical calibration
    transfer         Phase 2: calibration-transfer model (novel loci)
    global           Phase 3: stacked quantile model (shipped predictor)
    per-allele       evaluate the opt-in per-allele estimate
    haplotypes       per-haplotype lengths: allele catalog + allele-gap evidence
    haplotypes-snp   per-haplotype lengths + phase from flanking SNP haplotypes (+ spanning reads)
    evaluate         score predictions by quadrant vs naive baseline
    all-model        split -> assemble -> calibrate -> transfer -> global -> per-allele -> evaluate
"""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> None:
    argv = argv if argv is not None else sys.argv[1:]
    ap = argparse.ArgumentParser(prog="pipeline.cli", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command")
    ap.add_argument("--config", default=None)
    ap.add_argument("--sample", default=None, help="extract-one: sample_id")
    ap.add_argument("--only", nargs="+", default=None, help="extract: restrict sample_ids")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-reads", action="store_true", help="haplotypes-snp: skip the spanning-read term")
    args = ap.parse_args(argv)
    c = args.command
    cfgarg = ["--config", args.config] if args.config else []

    if c == "ground-truth":
        from .ground_truth.parse_vamos import main as m
        m(cfgarg + ["--run"])
    elif c == "covariates":
        from .covariates.build_locus_covariates import main as m
        m(cfgarg)
    elif c == "list-vamos":
        from .ground_truth.parse_vamos import _iter_vcf_records
        from .paths import load_config
        cfg = load_config(args.config)
        for rec in _iter_vcf_records(cfg.path("inputs.vamos_vcf")):
            if rec[0] == "__SAMPLES__":
                print("\n".join(rec[1]))
                break
    elif c == "extract":
        from .run.run_all import main as m
        extra = (["--only"] + args.only if args.only else []) + \
                (["--limit", str(args.limit)] if args.limit else [])
        m(cfgarg + extra)
    elif c == "extract-one":
        from .run.process_sample import main as m
        if not args.sample:
            ap.error("extract-one needs --sample")
        m(cfgarg + ["--sample", args.sample])
    elif c == "split":
        from .model.splits import main as m
        m(cfgarg)
    elif c == "assemble":
        from .model.assemble import main as m
        m(cfgarg)
    elif c == "calibrate":
        from .model.per_locus_calibration import main as m
        m(cfgarg)
    elif c == "transfer":
        from .model.calibration_transfer import main as m
        m(cfgarg)
    elif c == "global":
        from .model.global_residual import main as m
        m(cfgarg)
    elif c == "per-allele":
        from .model.per_allele import main as m
        m(cfgarg)
    elif c == "haplotypes":
        from .model.haplotypes import main as m
        m(cfgarg)
    elif c == "haplotypes-snp":
        from .model.haplotypes_snp import main as m
        m(cfgarg + (["--no-reads"] if args.no_reads else []))
    elif c == "evaluate":
        from .model.evaluate import main as m
        m(cfgarg)
    elif c == "all-model":
        for step in ("split", "assemble", "calibrate", "transfer", "global",
                     "per-allele", "evaluate"):
            print(f"\n########## {step} ##########")
            main([step] + cfgarg)
    else:
        ap.error(f"unknown command: {c}")


if __name__ == "__main__":
    main()
