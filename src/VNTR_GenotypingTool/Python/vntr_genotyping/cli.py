"""
cli.py — Command-line entry point for vntr_genotyping (vntr-count).
"""

import argparse
import sys

from ._io import DEFAULT_PSL
from ._counting import count_vntrs

_USAGE = """
Count reads overlapping VNTR regions in SAM/BAM/CRAM files.
Outputs density ratios (VNTR / gene) when --gtf is provided, otherwise raw read counts.

Region selection (at least one required):
    --default               Use all VNTRs in the built-in BED file
    --gene  GENE [GENE ...] Use VNTRs from the built-in BED matching these gene(s)
    --vntr  ID   [ID ...]   Use VNTRs from the built-in BED matching these VNTR name(s)
    -r / --regions FILE     Use a custom BED file (combinable with any of the above)

    --default and --gene/--vntr are mutually exclusive.
    --regions can be combined with any other option.

Examples:
    # All built-in VNTRs
    vntr-count --default -i s1.cram -g annotation.gtf.gz -T ref.fa -o out.csv

    # Subset by gene
    vntr-count --gene SLC6A3 BRCA1 -i s1.cram -g annotation.gtf.gz -T ref.fa -o out.csv

    # Subset by VNTR name
    vntr-count --vntr TR30 TR22 -i s1.cram -g annotation.gtf.gz -T ref.fa -o out.csv

    # Custom BED only
    vntr-count -r my_vntrs.bed -i s1.cram -g annotation.gtf.gz -T ref.fa -o out.csv

Output:
    CSV with columns: sample, metric, <vntr_1>, <vntr_2>, ...
    metric = 'density_ratio' when --gtf is supplied, 'read_count' otherwise.
"""


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Count reads overlapping VNTR regions in SAM/BAM/CRAM files. "
            "Outputs density ratios (VNTR / gene) when --gtf is provided, "
            "otherwise raw read counts."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_USAGE,
    )

    # --- Region selection ---
    region_group = parser.add_argument_group("region selection (at least one required)")
    region_group.add_argument(
        "--default", action="store_true",
        help="Use all VNTRs in the built-in BED file",
    )
    region_group.add_argument(
        "--gene", nargs="+", metavar="GENE",
        help="Use VNTRs from the built-in BED matching these gene name(s)",
    )
    region_group.add_argument(
        "--vntr", nargs="+", metavar="ID",
        help="Use VNTRs from the built-in BED matching these VNTR name(s)",
    )
    region_group.add_argument(
        "-r", "--regions", metavar="BED",
        help="Custom BED file (combinable with --default / --gene / --vntr)",
    )

    # --- Required inputs ---
    parser.add_argument(
        "-i", "--input", nargs="+", required=True, metavar="FILE",
        help="Input SAM/BAM/CRAM file(s)",
    )

    # --- Optional inputs ---
    parser.add_argument(
        "-g", "--gtf", default=None, metavar="GTF",
        help="Gene annotation GTF or GTF.gz (enables density ratio normalization)",
    )
    parser.add_argument(
        "-o", "--output", default="read_counts.csv", metavar="CSV",
        help="Output CSV file (default: read_counts.csv)",
    )
    parser.add_argument(
        "-T", "--reference", default=None, metavar="FASTA",
        help="Reference FASTA (required for CRAM without embedded reference)",
    )
    parser.add_argument(
        "--psl", default=DEFAULT_PSL, metavar="PSL",
        help=(
            "UCSC altSeqLiftOverPsl.txt(.gz) for alt-contig gene normalization. "
            f"Default: {DEFAULT_PSL}"
        ),
    )
    parser.add_argument(
        "--no-alt-contigs", action="store_true",
        help=(
            "Disable alt-contig gene normalization. Use this when your reference "
            "genome contains only primary chromosome sequences and no alt scaffolds."
        ),
    )
    parser.add_argument(
        "--norm-gene-types", nargs="+", default=["protein_coding"], metavar="TYPE",
        help=(
            "Gene biotype(s) eligible as the normalization gene when a VNTR has no "
            "overlapping gene. Default: protein_coding. "
            "Use 'any' to allow all biotypes. "
            "Examples: --norm-gene-types protein_coding lncRNA"
        ),
    )

    args = parser.parse_args()

    # Validate: at least one region source
    if not args.default and not args.gene and not args.vntr and not args.regions:
        parser.error(
            "No regions specified. Provide at least one of: "
            "--default, --gene, --vntr, or --regions."
        )

    # Validate: --default is incompatible with --gene/--vntr
    if args.default and (args.gene or args.vntr):
        parser.error(
            "--default cannot be combined with --gene or --vntr. "
            "Use --default to select all built-in regions, or use "
            "--gene/--vntr to select a subset."
        )

    count_vntrs(
        args.input,
        default=args.default,
        gene=args.gene,
        vntr=args.vntr,
        regions=args.regions,
        gtf=args.gtf,
        reference=args.reference,
        psl=args.psl,
        no_alt_contigs=args.no_alt_contigs,
        norm_gene_types=args.norm_gene_types,
        output_csv=args.output,
    )


if __name__ == "__main__":
    main()
