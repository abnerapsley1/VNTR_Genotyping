#!/usr/bin/env python3
"""
count_reads.py — Legacy entry point. Delegates to vntr_genotyping.cli.

Prefer using the installed 'vntr-count' command (pip install ./Python)
or importing the package directly:

    import vntr_genotyping as vg
    df = vg.count_vntrs("sample.cram", default=True, gtf="annotation.gtf.gz")
"""

import sys
import os

# Allow running this script directly without installing the package,
# by adding the package directory to the path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vntr_genotyping.cli import main

if __name__ == "__main__":
    main()
