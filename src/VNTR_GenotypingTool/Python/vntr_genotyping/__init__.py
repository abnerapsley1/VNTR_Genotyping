"""
vntr_genotyping — VNTR read-count and density-ratio genotyping tool.

Quick start
-----------
>>> import vntr_genotyping as vg
>>> df = vg.count_vntrs(
...     "sample.cram",
...     default=True,
...     reference="hg38.fa",
...     output_csv="results.csv",
... )
"""

from ._counting import count_vntrs
from ._io import DEFAULT_BED, DEFAULT_GTF, DEFAULT_PSL

__all__ = ["count_vntrs", "DEFAULT_BED", "DEFAULT_GTF", "DEFAULT_PSL"]
__version__ = "0.1.0"
