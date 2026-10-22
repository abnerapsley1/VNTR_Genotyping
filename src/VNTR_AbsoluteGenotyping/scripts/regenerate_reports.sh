#!/bin/bash
# Regenerate every reported result from the current default model (work/model).
set -uo pipefail
R=/data/GitHub_Repos/VNTR_Genotyping
cd $R/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
export VNTR_GBM=xgboost
run() { echo "$(date) >> $*"; "$@" || echo "$(date) FAILED: $*"; }
run python -u -m pipeline.cli per-allele
run python -u scripts/feature_inventory.py
cp work/model/feature_inventory.csv docs/feature_inventory.csv
run python -u work/figures/make_bucket_statistics_csv.py
run python -u work/figures/make_figures_v3.py
run python -u work/figures/make_bp_density_figures.py
cd $R
run python -u src/VNTR_Benchmark/compare_tools.py
run python -u src/VNTR_Benchmark/make_compare_figures.py
run python -u src/VNTR_Benchmark/stats_ours_vs_danbing.py
run python -u src/VNTR_Benchmark/stats_supplement.py
echo "$(date) quick reports done"
cd $R/src/VNTR_AbsoluteGenotyping
run python -u -m pipeline.cli haplotypes
run python -u -m pipeline.cli haplotypes-snp
run python -u scripts/family_check.py
echo "$(date) all reports regenerated"
