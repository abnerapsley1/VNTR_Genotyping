#!/bin/bash
# After danbing-tk genotyping finishes: fit its length model on the training
# people, then rebuild the three-tool comparison table and figures.
set -uo pipefail
cd /data/GitHub_Repos/VNTR_Genotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
L=src/VNTR_Benchmark/work/logs/danbing_290.log
until grep -q "all danbing-tk alignments finished" $L; do sleep 600; done
echo "$(date) danbing-tk genotyping finished: $(tail -1 $L)"
python -u src/VNTR_Benchmark/fit_danbing_length.py || { echo "$(date) length model FAILED"; exit 1; }
python -u src/VNTR_Benchmark/compare_tools.py || { echo "$(date) comparison FAILED"; exit 1; }
python -u src/VNTR_Benchmark/make_compare_figures.py || { echo "$(date) figures FAILED"; exit 1; }
echo "$(date) three-tool comparison done"
