#!/bin/bash
# systemd user service: danbing-tk genotyping (waits for RPGG_290, then test + training people)
cd /data/GitHub_Repos/VNTR_Genotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
exec python -u src/VNTR_Benchmark/run_danbing_290.py >> src/VNTR_Benchmark/work/logs/danbing_290.log 2>&1
