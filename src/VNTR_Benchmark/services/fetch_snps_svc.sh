#!/bin/bash
cd /data/GitHub_Repos/VNTR_Genotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
exec python -u src/VNTR_Benchmark/fetch_phased_snps.py --jobs 4 --use-2020 chr6 chr9 >> src/VNTR_Benchmark/work/logs/fetch_snps.log 2>&1
