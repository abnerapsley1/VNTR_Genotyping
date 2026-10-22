#!/bin/bash
# systemd user service: download training-individual assemblies for the RPGG
cd /data/GitHub_Repos/VNTR_Genotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
export PATH=$PATH:/data/miniconda3/envs/vntr_benchmark/bin
exec python -u src/VNTR_Benchmark/fetch_assemblies.py
