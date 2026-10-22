#!/bin/bash
# systemd user service: feature extraction for all pending manifest samples
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
exec python -u -m pipeline.cli extract
