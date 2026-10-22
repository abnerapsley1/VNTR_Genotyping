#!/bin/bash
# spanning-read allele lengths: training people with local CRAMs (offsets), test (local), validation (remote)
cd /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
python -u -m pipeline.features.spanning_alleles --fold train --local-only --jobs 20
python -u -m pipeline.features.spanning_alleles --fold test --jobs 20
python -u -m pipeline.features.spanning_alleles --fold val --jobs 24
echo "spanning extraction done"
