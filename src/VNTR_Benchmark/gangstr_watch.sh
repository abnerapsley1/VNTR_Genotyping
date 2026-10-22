#!/bin/bash
# Run GangSTR on every TEST individual (from the current sample split) as soon as
# a verified CRAM for them is in data/Benchmark/cram. Exits when all are done.
cd "$(dirname "$0")/../.."
source /data/miniconda3/etc/profile.d/conda.sh && conda activate vntr_analysis
while true; do
  TEST=$(python -c "import polars as pl; print(' '.join(sorted(pl.read_parquet('src/VNTR_AbsoluteGenotyping/work/model/sample_split.parquet').filter(pl.col('fold')=='test')['sample_id'])))")
  left=0
  for s in $TEST; do
    [ -s src/VNTR_Benchmark/work/gangstr/calls/$s.tsv ] && continue
    left=$((left+1))
    [ -f data/Benchmark/cram/$s.cram ] && [ -f data/Benchmark/cram/$s.cram.crai ] && python -u src/VNTR_Benchmark/run_gangstr.py $s --jobs 8
  done
  [ $left -eq 0 ] && { echo "all GangSTR runs done ($(echo $TEST | wc -w) test individuals)"; exit 0; }
  sleep 120
done
