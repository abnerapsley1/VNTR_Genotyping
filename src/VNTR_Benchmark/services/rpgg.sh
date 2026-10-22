#!/bin/bash
# systemd user service: RPGG_290 build (batched lift, then joint annotation + graph)
exec /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_Benchmark/build_rpgg_290.sh >> /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_Benchmark/work/logs/build_rpgg_290.log 2>&1
