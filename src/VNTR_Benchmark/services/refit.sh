#!/bin/bash
# systemd user service: 290-sample model refit (GPU) + one-time CPU check, after extraction
exec /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping/scripts/refit_290.sh >> /data/GitHub_Repos/VNTR_Genotyping/src/VNTR_AbsoluteGenotyping/work/logs/refit_290.log 2>&1
