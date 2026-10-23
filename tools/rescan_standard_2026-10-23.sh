#!/bin/zsh
# Run the 2026-10-23 standard-model rescans (food_cap=1) one after another, each with its own
# run.log, from the repo root of a checkout of branch food-cap-standard. Meant for torque:
#   nohup zsh tools/rescan_standard_2026-10-23.sh .venv/bin/python > rescan.log 2>&1 &! 
PY=${1:-python3}
for d in Resource-driving/Scans/2026-10-23_R1_mscale_bracket_standard \
         Activity-metrics/Scans/2026-10-23_pure_evo_joint_standard \
         Activity-metrics/Scans/2026-10-23_pure_evo_lut_only_standard \
         Lineage-co-evolution/Scans/2026-10-23_pure_evo_egene_on_winners_standard \
         Activity-metrics/Scans/2026-10-23_evo_focus_standard \
         Lineage-co-evolution/Scans/2026-10-23_3c_coevolution_substrate_standard \
         Activity-metrics/Scans/2026-10-23_causal_control_standard \
         Activity-metrics/Scans/2026-10-23_causal_control_v2_standard; do
  drv=$(ls Research/$d/*.py | head -1)
  echo "=== $(date '+%H:%M:%S') start $d"
  (cd Research/$d && N_WORKERS=30 $PY $(basename $drv) > run.log 2>&1); echo "=== $(date '+%H:%M:%S') exit $? $d"
done
echo "=== ALL DONE"
