#!/usr/bin/env bash
# F3 scaling re-measurement of the manuscript's Section 7 benchmark.
#
# 13,216-atom two-electrode Cu cell (in.cucell_bench, 500 steps) on Perlmutter CPU nodes.
# Ladder is deliberately split so intra-node and inter-node scaling can be told apart:
#   16, 32, 64, 128 ranks all on ONE node (128 physical cores)
#   256 on two nodes, 512 on four   -> the first points that cross Slingshot
#
# -c2 pins one physical core per rank (128 phys / 256 logical per node).
set -u
JID=57633318
LMP=$HOME/codes/lammps/bin/lmp_cpu
cd /pscratch/sd/t/tpascal/alphasigma/f3_scaling || exit 1

module load cpu                                  >/dev/null 2>&1
module unload craype-accel-nvidia80 gpu          >/dev/null 2>&1
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1

echo "host=$(hostname)  jobid=$JID  binary=$(readlink -f $LMP)"
echo "start $(date)"

# "nodes:ranks"
IONF=${IONF:-off}
for pair in 1:16 1:32 1:64 1:128 2:256 4:512; do
  N=${pair%%:*}; NP=${pair##*:}
  LOG=log.np${NP}_N${N}.${IONF}
  echo "=== N=$N np=$NP ionfield=$IONF -> $LOG ==="
  srun --jobid=$JID -N$N -n$NP -c2 --cpu-bind=cores \
       $LMP -in in.cucell_bench_f3 -var NB 500 -var GRID none -var IONF $IONF -log $LOG > out.np${NP}_N${N}.${IONF} 2>&1
  rc=$?
  perf=$(grep -h "Performance:" $LOG 2>/dev/null | tail -1)
  loop=$(grep -h "^Loop time" $LOG 2>/dev/null | tail -1)
  echo "    rc=$rc"
  echo "    $loop"
  echo "    $perf"
done

echo "end $(date)"
echo "DONE_F3_SCALING"
