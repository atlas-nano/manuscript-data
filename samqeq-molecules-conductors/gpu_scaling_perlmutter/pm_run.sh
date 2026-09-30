#!/bin/bash
# samQEq alpha-sigma CPC-1: Perlmutter GPU benchmark runner.
#   pm_run.sh equil                      -> builds base.data (host, 1 rank)
#   pm_run.sh bench REP NNODE NSTEP TAG  -> timed run on NNODE*4 A100s, device solve
set -uo pipefail
cd "$HOME/alphasigma_pm_bench" || exit 2
LMP=/global/u1/t/tpascal/codes/lammps/bin/lmp_gpu
module load gpu cudatoolkit >/dev/null 2>&1
export OMP_NUM_THREADS=1

MODE=$1
if [ "$MODE" = equil ]; then
  export MPICH_GPU_SUPPORT_ENABLED=0
  srun --jobid="$JID" -n 1 -c 32 --gpus-per-node=1 --gpu-bind=none \
       "$LMP" -in in.equil -log log.equil > out.equil 2>&1
  rc=$?
  echo "equil rc=$rc  base.data=$(ls -l base.data 2>/dev/null | awk '{print $5}') bytes"
  grep -E "^ERROR|Total wall" out.equil log.equil 2>/dev/null | head -3
  exit $rc
fi

REP=$2; NNODE=$3; NSTEP=$4; TAG=$5
NG=$((NNODE*4))
# GPU-aware MPI is REQUIRED for multi-GPU KOKKOS here; --gpus-per-task=1 hangs in GTL IPC.
export MPICH_GPU_SUPPORT_ENABLED=1
srun --jobid="$JID" -N "$NNODE" -n "$NG" -c 32 --gpus-per-node=4 --gpu-bind=none \
     "$LMP" -k on g 4 -sf kk -pk kokkos gpu/aware on newton on neigh half \
     -in in.bench -log "log.$TAG" \
     -var REP "$REP" -var NSTEP "$NSTEP" -var DEV 1 -var ASPC 0 -var TAG "$TAG" \
     > "out.$TAG" 2>&1
rc=$?
echo "=== $TAG rc=$rc nodes=$NNODE gpus=$NG rep=$REP"
grep -E "BENCH_DONE" "out.$TAG" | head -1
# LAMMPS prints "Performance: ... timesteps/s"
grep -E "Performance:" "log.$TAG" | head -1
grep -E "^ERROR|Segmentation|not found" "out.$TAG" | head -3
