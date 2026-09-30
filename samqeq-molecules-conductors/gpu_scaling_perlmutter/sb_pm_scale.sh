#!/bin/bash
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -A m5400
#SBATCH -t 00:30:00
#SBATCH -o %x_%j.out
# CPC-1: samQEq multi-node GPU strong scaling on Perlmutter A100.
# Submit as: sbatch -N <nodes> -J pmscale_n<nodes> sb_pm_scale.sh
set -uo pipefail
cd "$HOME/alphasigma_pm_bench" || exit 2
LMP=/global/u1/t/tpascal/codes/lammps/bin/lmp_gpu
module load gpu cudatoolkit
export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1     # required: multi-GPU KOKKOS halos

NN=$SLURM_JOB_NUM_NODES
NG=$((NN*4))
# Guard: a binary that cannot resolve its libraries must not be timed.
MISS=$(ldd "$LMP" 2>/dev/null | grep -c "not found")
echo "GUARD nodes=$NN gpus=$NG missing_libs=$MISS"
[ "$MISS" -ne 0 ] && { echo "FAIL unresolved libraries"; exit 3; }

OK=0; TOT=0
for REP in 6 8; do
  TAG="rep${REP}_n${NN}"
  TOT=$((TOT+1))
  srun -N "$NN" -n "$NG" -c 32 --gpus-per-node=4 --gpu-bind=none \
       "$LMP" -k on g 4 -sf kk -pk kokkos gpu/aware on newton on neigh half \
       -in in.bench -log "log.$TAG" \
       -var REP "$REP" -var NSTEP 200 -var DEV 1 -var ASPC 0 -var TAG "$TAG" \
       > "out.$TAG" 2>&1
  P=$(grep -E "Performance:" "log.$TAG" | head -1)
  N=$(grep -oE "natoms=[0-9]+" "out.$TAG" | head -1)
  if [ -n "$P" ]; then
    OK=$((OK+1))
    echo "RESULT $TAG nodes=$NN gpus=$NG $N | $P"
  else
    echo "RESULT $TAG FAILED"
    grep -E "^ERROR|Segmentation|error" "out.$TAG" | head -2
  fi
done
# Markers are earned: only a run that produced a Performance line counts.
if [ "$OK" -eq "$TOT" ]; then echo "$OK/$TOT" > "SCALE_n${NN}.DONE"; else echo "$OK/$TOT" > "SCALE_n${NN}.FAIL"; fi
echo "ARMS_OK $OK/$TOT"
