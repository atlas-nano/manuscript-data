#!/bin/bash
# Two replicate campaigns for the alpha-sigma manuscript, packed onto one CPU node.
#
#   A  NVE shield-cutoff replicate  : 3 seeds x {5 A, 8 A}, 64 SPC-FQ waters, 200 ps, dt 1 fs
#   B  tip4pfq pair-scope replicate : 3 seeds x {intra, all}, 256 TIP4P-FQ, 10 ps, dt 0.25 fs
#
# Core budget on a 256-CPU node:  A 6 x (8 x c2) = 96   +   B 6 x (12 x c2) = 144   = 240.
# --exact so slurm packs each step onto its own cores; NEVER --overlap (that shares them).
set -u
cd "$(dirname "$0")"
JID=${JID:?set JID}
LMP=$HOME/codes/lammps/bin/lmp_cpu
SEEDS="731 4127 9013"

module load cpu 2>/dev/null
module unload craype-accel-nvidia80 gpu 2>/dev/null
export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1

R(){ srun --jobid=$JID --exact -N1 "$@"; }

echo "=== stage 1: three independent starting states for campaign A ==="
for S in $SEEDS; do
  R -n8 -c2 --cpu-bind=cores $LMP -in in.seed_eq -var seed $S -log log.eqS.$S > out.eqS.$S 2>&1 &
  sleep 2
done
wait
for S in $SEEDS; do
  [ -s eqS_$S.data ] || { echo "FAIL: eqS_$S.data missing"; echo FAIL > PM_MARKER; exit 1; }
done
echo "stage 1 ok"

echo "=== stage 2: 12 production steps ==="
for S in $SEEDS; do
  for R5 in 5.0 8.0; do
    RT=$(echo $R5 | tr -d '.')
    R -n8 -c2 --cpu-bind=cores $LMP -in in.nve_seed -var seed $S -var rsh $R5 \
        -var tag A_s${S}_r${RT} -log log.A.s${S}_r${RT} > out.A.s${S}_r${RT} 2>&1 &
    sleep 1
  done
  for SC in base all; do
    R -n12 -c2 --cpu-bind=cores $LMP -in in.mg_$SC -var seed $S \
        -log log.B.s${S}_$SC > out.B.s${S}_$SC 2>&1 &
    sleep 1
  done
done
wait

# Markers must be EARNED: `wait` returns 0 even if every srun died instantly, so check
# that the outputs actually carry the numbers we came for.
{
  echo "=== $(date -u +%FT%TZ) ==="
  ok=0; tot=0
  for S in $SEEDS; do
    for RT in 50 80; do
      tot=$((tot+1)); f=out.A.s${S}_r${RT}
      n=$(grep -cE '^ *[0-9]+ +[0-9.]+ +-?[0-9.]+ +-?[0-9.]+ +-?[0-9.]+' $f 2>/dev/null || echo 0)
      if [ "$n" -ge 1900 ]; then ok=$((ok+1)); st=OK; else st=FAIL; fi
      echo "  A s$S r$RT : $st rows=$n"
    done
    for SC in base all; do
      tot=$((tot+1)); f=out.B.s${S}_$SC
      m=$(grep -c "LIQUID_TAIL" $f 2>/dev/null || echo 0)
      if [ "$m" -ge 1 ]; then ok=$((ok+1)); st=OK; else st=FAIL; fi
      echo "  B s$S $SC : $st tail=$m  $(grep -h LIQUID_TAIL $f 2>/dev/null)"
    done
  done
  echo "complete: $ok / $tot"
  [ "$ok" -eq "$tot" ] && echo PM_OK || echo PM_INCOMPLETE
} > PM_STATUS 2>&1
cp PM_STATUS PM_MARKER
