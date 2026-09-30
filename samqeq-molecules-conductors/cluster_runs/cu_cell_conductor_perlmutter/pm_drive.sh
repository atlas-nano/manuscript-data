#!/bin/bash
# Perlmutter: inert (published deck) and conductor-mode decks at 16/32/64/128 ranks, run one at a time on one
# interactive CPU node so the timings do not interfere. Usage: bash pm_drive.sh <jobid>
JID=$1; LMP=$HOME/codes/lammps/bin/lmp_cpu
module load cpu; module unload craype-accel-nvidia80 gpu; export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
cd "$(dirname "$0")"; md5sum $LMP > binary.md5
for arm in conductor inert; do
  for np in 16 32 64 128; do
    d=runs/${arm}_np${np}; mkdir -p $d; cp data.cell_L0*.data cucell.param CuZn.eam.alloy ffield.RexPoN in.cucell_$arm $d/
    ( cd $d && srun --jobid=$JID -N1 -n$np -c2 --cpu-bind=cores $LMP -in in.cucell_$arm -var NB 500 -log log.lammps -screen none > run.out 2>&1; echo "rc=$?" > rc.txt )
  done
done
echo ALL_DONE > all_done.txt
