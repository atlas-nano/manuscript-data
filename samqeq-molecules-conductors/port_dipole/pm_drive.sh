#!/bin/bash
# Perlmutter: 8 replicas x 16 ranks packed on one interactive CPU node. Usage: bash pm_drive.sh <jobid>
JID=$1; LMP=$HOME/codes/lammps/bin/lmp_cpu
module load cpu; module unload craype-accel-nvidia80 gpu; export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
cd "$(dirname "$0")"
md5sum $LMP > binary.md5
for d in runs/*/; do
  ( cd $d && srun --jobid=$JID --exact -N1 -n16 -c2 --cpu-bind=cores $LMP -in in.lammps -log log.lammps -screen none > run.out 2>&1; echo "rc=$?" > rc.txt ) &
  sleep 2
done
wait
echo ALL_DONE > all_done.txt
