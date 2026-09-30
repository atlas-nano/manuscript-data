#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -A m5400
#SBATCH -t 00:30:00
#SBATCH -J cucell_cond
#SBATCH -o sb_%j.out
# Inert (published deck) and conductor-mode decks at 16/32/64/128 ranks, one at a time on one node.
LMP=$HOME/codes/lammps/bin/lmp_cpu
module load cpu; module unload craype-accel-nvidia80 gpu; export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR; md5sum $LMP > binary.md5
# smoke: 5 steps of each deck at 16 ranks; stop if either fails
for arm in conductor inert; do
  d=runs/smoke_$arm; mkdir -p $d; cp data.cell_L0*.data cucell.param CuZn.eam.alloy ffield.RexPoN in.cucell_$arm $d/
  ( cd $d && srun -N1 -n16 -c2 --cpu-bind=cores $LMP -in in.cucell_$arm -var NB 5 -var GRID none -log log.lammps -screen none > run.out 2>&1 ) || { echo "SMOKE FAIL $arm" > all_done.txt; exit 1; }
done
for arm in conductor inert; do
  for np in 16 32 64 128; do
    d=runs/${arm}_np${np}; mkdir -p $d; cp data.cell_L0*.data cucell.param CuZn.eam.alloy ffield.RexPoN in.cucell_$arm $d/
    ( cd $d && srun -N1 -n$np -c2 --cpu-bind=cores $LMP -in in.cucell_$arm -var NB 500 -var GRID none -log log.lammps -screen none > run.out 2>&1; echo "rc=$?" > rc.txt )
  done
done
echo ALL_DONE > all_done.txt
