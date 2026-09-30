#!/bin/bash
# Mobile-electrode conductor deck on one interactive CPU node: smoke (5 steps), throughput at 16/32/64/128 ranks
# (500 steps each), then a 5 ps stability run (20000 steps) at 64 ranks. Usage: bash pm_mobile.sh <jobid>
JID=$1; LMP=$HOME/codes/lammps/bin/lmp_cpu
module load cpu; module unload craype-accel-nvidia80 gpu; export MPICH_GPU_SUPPORT_ENABLED=0 OMP_NUM_THREADS=1
cd "$(dirname "$0")"; md5sum $LMP > binary_mobile.md5
IN="data.cell_L0_conductor.data cucell.param Mendelev_Cu2_2012.eam.fs.real ffield.RexPoN in.cucell_mobile"
run() { d=$1; np=$2; nb=$3; mkdir -p $d; cp $IN $d/; ( cd $d && srun --jobid=$JID -N1 -n$np -c2 --cpu-bind=cores $LMP -in in.cucell_mobile -var NB $nb -var GRID none -log log.lammps -screen none > run.out 2>&1; echo "rc=$?" > rc.txt ); }
run runs_mobile/smoke 16 5; grep -q "rc=0" runs_mobile/smoke/rc.txt || { echo "SMOKE FAIL" > mobile_done.txt; exit 1; }
for np in 16 32 64 128; do run runs_mobile/np$np $np 500; done
run runs_mobile/stab_np64 64 20000
echo DONE > mobile_done.txt
