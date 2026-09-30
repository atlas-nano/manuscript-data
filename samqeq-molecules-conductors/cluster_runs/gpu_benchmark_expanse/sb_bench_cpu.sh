#!/bin/bash
#SBATCH -J asg_bcpu
#SBATCH -A csd799
#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -t 01:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# D3 GPU benchmark, CPU arm: one full Expanse CPU node (2x AMD EPYC 7742, 128 cores), 128 MPI ranks,
# host samQEq solve. Sizes and step counts come from the environment (ARMS="REP:NSTEP:ASPC ...").
# One srun step at a time (no packing), so no --mem slice is needed.
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
LMP="$lmp_exec"
SO=$HOME/codes/lammps/install/cpu/lib64/liblammps_expanse_cpu.so.0
echo "host $(hostname)  ldd not-found (.so): $(ldd $SO | grep -c 'not found')"
ARMS=${ARMS:-"2:500:0 4:200:0 6:100:0 8:50:0 2:500:1 4:200:1 6:100:1 8:50:1"}
n_ok=0; n_arm=0
for a in $ARMS; do
  IFS=: read REP NSTEP ASPC <<< "$a"
  tag=cpu128_r${REP}_a${ASPC}
  n_arm=$((n_arm+1))
  SAMQEQ_HOST_TIME=1 srun --mpi=pmi2 -n 128 -c 1 --cpu-bind=cores "$LMP" -in in.bench \
      -var REP $REP -var NSTEP $NSTEP -var ASPC $ASPC -log log.$tag > out.$tag 2>&1
  echo "$tag rc=$?"
  grep -h 'Loop time\|Performance:' log.$tag
  grep -q '^BENCH_DONE' log.$tag && grep -q '^Loop time' log.$tag && n_ok=$((n_ok+1))
done
if [ $n_ok -eq $n_arm ]; then echo DONE > BCPU.DONE; else echo "FAIL $n_ok/$n_arm" > BCPU.FAIL; fi
