#!/bin/bash
#SBATCH -J asg_equil
#SBATCH -A csd799
#SBATCH -p debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# D3 GPU benchmark, step 1: equilibrate the 512-molecule native-water base box on CPU and record
# the CPU binary's style list. Writes base.data. Exit 1 unless base.data exists AND log says EQUIL_DONE.
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
LMP="$lmp_exec"
SO=$HOME/codes/lammps/install/cpu/lib64/liblammps_expanse_cpu.so.0
echo "exe: $(readlink -f $LMP)"
echo "ldd not-found (exe): $(ldd "$(readlink -f $LMP)" | grep -c 'not found')"
echo "ldd not-found (.so): $(ldd $SO | grep -c 'not found')"
ldd $SO | grep 'not found'
srun --mpi=pmi2 -n 1 -c 1 "$LMP" -h > styles_cpu.txt 2>&1
echo "style check rc=$? lines=$(wc -l < styles_cpu.txt) qeq/sam=$(grep -c 'qeq/sam' styles_cpu.txt) pppm/samqeq=$(grep -c 'pppm/samqeq' styles_cpu.txt) shake=$(grep -c 'shake' styles_cpu.txt)"
srun --mpi=pmi2 -n 32 -c 1 --cpu-bind=cores "$LMP" -in in.equil -var NEQ 10000 -log log.equil > out.equil 2>&1
rc=$?
echo "equil rc=$rc"
if [ $rc -eq 0 ] && grep -q '^EQUIL_DONE' log.equil && [ -s base.data ] && grep -q ' atoms$' base.data; then
  echo DONE > EQUIL.DONE; exit 0
else
  echo "FAIL rc=$rc" > EQUIL.FAIL; exit 1
fi
