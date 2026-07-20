#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=shared
#SBATCH -t 48:00:00
#SBATCH -J H2O_SHG
#SBATCH --account=csd816
#SBATCH --tasks-per-node=16
#cd $SCRATCH
#mkdir TiSe21E12_T300
#cd TiSe21E12_T300

#date

#cp $SLURM_SUBMIT_DIR/*xml* .

#OpenMP settings:
ulimit -s unlimited
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true
export SLURM_CPU_BIND="cores"

#run the application:
#srun  /global/homes/s/sjamnuch/Perlmutter/Siesta_LAPACK/siesta-v4.1.5/Obj/siesta <  input.fdf > $SLURM_SUBMIT_DIR/OUT
srun -n 16 /home/sjamnuch/Siesta/VG-RTTDDFT/Src/siesta < input.fdf > $SLURM_SUBMIT_DIR/OUT
#srun /expanse/lustre/projects/csd626/sjamnuch/Comet/home/Siesta/VG-TDDFT/VG-TDDFT-master/tddft/code/siesta-2.0.1/Src/siesta < input.fdf > $SLURM_SUBMIT_DIR/OUT
#cp *EIG $SLURM_SUBMIT_DIR
