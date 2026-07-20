#!/bin/bash
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -J na_conductor
#SBATCH -A csd799
#SBATCH --export=ALL
#SBATCH --mem=32G
#Wall clock limit:
#SBATCH --time=12:00:00
#
##SBATCH --mail-type=ALL
#
module purge
module load cpu/0.15.4  gcc/10.2.0  mvapich2/2.3.6
module load amber/20.21
export nprocs=$SLURM_NTASKS_PER_NODE
export temp_dir=/expanse/lustre/scratch/xiruan/temp_project/$SLURM_JOB_ID/
echo "Running in ${temp_dir}"
mkdir $temp_dir
cp *.mol2 $temp_dir
cd ${temp_dir}
~/ATLAS-toolkit/scripts/amberAutoType.pl -m resp_full.mol2 -s bonded_forces -c -4 -t gaff
cp * $SLURM_SUBMIT_DIR

rm -r $temp_dir