#!/bin/bash
#Job name:
#SBATCH --job-name=glyme_tetramer_32_4_small_box
#
#Partition:
#SBATCH --partition=shared
#
#Account
#SBATCH --account=csd799
#
#Nodes
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks-per-node=32
#
#Wall clock limit:
#SBATCH --time=48:00:00

project=na_conductor
prefix=lammps
name=glyme_tetramer_8_4_small_box
rtemp=298
press=1
# rep_id=__REPID__
# restart_id=__RESTARTID__
# i=__i__
curr_dir=$SLURM_SUBMIT_DIR
temp_dir=/expanse/lustre/scratch/$USER/temp_project/${project}/$SLURM_JOB_ID/

lmp_equil_file=in.${prefix}
lmp_data_file=data.${prefix}
# lmp_data_file=lammps.298K.6000000.restart # 0.5M 

module purge
# module load cpu/0.15.4 slurm gcc/10.2.0 openmpi cmake gsl intel-mkl amdfftw
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 cmake/3.21.4/teqow32 openmpi/4.1.3/oq3qvsv fftw/3.3.10/bmvqsbr python/3.8.12/7zdjza7 intel-mkl/2020.4.304/ghfk3mu gsl/2.7/wtlsmyy gnuplot/5.4.2/mfinpvw sdsc slurm

nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES / 2 ))
PARALLEL="mpirun -n $nprocs -mca btl vader,self"

# LMP="/expanse/lustre/projects/csd626/tpascal/codes/bin/lmp_expanse -screen none -var rtemp $rtemp -var press $press -var i $i -var rep_id $rep_id -var restart_id $restart_id"
# mkdir ${CURR_DIR}/results
LMP="/expanse/lustre/projects/csd626/tpascal/codes/bin/lmp_expanse -screen none -var rtemp $rtemp -var press $press"

mkdir -p $temp_dir/analysis ${curr_dir}/results ${curr_dir}/restart
cd $temp_dir
for i in $lmp_data_file $lmp_equil_file 
do
	cp ${curr_dir}/${i} ./
done

echo "LAMMPS dynamics of ${prefix} at ${rtemp}K"
echo "running in $temp_dir"
$PARALLEL $LMP -in ${lmp_equil_file} -log ${prefix}.${rtemp}K.equil.lammps.log
cp *.log *.lammpstrj ${curr_dir}/results
echo "LAMMPS dynamics of ${name} lammps at ${rtemp}K"
cp *.log *.lammpstrj ${curr_dir}/results
cp *.restart ${curr_dir}/restart



