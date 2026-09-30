#!/bin/bash
#SBATCH -J asg_gaprobe
#SBATCH -A sdp119
#SBATCH -p gpu-debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus=4
#SBATCH --mem=180G
#SBATCH -t 00:25:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Scoping probe: is the linked OpenMPI/UCX GPU-aware on a V100 node, which transports can it use, and what
# does LAMMPS's GPU-aware path buy? No build. Writes GAPROBE.DONE only if all three runs finish.
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
LMPG=${LAMMPS_GPU_ROOT}/bin/lmp_expanse_gpu_v100
echo "== node $(hostname)"; nvidia-smi topo -m 2>&1 | head -12
echo "== kernel modules"; lsmod | grep -E "gdrdrv|nvidia_peermem|nv_peer_mem|nvidia_uvm" ; ls -l /dev/gdrdrv 2>&1
echo "== mpi"; which mpirun; ompi_info --parsable --all 2>&1 | grep -E "mpi_built_with_cuda_support:value|mca:pml:ucx:version|mca:btl:smcuda" | head -5
echo "== ucx"; ucx_info -v 2>&1 | head -2; ucx_info -d 2>&1 | grep -E "Transport: (cuda|gdr|rc|dc|ud|sysv|posix)" | sort | uniq -c
cat > mpix_probe.c <<'C'
#include <mpi.h>
#include <stdio.h>
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
#if defined(MPIX_CUDA_AWARE_SUPPORT)
  printf("compile-time MPIX_CUDA_AWARE_SUPPORT=%d runtime MPIX_Query_cuda_support()=%d\n", MPIX_CUDA_AWARE_SUPPORT, MPIX_Query_cuda_support());
#else
  printf("MPIX_CUDA_AWARE_SUPPORT not defined\n");
#endif
  MPI_Finalize(); return 0;
}
C
mpicc -o mpix_probe mpix_probe.c && mpirun -n 1 ./mpix_probe
KK="-k on g 4 -sf kk -pk kokkos newton on neigh half"
n=0
for mode in off on; do
  mpirun -n 4 $LMPG $KK gpu/aware $mode -in in.bench -var REP 6 -var NSTEP 50 -var ASPC 0 -var DEV 1 -log log.ga_$mode > out.ga_$mode 2>&1
  echo "gpu/aware $mode rc=$?"; grep -h "Performance:\|aware" log.ga_$mode out.ga_$mode | sort -u | head -4
  grep -q '^BENCH_DONE' log.ga_$mode && n=$((n+1))
done
UCX_TLS=self,sm,cuda_copy mpirun -x UCX_TLS -n 4 $LMPG $KK gpu/aware on -in in.bench -var REP 6 -var NSTEP 50 -var ASPC 0 -var DEV 1 -log log.ga_on_noipc > out.ga_on_noipc 2>&1
echo "gpu/aware on, UCX_TLS without cuda_ipc rc=$?"; grep -h "Performance:" log.ga_on_noipc | head -2
grep -q '^BENCH_DONE' log.ga_on_noipc && n=$((n+1))
[ $n -eq 3 ] && echo DONE > GAPROBE.DONE || echo "FAIL $n/3" > GAPROBE.FAIL
