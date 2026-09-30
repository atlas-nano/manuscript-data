#!/bin/bash
#SBATCH -J asg_bgpu
#SBATCH -A sdp119
#SBATCH -p gpu-debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --gpus=4
#SBATCH --mem=180G
#SBATCH -t 00:30:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# D3 GPU benchmark, GPU arms: Expanse V100 node (4x V100 32 GB), -sf kk with `devsolve on`.
# ARMS="NGPU:REP:NSTEP:ASPC:DEV ..." from the environment; launched with mpirun (never srun).
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
LMPG=${LAMMPS_GPU_ROOT}/bin/lmp_expanse_gpu_v100
SO=${LAMMPS_GPU_ROOT}/lib64/liblammps_expanse_gpu_v100.so.0
echo "host $(hostname)  CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES  ldd not-found (.so): $(ldd $SO | grep -c 'not found')"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
ARMS=${ARMS:?set ARMS}
n_ok=0; n_arm=0
for a in $ARMS; do
  IFS=: read NG REP NSTEP ASPC DEV <<< "$a"
  tag=gpu${NG}_r${REP}_a${ASPC}_d${DEV}
  n_arm=$((n_arm+1))
  ( while true; do nvidia-smi --query-gpu=index,memory.used --format=csv,noheader >> mem.$tag; sleep 5; done ) &
  MON=$!
  mpirun -n $NG $LMPG -k on g $NG -sf kk -pk kokkos newton on neigh half -in in.bench \
      -var REP $REP -var NSTEP $NSTEP -var ASPC $ASPC -var DEV $DEV -log log.$tag > out.$tag 2>&1
  echo "$tag rc=$?"
  kill $MON 2>/dev/null
  echo "   peak GPU memory (MiB): $(awk -F', ' '{gsub(/ MiB/,"",$2); if ($2+0>m) m=$2+0} END {print m}' mem.$tag)"
  grep -h 'Loop time\|Performance:\|stored-H\|ERROR' log.$tag out.$tag | sort -u
  grep -q '^BENCH_DONE' log.$tag && grep -q '^Loop time' log.$tag && n_ok=$((n_ok+1))
done
if [ $n_ok -eq $n_arm ]; then echo DONE > BGPU_${SLURM_JOB_ID}.DONE; else echo "FAIL $n_ok/$n_arm" > BGPU_${SLURM_JOB_ID}.FAIL; fi
