#!/bin/bash
#SBATCH -J asg_gate
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
# D3 GPU benchmark, step 2 (needs base.data from sb_equil_cpu.sh):
#   (a) V100 binary: ldd + style list inside the job (mpirun, never srun)
#   (b) force gate: host vs -sf kk reciprocal forces at fixed charges (in.force_bisect V=4 and V=5)
#   (c) charge agreement on an identical short run: REP 2, 10 steps, host (8 ranks) vs kk devsolve (1 GPU)
#   (d) pilot timing, REP 4, 50 steps: kk devsolve on 1 GPU with and without the phase timers,
#       kk host-solve on 1 GPU, kk devsolve on 4 GPUs
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
LMPG=${LAMMPS_GPU_ROOT}/bin/lmp_expanse_gpu_v100
SO=${LAMMPS_GPU_ROOT}/lib64/liblammps_expanse_gpu_v100.so.0
echo "host $(hostname)  CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES  SLURM_GPUS=$SLURM_GPUS"
nvidia-smi --query-gpu=name,memory.total --format=csv
echo "ldd not-found (exe): $(ldd $LMPG | grep -c 'not found')"
echo "ldd not-found (.so): $(ldd $SO | grep -c 'not found')"
ldd $SO | grep 'not found'
mpirun -n 1 $LMPG -h > styles_gpu.txt 2>&1
echo "style check rc=$? lines=$(wc -l < styles_gpu.txt) qeq/sam/kk=$(grep -c 'qeq/sam/kk' styles_gpu.txt) pppm/samqeq/kk=$(grep -c 'pppm/samqeq/kk' styles_gpu.txt) shake/kk=$(grep -c 'shake/kk' styles_gpu.txt) coul/shield/intra/kk=$(grep -c 'coul/shield/intra/kk' styles_gpu.txt)"
KK1="-k on g 1 -sf kk -pk kokkos newton on neigh half"
KK4="-k on g 4 -sf kk -pk kokkos newton on neigh half"

# (b) force gate
for V in 4 5; do
  mpirun -n 1 $LMPG -in in.force_bisect -var V $V -var TAG gh$V -log log.fgate_h$V > out.fgate_h$V 2>&1; echo "fgate host V$V rc=$?"
  mpirun -n 1 $LMPG $KK1 -in in.force_bisect -var V $V -var TAG gk$V -log log.fgate_k$V > out.fgate_k$V 2>&1; echo "fgate kk   V$V rc=$?"
  python3 compare_dumps.py fb_gh$V.dump fb_gk$V.dump --faces 24.8
done

# (c) charge agreement, identical 10-step run
mpirun -n 8 --bind-to core $LMPG -in in.bench -var REP 2 -var NSTEP 10 -var DUMP 1 -var TAG agree_host -log log.agree_host > out.agree_host 2>&1; echo "agree host rc=$?"
mpirun -n 1 $LMPG $KK1 -in in.bench -var REP 2 -var NSTEP 10 -var DEV 1 -var DUMP 1 -var TAG agree_dev -log log.agree_dev > out.agree_dev 2>&1; echo "agree dev rc=$?"
python3 compare_dumps.py q_agree_host.dump q_agree_dev.dump
grep -h 'samqeq/kk S[789]\|stored-H' log.agree_dev

# (d) pilot timing
mpirun -n 1 $LMPG $KK1 -in in.bench -var REP 4 -var NSTEP 50 -var DEV 1 -log log.pilot_dev1 > out.pilot_dev1 2>&1; echo "pilot dev1 rc=$?"
SAMQEQ_KK_TIME=1 SAMQEQ_HOST_TIME=1 mpirun -n 1 -x SAMQEQ_KK_TIME -x SAMQEQ_HOST_TIME $LMPG $KK1 -in in.bench -var REP 4 -var NSTEP 50 -var DEV 1 -log log.pilot_dev1_timers > out.pilot_dev1_timers 2>&1; echo "pilot dev1 timers rc=$?"
mpirun -n 1 $LMPG $KK1 -in in.bench -var REP 4 -var NSTEP 50 -var DEV 0 -log log.pilot_kkhost1 > out.pilot_kkhost1 2>&1; echo "pilot kkhost1 rc=$?"
mpirun -n 4 $LMPG $KK4 -in in.bench -var REP 4 -var NSTEP 50 -var DEV 1 -log log.pilot_dev4 > out.pilot_dev4 2>&1; echo "pilot dev4 rc=$?"
for f in log.pilot_dev1 log.pilot_dev1_timers log.pilot_kkhost1 log.pilot_dev4; do
  echo "== $f"; grep -h 'Loop time\|Performance:\|^Modify\|^Kspace\|^Pair\|^Comm\|BENCH_DONE\|stored-H\|WARNING' $f
done
grep -h -A8 'HOST TIMING' log.pilot_dev1_timers; grep -h -A7 'TIMING over' log.pilot_dev1_timers

n_ok=$(grep -l '^BENCH_DONE' log.agree_host log.agree_dev log.pilot_dev1 log.pilot_dev1_timers log.pilot_kkhost1 log.pilot_dev4 2>/dev/null | wc -l)
if [ -s fb_gk5.dump ] && [ "$n_ok" -eq 6 ]; then echo DONE > GATE.DONE; else echo "FAIL n_ok=$n_ok" > GATE.FAIL; fi
