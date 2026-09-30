#!/bin/bash
#SBATCH -J asg_fgate_fix
#SBATCH -A sdp119
#SBATCH -p gpu-debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --gpus=2
#SBATCH --mem=90G
#SBATCH -t 00:20:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Post-fix verification of the pppm/samqeq/kk grid-comm flag fix (2026-09-16):
#   host vs -sf kk forces at fixed charges (in.force_bisect V=5 and V=1) on 1 rank/1 GPU and 2 ranks/2 GPUs,
#   and host vs kk devsolve charges after an identical 10-step run (pre-fix they drifted by 3e-3 e).
# FGATE_FIX.DONE is written only if every force comparison is below 1e-8 eV/A; otherwise FGATE_FIX.FAIL.
set +u
source ~/scripts/lammps_cluster_setup_files/expanse.sh
cd "$SLURM_SUBMIT_DIR"
rm -f FGATE_FIX.DONE FGATE_FIX.FAIL
LMPG=${LAMMPS_GPU_ROOT}/bin/lmp_expanse_gpu_v100
SO=${LAMMPS_GPU_ROOT}/lib64/liblammps_expanse_gpu_v100.so.0
echo "host $(hostname)  CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
ls -la --time-style=+%F_%T $SO
echo "ldd not-found (.so): $(ldd $SO | grep -c 'not found')"
KK1="-k on g 1 -sf kk -pk kokkos newton on neigh half"
KK2="-k on g 2 -sf kk -pk kokkos newton on neigh half"
: > fgate_fix_compare.txt
for V in 5 1; do
  for NP in 1 2; do
    KK=$KK1; [ $NP = 2 ] && KK=$KK2
    mpirun -n $NP $LMPG -in in.force_bisect -var V $V -var TAG fh${V}n$NP -log log.ffix_h${V}n$NP > out.ffix_h${V}n$NP 2>&1; echo "host V$V np$NP rc=$?"
    mpirun -n $NP $LMPG $KK -in in.force_bisect -var V $V -var TAG fk${V}n$NP -log log.ffix_k${V}n$NP > out.ffix_k${V}n$NP 2>&1; echo "kk   V$V np$NP rc=$?"
    grep -m1 'PPPM Kokkos initialization' out.ffix_k${V}n$NP || echo "WARNING: no Kokkos PPPM banner for V$V np$NP"
    echo "== V$V np$NP" >> fgate_fix_compare.txt
    python3 compare_dumps.py fb_fh${V}n$NP.dump fb_fk${V}n$NP.dump --faces 24.8 >> fgate_fix_compare.txt 2>&1
  done
done
cat fgate_fix_compare.txt

mpirun -n 2 --bind-to core $LMPG -in in.bench -var REP 2 -var NSTEP 10 -var DUMP 1 -var TAG fix_agree_host -log log.fix_agree_host > out.fix_agree_host 2>&1; echo "agree host rc=$?"
mpirun -n 1 $LMPG $KK1 -in in.bench -var REP 2 -var NSTEP 10 -var DEV 1 -var DUMP 1 -var TAG fix_agree_dev -log log.fix_agree_dev > out.fix_agree_dev 2>&1; echo "agree dev rc=$?"
python3 compare_dumps.py q_fix_agree_host.dump q_fix_agree_dev.dump | tee fgate_fix_agree.txt

python3 - <<'EOF'
import re, sys
t = open('fgate_fix_compare.txt').read()
blocks = re.findall(r'== (V\d np\d)\n(step .*)', t)
bad = []
if len(blocks) != 4:
    bad.append(f'expected 4 comparisons, found {len(blocks)}')
for tag, line in blocks:
    m = re.search(r'fx ([0-9.e+-]+)\s+fy ([0-9.e+-]+)\s+fz ([0-9.e+-]+)', line)
    if not m:
        bad.append(f'{tag}: unparsed'); continue
    worst = max(float(x) for x in m.groups())
    print(f'{tag}: max force diff {worst:.3e}')
    if worst > 1e-8:
        bad.append(f'{tag}: {worst:.3e}')
open('FGATE_FIX.FAIL' if bad else 'FGATE_FIX.DONE', 'w').write('\n'.join(bad) + '\n')
print('FAIL' if bad else 'PASS', bad)
EOF
