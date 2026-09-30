#!/bin/bash
# A5 rerun (2026-09-16): tolerance scan on the NPT snapshot with the production pin, g_ewald and mesh.
set -u
cd "$(dirname "$0")"
LMP=/home/tpascal/codes/lammps/lammps-dev2/build/lmp_parallel
export TF_CPP_MIN_LOG_LEVEL=3
TOLS="1.0e-1 1.0e-2 1.0e-3 1.0e-4 1.0e-5 1.0e-6 1.0e-7 1.0e-8 1.0e-9 1.0e-10 1.0e-11"
for t in $TOLS; do
  mpirun -np 4 --bind-to none $LMP -in in.s8 -var start npt -var tol $t -log log.npt.$t > out.npt.$t 2>&1
done
: > scan_npt.txt
for t in $TOLS; do
  it=$(grep -o 'CG [0-9]* iters' out.npt.$t | head -1 | awk '{print $2}')
  rb=$(grep -o 'resid/b=[0-9.eE+-]*' out.npt.$t | head -1 | cut -d= -f2)
  mu=$(grep -o 'dipole=[0-9.]*' out.npt.$t | cut -d= -f2)
  ge=$(grep -m1 'G vector' out.npt.$t | awk '{print $NF}'); gr=$(grep -m1 'grid = ' out.npt.$t | awk '{print $3"x"$4"x"$5}')
  echo "$t ${it:-NA} ${mu:-NA} ${rb:-NA} $ge $gr" >> scan_npt.txt
done
if [ "$(awk '$2!="NA" && $3!="NA"' scan_npt.txt | wc -l)" -eq 11 ]; then echo OK > TOLSCAN_NPT.DONE; else echo FAIL > TOLSCAN_NPT.DONE; fi
cat scan_npt.txt
