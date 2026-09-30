#!/bin/bash
# Diagnostic only: smallest eigenvalue of the projected charge operator on each sweep arm, from the
# code's own m-step Lanczos estimate. `ridge eig <floor> <m>` with a floor far above any eigenvalue
# forces a ridge so that the estimate is printed; the charges of these runs are NOT used.
set -e
LMP=${LMP:-$HOME/codes/lammps/lammps-dev2/build/lmp_parallel}
HERE=$(cd "$(dirname "$0")" && pwd)
for d in $HERE/runs/eta*_gself*; do
  e=$d/lmin; mkdir -p $e
  cp $d/data.au_ortho $d/au_pqeq.param $e/
  sed 's/^run             0/fix_modify      chg ridge eig 1000.0 200\nrun             0/' $d/in.lammps > $e/in.lammps
  grep -q 'ridge eig' $e/in.lammps || { echo "ridge insert failed"; exit 1; }
  (cd $e && mpirun -np 1 $LMP -in in.lammps -log log.lammps -screen none > run.out 2>&1)
  l=$(grep -o 'Lanczos lambda_min=[^ ]*' $e/log.lammps | head -1)
  [ -n "$l" ] || { echo "no lambda_min in $e"; exit 1; }
  echo "$(basename $d) $l"
done
