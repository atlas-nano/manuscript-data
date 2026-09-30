#!/bin/bash
# gself at its DEFAULT width (the pair kernel's own Gaussian, R_c*sqrt(lambda)-derived, 2.38 A for Au),
# eta = 5.172 and 0.3: layer charges plus the Lanczos lambda_min (diagnostic ridge run, charges unused).
set -e
LMP=${LMP:-$HOME/codes/lammps/lammps-dev2/build/lmp_parallel}
HERE=$(cd "$(dirname "$0")" && pwd)
for eta in 5.1720 0.3000; do
  d=$HERE/runs/eta${eta}_gselfdefault; mkdir -p $d
  cp $HERE/data.au_ortho $d/
  sed "s/5\.1720/${eta}/" $HERE/au_pqeq.param > $d/au_pqeq.param
  sed 's/^fix_modify      chg gself on width 0.50/fix_modify      chg gself on/' $HERE/in.lammps > $d/in.lammps
  grep -q '^fix_modify      chg gself on$' $d/in.lammps || { echo "default-width substitution failed"; exit 1; }
  (cd $d && rm -f charges.dump pe.txt && mpirun -np 1 $LMP -in in.lammps -log log.lammps -screen none > run.out 2>&1)
  e=$d/lmin; mkdir -p $e; cp $d/data.au_ortho $d/au_pqeq.param $e/
  sed 's/^run             0/fix_modify      chg ridge eig 1000.0 200\nrun             0/' $d/in.lammps > $e/in.lammps
  (cd $e && mpirun -np 1 $LMP -in in.lammps -log log.lammps -screen none > run.out 2>&1)
  echo "$(basename $d) $(grep -o 'Lanczos lambda_min=[^ ]*' $e/log.lammps | head -1) $(grep -o 'gself per type[^\n]*' $d/log.lammps | head -1)"
done
