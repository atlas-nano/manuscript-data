#!/bin/bash
# eta sweep on the 240-atom Au(111) slab (regression case metal_slab_gself), with and
# without the finite-width self-energy. Single point, one rank, current local binary.
# Usage: bash run_sweep.sh        Harvest: python3 analyze.py
set -e
LMP=${LMP:-$HOME/codes/lammps/lammps-dev2/build/lmp_parallel}
HERE=$(cd "$(dirname "$0")" && pwd)
for eta in 5.1720 2.5000 1.0000 0.3000; do
  for g in on off; do
    d=$HERE/runs/eta${eta}_gself${g}
    mkdir -p $d
    cp $HERE/data.au_ortho $d/
    sed "s/5\.1720/${eta}/" $HERE/au_pqeq.param > $d/au_pqeq.param
    grep -q " ${eta} " $d/au_pqeq.param || { echo "eta substitution failed for $eta"; exit 1; }
    if [ $g = on ]; then
      cp $HERE/in.lammps $d/in.lammps
    else
      sed 's/^fix_modify      chg gself on width 0.50/fix_modify      chg gself off/' $HERE/in.lammps > $d/in.lammps
      grep -q 'gself off' $d/in.lammps || { echo "gself off substitution failed"; exit 1; }
    fi
    (cd $d && rm -f charges.dump pe.txt log.lammps && \
     mpirun -np 1 $LMP -in in.lammps -log log.lammps -screen none > run.out 2>&1)
    [ -s $d/charges.dump ] && grep -q 'ITEM: ATOMS' $d/charges.dump || { echo "no charges for $d"; exit 1; }
    echo "done $d"
  done
done
