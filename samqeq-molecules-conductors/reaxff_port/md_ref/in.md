# Reference: stock ReaxFF water (LAMMPS examples/reaxff/water, Achtyl et al. 2015 QEq water),
# single point, charges from fix qeq/reaxff at a tight tolerance.
boundary        p p p
units           real
atom_style      charge
read_data       ../data.water
pair_style      reaxff NULL safezone 3.0 mincap 150
pair_coeff      * * ../qeq_ff.water O H
neighbor        0.5 bin
neigh_modify    every 1 delay 0 check yes
fix             1 all qeq/reaxff 1 0.0 10.0 1.0e-10 reaxff maxiter 2000
compute         reax all pair reaxff
variable        eqeq equal c_reax[14]
thermo_style    custom step temp pe etotal press
thermo_modify   format float %.10f
velocity        all create 300.0 4928459 rot yes dist gaussian
fix             nvt all nvt temp 300 300 50.0
timestep        0.5
thermo          10
run             20
variable        ecv equal ecoul
