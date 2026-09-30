# samQEq on the same ReaxFF water: cube-root kernel, global neutrality (atom_style charge, no molecule IDs),
# taper-cutoff operator to 10 A with no reciprocal (fix_modify cutoff on), tolerance 1e-10.
# pair reaxff keeps supplying ReaxFF's own shielded Coulomb from the solved charges (checkqeq no).
# coul/long + pppm/samqeq are present only because the lr_ewald=2 path reads g_ewald from kspace; the
# comparison uses the reaxff sub-style energy, not pe.
boundary        p p p
units           real
atom_style      charge
read_data       ../data.water
pair_style      hybrid/overlay reaxff NULL checkqeq no safezone 3.0 mincap 150 coul/long 10.0
pair_coeff      * * reaxff ../qeq_ff.water O H
pair_coeff      * * coul/long
kspace_style    pppm/samqeq 1.0e-5
neighbor        0.5 bin
neigh_modify    every 1 delay 0 check yes
fix             chg all qeq/sam 1 0.0 10.0 1.0e-10 ../sam_opmatch.param
fix_modify      chg cutoff on
fix_modify      chg energy no          # reaxff already carries the polarization energy
pair_modify     pair coul/long compute no
kspace_modify   compute no
compute         reax all pair reaxff
variable        ereax equal c_reax
variable        eqeq equal c_reax[14]
thermo_style    custom step temp pe etotal press
thermo_modify   format float %.10f
velocity        all create 300.0 4928459 rot yes dist gaussian
fix             nvt all nvt temp 300 300 50.0
timestep        0.5
thermo          10
run             20
