Quantum vs classical weighting of the solid-like DoS component.
TIP4P/2005, 220 K isotherm, 100 ns trajectories (the R1 comment 7 campaign).
Source of response Table R2.

  Sq-thermo.csv          liquid entropy, quantum weighting   (J/mol/K, divide by N*R with N=216)
  Sc-thermo.csv          liquid entropy, classical weighting (same units)
  sex_tot-quantum.csv    excess entropy, quantum weighting   (k_B per molecule)
  sex_tot_classical.csv  excess entropy, classical weighting (k_B per molecule)

All four files carry columns in the order their header states, (Density, Temperature, thermoValue),
and can be merged on the column names.

Note for anyone comparing against an earlier copy of these files: as originally produced, the two
quantum files stored (Temperature, Density) beneath a header reading "Density, Temperature", while
the two classical files stored the order the header states. They were normalised for this deposit.
Only the column order changed; every value is bit-identical, and the two orders are trivially
distinguishable because temperature is 220 throughout while density runs 0.85-1.30.

Classical weighting is bounded above by quantum at every frequency, so S_ex(classical) must be more
negative than S_ex(quantum) at every state point. Here the separation is 1.698-1.773 k_B. A result
in which the classical value is LESS negative indicates the distinguishable-molecule column (Sd,
onsager's wspd) was read instead of the standard indistinguishable one (wsp); the two differ by
ln N - 1. That error produced the original, retracted version of Table R2.
