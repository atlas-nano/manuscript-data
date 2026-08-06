Quantum vs classical weighting of the solid-like DoS component.
TIP4P/2005, 220 K isotherm, 100 ns trajectories (the R1 comment 7 campaign).
Source of response Table R2.

  Sq-thermo.csv          liquid entropy, quantum weighting   (J/mol/K, divide by N*R with N=216)
  Sc-thermo.csv          liquid entropy, classical weighting (same units)
  sex_tot-quantum.csv    excess entropy, quantum weighting   (k_B per molecule)
  sex_tot_classical.csv  excess entropy, classical weighting (k_B per molecule)

WARNING - the column order is not the same in all four files. The two quantum files store
(Temperature, Density) beneath a header that reads "Density, Temperature"; the two classical files
store (Density, Temperature) as the header says. Merging them on the header names produces nonsense.
Temperature is 220 throughout and density runs 0.85-1.30, so the two are easy to tell apart.

Classical weighting is bounded above by quantum at every frequency, so S_ex(classical) must be more
negative than S_ex(quantum) at every state point. Here the separation is 1.698-1.773 k_B. A result
in which the classical value is LESS negative indicates the distinguishable-molecule column (Sd,
onsager's wspd) was read instead of the standard indistinguishable one (wsp); the two differ by
ln N - 1. That error produced the original, retracted version of Table R2.
