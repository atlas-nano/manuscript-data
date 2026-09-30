# The temporary operator-dump diagnostic (referee item 2)

`dense_exp1.py` reads `runs/exp1/*/opdump.txt`, written by a diagnostic that is NOT in the release: inserted into
`lammps-release/src/SAMQEQ/fix_qeq_sam_lr.cpp` at the top of `FixQEqSam::qeq_solve()` (after `build_mol_blocks()`),
built, run once per point with `SAMQEQ_OPDUMP=opdump.txt`, then reverted and rebuilt (2026-09-27). It writes the
columns of P·H by applying the code's own `qeq_matvec` to each unit vector (single rank), then stops the run.

```cpp
  if (const char *od = getenv("SAMQEQ_OPDUMP")) {
    if (comm->nprocs != 1) error->all(FLERR, "SAMQEQ_OPDUMP is single-rank only");
    std::vector<int> idx;
    for (int ii = 0; ii < nn; ii++) { int i = ilist[ii]; if (mask[i] & groupbit) idx.push_back(i); }
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return atom->tag[a] < atom->tag[b]; });
    std::vector<double> x(atom->nmax, 0.0), y(atom->nmax, 0.0);
    FILE *fp = fopen(od, "w");
    fprintf(fp, "%d\n", (int) idx.size());
    for (int i : idx) fprintf(fp, "%d ", (int) atom->tag[i]);
    fprintf(fp, "\n");
    for (int j : idx) {
      std::fill(x.begin(), x.end(), 0.0); x[j] = 1.0;
      qeq_matvec(x.data(), y.data());
      for (int i : idx) fprintf(fp, "%.15e ", y[i]);
      fprintf(fp, "\n");
    }
    fclose(fp);
    error->all(FLERR, "SAMQEQ_OPDUMP written to {} ({} columns)", od, idx.size());
  }
```
Validation: the dumped operator is symmetric to 1e-15 and reproduces the code's solved charges to ≤ 2.1e-7 e at
every point (`dense_exp1.csv`, column dq).
