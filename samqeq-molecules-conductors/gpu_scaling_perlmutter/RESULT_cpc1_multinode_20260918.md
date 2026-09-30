# CPC-1: multi-node GPU strong scaling on Perlmutter, 2026-09-18

Placeholder CPC-1 of the split plan. It was in the drop list of `AUDIT_split_plan_20260918.md` because
Expanse cannot do it (no `nvidia_peermem`, so no inter-node GPUDirect RDMA, row 117). Perlmutter came
back early on 2026-09-18 with rebuilt binaries, so it was measured instead.

## Binary verification, before any timing

| check | result |
|---|---|
| `~/codes/bin/lmp_perlmutter_gpu` → `~/codes/lammps/bin/lmp_gpu` | rebuilt **2026-09-18 06:14** |
| `~/codes/bin/lmp_perlmutter_cpu` → `~/codes/lammps/bin/lmp_cpu` | rebuilt **2026-09-18 06:09** |
| source `/pscratch/.../lammps-dev2/src/SAMQEQ` | synced 06:05 |
| `ldd lmp_gpu \| grep -c "not found"` | **0** |
| CUDA runtime | **`libcudart.so.13` only** — no regression to 12.9 (that has happened twice) |
| marker `a50d0e1` (`KSpace::FORWARD_AD`) | present, 3 occurrences |
| stale file-level `enum ... FORWARD_IK` | **0** — the Kokkos force fix is in |
| marker `fd1a75f` (`redox_dp`, chicn s-basis) | present, 13 occurrences |

## Protocol

`in.bench` from `gpu_bench_expanse/`, unchanged. Native SPC-FQ water, `base.data` (512 molecules)
regenerated on Perlmutter with `in.equil` at its shipped `NEQ 5000`, replicated REP³. 200 timesteps,
full solve every step (`DEV 1`, `ASPC 0`), Nosé--Hoover, PPPM 1e-6, solve tolerance 1e-5.

Launch, one rank per A100, GPU-aware MPI (required; `--gpus-per-task=1` hangs in GTL IPC):

    export MPICH_GPU_SUPPORT_ENABLED=1
    srun -N $NN -n $((NN*4)) -c 32 --gpus-per-node=4 --gpu-bind=none \
         lmp_gpu -k on g 4 -sf kk -pk kokkos gpu/aware on newton on neigh half -in in.bench ...

Jobs 58538410 (N=1), 58538411 (N=2), 58538412 (N=4), account m5400. Markers earned: a run counts only
if its log carries a `Performance:` line. All 3 wrote `SCALE_n*.DONE` at 2/2 arms.

## Result

| atoms | 4 A100 (1 node) | 8 (2 nodes) | 16 (4 nodes) | speedup 4→16 | efficiency |
|---:|---:|---:|---:|---:|---:|
| 331,776 | **3.412** | 3.381 | **4.472** | 1.31× | 33% |
| 786,432 | **1.472** | 1.385 | **1.815** | 1.23× | 31% |

timesteps s⁻¹. **Two nodes are marginally slower than one at both sizes.** Four nodes return 1.2–1.3×.

Against the existing Expanse table, one A100 node is about twice one V100 node on the same deck
(3.412 vs 1.643 at 331,776; 1.472 vs 0.568 at 786,432). The single-node claim in the abstract is
unaffected.

## Mechanism, measured

| device count | Modify (charge solve) | Kspace (force-time mesh) | LAMMPS Comm | CG iters/step |
|---:|---:|---:|---:|---:|
| 4 | 91.4% | 6.1% | 0.58% | 15 |
| 8 | 91.0% | 7.7% | 0.30% | 15 |
| 16 | 92.1% | 6.9% | 0.19% | 15 |

(786,432 atoms; the 331,776 rows behave the same, 90.6–91.6% Modify.)

Three things follow, and together they identify the bottleneck:

1. **The solve is 91–92% of the timestep everywhere.** Whatever limits scaling is inside it.
2. **The iteration count is constant at 15** across all six runs, so the solve is not doing more work
   on more nodes. The loss is per-iteration cost, not extra iterations.
3. **LAMMPS's own `Comm` falls** (0.58 → 0.19%) as atoms per device drop, which is the normal MD
   behaviour. The standard task breakdown therefore cannot see the problem: the solve's own
   collectives are charged to the fix, not to `Comm`.

Each CG iteration carries a global reduction for its inner products and a distributed transform for the
reciprocal term, so a 15-iteration step crosses the interconnect roughly thirty times. That is
latency-bound once it leaves the node.

## What went into the paper

- §9 Performance: three paragraphs and `Figure~\ref{fig:gpuscale}` (two panels, from `scaling.csv`).
- §8.2 Current limitations: the device path does not strong-scale across nodes; single-node device runs
  are the supported case.
- The abstract is unchanged: its claim is single-node and still holds.

## Honest scope

- Strong scaling only, two sizes, one repeat per point. No weak-scaling series and no error bars; the
  200-step timings are stable to the third digit but a single run each.
- Not tried: the predictor--corrector under multi-node, which would cut the collective count per step
  and is the obvious next measurement if a reviewer asks for one.
- The recommendation in §9 (communication-avoiding Krylov, or accepting the ASPC criterion) is stated
  as a direction, not a measured claim.
