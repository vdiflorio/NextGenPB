# Reproducing our results (nonlinear PBE project)

Quick guide to rerun everything we tested for the course project, on top
of the Docker setup described in the main `README.md`. All commands
below assume you're inside the dev container with the repo mounted at
`/App/NextGenPB` and already built (`cd src && make distclean && make -j2`).

Commit these results were produced from: `8ed1cc9`
(branch `nonlinear-solver`).

## 1. Linear vs. nonlinear on a real molecule (1CCM)

We use `data/1CCM.pqr` (a real 642-atom protein structure, not a sphere)
with Coulombic boundary conditions.

```bash
cd /App/NextGenPB
mkdir -p test_runs/repro_linear test_runs/repro_nonlinear
cp data/options.prm test_runs/repro_linear/options.prm
cp data/options.prm test_runs/repro_nonlinear/options.prm

for f in test_runs/repro_linear/options.prm test_runs/repro_nonlinear/options.prm; do
  sed -i -E \
    -e 's#^[[:space:]]*filename[[:space:]]*=.*#filename = ../../data/1CCM.pqr#' \
    -e 's#^[[:space:]]*bc_type[[:space:]]*=.*#bc_type = 2#' \
    -e 's#^[[:space:]]*linear_solver[[:space:]]*=.*#linear_solver = lis#' \
    -e 's#^[[:space:]]*calc_energy[[:space:]]*=.*#calc_energy = 0#' \
    -e 's#^[[:space:]]*potential_map[[:space:]]*=.*#potential_map = 1#' \
    -e 's#^[[:space:]]*map_type[[:space:]]*=.*#map_type = vtu#' \
    -e 's#^[[:space:]]*scale[[:space:]]*=.*#scale = 1.0#' \
    "$f"
done
sed -i -E 's#^[[:space:]]*linearized[[:space:]]*=.*#linearized = 1#' test_runs/repro_linear/options.prm
sed -i -E 's#^[[:space:]]*linearized[[:space:]]*=.*#linearized = 0#' test_runs/repro_nonlinear/options.prm

cd test_runs/repro_linear    && /App/NextGenPB/src/ngpb --prmfile options.prm 2>&1 | tee run.log
cd ../repro_nonlinear        && /App/NextGenPB/src/ngpb --prmfile options.prm 2>&1 | tee run.log
```

**What you should see:**
- Linear run: solves once, `linear solver status : normal end`.
- Nonlinear run: Newton loop converges in **6 iterations**, with
  `‖du‖_inf,solvent` starting around **6.7618** (iteration 0, which
  reproduces the linear solve exactly by construction) and the final
  solvent potential ending around **5.2340**.
- No NaN/Inf/segfaults in either log.

Each run folder gets a `potential_map_0000.vtu` you can open in ParaView
to compare the linear and nonlinear potential fields (color by `phi`,
add a Slice filter, use the same Origin/Normal/color range for both so
they're visually comparable).

## 2. Energy is intentionally disabled for the nonlinear model

If you set `calc_energy = 1` (or higher) together with `linearized = 0`,
the code will **not** silently compute a wrong number — it prints a
warning and skips the energy calculation instead:

```
[WARNING] Energy and potential/field post-processing use a linear-response
          functional and are NOT implemented for the nonlinear model.
          Skipping (set linearized=1 to compute them for the linear PBE).
```

This is by design (see `src/poisson_boltzmann.cpp` around line 271):
the energy functional we have is a linear-response quantity and isn't
meaningful applied to a nonlinear `phi`. Our test configs above already
set `calc_energy = 0` for both runs.

## 3. Globalization (clamping) actually engages

With the default `maxdu = 2.0` (in `newton_solve()`,
`src/pb_class.cpp`), the clamp never triggers on 1CCM — full Newton is
already stable there. To see it actually activate:

1. Open `src/pb_class.cpp`, find `const double maxdu = 2.0;` inside
   `newton_solve()`, and temporarily change it to `0.5`.
2. Rebuild: `cd src && make -j2`.
3. Rerun the nonlinear case from step 1.
4. You should see `[Newton] CLAMP active: scale = ...` lines in the
   output starting at iteration 1, one extra Newton iteration (7
   instead of 6), and the **same final converged potential** (~5.2340)
   as the unclamped run.
5. Set `maxdu` back to `2.0` and rebuild before doing anything else —
   don't leave the test value committed.

## 4. Weak / strong scalability

These need MPI ranks, so run with `mpirun`. Note: OpenMPI's default
slot count equals your number of **physical** cores — if you want to
oversubscribe (e.g. run 8 ranks on a 4-core machine), add
`--use-hwthread-cpus` to the `mpirun` command.

**Strong scaling** — fixed problem, vary rank count:

```bash
cd /App/NextGenPB
cp data/options.prm strong_test.prm
sed -i -E \
  -e 's#^[[:space:]]*filename[[:space:]]*=.*#filename = /App/NextGenPB/data/1CCM.pqr#' \
  -e 's#^[[:space:]]*bc_type[[:space:]]*=.*#bc_type = 2#' \
  -e 's#^[[:space:]]*scale[[:space:]]*=.*#scale = 1.0#' \
  strong_test.prm
for np in 1 2 4 8; do
  mpirun --allow-run-as-root -np $np src/ngpb --prmfile strong_test.prm 2>&1 | tee strong_np${np}.log
done
```

Look at the `Compute numerical solution` line in the `Timing Report:`
at the end of each log — that's the actual parallel solve time to
compare across rank counts.

**Weak scaling** — problem size grows with rank count, so DOF-per-rank
stays roughly constant. We used these `scale` values to get
approximately 90k DOF per rank at 1/2/4/8 ranks (93k / 176k / 326k /
635k total DOF — not exact, the mesh generator doesn't give an exact
linear knob, but close enough):

| ranks | scale |
|---|---|
| 1 | 1.0 |
| 2 | 1.26 |
| 4 | 1.587 |
| 8 | 2.0 |

Same as above but set `scale` per rank count before each run, and check
the DOF count actually produced by grepping `global nodes` in the log
if you want to double check the calibration still holds on your
machine.

**Heads up:** run each configuration a few times (we did 3 reps) and
take the median — timing on a shared/virtualized machine (e.g. a
laptop running Docker/WSL2) is noisy, especially once you oversubscribe
physical cores. Don't trust a single run.

**What we actually found:** at these problem sizes, adding MPI ranks
does **not** give a real speedup (strong scaling is roughly flat), and
weak-scaling efficiency drops off fast as both rank count and problem
size grow together. This is a real result, not a bug — we checked it's
not due to a crash or wrong answer, the solver just doesn't benefit
from parallelism here. See the report for the numbers and discussion.

## Notes

- All `test_runs/`, `*.log`, and generated `.prm`/`.vtu` files are
  gitignored on purpose — don't force-add them, regenerate instead.
- If you change `maxdu` or anything else in `src/` to test something,
  always check `git status` afterwards and revert before committing.
