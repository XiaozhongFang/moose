# MAS1998 Quantitative Reference and Reproduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a converged MAS1998 reference along exact solid-body characteristics, compare all six Table 3 runs against it with reproducible error norms, and provide restart-safe manual commands for computations that exceed an interactive session.

**Architecture:** Keep the six Eulerian Type I/II production solvers unchanged except for a shared, coordinate-correct local-solar-time helper. Add a separate characteristic reference solver which samples the same diagonal cells as each production grid, removes horizontal discretization by following the analytic solid-body trajectory, and integrates only the existing 15-layer chemistry-diffusion system. A Python comparison tool validates coordinate alignment and writes normalized L1/L2/Linf errors for the four species plotted in the paper.

**Tech Stack:** MOOSE C++17 application, existing KPP `Fun_SPLIT` mechanism, existing Euler/BDF2 Gauss-Seidel column integrator, GoogleTest, MOOSE TestHarness, Python 3 standard library, CSV diagnostics.

**Execution status (2026-08-17):** The implementation, short verification, long-run inputs, comparison
tools, and manual runner are complete. Internet Archive replay did not yield the author payload. The
14-day jobs and data-based convergence acceptance remain intentionally unchecked and must be run manually.

---

## Evidence Policy

The author-generated reference and an independently regenerated reference are related but not identical evidence:

1. Prefer the archived CWI data if the author page links and data payload can be recovered. Preserve its original bytes and SHA-256; never rewrite it as generated output.
2. If the payload remains unavailable, generate an **independent characteristic reference** using the mathematical construction in paper Section 4.5: exact horizontal characteristics plus the published semi-discrete 15-layer vertical operator.
3. Call the independent result converged only after step-halving and Gauss-Seidel convergence meet the numerical criteria in Task 5.
4. Do not call any new output the original CWI reference or claim bitwise agreement without the archived payload.
5. Quantitative reproduction means repeatable error tables against the best available reference, not visual agreement alone.

Internet Archive CDX currently lists author-page captures at:

```text
19970410111732
19980125163638
20000902233637
20010831050617
```

for `http://www.cwi.nl/ftp/edwins/Ref_Sol_Benchmark_Global.html`. At planning time the replay service returns a temporary-offline response, so data recovery is an opportunistic first task, not a blocker for the independent generator.

## File Structure

Standalone app additions and modifications:

- Modify `include/utils/MAS1998BenchmarkUtils.h`, `src/utils/MAS1998BenchmarkUtils.C`: analytic trajectory and coordinate-correct local solar time.
- Modify `src/utils/MAS1998BenchmarkSolver.C`, `include/utils/MAS1998BenchmarkSolver.h`: reuse the shared local-solar-time helper.
- Create `include/utils/MAS1998CharacteristicReferenceSolver.h`, `src/utils/MAS1998CharacteristicReferenceSolver.C`: diagonal sampling and exact-characteristic column integration.
- Create `include/executioners/MAS1998ReferenceExecutioner.h`, `src/executioners/MAS1998ReferenceExecutioner.C`: MOOSE input surface for reference jobs.
- Create `test/tests/mas1998/reference/common_characteristic.i` and three grid-specific inputs.
- Create `test/tests/mas1998/mas1998_characteristic_reference_smoke.i`: short full-stack regression.
- Modify `test/tests/mas1998/tests`: register the reference smoke and CSV check.
- Modify `unit/src/MAS1998BenchmarkUtilsTest.C`: trajectory and solar-time tests.
- Create `scripts/compare_diagonal.py`, `scripts/tests/test_compare_diagonal.py`: quantitative error calculation and synthetic regression.
- Create `scripts/run_quantitative_reproduction.sh`: bounded-parallel manual production/reference runner.
- Modify `README.md`, `test/tests/mas1998/README.md`: exact commands, timing, evidence labels, and outputs.

MOOSE documentation changes:

- Modify `.codex/docs/numerical_mapping.md`: replace the current regular-grid reference description with the characteristic-reference hierarchy and quantitative acceptance criteria.
- Modify `.codex/docs/architecture.md`: add the reference solver without moving benchmark-specific code into Atmospheric Chemistry.

### Task 1: Recover and Classify Historical CWI Artifacts

**Files:**
- Create if recovered: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/cwi/PROVENANCE.md`
- Create if recovered: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/cwi/SHA256SUMS`
- Do not commit a downloaded data payload until its redistribution status and size are reviewed.

- [x] **Step 1: Query the archived author page**

Run:

```bash
curl -L --fail --silent --show-error \
  'https://web.archive.org/web/19980125163638id_/http://www.cwi.nl:80/ftp/edwins/Ref_Sol_Benchmark_Global.html'
```

Expected when replay is available: HTML from the CWI page, not an Internet Archive offline page.

- [ ] **Step 2: Enumerate page links without guessing filenames**

Use Python's `html.parser.HTMLParser` to print every `href`, resolve it against the captured CWI URL, and query each candidate through the same `19980125163638id_` replay prefix. Select files described by the page as reference solution or benchmark data; do not select figures merely because their names contain `Ref`.

- [ ] **Step 3: Record provenance before parsing**

For each recovered artifact record original URL, archive timestamp, retrieval UTC time, byte length, media type, and SHA-256 in `PROVENANCE.md` and `SHA256SUMS`.

- [x] **Step 4: Classify the outcome**

If numeric author data is recovered, label it `CWI_AUTHOR_REFERENCE`. If only the page or source is recovered, label it `CWI_DOCUMENTATION_ONLY`. If replay remains unavailable, record the four CDX timestamps in the documentation and proceed with Task 2.

### Task 2: Correct Shared Longitude/Time Semantics and Add Exact Trajectories

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998BenchmarkUtils.h`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998BenchmarkUtils.C`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998BenchmarkSolver.h`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998BenchmarkSolver.C`
- Test: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998BenchmarkUtilsTest.C`

- [x] **Step 1: Write failing coordinate tests**

Add tests which establish the paper's coordinate definition `lambda' = lambda - pi`:

```cpp
TEST(MAS1998BenchmarkUtils, convertsPrimeLongitudeToLocalSolarTime)
{
  EXPECT_DOUBLE_EQ(MAS1998::localSolarTime(0.0, -180.0), 0.0);
  EXPECT_DOUBLE_EQ(MAS1998::localSolarTime(0.0, 0.0), 43200.0);
  EXPECT_DOUBLE_EQ(MAS1998::localSolarTime(0.0, 90.0), 64800.0);
}
```

The second argument is `lambda'` in degrees. The first assertion represents Greenwich at `lambda=0`, hence `lambda'=-180 deg`.

- [x] **Step 2: Write failing trajectory tests**

Define a public POD result:

```cpp
struct SphericalLocation
{
  Real longitude_prime_degrees;
  Real latitude_degrees;
};
```

Test that `solidBodyTrajectory(initial, 14 days)` returns to the initial point, forward then backward rotations cancel, and a centered finite difference of the trajectory agrees with `solidBodyWind()` in both eastward and northward components at a non-polar point.

- [ ] **Step 3: Run unit tests and confirm failure**

Run from the app:

```bash
conda run -n moose --no-capture-output make -C unit -j4
conda run -n moose --no-capture-output ./unit/mas1998_benchmark-unit-opt \
  --gtest_filter='MAS1998BenchmarkUtils.*'
```

Expected: compile failure because `localSolarTime`, `SphericalLocation`, and `solidBodyTrajectory` do not exist.

- [x] **Step 4: Implement coordinate-correct local solar time**

Implement:

```text
lambda = lambda_prime + 180 deg
local_seconds = (gmt_seconds + lambda/360 deg * 86400 s) mod 86400 s
```

Normalize the result into `[0,86400)`. Replace `MAS1998BenchmarkSolver::localSolarTime()` with this shared helper and remove the private method. This intentionally fixes the current 12-hour offset caused by treating `lambda'` as geographic `lambda`.

- [x] **Step 5: Implement the analytic solid-body rotation**

Convert `(lambda'+pi, phi)` to a Cartesian unit vector, rotate it by Rodrigues' formula through
`theta=2*pi*elapsed/(14 days)` about

```text
n = (-sin(beta), 0, cos(beta)), beta=45 deg
```

and convert back to `(lambda',phi)`. Normalize `lambda'` to `[-180,180)` and clamp the Cartesian `z` component to `[-1,1]` before `asin`.

- [x] **Step 6: Pass focused unit tests**

Run the command from Step 3. Expected: all `MAS1998BenchmarkUtils` tests pass, including velocity consistency within a tolerance justified by the finite-difference interval.

- [x] **Step 7: Commit the coordinate correction**

```bash
git add include/utils/MAS1998BenchmarkUtils.h src/utils/MAS1998BenchmarkUtils.C \
  include/utils/MAS1998BenchmarkSolver.h src/utils/MAS1998BenchmarkSolver.C \
  unit/src/MAS1998BenchmarkUtilsTest.C
git commit -m 'Correct MAS1998 solar time and add exact trajectories'
```

### Task 3: Implement the Characteristic Reference Solver

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998CharacteristicReferenceSolver.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998CharacteristicReferenceSolver.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/executioners/MAS1998ReferenceExecutioner.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/executioners/MAS1998ReferenceExecutioner.C`

- [x] **Step 1: Define reference options**

Use this minimal options surface:

```cpp
struct Options
{
  unsigned int n_longitude;
  unsigned int n_latitude;
  bool reduced_grid;
  Real start_time;
  Real end_time;
  Real chemistry_step;
  unsigned int gauss_seidel_iterations;
  Real plume_radius_degrees;
  std::string mechanism_library;
  std::string output_file;
};
```

Reject non-positive time intervals/steps, zero sweeps, empty mechanism path, and intervals not exactly divisible except for the integrator's supported final short step.

- [x] **Step 2: Select exactly the production diagnostic locations**

Construct `MAS1998ReducedGrid(options.n_longitude, options.n_latitude, options.reduced_grid)`. For every latitude row, use the same rule as production diagnostics:

```cpp
target_lambda_prime = 0.5 * ring.latitude_center;
cell = grid.cellAtLongitude(row, target_lambda_prime);
sample_lambda_prime = grid.longitudeCenter(row, longitude_index_of_cell);
sample_phi_prime = ring.latitude_center;
```

Store the actual cell-center longitude so reference and production rows align exactly without interpolation.

- [x] **Step 3: Trace the final sample backward to its initial foot**

For duration `T=end_time-start_time`, compute the initial foot with
`solidBodyTrajectory(final_location, -T)`. Initialize all 15x17 values from
`MAS1998::initialConcentration()` at that foot. At rate-evaluation time `t`, evaluate chemistry at
`solidBodyTrajectory(initial_foot, t-start_time)`.

- [x] **Step 4: Integrate the coupled column without horizontal splitting**

Call the existing `MAS1998ChemistryIntegrator::integrate()` once per sampled row with
`chemistry=true`, `diffusion=true`, the reference step and sweep count, the same lowest-layer NO source,
and the same KPP fixed species/temperature/pressure setup as production. This is the Section 4.5
characteristic construction: horizontal advection is exact and only the published vertical spatial
discretization remains.

- [x] **Step 5: Write production-compatible output**

Write final-time ground rows using the existing column names:

```text
n_longitude,n_latitude,reduced_grid,splitting,time,row,cell,
lambda_prime,phi_prime,z_km,<17 species>,NOx
```

Set `splitting=characteristic_reference`. Use `setprecision(17)`. Also write metadata comments or a sidecar summary containing `chemistry_step`, sweep count, mechanism path, and evidence label `INDEPENDENT_CHARACTERISTIC_REFERENCE`.

- [x] **Step 6: Register `MAS1998ReferenceExecutioner`**

Expose the exact options from Step 1. Require one MPI rank for each reference process, matching KPP's process-local global state. The executable remains the same `mas1998_benchmark-opt`; the input selects the new executioner type.

- [x] **Step 7: Build the app**

```bash
conda run -n moose --no-capture-output make -j4
```

Expected: `mas1998_benchmark-opt` links with both executioners.

- [x] **Step 8: Commit the reference solver**

```bash
git add include/utils/MAS1998CharacteristicReferenceSolver.h \
  src/utils/MAS1998CharacteristicReferenceSolver.C \
  include/executioners/MAS1998ReferenceExecutioner.h \
  src/executioners/MAS1998ReferenceExecutioner.C
git commit -m 'Add MAS1998 characteristic reference solver'
```

### Task 4: Add Reference Inputs and Full-Stack Regression

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/common_characteristic.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/characteristic_64x32.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/characteristic_128x64.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/characteristic_256x128.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/mas1998_characteristic_reference_smoke.i`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/tests`

- [x] **Step 1: Add a 600 s regular-grid smoke input**

Use `8x4`, `reduced_grid=false`, `chemistry_step=300`, two sweeps, and output
`mas1998_characteristic_reference_smoke.csv`. The short interval exercises backward foot tracing rather
than relying on 14-day closure.

- [x] **Step 2: Register the smoke after KPP generation**

Add a heavy `RunApp` test with `prereq=build_mas1998_methane_kpp`, followed by the existing
`check_summary.py` or a reference-specific CSV checker requiring finite, non-negative values and four
rows.

- [x] **Step 3: Add the converged-reference candidate inputs**

Use the three reduced production sampling grids, 14 days, `chemistry_step=75 s`, and eight sweeps. Name
outputs `characteristic_reference_{64x32,128x64,256x128}.csv`. These are candidates until Task 5 passes;
the filenames must not include `CWI`.

- [x] **Step 4: Add convergence-check variants**

For `64x32`, add inputs at `(150 s, 4 sweeps)`, `(75 s, 8 sweeps)`, and `(37.5 s, 12 sweeps)`. They use
the same sample coordinates and mechanism. The finest run is the convergence oracle, not automatically
the production reference.

- [x] **Step 5: Run the short reference test**

```bash
conda run -n moose --no-capture-output ./run_tests \
  --re=mas1998_characteristic_reference -j1
```

Expected: KPP build, reference smoke, and CSV validation pass.

- [x] **Step 6: Commit inputs and regression**

```bash
git add test/tests/mas1998/reference test/tests/mas1998/mas1998_characteristic_reference_smoke.i \
  test/tests/mas1998/tests
git commit -m 'Add MAS1998 characteristic reference inputs'
```

### Task 5: Define and Verify Reference Convergence

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/compare_diagonal.py`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/tests/test_compare_diagonal.py`

- [x] **Step 1: Write synthetic failing tests**

Using `unittest`, create two tiny CSVs with exact metadata/coordinate matches and known values. Assert
the script computes, for each species,

```text
relative_L1   = sum(abs(case-ref)) / sum(abs(ref))
relative_L2   = sqrt(sum((case-ref)^2) / sum(ref^2))
relative_Linf = max(abs(case-ref)) / max(abs(ref))
absolute_Linf = max(abs(case-ref))
```

for `O3`, `NOx`, `HNO3`, and `HO2NO2`.

- [x] **Step 2: Test rejection paths**

Assert failure for mismatched grid metadata, row count, row index, coordinates, final time, NaN/Inf, and a
zero reference denominator. Coordinates must match within `1e-12 deg`; no hidden interpolation is allowed.

- [ ] **Step 3: Run tests and confirm failure**

```bash
conda run -n moose --no-capture-output python -m unittest \
  scripts.tests.test_compare_diagonal -v
```

Expected: import failure because `compare_diagonal.py` does not exist.

- [x] **Step 4: Implement the comparison CLI**

Interface:

```bash
python scripts/compare_diagonal.py \
  --reference characteristic_reference_64x32.csv \
  --output errors_64x32.csv \
  mas1998_type_i_64x32.csv mas1998_type_ii_64x32.csv
```

Read only each file's final-time ground rows, sort by row, validate alignment, and write one row per
`case,species` with the four norms and number of points. Preserve each case's splitting and grid metadata.

- [x] **Step 5: Pass Python tests**

Run Step 3. Expected: all tests pass.

- [x] **Step 6: Establish a convergence acceptance rule**

Compare `(150 s,4)` and `(75 s,8)` against `(37.5 s,12)`. For each of the four diagnostic species require:

```text
relative_L2(75s,8 vs 37.5s,12) <= 1.0e-3
relative_Linf(75s,8 vs 37.5s,12) <= 2.0e-3
```

Also require the 75 s errors to be smaller than the 150 s errors. If any condition fails, promote
`(37.5 s,12)` to the candidate and compare it against a new `(18.75 s,16)` run using the same thresholds.

- [x] **Step 7: Commit the comparison tool**

```bash
git add scripts/compare_diagonal.py scripts/tests/test_compare_diagonal.py
git commit -m 'Add MAS1998 quantitative diagonal comparison'
```

### Task 6: Add a Bounded-Parallel Manual Runner

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/run_quantitative_reproduction.sh`

- [x] **Step 1: Implement strict preflight**

The script must use `set -euo pipefail`, accept `--jobs N` and `--output-dir DIR`, require the executable
and KPP library, reject an existing non-empty output directory, and print all selected inputs before
launching. Default `N=4`; reject `N<1`.

- [x] **Step 2: Run each input in a separate process**

Launch the six production inputs and three characteristic references through:

```text
conda run -n moose --no-capture-output ./mas1998_benchmark-opt -i INPUT \
  Executioner/output_file=ABSOLUTE_OUTPUT
```

Use a FIFO/job-count loop or `xargs -P` to cap concurrency. Write one log and one `/usr/bin/time -v`
report per input. Do not use MPI ranks because both executioners require one rank.

- [x] **Step 3: Fail on partial output**

After all jobs finish, require each production CSV to contain final time `1209600`, each reference CSV to
contain final time `1209600`, and each summary to pass its checker. Preserve logs on failure and return
nonzero.

- [x] **Step 4: Generate three error tables and one plot**

Pair Type I/II with the matching characteristic reference for each grid, call
`compare_diagonal.py`, and call `plot_diagonal.py` to overlay the six final curves and references.

- [x] **Step 5: Add a dry-run mode and verify it**

`--dry-run` prints nine commands and creates no output directory. Verify:

```bash
conda run -n moose --no-capture-output \
  scripts/run_quantitative_reproduction.sh --dry-run --jobs 4
```

Expected: nine quoted commands and exit 0.

- [x] **Step 6: Commit the runner**

```bash
git add scripts/run_quantitative_reproduction.sh
git commit -m 'Add MAS1998 quantitative reproduction runner'
```

### Task 7: Document Runtime and Quantitative Evidence

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/README.md`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/README.md`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/numerical_mapping.md`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/architecture.md`

- [x] **Step 1: Record measured timing rather than paper C90 timing**

Record the 2026-08-17 local measurements on one process:

```text
64x32 Type I, one 4800 s split: 28.35 s execution, about 150 MB RSS
64x32 Type II, one 4800 s split: 29.90 s execution, about 150 MB RSS
```

State the resulting rough serial estimates: coarse 2.0-2.1 h per case, medium 7.5-7.9 h, fine
31-33 h, approximately 84 h for all six sequentially. Label these extrapolations, not measurements of
the complete runs.

- [x] **Step 2: Document the one-command manual workflow**

```bash
cd /home/fangxiaozhong/git_repo/mas1998_benchmark
nohup conda run -n moose --no-capture-output \
  scripts/run_quantitative_reproduction.sh \
  --jobs 4 --output-dir runs/mas1998_quantitative_2026-08-17 \
  > runs/mas1998_quantitative_2026-08-17.launch.log 2>&1 &
```

Explain that the custom executioner has no mid-run checkpoint/restart: an interrupted case restarts from
`t=0`. The runner preserves completed outputs and logs but must not silently accept missing final rows.

- [x] **Step 3: Update the evidence hierarchy**

In `numerical_mapping.md`, distinguish `CWI_AUTHOR_REFERENCE` if recovered, converged
`INDEPENDENT_CHARACTERISTIC_REFERENCE`, the older regular-grid refined experiment, and six production
results. Remove language implying the regular-grid input is the best replacement once characteristic
results exist.

- [x] **Step 4: Commit app and MOOSE documentation separately**

```bash
git -C /home/fangxiaozhong/git_repo/mas1998_benchmark add README.md test/tests/mas1998/README.md
git -C /home/fangxiaozhong/git_repo/mas1998_benchmark commit -m 'Document MAS1998 quantitative runs'

git -C /home/fangxiaozhong/git_repo/moose add \
  .codex/docs/numerical_mapping.md .codex/docs/architecture.md \
  .codex/docs/superpowers/plans/2026-08-17-mas1998-quantitative-reference.md
git -C /home/fangxiaozhong/git_repo/moose commit -m 'Document MAS1998 quantitative reference workflow'
```

### Task 8: Focused Verification and Long-Run Handoff

**Files:** all files above; preserve unrelated
`/home/fangxiaozhong/git_repo/moose/HETerogeneous-vectorized-or-Parallel/`.

- [x] **Step 1: Build in the established conda environment**

```bash
conda run -n moose --no-capture-output make -j4
conda run -n moose --no-capture-output make -C unit -j4
```

Run from `/home/fangxiaozhong/git_repo/mas1998_benchmark`. Expected: both targets are up to date or build successfully.

- [x] **Step 2: Run focused numerical and script tests**

```bash
conda run -n moose --no-capture-output ./unit/mas1998_benchmark-unit-opt \
  --gtest_filter='MAS1998BenchmarkUtils.*:MAS1998ChemistryIntegrator.*:MAS1998VerticalOperator.*'
conda run -n moose --no-capture-output python -m unittest \
  scripts.tests.test_compare_diagonal -v
conda run -n moose --no-capture-output ./run_tests \
  --re='mas1998_complete_type|mas1998_characteristic_reference' --heavy -j1
```

Expected: all selected checks pass.

- [x] **Step 3: Syntax-check all long inputs**

For each of the six production and three characteristic inputs run:

```bash
conda run -n moose --no-capture-output ./mas1998_benchmark-opt -i INPUT --check-input
```

Expected: `Syntax OK` nine times.

- [x] **Step 4: Check both repository diffs**

```bash
git -C /home/fangxiaozhong/git_repo/mas1998_benchmark diff --check
git -C /home/fangxiaozhong/git_repo/moose diff --check
```

- [x] **Step 5: Do not run the nine long jobs interactively**

Return the exact `nohup` command from Task 7, the measured timing basis, expected output filenames, and
the commands to inspect progress:

```bash
tail -f runs/mas1998_quantitative_2026-08-17.launch.log
find runs/mas1998_quantitative_2026-08-17 -name '*.csv' -o -name '*.log'
```

- [ ] **Step 6: Complete quantitative acceptance after manual jobs finish**

The reproduction is quantitatively complete only when all six final-time rows exist, the chosen
characteristic reference passes Task 5 convergence, the three error CSVs are generated, and the report
records all norms without replacing failed/negative/NaN results.

## Self-Review

- Spec coverage: historical author data, independent same-method reference generation, coordinate correctness, convergence evidence, all six quantitative comparisons, runtime measurement, and manual long-run commands are covered.
- Scope: all benchmark policy remains in the standalone app; Atmospheric Chemistry receives no MAS1998-specific code.
- Placeholder scan: no deferred implementation placeholder is present. The archive outage and convergence branch have explicit outcomes and commands.
- Type consistency: reference and production share the same 17-species order, vertical operator, KPP production/loss contract, NO source, and output coordinates.
- Numerical validity: the reference removes horizontal discretization rather than merely refining it, and it cannot be accepted without a documented step/sweep convergence check.
