# MAS1998 Complete Reproduction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the 1998 Spee et al. global transport-chemistry benchmark in the standalone `mas1998_benchmark` MOOSE application, including its reduced spherical grid, 17-species chemistry, vertical diffusion, Type I and Type II splitting, paper time integrators, nitrogen correction, 14-day production cases, and reproducible diagnostics.

**Architecture:** A MOOSE `MAS1998Executioner` owns a benchmark-specific array solver because libMesh's ordinary mesh topology cannot exactly represent the paper's latitude-ring reduction and coarse/fine flux summation. Focused classes implement the reduced grid, conservative horizontal transport, the 15-layer density-weighted vertical operator, and the Euler/BDF2 Gauss-Seidel chemistry integrator. The Atmospheric Chemistry module remains benchmark-neutral and only gains generic KPP `Fun_SPLIT` and fixed-species APIs.

**Tech Stack:** MOOSE C++17 application framework, KPP-generated C mechanism loaded through `KPPGeneratedMechanism`, GoogleTest, MOOSE TestHarness, CSV diagnostics.

---

## Reproduction Scope and Evidence Boundaries

The normative benchmark definition is Spee et al. (1998),
`/mnt/d/ZoteroData/storage/D7JMCSVH/1-s2.0-S0378475498001554-main.pdf`. CWI report
NM-R9505, `/mnt/d/BaiduSyncdisk/Zotero/ZoteroData/storage/DQYC9G5U/05028D.pdf`, is used to
check detailed chemistry formulas. It describes the 46-reaction/19-species predecessor, so the final
45-reaction/17-variable mechanism already recovered in the app remains authoritative.

The implementation covers the published model and numerical method, but does not claim to reproduce
Cray C90 hardware Mflop rates, compiler vectorization, or autotasking. Three source limitations must
remain visible in code and documentation:

1. The original `Ref_Sol_Benchmark_Global.text` is unavailable. Any newly generated refined solution
   is labelled an independent replacement reference, never the original CWI data.
2. Figure 1 and Table 2 determine the `64x32` ring schedule directly. For `128x64` and `256x128`, the
   schedules are reconstructed as symmetric, monotone powers-of-two grids with one polar ring at each
   coarser resolution and the half-resolution band width selected to reproduce Table 2 exactly.
3. The paper names a standard variable-step local-error controller but does not publish all controller
   constants. Paper runs therefore use the documented 300 s minimum as a deterministic constant
   chemistry step; a refined replacement-reference input uses smaller steps and is explicitly not a
   bitwise reconstruction of the unavailable CWI solution.

## File Structure

Atmospheric Chemistry module changes:

- `modules/atmospheric_chemistry/include/utils/KPPGeneratedMechanism.h`: generic fixed-species and
  KPP production/loss coefficient API.
- `modules/atmospheric_chemistry/src/utils/KPPGeneratedMechanism.C`: resolve `Fun_SPLIT`, preserve
  KPP pressure units, apply fixed-species overrides, and evaluate production/loss.

Standalone app numerical core:

- `include/utils/MAS1998BenchmarkUtils.h`, `src/utils/MAS1998BenchmarkUtils.C`: published constants,
  species metadata, vertical atmosphere values, and initial profiles.
- `include/utils/MAS1998ReducedGrid.h`, `src/utils/MAS1998ReducedGrid.C`: latitude-ring layout and
  longitude mapping at equal/coarse/fine interfaces.
- `include/utils/MAS1998TransportOperator.h`, `src/utils/MAS1998TransportOperator.C`: conservative
  limited third-order fluxes and two-stage explicit trapezoidal RK.
- `include/utils/MAS1998VerticalOperator.h`, `src/utils/MAS1998VerticalOperator.C`: Eq. (3.16)
  density-weighted tridiagonal operator and zero-flux boundaries.
- `include/utils/MAS1998ChemistryIntegrator.h`, `src/utils/MAS1998ChemistryIntegrator.C`: implicit
  Euler, variable-step BDF2 coefficients, fixed two-pass scalar Gauss-Seidel, and nitrogen correction.
- `include/utils/MAS1998BenchmarkSolver.h`, `src/utils/MAS1998BenchmarkSolver.C`: state layout,
  Type I/II Strang sequences, emissions, 14-day loop, and CSV diagnostics.
- `include/executioners/MAS1998Executioner.h`, `src/executioners/MAS1998Executioner.C`: MOOSE input
  interface and solver invocation.

Verification and user-facing artifacts:

- `unit/src/MAS1998ReducedGridTest.C`, `unit/src/MAS1998TransportTest.C`,
  `unit/src/MAS1998VerticalOperatorTest.C`, `unit/src/MAS1998ChemistryTest.C`: focused numerical tests.
- `test/tests/mas1998/mas1998_complete_smoke.i`: short full-stack Type I/II run.
- `test/tests/mas1998/production/type_i_64x32.i`, `type_ii_64x32.i`: published coarse-grid 14-day cases.
- `test/tests/mas1998/reference/type_ii_refined.i`: independently generated refined reference recipe.
- `scripts/plot_diagonal.py`: Figure 4-style diagonal plots from solver CSV.
- `README.md`, `test/tests/mas1998/README.md`: build, run, provenance, and limitations.
- `.codex/docs/architecture.md`, `.codex/docs/numerical_mapping.md`: ownership and equation mapping
  after migration to the standalone app.

### Task 1: Extend the Generic KPP Runtime Contract

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/moose/modules/atmospheric_chemistry/include/utils/KPPGeneratedMechanism.h`
- Modify: `/home/fangxiaozhong/git_repo/moose/modules/atmospheric_chemistry/src/utils/KPPGeneratedMechanism.C`
- Exercise: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/mas1998_complete_smoke.i`

- [x] **Step 1: Add the generic KPP APIs**

Add these public methods and the matching private function pointer/state:

```cpp
void setFixedSpecies(const std::string & species, Real concentration);
void clearFixedSpecies();
void computeProductionLoss(Real t,
                           const std::vector<Real> & concentrations,
                           const PhysParams & params,
                           std::vector<Real> & production,
                           std::vector<Real> & loss_coefficient) const;

using KppFunSplitFn = void (*)(double[], double[], double[], double[], double[], double[]);
KppFunSplitFn _fun_split;
std::map<std::string, Real> _fixed_species;
```

`loss_coefficient[j]` is the non-negative `D_VAR[j]` in
`dc[j]/dt = production[j] - loss_coefficient[j] * c[j]`.

- [x] **Step 2: Resolve and evaluate `Fun_SPLIT`**

Resolve `Fun_SPLIT` with `dlsym` without making it mandatory for existing KPP libraries. The new method
must lock the existing recursive mutex, validate all vector sizes, copy concentrations into `VAR`, call
`updateParams`, call `Fun_SPLIT`, and copy `P_VAR` and `D_VAR` to the output vectors. If the symbol is
absent, report a targeted error naming `Fun_SPLIT`.

- [x] **Step 3: Keep pressure in the documented KPP unit**

`PhysParams::pressure` is mbar. Continue converting mbar to Pa only for the ideal-gas density formula,
but set the exported KPP `PRESS` global to mbar:

```cpp
const double pressure_mbar =
    params.pressure > 0.0
        ? static_cast<double>(params.pressure)
        : air_dens * 1.0e6 * 1.380649e-23 * static_cast<double>(params.temperature) / 100.0;
set_global("PRESS", pressure_mbar);
```

- [x] **Step 4: Apply named fixed-species overrides**

After the generic `M/AIR`, `O2`, `N2`, and `H2O` defaults, find each named fixed species in
`_all_species_names[_n_variable.._n_species)` and assign its `FIX` entry. Reject unknown or variable
species in `setFixedSpecies` so misspelled mechanism input cannot silently run.

- [x] **Step 5: Make `computeSpeciesRates` consistent with `IMechanism`**

Use `computeProductionLoss`; return production rates unchanged and convert destruction coefficients to
actual destruction rates with `rates.loss[j] = loss_coefficient[j] * C[j]`.

- [x] **Step 6: Verify formatting and later full-stack coverage**

Run after the app smoke input exists:

```bash
conda run -n moose --no-capture-output make -j4
conda run -n moose --no-capture-output ./run_tests --re=mas1998_complete_smoke -j2
```

Expected: the app compiles against the changed module and the full-stack test evaluates
`Fun_SPLIT` with height-dependent fixed species.

### Task 2: Complete Benchmark Constants and Species Metadata

**Files:**
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998BenchmarkUtils.h`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998BenchmarkUtils.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998BenchmarkUtilsTest.C`

- [x] **Step 1: Add published grid/species constants**

Expose immutable accessors for the 15 centers
`0.3,1.0,2.2,4.3,6.5,8.4,10.0,11.3,13.0,15.2,17.6,19.8,22.5,27.6,34.7`,
the 38.2 km top, the 17 KPP variable species in generated order, and nitrogen atom counts:

```text
O1D CH4 HNO2 H2O2 N2O5 HNO3 HO2NO2 CH3OOH HCHO CH3O2 NO3 O3P NO OH NO2 O3 HO2
0   0   1    0    2    1    1      0       0    0     1   0   1  0  1   0  0
```

- [x] **Step 2: Add standard-atmosphere temperature and pressure**

Reuse the existing 1976-standard-atmosphere layer calculation so `temperature(z)` and
`pressureMbar(z)` are internally consistent with `airNumberDensity(z)`. Do not duplicate a second
atmosphere table.

- [x] **Step 3: Centralize the cylinder profile**

Add `initialConcentration(species, longitude_degrees, latitude_degrees, z_km,
plume_radius_degrees)`. It uses the Table 1 ground value, the HNO3/NO cylinder centered at `(0,0)`, and
the density ratio so each initial vertical mixing ratio is constant.

- [x] **Step 4: Test exact metadata and atmosphere consistency**

Verify 17 names, nitrogen weights, 15 heights, 38.2 km top, sea-level pressure near 1013.25 mbar,
temperature 288.15 K, and that all initial concentrations divided by density are layer-independent.

### Task 3: Implement the Paper Reduced Grid

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998ReducedGrid.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998ReducedGrid.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998ReducedGridTest.C`

- [x] **Step 1: Define ring and cell geometry**

Store south-to-north rings with latitude bounds/center, longitude count/spacing, and flat offset. Provide
periodic longitude indexing plus `cellAtLongitude(row, lambda)` for piecewise-constant virtual values.

- [x] **Step 2: Encode the three paper schedules**

Use symmetric hemispheres with the following pole-to-equator ring counts:

```text
64x32:   4,8,16,32, then 12 rings of 64
128x64:  4,8,16,32, then 8 rings of 64, then 20 rings of 128
256x128: 4,8,16,32,64, then 17 rings of 128, then 42 rings of 256
```

Also support unreduced regular even grids for small CI and refined-reference experiments.

- [x] **Step 3: Represent coarse/fine latitude interfaces**

For each adjacent ring pair, divide the shared boundary into `max(n_south,n_north)` equal segments.
Each segment stores the south and north real-cell indices selected by its longitude midpoint. This is
piecewise-constant interpolation; summing segment fluxes gives the coarse-cell flux required by the CWI
reduced-grid report.

- [x] **Step 4: Verify Table 2 exactly**

Tests must assert reduced 3D cell counts `24840`, `93960`, `391560` and 17-species unknown counts
`422280`, `1597320`, `6656520`. Also assert periodic seam mapping, north/south symmetry, and that every
interface segment maps to valid cells.

### Task 4: Implement Conservative Spherical Advection

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998TransportOperator.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998TransportOperator.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998TransportTest.C`

- [x] **Step 1: Implement Eq. (3.7)-(3.9) reconstruction**

Use the generic upwind triplet `(far_upwind, upwind, downwind)` and

```cpp
theta = (upwind - far_upwind) / (downwind - upwind);
limiter = std::max(0.0, std::min({1.0, theta, 1.0 / 3.0 + theta / 6.0}));
face = upwind + limiter * (downwind - upwind);
```

Use the upwind value when the denominator is at most `1e-30`; the tolerance matches the existing app
FV kernel's roundoff guard.

- [x] **Step 2: Accumulate each physical face flux once**

Longitude faces are periodic within each ring. Latitude faces use the interface segments from Task 3.
Accumulate equal and opposite integrated fluxes into adjacent cells, including `cos(phi_face)` on
meridional fluxes, then divide by the cell's `a^2 cos(phi_center) d_lambda d_phi` measure. Polar boundary
fluxes are zero because `cos(phi_face)=0`.

- [x] **Step 3: Implement Eq. (3.11) explicit trapezoidal RK**

For every species/layer horizontal field, compute
`w = c + dt*f(c)` and `c_new = c + 0.5*dt*(f(c)+f(w))`. Clamp only negative roundoff below a documented
absolute threshold; do not hide a material positivity failure.

- [x] **Step 4: Test conservation, constants, seam, reduction, and CFL**

Assert a constant field has zero RHS, the volume-weighted sum of RHS is roundoff zero, a cylinder crossing
the longitude seam remains conservative/non-negative, and the paper steps `2400/1200/600 s` satisfy the
published `max(nu_lambda+nu_phi) <= 2/3` on their respective grids.

### Task 5: Implement the 15-Layer Vertical Operator

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998VerticalOperator.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998VerticalOperator.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998VerticalOperatorTest.C`

- [x] **Step 1: Assemble Eq. (3.16)**

Build lower/diagonal/upper coefficients for
`d_z[rho K d_z(c/rho)]`. Evaluate `rho*K` halfway between centers, use meters in all differences, ghost
centers reflected across 0 and 38.2 km, and set the exterior interface coefficients to zero for Eq. (2.5).

- [x] **Step 2: Provide operator application and tridiagonal solve**

Expose `apply(c)` and `solveImplicit(tau, loss, rhs)`, where the latter solves
`(I - tau*A + tau*diag(loss))*c = rhs` by the Thomas algorithm. Validate all diagonals remain positive.

- [x] **Step 3: Verify the vertical invariants**

Assert constant mixing ratio `c/rho` has zero diffusion, the layer-integrated tendency is conservative,
zero `K` is identity under implicit solve, and a nonuniform mixing-ratio pulse is smoothed without
creating negative concentrations.

### Task 6: Implement the Paper Chemistry/Diffusion Integrator

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998ChemistryIntegrator.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998ChemistryIntegrator.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/unit/src/MAS1998ChemistryTest.C`

- [x] **Step 1: Implement Euler and variable-step BDF2 coefficients**

At each split-subproblem start, use implicit Euler (`alpha=1`, `C=c_n`). For subsequent steps use Eq.
(3.21): `q=h_new/h_old`, `alpha=(1+q)/(1+2q)`,
`C=((1+q)^2*c_n-q^2*c_nm1)/(1+2q)`, and `tau=alpha*h_new`.

- [x] **Step 2: Implement exactly two ordered Gauss-Seidel passes**

Start with Eq. (3.24)'s non-negative extrapolation. For species in KPP generated order, evaluate the
current mixed iterate with `computeProductionLoss`, add `1e4 molecules/(cm^3 s)` to lowest-layer NO
production, and solve Eq. (3.23) across the whole column. Chemistry-only mode uses the same scalar
formula with `A=0`; diffusion-only mode skips KPP and Gauss-Seidel.

- [x] **Step 3: Apply nitrogen correction after each pass**

Before a pass, record nitrogen mass. Afterwards scale all nitrogen compounds by `(mass_before / mass_after)`.
Use each cell independently for Type I chemistry and layer-width-weighted column totals for Type II.
Reject a non-positive denominator instead of silently skipping correction.

- [x] **Step 4: Test equations independent of the production run**

Use a deterministic mock production/loss evaluator to verify Euler, equal/unequal-step BDF2 coefficients,
ordered species updates, exactly two passes, and nitrogen restoration. The KPP-backed smoke test in Task 8
verifies the real mechanism contract.

### Task 7: Assemble Type I and Type II Benchmark Sequences

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/utils/MAS1998BenchmarkSolver.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/utils/MAS1998BenchmarkSolver.C`

- [x] **Step 1: Allocate and initialize the 4D state**

Use flat storage ordered as horizontal cell, vertical layer, species. Initialize all 17 variables with
Task 2's profile. Set KPP fixed `M`, `H2O`, `CO`, and `O2` at each height using constant mixing ratio;
set `MAS1998_PHI` from latitude, `MAS1998_DELTA=0` for the equinox benchmark, and pass local solar time
`GMT + longitude/360*86400` modulo one day.

- [x] **Step 2: Implement Type I Strang splitting**

For each `dt_split=2*dt_adv`, execute
`A(dt/2) D(dt/2) C(dt) D(dt/2) A(dt/2)`. Use one RK advection step per half, diffusion step
`dt_dif=dt_adv/2`, and 300 s chemistry steps.

- [x] **Step 3: Implement Type II Strang splitting**

Execute `A(dt/2) [C+D](dt) A(dt/2)` with one RK step per advection half and 300 s coupled steps.

- [x] **Step 4: Emit reproducible diagnostics**

Write metadata followed by CSV rows sampled along `lambda'=phi'/2`, including time, grid indices,
coordinates, height, all 17 species, and derived `NOx=NO+NO2`. Also write global volume-weighted species
totals and minimum concentrations so conservation/positivity can be audited.

### Task 8: Add the MOOSE Executioner and Input Decks

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/include/executioners/MAS1998Executioner.h`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/src/executioners/MAS1998Executioner.C`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/mas1998_complete_smoke.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/production/type_i_64x32.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/production/type_ii_64x32.i`
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/reference/type_ii_refined.i`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/tests`

- [x] **Step 1: Expose only benchmark controls**

Required/validated parameters are mechanism library, grid size/reduction, split type, end time,
`dt_adv`, chemistry step, Gauss-Seidel passes, cylinder radius, output file, and output interval. Defaults
match the `64x32x15` paper case: Type II, 14 days, `dt_adv=2400`, chemistry step 300 s, two passes.

- [x] **Step 2: Create a short full-stack regression**

Use a regular `8x4x15` grid, 600 s end time, 300 s advection, 600 s split, and 300 s chemistry. Run both
Type I and Type II, assert CSV creation and finite/non-negative summary values, and mark the KPP build as
a prerequisite.

- [x] **Step 3: Create the two paper production inputs**

Both inputs use reduced `64x32x15`, 14 days, `dt_adv=2400`, `dt_split=4800`, 300 s chemistry minimum,
two Gauss-Seidel passes, NO emission, and the recovered 45-reaction mechanism. They differ only in the
published Type I/II sequence.

- [x] **Step 4: Create the replacement-reference recipe**

Use Type II with regular `256x128x15`, a smaller split/advection step, smaller chemistry step, and more
Gauss-Seidel passes. Name its output `independent_refined_reference.csv` and state in the input comments
that it is not the unavailable CWI file.

### Task 9: Add Plotting and Documentation

**Files:**
- Create: `/home/fangxiaozhong/git_repo/mas1998_benchmark/scripts/plot_diagonal.py`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/README.md`
- Modify: `/home/fangxiaozhong/git_repo/mas1998_benchmark/test/tests/mas1998/README.md`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/architecture.md`
- Modify: `/home/fangxiaozhong/git_repo/moose/.codex/docs/numerical_mapping.md`

- [x] **Step 1: Plot the paper diagnostics**

Read solver CSV with Python's standard `csv` module and plot ground-level O3, NOx, HNO3, and HO2NO2
against latitude along `lambda'=phi'/2`. Accept multiple CSVs so Type I, Type II, and refined reference
can be overlaid. Require matplotlib and fail with a direct installation message if absent.

- [x] **Step 2: Replace the old app scope statement**

Document complete model components, exact production commands, expected Table 2 sizes, output schema,
the NM-R9505 chemistry provenance, and the three evidence limitations at the top of this plan.

- [x] **Step 3: Rewrite MOOSE ownership documentation**

`architecture.md` must show MAS1998-specific code under the standalone app and only generic KPP/transport
facilities in the module. `numerical_mapping.md` must map Eqs. (3.5)-(3.24), (4.1)-(4.2), both split
sequences, Table 1, Table 2, and Table 3 time steps to exact app classes/tests.

### Task 10: Build, Verify, Commit, and Rebase

**Files:** all files above; preserve unrelated
`/home/fangxiaozhong/git_repo/moose/HETerogeneous-vectorized-or-Parallel/`.

- [x] **Step 1: Build the changed module and app in the established environment**

```bash
conda run -n moose --no-capture-output make -j4
```

Run from `modules/atmospheric_chemistry`, then from `mas1998_benchmark`.

- [x] **Step 2: Run focused unit and TestHarness checks**

```bash
conda run -n moose --no-capture-output make -C unit -j4
conda run -n moose --no-capture-output ./unit/run_tests
conda run -n moose --no-capture-output ./run_tests -j4
conda run -n moose --no-capture-output ./run_tests --heavy -j2
```

Expected: all numerical unit tests, migrated smoke tests, KPP build/box tests, and complete Type I/II smoke
tests pass. The 14-day and refined cases are syntax-checked but not run in CI due their cost.

- [x] **Step 3: Check both diffs and input syntax**

```bash
git -C /home/fangxiaozhong/git_repo/moose diff --check
git -C /home/fangxiaozhong/git_repo/mas1998_benchmark diff --check
conda run -n moose --no-capture-output ./mas1998_benchmark-opt \
  -i test/tests/mas1998/production/type_i_64x32.i --check-input
conda run -n moose --no-capture-output ./mas1998_benchmark-opt \
  -i test/tests/mas1998/production/type_ii_64x32.i --check-input
conda run -n moose --no-capture-output ./mas1998_benchmark-opt \
  -i test/tests/mas1998/reference/type_ii_refined.i --check-input
```

- [x] **Step 4: Commit each repository with scoped paths**

Create one MOOSE commit for the generic KPP extension and migrated architecture docs, and one app commit
for the benchmark solver/tests/docs. Do not add generated KPP products, output CSVs, plots, or unrelated
untracked paths.

- [x] **Step 5: Rebase each branch onto its latest configured remote**

Fetch first. Rebase the MOOSE branch onto `upstream/devel` and the app branch onto `origin/main`. If a
remote advanced into overlapping files, stop at the conflict, preserve both intended changes, rerun the
focused verification, then report the resulting commit IDs.

## Self-Review

- Spec coverage: tasks cover the 17-variable/45-reaction model, fixed species, emissions, initial data,
  reduced grids, spherical fluxes, RK2, 15-layer diffusion, Euler/BDF2, two Gauss-Seidel passes, nitrogen
  correction, both splittings, 14-day decks, diagnostics, replacement reference, and migrated docs.
- Placeholder scan: no deferred implementation placeholder is present. Unavailable historical data and
  incompletely published controller parameters are handled as explicit evidence boundaries rather than
  guessed gold data.
- Type consistency: KPP production/loss names and units are consistent across the module interface,
  chemistry integrator, and tests; state order follows generated KPP order everywhere.
- Scope control: generic runtime support stays in Atmospheric Chemistry; every MAS1998 numerical policy,
  object, input, and diagnostic remains in the standalone app.
