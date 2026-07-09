# Atmospheric Chemistry Module

The atmospheric chemistry module provides tools for 0-D box model and spatially-resolved
atmospheric chemistry simulations using Master Chemical Mechanism (MCM) chemical mechanisms.

## Overview

The module uses dedicated blocks under `[AtmosphericChemistry]`:

1. **Box mode** (`[Box]`) — 0-D ODE using ScalarVariable + ChemistryODEKernel + MCMBoxModel.
   Suitable for large mechanisms (up to full MCM ~5832 species).
2. **Coupled mode** (`[Coupled]`) — FEM transport + chemistry using MooseVariableFE +
   ChemicalSourceKernel + MCMRatesMaterial. Suitable for spatially-resolved simulations (5--50 species).

## Quick Start

```moose
# Box mode (0-D ODE, ScalarVariable)
[AtmosphericChemistry]
  [Box]
    mechanism_file = 'mechanism.fac'
    temperature = 298
  []
[]

# Coupled mode (FEM transport + chemistry)
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'mechanism.fac'
    temperature = 298
    include_transport = true
  []
[]
```

## Concentration Units

The `units` parameter (default: `molec_cm3`) controls input/output units:

```moose
[AtmosphericChemistry]
  [Box]
    units = ppb           # input ICs in ppb, output in ppb
    air_density = 2.46e19
    ...
  []
[]
```

When `units = ppb`:
- Species ICs are specified in ppb (auto-converted to molec/cm³ internally)
- The solver works in molec/cm³; results are converted back to ppb for output
- Conversion factor: `M / 1e9` where `M` = air number density (molec/cm³)
- For box mode: conversion is dynamic (updates with T/P changes)
- For coupled mode: conversion is evaluated at each quadrature point

Default `molec_cm3` (no conversion) preserves backward compatibility.
See `MCMBoxModel` and `ChemistryODEKernel` documentation for implementation details.

## Photolysis Schemes

Three photolysis calculation methods are supported via the `photolysis_scheme` parameter:

| Scheme | Input | Dependencies | Use Case |
|--------|-------|-------------|----------|
| `MCM_SZA` (default) | SZA only (lat/lon/day/time) | MCM photolysis-rates file | AtChem2-compatible, fast |
| `HYBRID` | SZA + albedo + O3col + altitude | F0AM TUV lookup tables | O3 column / albedo dependency |
| `BOTTOMUP` | Lamp flux + CS/QY data | Cross-section + quantum yield files | Lab chamber experiments |

```moose
# Hybrid scheme example
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'mechanism.fac'
    photolysis_scheme = HYBRID
    hybrid_table_dir = 'tuv_tables'
    albedo = 0.1
  []
[]

# BottomUp (chamber) scheme example
[AtmosphericChemistry]
  [Box]
    mechanism_file = 'mechanism.fac'
    temperature = 298.0
    press = 1013.0
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = 'database/photolysis/bottomup'
  []
[]
```

The BottomUp photolysis engine is implemented natively in C++ (`BottomUpJIntegrator`),
reproducing F0AM's photolysis computation without external tools. The workflow is:

```
bottomup_jmap.dat  (reaction map, cs_type/qy_type=10 for built-in formulas)
         │
BottomUpJIntegrator::computeJ(jname, T, P)
  → computeCS_builtin(species, T, P)     // 18 CS formulas ported from F0AM
  → computeQY_builtin(species, T, P)     // 13 QY formulas ported from F0AM
  → smear() + trapz()                     // IntegrateJ.m numerical integration
  → J-value
```

The `database/photolysis/bottomup/` directory contains:

| File type | Purpose | Source |
|-----------|---------|--------|
| `bottomup_jmap.dat` | Reaction mapping (J-name → CS/QY formula) | Built-in, `cs_type/qy_type=10` |
| `CrossSections/*.csv` | Raw cross-section data (wavelength, sigma values) | Copied from F0AM `Chem/Photolysis/CrossSections/` |
| `QuantumYields/*.csv` | Raw quantum-yield data (wavelength, QY values) | Copied from F0AM `Chem/Photolysis/QuantumYields/` |
| `ExampleLightFlux.txt` | Lamp spectrum for chamber simulations | Copied from F0AM `Setups/Examples/` |

The raw CSV files serve as the base data; temperature and pressure corrections
are applied at runtime by the C++ built-in formulas (e.g., linear T-correction
for HCHO, Arrhenius for N2O5, 4-regime spectroscopic for O3, etc.).
No Python or MATLAB precomputation is needed.

The chamber validation path follows F0AM's `J_BottomUp` implementation closely:
HNO3 and NO3 cross-section scaling, O3 O(1D)/O(3P) quantum yields, H2O2
temperature-dependent cross sections, acetone branching yields, and glyoxal /
methylglyoxal quantum yields are evaluated from the same spectral data and
formula conventions used by F0AM. This matters for `MCMv331_Inorg_Isoprene`,
where incorrect J-values propagate rapidly into NO3, HNO3, and isoprene
oxidation products.

### Internal: cs_type / qy_type system

The `bottomup_jmap.dat` file assigns each reaction a `cs_type` and `qy_type`
that tells `BottomUpJIntegrator` how to compute the cross-section and quantum yield.
These are **internal implementation details**, not user-facing parameters:

| Type | Meaning | Data source |
|------|---------|-------------|
| `0` | Scalar constant | Value embedded in jmap (e.g., `0.158`) |
| `1` | 2-column CSV | `CrossSections/*.csv` or `QuantumYields/*.csv` |
| `2` | 3-column CSV with T-interpolation | File with columns `[wl, val@T1, val@T2]` |
| `3` | TXT file | Space/tab-delimited `[wl, val]` |
| `10` | **Built-in C++ formula** | C++ `computeCS_builtin` / `computeQY_builtin` |

Switching between types would occur in these scenarios:

1. **New photolysis reaction** added to F0AM that doesn't have a C++ built-in formula yet → use `cs_type=1` with the raw CSV as temporary fallback, then add a C++ formula later.
2. **Override built-in formula** with custom measured data → change to `cs_type=1` and point to a user-supplied CSV file.
3. **TXT-formatted data** from external sources (e.g., JPL-recommended spectral files) → `cs_type=3`.

The current jmap uses `type=10` for all reactions where F0AM's MATLAB formulas
have been ported to C++ (18 cross-sections, 13 quantum yields), and `type=1` or
`type=3` for species that only need a simple CSV read (e.g., MACR, MEK) or
halogen text files from JPL/IUPAC references.

## RO2 Species Extraction

The parser detects peroxy radical (RO₂) species using explicit family declarations
when available: `RO2 = CH3O2 + ... ;` in `.fac` files or the corresponding
`RO2 = & ... )` section in `.kpp` files. If no explicit family declaration exists,
it falls back to species-name detection (`O2` suffix or names containing `RO2`,
excluding known false positives such as HO2, NO2, SO2, and H2O2). This fallback
covers simplified mechanisms that declare `RO2` directly as a species.

Regenerate the RO₂ CSVDiff gold directly from a `.fac` mechanism's explicit
`RO2 = ... ;` declaration:

```bash
python3 scripts/check_ro2.py mechanism.fac --gold-csv gold/test_ro2_detection_ro2_list_0001.csv
```

To validate the parser output against that explicit declaration, add `--run-app`.
For simplified mechanisms without an explicit declaration, `--run-app` reports the
fallback-detected species names.

## Chemical Mechanism Files

Mechanism files in FACSIMILE (`.fac`) format are supported directly. For MCM website-exported `.kpp`
files, use the included converter:

```bash
python3 modules/atmospheric_chemistry/scripts/kpp_to_fac.py mechanism.kpp
```

## Validation

### AtChem2 Validation (MCM_SZA photolysis)

The MCM_SZA photolysis scheme (default) has been validated against AtChem2 (v1.2.1)
for the MCM v3.3.1 inorganic subset (29 species, 71 reactions, 12 h diurnal cycle).
Both box and coupled modes produce identical photolysis rates and match AtChem2's
CVODE reference:

| Mode | Steps | Converged | Photolysis |
|------|-------|-----------|------------|
| `vs_AtChem2_inorg_box` | 480 | 480/480 | MCM_SZA |
| `vs_AtChem2_inorg_coupled` | 480 | 480/480 | MCM_SZA |

O$_3$, NO$_2$, and NO concentrations match to within 0.1--0.3 % relative error
during daytime and < 0.05 % at night. Solar parameters (declination, zenith angle)
follow the Madronich (1993) formulation identically to AtChem2.

To regenerate AtChem2 validation data:

```bash
python3 scripts/generate_atchem2_gold.py \
    --atchem2 <path-to-AtChem2>/model/output \
    --output test/tests/actions/gold/vs_AtChem2_inorg_box.csv

python3 scripts/plot_atchem2.py \
    --moose test/tests/actions/vs_AtChem2_inorg_box.csv \
    --atchem2 <path-to-AtChem2>/model/output
```

### F0AM Chamber Validation (610 species, BOTTOMUP photolysis)

The `MCMv331_Inorg_Isoprene` mechanism (610 species, 1974 reactions) is validated
against F0AM v4.4's `ExampleSetup_Chamber.m` using the BOTTOMUP photolysis scheme
(cross-section × quantum-yield × lamp-flux integration, ported from F0AM's
MATLAB functions to C++ `BottomUpJIntegrator`). Four scenarios matching F0AM's
3-step + 1-event simulation are included:

| Scenario | NO₂ (ppb) | jcorr | Duration | F0AM source |
|----------|-----------|-------|----------|-------------|
| S1_box   | 0.1       | 1     | 3 h      | `S.StepIndex == 1` |
| S2_box   | 1         | 1     | 3 h      | `S.StepIndex == 2` |
| S3_box   | 10        | 1     | 3 h      | `S.StepIndex == 3` |
| S2b_box  | 1         | 10    | 1 h      | Restart from S2 end state |

**Gold CSV generation** — run the F0AM chamber export script from MATLAB:

```bash
# From the MOOSE repository root. Set F0AM_ROOT if the default .reasonix path
# is not available on your machine.
matlab -batch "addpath('modules/atmospheric_chemistry/scripts'); export_chamber_gold"
```

The export script:
1. Reads air number density from `S.Met.M` (dynamic, not hardcoded)
2. Extracts S1/S2/S3 by `StepIndex`, then runs S2b as a restart from S2
3. Removes F0AM pseudo-species ONE and CH3ONO (not in MOOSE mechanism);
   keeps RO2 (now output as a diagnostic ScalarVariable in box mode)
4. Reorders columns to match MOOSE output order (611 columns: 610 species + RO2)
5. Writes to `gold/vs_F0AM_chamber_{S1,S2,S3,S2b}_box.csv`
6. Prints the F0AM compute time for the current four-test chamber workflow
   (`S1 + S2 + S2b + S3`)

The MOOSE chamber tests use PETSc TS BDF with `chem_solver_rtol = 5e-4` and
`chem_solver_atol = 5e-8`. The CSVDiff comparison keeps a 5% relative tolerance
and uses `abs_zero = 1.15e-6` to ignore roundoff-scale species concentrations
near zero. S2 and S3 stop at the final F0AM output times in their time sequences
(`9957 s` and `9832 s`) instead of a nominal 3-hour endpoint.

### Running Tests

All chamber tests are marked `heavy = true` (610-species ODE with LU factorization,
~30-180s per scenario for PETSc TS BDF on a laptop-class CPU). Run the full suite:

```bash
cd /path/to/moose
./run_tests --heavy -j1 modules/atmospheric_chemistry
```

Or run the chamber validation subset:

```bash
cd modules/atmospheric_chemistry
./run_tests -C test/tests/chamber --heavy --re 'vs_F0AM_chamber_(S1|S2|S2b|S3)_box'
```

Or run a single scenario directly:

```bash
cd modules/atmospheric_chemistry/test/tests/chamber
../../../atmospheric_chemistry-opt -i vs_F0AM_chamber_S1_box.i
```

### Solver Runtime Benchmark

F0AM chamber timing should be compared against the same four work items used by
the current MOOSE chamber tests: `S1`, `S2`, `S2b` restart from the S2 checkpoint,
and `S3`. The reference F0AM timing of about 6 seconds is therefore the total
time for those four items, not a per-scenario value.

The benchmark helper generates solver-specific copies of the current chamber
inputs under `test/tests/chamber/benchmark_runs/`, runs each scenario, and writes
both per-scenario and total timing CSV files:

```bash
cd modules/atmospheric_chemistry

# List solver selections covered by the benchmark.
python3 scripts/benchmark_chamber_solvers.py --list-solvers

# Run all directly supported FAC solvers and generated KPP solvers.
python3 scripts/benchmark_chamber_solvers.py --solvers all --f0am-seconds 6

# Run only the current production solver.
python3 scripts/benchmark_chamber_solvers.py --solvers petsc_bdf --f0am-seconds 6

# Run the generated KPP Rosenbrock mechanism with MATLAB/F0AM default tolerances.
python3 scripts/benchmark_chamber_solvers.py \
    --solvers kpp_rosenbrock --rtol 1e-3 --atol 1e-6 --f0am-seconds 6
```

The full solver sweep is a performance experiment and can run for many minutes
or hours depending on the host. Use `--timeout <seconds>` to cap each scenario.
The script writes `chamber_solver_timing.csv` and `chamber_solver_summary.csv`
incrementally, so partial results are preserved if the run is interrupted.

For `.fac` chamber inputs, the directly runnable solver set is MOOSE implicit,
SUNDIALS CVODE, and PETSc TS with each exposed `chem_solver_type` (`bdf`,
`arkimex`, `eimex`, `rosw`, `mimex`, `beuler`, `cn`, `rk`, `theta`, `ssp`, and
`sundials`). KPP selections first convert `MCMv331_Inorg_Isoprene.fac` to a
KPP mechanism with `scripts/fac_to_kpp.py`, build the selected KPP integrator
shared library under `test/tests/chamber/kpp_chamber/`, then run the same S1,
S2, S2b, and S3 inputs against the generated `.kpp` mechanism. These generated
KPP files are ignored by git and can be recreated by the benchmark script.

`fac_to_kpp.py` filters KPP dummy species and FAC species that are declared but
not present in the reaction graph, keeping the generated `#DEFVAR` list aligned
with KPP's generated `NVAR` and `SPC_NAMES`.

On the reference workstation used during development, `kpp_rosenbrock` completed
the four chamber work items in 9.38 s with `rtol=1e-3` and `atol=1e-6`
(1.58x the 6 s F0AM total). With the stricter chamber-test tolerances
`rtol=5e-4` and `atol=5e-8`, the same KPP run completed in 11.94 s.

The lightweight test `chamber_solver_timing_inputs` checks that benchmark inputs
can be generated for the current four chamber test files without running the full
solver sweep:

```bash
cd modules/atmospheric_chemistry
./run_tests -C test/tests/chamber --re chamber_solver_timing_inputs
```

### Comparing Results

**Numerical comparison** (per-species ratio report):

```bash
# Checked-in comparison scenarios (S1, S2, S3)
python3 scripts/compare_chamber_f0am.py

# Single scenario
python3 scripts/compare_chamber_f0am.py \
    --moose test/tests/chamber/vs_F0AM_chamber_S1_box.csv \
    --gold test/tests/chamber/gold/vs_F0AM_chamber_S1_box.csv
```

**Visual comparison** (2x2 panel PDFs):

```bash
# Check files first
python3 scripts/plot_chamber_comparison.py --check

# Plot checked-in comparison scenarios (S1, S2, S3)
python3 scripts/plot_chamber_comparison.py

# Single scenario, custom output
python3 scripts/plot_chamber_comparison.py \
    --moose test/tests/chamber/vs_F0AM_chamber_S1_box.csv \
    --gold test/tests/chamber/gold/vs_F0AM_chamber_S1_box.csv \
    --save chamber_S1.pdf
```

The comparison scripts convert MOOSE's molec/cm³ output to ppb using
`M = 2.46e19 molecules/cm³` (standard conditions at 298 K, 1013 mbar).
Use the single-scenario options for S2b after generating an S2b gold CSV.

### Reproducing F0AM Figures

`ExampleSetup_Chamber.m` also creates F0AM-only figures with `PlotConc`,
`PlotConcGroup`, `PlotRates`, `PlotRatesAvg`, `PlotReactivity`, and `PlotYield`.
Use the MOOSE-side MATLAB helper to rerun the same F0AM chamber setup, save those
figures, and record the F0AM compute time for `S1 + S2 + S2b + S3`:

```bash
matlab -batch "addpath('modules/atmospheric_chemistry/scripts'); plot_f0am_chamber_figures"
```

The helper writes PNG/FIG outputs and `f0am_chamber_timing.csv` under
`test/tests/chamber/f0am_figures/` by default. Pass explicit arguments to use a
different F0AM checkout or output directory:

```matlab
plot_f0am_chamber_figures('/path/to/F0AM', '/path/to/output')
```

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `scripts/kpp_to_fac.py` | Convert `.kpp` → `.fac` |
| `scripts/fac_to_kpp.py` | Convert FACSIMILE chamber mechanisms to KPP C mechanisms for generated KPP solver benchmarks |
| `scripts/check_ro2.py` | Extract RO2 species names from a `.fac` mechanism's explicit `RO2 = ... ;` declaration, write CSVDiff gold with `--gold-csv`, and optionally compare parser output with `--run-app`. |
| `scripts/generate_atchem2_gold.py` | Generate gold CSV from AtChem2 reference output |
| `scripts/plot_atchem2.py` | Plot MOOSE vs AtChem2 comparison (3×3 grid PDF) |
| `scripts/plot_vs_f0am.py` | Plot MOOSE vs analytical solution for F0AM tutorial |
| `scripts/gen_tutorial_gold.py` | Generate F0AM tutorial gold CSV from scipy ODE |

| `scripts/convert_soas_mat.py` | Convert SOAS campaign .mat observation files to CSV |
| `scripts/export_chamber_gold.m` | MATLAB script to run F0AM chamber scenarios and export chamber gold CSVs |
| `scripts/compare_chamber_f0am.py` | Numerical comparison (per-species ratio) of MOOSE vs F0AM gold CSV |
| `scripts/plot_chamber_comparison.py` | Visual comparison (2x2 panel PDFs) of MOOSE vs F0AM gold CSV |
| `scripts/benchmark_chamber_solvers.py` | Generate/run chamber solver runtime comparisons for S1, S2, S2b, and S3 |
| `scripts/run_ts_sweep.sh` | Wrapper for `benchmark_chamber_solvers.py` |
| `scripts/plot_f0am_chamber_figures.m` | MATLAB script to reproduce the F0AM ExampleSetup_Chamber figures and timing |

## Database

Pre-converted mechanism files are available in `doc/content/modules/atmospheric_chemistry/database/`:

| File | Species | Reactions |
|------|---------|-----------|
| `tutorial_5sp.fac` | 5 | 6 |
| `atchem2_example.fac` | 29 | 71 |
| `mcm_export.fac` | 20 | 48 |
| `mcm_export_all.fac` | 5832 | 17224 |
| `MCMv331_Inorg_Isoprene.fac` | 610 | 1974 |
| `MCMv331_DielExampleChemistry.fac` | 2908 | 8797 |
| `SOAS_DielCycle.csv` | — | 24h × 61 fields (SOAS campaign observations) |

## Test Suite

| Test | Type | Species | Description |
|------|------|---------|-------------|
| `vs_F0AM_tutorial5_box` | CSVDiff | 5 | Box mode vs analytical solution |
| `vs_F0AM_tutorial5_coupled` | CSVDiff | 5 | Coupled FEM vs F0AM reference |
| `vs_AtChem2_inorg_box` | CSVDiff | 29 | Box mode vs AtChem2 CVODE |
| `vs_AtChem2_inorg_coupled` | CSVDiff | 29 | Coupled FEM vs AtChem2 |
| `vs_F0AM_tutorial5_ns_fe` | Exodiff | 5 | NS + chemistry fully coupled |
| `vs_F0AM_chamber_S1_box` | CSVDiff | 610 | S1 (NO₂=0.1ppb): BOTTOMUP photolysis vs F0AM |
| `vs_F0AM_chamber_S2_box` | CSVDiff | 610 | S2 (NO₂=1ppb): BOTTOMUP photolysis vs F0AM |
| `vs_F0AM_chamber_S3_box` | CSVDiff | 610 | S3 (NO₂=10ppb): BOTTOMUP photolysis vs F0AM |
| `vs_F0AM_chamber_S2b_box` | RunApp | 610 | S2b (restart, jcorr=10): BOTTOMUP continuation from S2 |
| `vs_F0AM_dielcycle_box` | RunApp | 2908 | SOAS diel cycle: MCM_SZA |
| `vs_atchem2_transport_building` | Exodiff | 29 | NS + chemistry + building emission |
| `test_family_conservation` | RunApp | 29 | Family conservation (NOx = NO2 + NO) |
| `box_chamber` | — | 610 | F0AM chamber example (3h, MCM_SZA) |
| `box_dielcycle` | — | 2908 | F0AM diel cycle example (24h) |

## Performance

### PETSc TS Standalone Integrator

An optional standalone ODE integrator using PETSc's TS module is available for
box-mode simulations, achieving dramatic speedups over the default MOOSE LU solver.

| Feature | Description |
|---------|-------------|
| Standalone ODE integrator | PETSc TS (BDF / ARKIMEX / SUNDIALS), adaptive step size, 28-32× speedup |
| Configuration | `box_solver = true` with `box_solver_type`, `box_solver_rtol`, `box_solver_atol` |
| Limitation | Box mode only; coupled mode with transport triggers input parsing error |

For the S1 chamber test case (610 species, 1974 reactions, 3h simulation):

| Integrator | Wall Time | vs LU Baseline | C5H8 Error |
|------------|-----------|---------------|------------|
| MOOSE LU (baseline) | 477.8s | 1× | — |
| PETSc TS BDF, rtol=1e-2 | 17.3s | **28×** | 1.5% |
| PETSc TS BDF, rtol=1e-1 | 14.8s | **32×** | 1.5% |

## Objects

- [`AtmosphericChemistryAction`](source/actions/AtmosphericChemistryAction.md) — Dedicated `[Box]` and `[Coupled]` action blocks
- [`MechanismLoader`](source/utils/MechanismLoader.md) — Standalone mechanism loading utility (path resolution + .fac parsing + photolysis set loading)
- [`BoxIntegrator`](source/utils/BoxIntegrator.md) — Strategy interface for box-model integration (MOOSE implicit / PETSc TS / KPP)
- [`MCMBoxModel`](source/userobjects/MCMBoxModel.md) — 0-D chemical ODE engine with caching interface, pluggable sparse matrix backends (CSR/COO/DENSE/CSC), and optional F0AM-style limiting-reagent (LR) formulation for RO₂ termination
- [`ChemistryODEKernel`](source/kernels/ChemistryODEKernel.md) — Box mode ScalarKernel bridge to BoxIntegrator strategy
- [`ChemicalSourceKernel`](source/kernels/ChemicalSourceKernel.md) — FEM chemical source with analytical Jacobian + CSR sparse reactant matrix
- [`MCMRatesMaterial`](source/materials/MCMRatesMaterial.md) — Runtime rate evaluation (coupled mode) with roof_open photolysis switch
- [`MCMConstraintKernel`](source/kernels/MCMConstraintKernel.md) — Fixed-species constraint
- [`MCMEmissionKernel`](source/kernels/MCMEmissionKernel.md) — Emission source
- [`MCMDepositionKernel`](source/kernels/MCMDepositionKernel.md) — Dry deposition
- [`MCMFacsimileParser`](source/userobjects/MCMFacsimileParser.md) — Standalone `.fac` parser
- [`HybridJTableReader`](source/utils/HybridJTableReader.md) — TUV photolysis 4D interpolation
- [`BottomUpJIntegrator`](source/utils/BottomUpJIntegrator.md) — F0AM BottomUp photolysis: cross-section × QY × lamp-flux integration
- [`MCMSolarPostprocessor`](source/postprocessors/MCMSolarPostprocessor.md) — Solar parameters (cosx, secx, lha, etc.)
- [`MCMPhotolysisPostprocessor`](source/postprocessors/MCMPhotolysisPostprocessor.md) — Individual photolysis J-values
- [`MCMRO2Kernel`](source/kernels/MCMRO2Kernel.md) — RO2 algebraic kernel (pins RO2 to sum of peroxy radicals)
- [`MCMRO2Postprocessor`](source/postprocessors/MCMRO2Postprocessor.md) — Total peroxy radical (RO2) concentration
- [`MCMRO2ListPostprocessor`](source/postprocessors/MCMRO2ListPostprocessor.md) — RO2 species count for CSVDiff validation
- [`MCMFamilyConstraint`](source/userobjects/MCMFamilyConstraint.md) — Chemical family conservation (DAE method)
- [`MCMFamilyScalarKernel`](source/kernels/MCMFamilyScalarKernel.md) — Family conservation ScalarKernel (DAE slack variable)
- [`MCMReactionRatePostprocessor`](source/postprocessors/MCMReactionRatePostprocessor.md) — Single reaction rate output
- [`MCMSpeciesLossRatePostprocessor`](source/postprocessors/MCMSpeciesLossRatePostprocessor.md) — Species loss rate diagnostic
- [`MCMSpeciesProductionRatePostprocessor`](source/postprocessors/MCMSpeciesProductionRatePostprocessor.md) — Species production rate diagnostic
- [`MCMLifetimePostprocessor`](source/postprocessors/MCMLifetimePostprocessor.md) — Chemical lifetime (tau = C / loss_rate)
- [`JCalibrator`](source/utils/JCalibrator.md) — JFAC auto-calibration from observed J values (F0AM jcorr)
- [`MCMInputInterpolator`](source/utils/MCMInputInterpolator.md) — Time-interpolation for input data (F0AM InputInterp)
- [`EKMAScanner`](source/utils/EKMAScanner.md) — EKMA sensitivity scan (VOC × NOx grid)

!syntax complete groups=AtmosphericChemistryApp level=3
