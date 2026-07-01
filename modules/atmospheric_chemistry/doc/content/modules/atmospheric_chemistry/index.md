# Atmospheric Chemistry Module

The atmospheric chemistry module provides tools for 0-D box model and spatially-resolved
atmospheric chemistry simulations using Master Chemical Mechanism (MCM) chemical mechanisms.

## Overview

The module uses a unified `[AtmosphericChemistry]` Action with two modes:

1. **Box mode** (`mode = box`) — 0-D ODE using ScalarVariable + ChemistryODEKernel + MCMBoxModel.
   Suitable for large mechanisms (up to full MCM ~5832 species).
2. **Coupled mode** (`mode = coupled`) — FEM transport + chemistry using MooseVariableFE +
   ChemicalSourceKernel + MCMRatesMaterial. Suitable for spatially-resolved simulations (5--50 species).

## Quick Start

```moose
# Box mode (0-D ODE, ScalarVariable)
[AtmosphericChemistry]
  mode = box
  mechanism_file = 'mechanism.fac'
  temperature = 298
[]

# Coupled mode (FEM transport + chemistry)
[AtmosphericChemistry]
  mode = coupled
  mechanism_file = 'mechanism.fac'
  temperature = 298
  include_transport = true
[]
```

## Concentration Units

The `units` parameter (default: `molec_cm3`) controls input/output units:

```moose
[AtmosphericChemistry]
  mode = box
  units = ppb           # input ICs in ppb, output in ppb
  air_density = 2.46e19
  ...
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

Two photolysis calculation methods are supported via the `photolysis_scheme` parameter:

| Scheme | Input | Dependencies | Use Case |
|--------|-------|-------------|----------|
| `MCM_SZA` (default) | SZA only (lat/lon/day/time) | MCM photolysis-rates file | AtChem2-compatible, fast |
| `HYBRID` | SZA + albedo + O3col + altitude | F0AM TUV lookup tables | O3 column / albedo dependency |
| `BOTTOMUP` | Lamp flux + CS/QY data | Cross-section + quantum yield files | Lab chamber experiments |

```moose
# Hybrid scheme example
[AtmosphericChemistry]
  mode = coupled
  mechanism_file = 'mechanism.fac'
  photolysis_scheme = HYBRID
  hybrid_table_dir = 'tuv_tables'
  albedo = 0.1
[]

# BottomUp (chamber) scheme example
[AtmosphericChemistry]
  mode = box
  mechanism_file = 'mechanism.fac'
  temperature = 298.0
  press = 1013.0
  photolysis_scheme = BOTTOMUP
  lamp_flux_file = 'ExampleLightFlux.txt'
  bottomup_data_dir = 'database/photolysis/bottomup'
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

The parser detects peroxy radical (RO₂) species using the same logic as AtChem2:
explicit `RO2 = CH3O2 + ...` section in `.kpp` files, or O₂-suffix heuristic (excluding
known false positives: HO2, NO2, SO2, H2O2 etc.).

Extract and validate the RO₂ list standalone:

```bash
# Extract RO2 from a .fac mechanism and compare against MCM reference
python3 scripts/check_ro2.py \
    doc/content/modules/atmospheric_chemistry/database/mcm_peroxy_radicals_v3.3.1.dat \
    --fac mechanism.fac \
    -o ro2_detected.txt
```

Output: `ro2_detected.txt` (one species per line) + comparison summary (detected / missing / extra).

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

**Gold CSV generation** — from F0AM `.mat` output (no MATLAB required):

```bash
# Prerequisites: F0AM ExampleSetup_Chamber.m must have been run to produce:
#   Runs/ChamberExampleOutput.mat
#   Runs/ChamberExampleHighLightsOutput.mat

python3 scripts/generate_gold_vs_F0AM_chamber.py \
    --mat /path/to/ChamberExampleOutput.mat \
    --mat-s2b /path/to/ChamberExampleHighLightsOutput.mat

# Custom output directory:
python3 scripts/generate_gold_vs_F0AM_chamber.py \
    --gold-dir test/tests/actions/gold
```

The script:
1. Reads air number density from `S.Met.M` (dynamic, not hardcoded)
2. Extracts S1/S2/S3 by `StepIndex`, S2b from the high-lights file
3. Removes F0AM pseudo-species ONE and CH3ONO (not in MOOSE mechanism);
   keeps RO2 (now output as a diagnostic ScalarVariable in box mode)
4. Reorders columns to match MOOSE output order (611 columns: 610 species + RO2)
5. Writes to `gold/vs_F0AM_chamber_{S1,S2,S3,S2b}_box.csv`

### Running Tests

All chamber tests are marked `heavy = true` (610-species ODE with LU factorization,
~30-180s per timestep). Run the full suite:

```bash
cd /path/to/moose
./run_tests --heavy -j1 modules/atmospheric_chemistry
```

Or run a single scenario:

```bash
cd modules/atmospheric_chemistry/test/tests/actions
../../../atmospheric_chemistry-opt -i vs_F0AM_chamber_S1_box.i
```

### Comparing Results

**Numerical comparison** (per-species ratio report):

```bash
# All 4 scenarios
python3 scripts/compare_chamber_f0am.py

# Single scenario
python3 scripts/compare_chamber_f0am.py \
    --moose test/tests/actions/vs_F0AM_chamber_S1_box.csv \
    --gold test/tests/actions/gold/vs_F0AM_chamber_S1_box.csv
```

**Visual comparison** (2x2 panel PDFs):

```bash
# Check files first
python3 scripts/plot_chamber_comparison.py --check

# Plot all 4 scenarios
python3 scripts/plot_chamber_comparison.py

# Single scenario, custom output
python3 scripts/plot_chamber_comparison.py \
    --moose test/tests/actions/vs_F0AM_chamber_S1_box.csv \
    --gold test/tests/actions/gold/vs_F0AM_chamber_S1_box.csv \
    --save chamber_S1.pdf
```

The comparison scripts convert MOOSE's molec/cm³ output to ppb using
`M = 2.46e19 molecules/cm³` (standard conditions at 298 K, 1013 mbar).

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `scripts/kpp_to_fac.py` | Convert `.kpp` → `.fac` |
| `scripts/check_ro2.py` | Extract / validate RO2 species list against MCM reference |
| `scripts/generate_atchem2_gold.py` | Generate gold CSV from AtChem2 reference output |
| `scripts/plot_atchem2.py` | Plot MOOSE vs AtChem2 comparison (3×3 grid PDF) |
| `scripts/plot_vs_f0am.py` | Plot MOOSE vs analytical solution for F0AM tutorial |
| `scripts/gen_tutorial_gold.py` | Generate F0AM tutorial gold CSV from scipy ODE |

| `scripts/convert_soas_mat.py` | Convert SOAS campaign .mat observation files to CSV |
| `scripts/generate_gold_vs_F0AM_chamber.py` | Generate gold CSVs for chamber validation from F0AM .mat output |
| `scripts/compare_chamber_f0am.py` | Numerical comparison (per-species ratio) of MOOSE vs F0AM gold CSV |
| `scripts/plot_chamber_comparison.py` | Visual comparison (2x2 panel PDFs) of MOOSE vs F0AM gold CSV |
| `scripts/export_f0am_all_scenarios.m` | MATLAB script to run F0AM and export all 4 chamber scenarios |

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
| `vs_F0AM_chamber_S2b_box` | CSVDiff | 610 | S2b (restart, jcorr=10): BOTTOMUP vs F0AM |
| `vs_F0AM_dielcycle_box` | RunApp | 2908 | SOAS diel cycle: MCM_SZA |
| `vs_atchem2_transport_building` | Exodiff | 29 | NS + chemistry + building emission |
| `test_family_conservation` | RunApp | 29 | Family conservation (NOx = NO2 + NO) |
| `box_chamber` | — | 610 | F0AM chamber example (3h, MCM_SZA) |
| `box_dielcycle` | — | 2908 | F0AM diel cycle example (24h) |

## Objects

- [`AtmosphericChemistryAction`](source/actions/AtmosphericChemistryAction.md) — Unified Action (box / coupled modes)
- [`MCMBoxModel`](source/userobjects/MCMBoxModel.md) — 0-D chemical ODE engine with caching interface, pluggable sparse matrix backends (CSR/COO/DENSE/CSC), and optional F0AM-style limiting-reagent (LR) formulation for RO₂ termination
- [`ChemistryODEKernel`](source/kernels/ChemistryODEKernel.md) — Box mode ScalarKernel bridge to MCMBoxModel
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
