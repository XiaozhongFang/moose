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

To generate the BottomUp data files from an F0AM installation:

```bash
python3 scripts/generate_bottomup_jmap.py \
    --f0am-photolysis /path/to/F0AM/Chem/Photolysis \
    --output-dir database/photolysis/bottomup \
    --temperature 298.0
```

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

The box mode has been validated against AtChem2 (v1.2.1) for the MCM v3.3.1 inorganic
subset (29 species, 71 reactions, 12 h diurnal cycle).  O$_3$, NO$_2$, and NO
concentrations match to within 0.1--0.3 % relative error during daytime and
< 0.05 % at night.  Solar parameters (declination, zenith angle) follow the
Madronich (1993) formulation identically to AtChem2.

To regenerate validation data:

```bash
python3 scripts/generate_atchem2_gold.py \
    --atchem2 <path-to-AtChem2>/model/output \
    --output test/tests/actions/gold/vs_AtChem2_inorg_box.csv

python3 scripts/plot_atchem2.py \
    --moose test/tests/actions/vs_AtChem2_inorg_box.csv \
    --atchem2 <path-to-AtChem2>/model/output
```

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `scripts/kpp_to_fac.py` | Convert `.kpp` → `.fac` |
| `scripts/check_ro2.py` | Extract / validate RO2 species list against MCM reference |
| `scripts/generate_atchem2_gold.py` | Generate gold CSV from AtChem2 reference output |
| `scripts/plot_atchem2.py` | Plot MOOSE vs AtChem2 comparison (3×3 grid PDF) |
| `scripts/plot_vs_f0am.py` | Plot MOOSE vs analytical solution for F0AM tutorial |
| `scripts/gen_tutorial_gold.py` | Generate F0AM tutorial gold CSV from scipy ODE |
| `scripts/generate_bottomup_jmap.py` | Pre-compute BottomUp photolysis CS/QY data from F0AM installation |
| `scripts/convert_soas_mat.py` | Convert SOAS campaign .mat observation files to CSV |
| `scripts/extract_mcm_k.py` | Extract MCM standard rate constants from MCMv331_K.m and inject into .fac files |

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
| `vs_F0AM_chamber_box` | RunApp | 610 | Chamber: BOTTOMUP photolysis |
| `vs_F0AM_dielcycle_box` | RunApp | 2908 | SOAS diel cycle: MCM_SZA |
| `vs_atchem2_transport_building` | Exodiff | 29 | NS + chemistry + building emission |

## Objects

- [`AtmosphericChemistryAction`](source/actions/AtmosphericChemistryAction.md) — Unified Action (box / coupled modes)
- [`MCMBoxModel`](source/userobjects/MCMBoxModel.md) — 0-D chemical ODE engine with caching interface and pluggable sparse matrix backends (CSR/COO/DENSE/CSC)
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
- [`MCMRO2Postprocessor`](source/postprocessors/MCMRO2Postprocessor.md) — Total peroxy radical (RO2) concentration

!syntax complete groups=AtmosphericChemistryApp level=3
