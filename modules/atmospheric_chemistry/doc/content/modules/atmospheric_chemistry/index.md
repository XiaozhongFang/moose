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

## Chemical Mechanism Files

Mechanism files in FACSIMILE (`.fac`) format are supported directly. For MCM website-exported `.kpp`
files, use the included converter:

```bash
python3 modules/atmospheric_chemistry/scripts/kpp_to_fac.py mechanism.kpp
```

## Utility Scripts

| Script | Purpose |
|--------|---------|
| `scripts/kpp_to_fac.py` | Convert `.kpp` → `.fac` |
| `scripts/check_ro2.py` | Compare RO2 list against AtChem2 reference |
| `scripts/generate_atchem2_gold.py` | Generate gold CSV from AtChem2 reference output |
| `scripts/plot_atchem2.py` | Plot MOOSE vs AtChem2 comparison (3×3 grid PDF) |
| `scripts/plot_vs_f0am.py` | Plot MOOSE vs analytical solution for F0AM tutorial |
| `scripts/gen_tutorial_gold.py` | Generate F0AM tutorial gold CSV from scipy ODE |

## Database

Pre-converted mechanism files are available in `doc/content/modules/atmospheric_chemistry/database/`:

| File | Species | Reactions |
|------|---------|-----------|
| `tutorial_5sp.fac` | 5 | 6 |
| `atchem2_example.fac` | 29 | 71 |
| `mcm_export.fac` | 20 | 48 |
| `mcm_export_all.fac` | 5832 | 17224 |

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
- [`MCMSolarPostprocessor`](source/postprocessors/MCMSolarPostprocessor.md) — Solar parameters (cosx, secx, lha, etc.)
- [`MCMPhotolysisPostprocessor`](source/postprocessors/MCMPhotolysisPostprocessor.md) — Individual photolysis J-values
- [`MCMRO2Postprocessor`](source/postprocessors/MCMRO2Postprocessor.md) — Total peroxy radical (RO2) concentration

!syntax complete groups=AtmosphericChemistryApp level=3
