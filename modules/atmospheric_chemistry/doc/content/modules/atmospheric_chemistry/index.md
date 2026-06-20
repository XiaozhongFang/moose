# Atmospheric Chemistry Module

The atmospheric chemistry module provides tools for 0-D box model and spatially-resolved
atmospheric chemistry simulations using Master Chemical Mechanism (MCM) chemical mechanisms.

## Overview

The module is built around two complementary interfaces:

1. **MCMFacsimileAction** — A high-level MOOSE Action that parses `.fac` mechanism files and
   automatically creates the complete ODE system (variables, kernels, materials).
2. **MCMBoxModel** — A standalone UserObject for direct ODE computation, providing an F0AM-compatible
   interface for dC/dt, analytical Jacobian, solar cycles, and dilution.

## Quick Start

```moose
[MCMFacsimileAction]
  mechanism_file = 'mechanism.fac'
  temperature = 298
  include_transport = true   # enable diffusion + reaction coupling
[]
```

## Chemical Mechanism Files

Mechanism files in FACSIMILE (`.fac`) format are supported directly. For MCM website-exported `.kpp`
files, use the included converter:

```bash
python3 modules/atmospheric_chemistry/scripts/kpp_to_fac.py mechanism.kpp
```

## Utility Scripts

| Script | Path | Purpose |
|--------|------|---------|
| KPP → FACSIMILE | `scripts/kpp_to_fac.py` | Convert `.kpp` → `.fac` |
| RO2 Validation | `scripts/check_ro2.py` | Compare RO2 list against AtChem2 reference |

## Database

Pre-converted mechanism files are available in `doc/content/modules/atmospheric_chemistry/database/`:

| File | Species | Reactions |
|------|---------|-----------|
| `tutorial_5sp.fac` | 5 | 6 |
| `atchem2_example.fac` | 29 | 71 |
| `mcm_export.fac` | 20 | 48 |
| `mcm_export_all.fac` | 5832 | 17224 |

## Objects

- [`MCMFacsimileAction`](source/actions/MCMFacsimileAction.md) — Parses `.fac` and builds MOOSE system
- [`MCMBoxModel`](source/userobjects/MCMBoxModel.md) — 0-D chemical ODE engine
- [`MCMRatesMaterial`](source/materials/MCMRatesMaterial.md) — Runtime rate evaluation
- [`ChemicalSourceKernel`](source/kernels/ChemicalSourceKernel.md) — ODE source term with analytical Jacobian
- [`MCMConstraintKernel`](source/kernels/MCMConstraintKernel.md) — Fixed-species constraint
- [`MCMEmissionKernel`](source/kernels/MCMEmissionKernel.md) — Emission source
- [`MCMDepositionKernel`](source/kernels/MCMDepositionKernel.md) — Dry deposition
- [`MCMFacsimileParser`](source/userobjects/MCMFacsimileParser.md) — Standalone `.fac` parser
- [`HybridJTableReader`](source/utils/HybridJTableReader.md) — TUV photolysis 4D interpolation

!syntax complete groups=AtmosphericChemistryApp level=3
