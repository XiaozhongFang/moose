# KPP Workflow for Atmospheric Chemistry Box Models

## Overview

The Kinetic Pre-Processor (KPP) generates optimized Fortran/C code for chemical
ODE systems from a high-level `.kpp` mechanism description.  This module uses
KPP-generated shared libraries (`.so`) as the runtime backend for box-mode
chemistry integration.

## Prerequisites

- **KPP tool** installed and available on `PATH`, or `KPP_HOME` environment
  variable set.
- **KPP_ENABLED=1** compile flag (enabled automatically when the configure
  script detects KPP).

## Workflow

```
.kpp mechanism file
        │
        ▼
  kpp/build/Makefile ──► kpp_build_<mech>/
        │                    ├── libkpp_<mech>.so
        │                    └── kpp_<mech>.json  (species metadata)
        ▼
  KppBoxIntegrator (dlopen .so, call kpp_integrate)
        │
        ▼
  KPPGeneratedMechanism (IMechanism backend for RHS/Jacobian)
```

### 1. Build the KPP shared library

```bash
# From the atmospheric_chemistry module root:
make -f kpp/build/Makefile MECH=path/to/mechanism.kpp
```

This produces:
- `kpp_build_<mech>/libkpp_<mech>.so` — shared library for runtime loading
- `kpp_build_<mech>/kpp_<mech>.json` — species metadata (species order,
  NSPEC, NVAR, NREACT)

### 2. Configure the input file

```bash
[AtmosphericChemistryBox]
  mechanism_file = 'path/to/mechanism.kpp'
  chem_solver = kpp_rosenbrock     # or kpp_sdirk, kpp_runge_kutta
  # ... other chemistry parameters ...
[]
```

The mechanism format is auto-detected from:
- `chem_solver = kpp_*` → KPP format
- File extension `.kpp` → KPP format

### 3. Set KPP_LIB (optional)

By default the runtime auto-discovers the `.so` path from the mechanism file
location:
- `<mech_dir>/kpp_build_<mech>/libkpp_<mech>.so`

Override with the `KPP_LIB` environment variable:
```bash
KPP_LIB=/custom/path/libkpp_mech.so ./atmospheric_chemistry-opt -i input.i
```

## Metadata JSON Format

The build generates `<mech>.json` alongside the `.so`:

```json
{
  "mechanism": "small_strato",
  "n_species": 14,
  "n_reactions": 32,
  "n_var": 14,
  "species": ["O3", "NO2", "NO", "OH", "CO", "CH4", ...],
  "jacobian_nnz": 48
}
```

This metadata enables runtime validation (species order mismatch) and is used
by `KPPGeneratedMechanism` to set up the IMechanism interface.

## Environment Variables

| Variable | Purpose | Example |
|---|---|---|
| `KPP_HOME` | KPP installation root (for `#MODEL` includes in `.kpp` files) | `/opt/KPP` |
| `KPP_LIB` | Override path to the KPP `.so` | `/path/to/libkpp_mech.so` |
| `KPP_ENABLED` | Compile-time flag (set by build system) | `1` |

## Test Files

- `test/tests/actions/kpp_small_strato.i` — KPP box validation test
- `test/tests/actions/kpp_small_strato/small_strato.kpp` — KPP mechanism

## Updating Gold Files

After modifying a KPP mechanism reference output:

```bash
cd modules/atmospheric_chemistry
make -f kpp/build/Makefile MECH=test/tests/actions/kpp_small_strato/small_strato.kpp
KPP_LIB=test/tests/actions/kpp_small_strato/kpp_build_small_strato/libkpp_small_strato.so \
  ../../atmospheric_chemistry-opt -i test/tests/actions/kpp_small_strato.i
cp test/tests/actions/kpp_small_strato_out.csv \
   test/tests/actions/gold/kpp_small_strato.csv
```

## Error Messages

### "no species found"
Ensure the `.kpp` file has a `#DEFVAR` section with species definitions.
Check that `#INCLUDE` and `#MODEL` directives resolve correctly.

### "cannot open KPP include file"
The parser searches:
1. Relative to the `.kpp` file directory
2. `$KPP_HOME/models/` (for KPP built-in models)
3. Current working directory

Set `KPP_HOME` if your mechanism uses `#MODEL` directives from the KPP
standard library.

### "KPP integration failed"
Check tolerances (`chem_solver_rtol`, `chem_solver_atol`).  Very stiff
systems may require tighter tolerances or a different solver type
(e.g., `kpp_rosenbrock` instead of `kpp_sdirk`).
