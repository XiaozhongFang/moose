# BottomUpJIntegrator

!syntax description /UserObjects/BottomUpJIntegrator

## Overview

`BottomUpJIntegrator` implements the F0AM "BottomUp" photolysis method
(Jmethod=1), which computes photolysis frequencies by directly integrating
the product of absorption cross-section, quantum yield, and lamp/actinic
flux over wavelength:

\[
J = \int QY(\lambda) \cdot CS(\lambda) \cdot F(\lambda) \, d\lambda
\]

This is the standard method for laboratory chamber experiments where
a known lamp spectrum replaces solar radiation.  The implementation
closely follows the F0AM `IntegrateJ.m` algorithm, including the
integral-based convolution ("smearing") from the TUV model `numer.f/interp2`.

## Supported Data Formats

| Type | Description | Example |
|------|-------------|---------|
| 0 | Scalar constant (quantum yield only) | `QY = 0.76` |
| 1 | 2-column CSV: `wl[nm], value` | `Cross_Section_HONO.csv` |
| 2 | 3-column CSV: `wl, val@T1, val@T2` (temperature-dependent, linear interp) | `Cross_Section_O3_JPL.csv` |
| 3 | TXT: space/tab-separated `wl value` | JPL halogen cross-sections |

Temperature-dependent `.m` function files from F0AM are pre-computed
at the target temperature using the `generate_bottomup_jmap.py` script
and saved as 2-column CSV.

## Data Files

The integrator expects the following directory layout under `bottomup_data_dir`:

```
bottomup_data_dir/
├── ExampleLightFlux.txt          # Lamp flux file (2-col: wl, flux)
├── bottomup_jmap.dat             # Reaction mapping file
├── CrossSections/                # Cross-section CSV/TXT files
│   ├── Cross_Section_O3_JPL_precomp.csv
│   ├── Cross_Section_NO2_precomp.csv
│   └── ...
└── QuantumYields/                # Quantum yield CSV/TXT files
    ├── Quantum_Yield_O3_O1D_JPL_precomp.csv
    └── ...
```

## Reaction Mapping File

`bottomup_jmap.dat` maps each J-value to its cross-section and quantum yield sources:

```
# Format: JNAME  CS_FILE  CS_TYPE  QY_FILE  QY_TYPE
J1   Cross_Section_O3_JPL_precomp.csv   1   Quantum_Yield_O3_O1D_JPL_precomp.csv   1
J4   Cross_Section_NO2_precomp.csv      1   Quantum_Yield_NO2_precomp.csv          1
J7   Cross_Section_HONO.csv             1   1                                       0
```

- `CS_TYPE`: 1=2-col CSV, 2=3-col CSV, 3=TXT
- `QY_TYPE`: 0=scalar, 1=2-col CSV, 2=3-col CSV, 3=TXT

## Generating Data Files

Use the `generate_bottomup_jmap.py` script to pre-compute all CS/QY data from
an F0AM installation:

```bash
python3 scripts/generate_bottomup_jmap.py \
    --f0am-photolysis /path/to/F0AM/Chem/Photolysis \
    --output-dir doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup \
    --temperature 298.0
```

## Input File Usage

```moose
[AtmosphericChemistry]
  mode = box
  mechanism_file = 'mechanism.fac'
  temperature = 298.0
  press = 1013.0
  photolysis_scheme = BOTTOMUP
  lamp_flux_file = 'ExampleLightFlux.txt'
  bottomup_data_dir = 'database/photolysis/bottomup'
  jfac = 1.0
[]
```

## Algorithm Details

For each registered J-value:

1. Load cross-section data as (wavelength, value) at the current T
2. Load quantum yield data similarly, or use scalar constant
3. Convolve ("smear") CS and QY onto the lamp flux wavelength grid using
   the TUV `interp2` trapezoidal integral convolution
4. Compute J = trapezoidal integration of QY_out × CS_out × Flux over wavelength

The smearing algorithm ensures correct handling of narrow spectral features
even when CS/QY data have finer wavelength resolution than the flux grid.
This replicates the behavior of the F0AM `IntegrateJ.m` subroutine exactly.

!syntax parameters /UserObjects/BottomUpJIntegrator

!syntax inputs /UserObjects/BottomUpJIntegrator

!syntax children /UserObjects/BottomUpJIntegrator
