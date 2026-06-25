# MCMRatesMaterial

!syntax description /Materials/MCMRatesMaterial

## Overview

`MCMRatesMaterial` evaluates MCM (Master Chemical Mechanism) reaction rates at
each quadrature point using sequential fparser evaluation. Rate coefficients are
computed in topological order via a shared parameter array that holds background
atmospheric variables (TEMP, M, O2, N2, H2O), all coefficient results, species
concentrations, and photolysis J-values.

Reaction rates $R_i = k_i \cdot \prod [C_j]^{\nu_{ij}}$ are then computed and
stored in the `reaction_rates` material property for consumption by
[ChemicalSourceKernel](/kernels/ChemicalSourceKernel.md).

## Photolysis Schemes

Two photolysis schemes are supported, selected via the `photolysis_scheme` parameter:

### 1. `MCM_SZA` (default)

Uses the MCM empirical solar-zenith-angle parameterization:

$$J_N = \mathrm{CL}_N \cdot \cos^{\mathrm{CMM}_N}(\mathrm{SZA}) \cdot \exp\left(-\mathrm{CNN}_N \cdot \sec(\mathrm{SZA})\right) \cdot \mathrm{JFAC}$$

CL/CMM/CNN coefficients are read from the MCM photolysis-rates parameter file
(e.g., `mcm_photolysis_rates_v3.3.1.dat`). This scheme is identical to AtChem2's
`zenith_data_mod` / `solar_functions_mod`.

Solar zenith angle is computed via the Madronich (1993) formula, using latitude,
longitude, day-of-year, and local time.

### 2. `HYBRID`

Uses pre-computed F0AM Hybrid J-value lookup tables (TUVv5.2 4D grid):

$$J_N = 10^{\,\log_{10}J_N(\mathrm{SZA}, \alpha, \Omega, z)} \cdot \mathrm{JFAC}$$

Where:
- $\mathrm{SZA}$ = solar zenith angle (degrees)
- $\alpha$ = surface albedo (0--1)
- $\Omega$ = O$_3$ column (Dobson Units)
- $z$ = altitude (meters)

The lookup table grid covers SZA=0:5:80,82:2:100, albedo=0:0.2:1,
O$_3$=100:50:600 DU, altitude=0:1:22 km. 4D linear interpolation is performed
in log$_{10}$(J) space via [HybridJTableReader](/utils/HybridJTableReader.md).

Requires `hybrid_table_dir` pointing to a directory containing F0AM table files
(`table_J<N>.dat`, `axis_*.dat`, `index.txt`).

!syntax parameters /Materials/MCMRatesMaterial

!syntax inputs /Materials/MCMRatesMaterial

!syntax children /Materials/MCMRatesMaterial
