# AtmosphericChemistryAction

!syntax description /AtmosphericChemistry/AtmosphericChemistryAction

`AtmosphericChemistryAction` is the compatibility wrapper for the new
`[AtmosphericChemistryBox]` and `[AtmosphericChemistryCoupled]` Actions.
New input files should use the dedicated blocks instead of `mode = box|coupled`.

The format is auto-detected: `chem_solver = kpp_*` or `.kpp` extension -> KPP format;
otherwise MCM_FACSIMILE is assumed.

## [AtmosphericChemistryBox]

Dedicated block for 0-D ODE box model simulations. Creates scalar variables,
an `MCMBoxModel` UserObject, and `ChemistryODEKernel` for each species.
Suitable for large mechanisms (tested up to full MCM v3.3.1, ~5832 species).

## [AtmosphericChemistryCoupled]

Dedicated block for FEM transport + chemistry simulations. Creates FE variables,
an `MCMRatesMaterial`, and `ChemicalSourceKernel` for each species.
Suitable for spatially-resolved simulations (5--50 species).

## Mechanism Formats

### `MCM_FACSIMILE` (default)

Standard MCM `.fac` format, parsed by [`MechanismLoader`](MechanismLoader.md).
Used with `chem_solver = moose_implicit` or `petsc_ts`.

### `KPP` (Kinetic Pre-Processor)

Native KPP `.kpp` format. In the box KPP workflow, the generated KPP shared
library provides the authoritative `#DEFVAR` species order through `SPC_NAMES`.
Use `chem_solver = kpp_rosenbrock` for the supported KPP box integrator.

The chemistry is integrated by a pre-compiled KPP shared library. By default,
the library is auto-discovered beside the `.kpp` file using
`kpp_build_<mechanism>/libkpp_<mechanism>.so`; `KPP_LIB` can override this path.
Build the `.so` with:

```bash
make -f modules/atmospheric_chemistry/kpp/build/Makefile MECH=path/to/mech.kpp
```

Example:
```moose
[AtmosphericChemistry]
  [Box]
    mechanism_file = 'mechanism.kpp'
    chem_solver = kpp_rosenbrock
  []
[]
```

See [`KPP Workflow for Atmospheric Chemistry Box Models`](/modules/atmospheric_chemistry/kpp_workflow.md)
for different mechanism files, user-defined mechanisms, and validation guidance.

Format is auto-detected:
- `chem_solver = kpp_*` -> KPP format, or
- File extension `.kpp` -> KPP format
- Otherwise -> `MCM_FACSIMILE` (`.fac`)

## Photolysis Schemes

Two photolysis schemes are available via the `photolysis_scheme` parameter:

### `MCM_SZA` (default)

MCM empirical solar-zenith-angle parameterization, identical to AtChem2:

$$J_N = \mathrm{CL}_N \cdot \cos^{\mathrm{CMM}_N}(\mathrm{SZA}) \cdot \exp(-\mathrm{CNN}_N \cdot \sec(\mathrm{SZA})) \cdot \mathrm{JFAC}$$

Requires `mcm_photolysis_file` pointing to a CL/CMM/CNN parameter file
(`mcm_photolysis_rates_v3.3.1.dat`).

### `HYBRID`

F0AM TUVv5.2 4D lookup table interpolation:

$$J_N = 10^{\,\log_{10}J_N(\mathrm{SZA}, \alpha, \Omega, z)} \cdot \mathrm{JFAC}$$

Where $\mathrm{SZA}$ is solar zenith angle, $\alpha$ is surface albedo,
$\Omega$ is O$_3$ column (Dobson Units), and $z$ is altitude (meters).

Requires `hybrid_table_dir` pointing to F0AM Hybrid table files. Uses
[`HybridJTableReader`](/utils/HybridJTableReader.md) for table loading
and 4D linear interpolation.

### Usage

```moose
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'mechanism.fac'
    photolysis_scheme = HYBRID
    hybrid_table_dir = 'tuv_tables'
    albedo = 0.1
    o3column = 350
    altitude = 500
  []
[]
```

## Limiting Reagent (RO₂ Termination)

The [!param](/AtmosphericChemistry/use_limiting_reagent) parameter (default
`false`) enables the F0AM-style limiting-reagent formulation for RO₂+RO₂
termination reactions. See the [`MCMBoxModel` documentation](MCMBoxModel.md#limiting-reagent-ro₂-termination)
for the mathematical formulation and usage guidance.

When enabled, the parser detects RO₂+RO₂ reactions during mechanism loading
and applies:
$$ R_r = k_r \cdot \min([\mathrm{RO_2^{\it i}}], [\mathrm{RO_2^{\it j}}]) $$
instead of the standard $k_r \cdot [\mathrm{RO_2^{\it i}}] \cdot [\mathrm{RO_2^{\it j}}]$.

This parameter is forwarded to the underlying `MCMBoxModel` UserObject in
box mode. In coupled mode, it is not applicable (the `MCMRatesMaterial`
always uses the standard formulation).

## Family Conservation (DAE Method)

The [!param](/AtmosphericChemistry/family_names), 
[!param](/AtmosphericChemistry/family_members), and
[!param](/AtmosphericChemistry/family_scaling) parameters enable F0AM-style
chemical family conservation via DAE slack variables. See the
[`MCMFamilyConstraint` documentation](MCMFamilyConstraint.md) for details.

### Chemical Solver (box mode only)

The `chem_solver` parameter selects the embedded chemical ODE solver inside
`MCMBoxModel`. This is separate from the outer `[Executioner]` time integration:

| What | Controls |
|------|----------|
| `[Executioner]` | Outer time step advancement + linear solver (`-pc_type lu`, etc.) |
| `chem_solver` | Chemical ODE integration within each time step |

Parameters:

- `chem_solver` (enum, default: `moose_implicit`):
  - `moose_implicit` — MOOSE Newton solver owns the integration (default)
  - `petsc_ts` — PETSc TS (BDF/ARKIMEX), 28-32× speedup for box mode
  - `sundials` — SUNDIALS CVODE (requires `--with-sundials`)
  - `kpp_rosenbrock` — KPP Rosenbrock through a generated shared library
  - `kpp_sdirk` — reserved KPP option; not wired to a supported adapter entry point
  - `kpp_runge_kutta` — reserved KPP option; not wired to a supported adapter entry point
- `chem_solver_type` (enum: bdf|arkimex|..., default: bdf): ODE solver type for `petsc_ts` only.
- `chem_solver_rtol` (real, default: 1e-6): Relative tolerance for the adaptive integrator.
- `chem_solver_atol` (real, default: 1e-10): Absolute tolerance for the adaptive integrator.

> **Note:** The old `box_solver`, `box_solver_type`, `box_solver_rtol`, and `box_solver_atol`
> parameters are **deprecated**. Use `chem_solver*` instead. Forward-compatible aliases
> remain available but emit a deprecation warning.

## Example Input File Syntax

### Box Mode

```moose
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = 'mechanism.fac'
    temperature = 298.15
    air_density = 2.46e19
    water_vapor = 2.46e17
    press = 1013.25
    mcm_photolysis_file = 'mcm_photolysis_rates_v3.3.1.dat'
    latitude = 51.51
    longitude = 0.13
    day = 21
    month = 6
    year = 2010
  []
[]

[ICs]
  [o3_ic]
    type = ScalarConstantIC
    variable = O3
    value = 5.2e11
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 900
  end_time = 43200
[]

[Outputs]
  csv = true
[]
```

### Coupled Mode

```moose
[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'mechanism.fac'
    temperature = 298
    include_transport = true
  []
[]
```

## Supported `.fac` Format

The underlying [`MCMFacsimileParser`](MCMFacsimileParser.md) handles:
- Fortran `D`/`d` notation ($1.0\mathrm{D}{-}31 \to 1.0\times 10^{-31}$)
- `@` and `**` exponentiation operators
- `EXP()`, `LOG10()` function names
- Multi-line statements terminated by `;`
- Automatic J\<N\> photolysis detection from coefficient expressions

!syntax parameters /AtmosphericChemistry/AtmosphericChemistryAction

!syntax inputs /AtmosphericChemistry/AtmosphericChemistryAction

!syntax children /AtmosphericChemistry/AtmosphericChemistryAction
