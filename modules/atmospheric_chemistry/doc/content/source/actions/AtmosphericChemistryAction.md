# AtmosphericChemistryAction

!syntax description /AtmosphericChemistry/AtmosphericChemistryAction

`AtmosphericChemistryAction` is the unified entry point for atmospheric chemistry
simulations. It parses MCM (Master Chemical Mechanism) Facsimile-format (`.fac`)
mechanism files and sets up the simulation system according to the selected `mode`.

## Modes

### Box Mode (`mode = box`)

Creates a 0-D ODE box model suitable for large mechanisms (tested up to full MCM
v3.3.1, ~5832 species). Mechanism loading is orchestrated by the
[`MechanismLoader`](MechanismLoader.md) utility during Action construction.

1. **`add_variable`** — Creates `MooseVariableScalar` (family=SCALAR) for each species
2. **`add_user_object`** — Creates an [`MCMBoxModel`](MCMBoxModel.md) UserObject
   that parses the `.fac` mechanism and provides cached dC/dt and Jacobian access.
   A [`BoxIntegrator`](BoxIntegrator.md) strategy (`MooseImplicitIntegrator` or
   `PetscTSIntegrator`) is created alongside the UO to encapsulate the integration mode.
3. **`add_scalar_kernel`** — For each species, creates:
   - [`ODETimeDerivative`](ODETimeDerivative.md) — contributes $du/dt$ (skipped in PETSc TS mode)
   - [`ChemistryODEKernel`](ChemistryODEKernel.md) — contributes $-dC/dt$ (chemical source),
     delegates to the `BoxIntegrator` strategy for all residual/Jacobian evaluation,
     with no mode-specific branching in the kernel code

### Coupled Mode (`mode = coupled`)

Creates a FEM transport + chemistry system for spatially-resolved simulations
(5--50 species):

1. **`add_variable`** — Creates `MooseVariableFE` (family=LAGRANGE) for each species
2. **`add_material`** — Creates an [`MCMRatesMaterial`](MCMRatesMaterial.md) for
   runtime rate evaluation via fparser
3. **`add_kernel`** — For each species, creates:
   - [`TimeDerivative`](TimeDerivative.md) — time derivative term
   - [`ChemicalSourceKernel`](ChemicalSourceKernel.md) — chemical source with analytical Jacobian
   - Optional [`Diffusion`](Diffusion.md) (when `include_transport = true`)

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
  mode = coupled
  mechanism_file = 'mechanism.fac'
  photolysis_scheme = HYBRID
  hybrid_table_dir = 'tuv_tables'
  albedo = 0.1
  o3column = 350
  altitude = 500
[]
```

## Limiting Reagent (RO₂ Termination)

The [!param](/AtmosphericChemistry/use_limiting_reagent) parameter (default
`false`) enables the F0AM-style limiting-reagent formulation for RO₂+RO₂
termination reactions. See the [`MCMBoxModel` documentation](MCMBoxModel.md#limiting-reagent-ro₂-termination)
for the mathematical formulation and usage guidance.

When enabled, the parser detects RO₂+RO₂ reactions during mechanism loading
and applies:
$$ R_r = k_r \cdot \min([\mathrm{RO_2^{\it i}}], [\mathrm{RO_2^{\it j}}])^2 $$
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
  - `kpp_rosenbrock` — KPP Rosenbrock (requires `--enable-kpp`)
  - `kpp_sdirk` — KPP SDIRK (requires `--enable-kpp`)
  - `kpp_runge_kutta` — KPP Runge-Kutta (requires `--enable-kpp`)
- `chem_solver_type` (enum: bdf|arkimex|..., default: bdf): ODE solver type for `petsc_ts` only.
- `chem_solver_rtol` (real, default: 1e-6): Relative tolerance for the adaptive integrator.
- `chem_solver_atol` (real, default: 1e-10): Absolute tolerance for the adaptive integrator.

> **Note:** The old `box_solver`, `box_solver_type`, `box_solver_rtol`, and `box_solver_atol`
> parameters are **deprecated**. Use `chem_solver*` instead. Forward-compatible aliases
> remain available but emit a deprecation warning.

### Mechanism Format

The `mechanism_format` parameter selects the RHS/Jacobian source:

- `MCM_FACSIMILE` (default): Runtime `.fac` parsing via `MCMFacsimileParser`
- `KPP`: Pre-generated KPP Fortran code (requires `--enable-kpp`)

Not all `(mechanism_format, chem_solver)` combinations are valid:

| mechanism_format | chem_solver | Valid |
|-----------------|-------------|-------|
| `MCM_FACSIMILE` | `moose_implicit` | ✅ |
| `MCM_FACSIMILE` | `petsc_ts` | ✅ |
| `MCM_FACSIMILE` | `sundials` | ✅ |
| `MCM_FACSIMILE` | `kpp_*` | ❌ |
| `KPP` | `moose_implicit` | ❌ |
| `KPP` | `petsc_ts` | ❌ |
| `KPP` | `sundials` | ❌ |
| `KPP` | `kpp_*` | ✅ |

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
  mode = box
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
  mode = coupled
  mechanism_file = 'mechanism.fac'
  temperature = 298
  include_transport = true
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
