# AtmosphericChemistryAction

!syntax description /AtmosphericChemistry/AtmosphericChemistryAction

`AtmosphericChemistryAction` is the unified entry point for atmospheric chemistry
simulations. It parses MCM (Master Chemical Mechanism) Facsimile-format (`.fac`)
mechanism files and sets up the simulation system according to the selected `mode`.

## Modes

### Box Mode (`mode = box`)

Creates a 0-D ODE box model suitable for large mechanisms (tested up to full MCM
v3.3.1, ~5832 species):

1. **`add_variable`** — Creates `MooseVariableScalar` (family=SCALAR) for each species
2. **`add_user_object`** — Creates an [`MCMBoxModel`](MCMBoxModel.md) UserObject
   that parses the `.fac` mechanism and provides cached dC/dt and Jacobian access
3. **`add_scalar_kernel`** — For each species, creates:
   - [`ODETimeDerivative`](ODETimeDerivative.md) — contributes $du/dt$
   - [`ChemistryODEKernel`](ChemistryODEKernel.md) — contributes $-dC/dt$ (chemical source)

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

### PETSc TS Standalone Integrator

New parameters available since implementation:

- `petsc_ts` (bool, default: false): Enable PETSc TS standalone integrator for box mode.
- `petsc_ts_type` (enum: bdf|arkimex|sundials, default: bdf): TS integrator type.
- `petsc_ts_rtol` (real, default: 1e-6): Relative tolerance.
- `petsc_ts_atol` (real, default: 1e-10): Absolute tolerance.

When `petsc_ts = true` and `mode = box`, the MCMBoxModel handles integration via
PETSc TS, achieving 28-32× speedup over the MOOSE LU solver for the S1 chamber case.
Coupled mode with `petsc_ts = true` triggers an error at input parsing time.

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
