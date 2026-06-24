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
