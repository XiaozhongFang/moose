# MCMFacsimileAction

!syntax description /Actions/MCMFacsimileAction

`MCMFacsimileAction` parses an MCM (Master Chemical Mechanism) Facsimile-format
(`.fac`) mechanism file and sets up a complete zero-dimensional atmospheric
chemistry box model ODE system in a single input file block.

The Action is registered for three MOOSE tasks and performs the following:

### `add_variable` task
Creates a nonlinear variable for each chemical species declared in the `.fac`
file's `VARIABLE` block. All variables use `LAGRANGE` first-order finite elements
on the (single-element) mesh.

### `add_material` task
Creates an [`MCMRatesMaterial`](MCMRatesMaterial.md) that evaluates all
rate coefficient expressions at runtime using `fparser`. The material:
- Evaluates intermediate coefficients in topological (Kahn's algorithm) order
- Computes photolysis rates via the MCM SZA parameterization:
  $J = C_L \cdot \cos^C_{MM}(\theta) \cdot \exp(-C_{NN}/\cos\theta) \cdot F_{JFAC}$
- Pre-computes per-reaction rate values from reactant concentrations

### `add_kernel` task
For each species, creates a [`TimeDerivative`](TimeDerivative.md) and a
[`ChemicalSourceKernel`](ChemicalSourceKernel.md). The source kernel computes
the ODE right-hand side from the stoichiometric matrix and reaction rates.

### Supported `.fac` Format

The parser handles:
- Fortran `D`/`d` notation ($1.0\mathrm{D}{-}31 \to 1.0\times 10^{-31}$)
- `@` and `**` exponentiation operators
- `EXP()`, `LOG10()` function names
- Multi-line statements terminated by `;`
- Automatic J<N> photolysis detection from coefficient expressions

## Example Input File Syntax

```moose
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[MCMFacsimileAction]
  mechanism_file = 'mcm_subset.fac'
  temperature = 298.15
  air_density = 2.46e19
  water_vapor = 2.46e17
  mcm_photolysis_file = 'mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [o3_ic]
    type = ConstantIC
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

!syntax parameters /Actions/MCMFacsimileAction

!syntax inputs /Actions/MCMFacsimileAction

!syntax children /Actions/MCMFacsimileAction
