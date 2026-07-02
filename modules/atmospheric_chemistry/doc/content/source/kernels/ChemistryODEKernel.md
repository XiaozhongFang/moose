# ChemistryODEKernel

!syntax description /ScalarKernels/ChemistryODEKernel

`ChemistryODEKernel` is a ScalarKernel (specifically an ODEKernel) that bridges
MOOSE's ScalarVariable system to the [`BoxIntegrator`](BoxIntegrator.md) strategy
interface for 0-D atmospheric chemistry box models.

One instance per chemical species. All residual and Jacobian evaluation is
delegated to the [`BoxIntegrator`](BoxIntegrator.md), which encapsulates the
integration strategy (MOOSE implicit or PETSc TS). The kernel has zero
mode-specific branching.

In MOOSE implicit mode, the first kernel instance triggers a full-system dC/dt
computation on the [`MCMBoxModel`](MCMBoxModel.md) via its caching interface,
and subsequent instances return cached results. In PETSc TS mode, the integrator
returns zero for all evaluations and the full-system solve is handled by
`MCMBoxModel::execute()`.

## Residual

$$R = -\frac{dC_s}{dt}$$

where $dC_s/dt$ is the chemical source term for species $s$ provided by the
`BoxIntegrator`.

## Jacobian

The analytical Jacobian is delegated to the `BoxIntegrator`:
- Diagonal: $\partial R_s / \partial C_s$
- Off-diagonal: $\partial R_s / \partial C_j$

## Caching (MOOSE implicit mode)

The kernel calls `BoxIntegrator::reinit()` at the start of each timestep. In
MOOSE implicit mode, this invalidates the `MCMBoxModel` cache via
`markDirty()`/`setCurrentTime()`. The first call to `computeResidual()` or
`computeJacobian*()` triggers a fresh full-system computation. Subsequent
calls within the same evaluation return cached values, avoiding redundant
$O(N^2)$ work.

In PETSc TS mode, `reinit()` is a no-op — the integrator handles all state.

## Example Input File Syntax

`ChemistryODEKernel` is created automatically by the
[`AtmosphericChemistryAction`](AtmosphericChemistryAction.md) when `mode = box`.
It can also be created manually:

```moose
[ScalarKernels]
  [chem_O3]
    type = ChemistryODEKernel
    variable = O3
    box_model = box_model
    species_index = 0
    species_variables = 'O3 NO NO2 ...'
  []
[]
```

!syntax parameters /ScalarKernels/ChemistryODEKernel

!syntax inputs /ScalarKernels/ChemistryODEKernel

!syntax children /ScalarKernels/ChemistryODEKernel
