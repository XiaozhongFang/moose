# ChemistryODEKernel

!syntax description /ScalarKernels/ChemistryODEKernel

`ChemistryODEKernel` is a ScalarKernel (specifically an ODEKernel) that bridges
MOOSE's ScalarVariable system to the [`MCMBoxModel`](MCMBoxModel.md) computation
engine for 0-D atmospheric chemistry box models.

One instance per chemical species. During each residual/Jacobian evaluation, the
first kernel instance triggers a full-system dC/dt computation on the MCMBoxModel
(via its caching interface), and subsequent instances return cached results.

## Residual

$$R = -\frac{dC_s}{dt}$$

where $dC_s/dt$ is the chemical source term for species $s$ computed by MCMBoxModel.

## Jacobian

The analytical Jacobian is provided by MCMBoxModel's sparse Jacobian cache:
- Diagonal: $\partial R_s / \partial C_s$
- Off-diagonal: $\partial R_s / \partial C_j$

## Caching

The kernel calls `MCMBoxModel::markDirty()` in `reinit()`, which invalidates the
internal cache. The first call to `getDCdt()` or `getJacobian*()` triggers a fresh
full-system computation. Subsequent calls within the same evaluation return cached
values, avoiding redundant $O(N^2)$ work.

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
