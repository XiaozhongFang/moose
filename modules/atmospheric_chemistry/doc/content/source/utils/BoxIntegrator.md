# BoxIntegrator

!syntax description /Utils/BoxIntegrator (abstract interface, not a direct MOOSE object)

`BoxIntegrator` is an abstract strategy interface that decouples the chemical
source-term evaluation from the integration strategy in box-mode simulations.
It has two concrete implementations:

- **`MooseImplicitIntegrator`** — Wraps `MCMBoxModel` for per-species residual/Jacobian
  evaluation in MOOSE-driven implicit mode. `computeResidual()` / `computeJacobian*()`
  forward to `MCMBoxModel::getDCdt()` and `getJacobian*()`. `solve()` is a no-op.

- **`PetscTSIntegrator`** — Returns zero for all residual/Jacobian evaluations
  (self-driven mode). The full-system integration is handled separately by
  `MCMBoxModel::execute()` running PETSc TS.

## Architecture

```
ChemistryODEKernel
  → BoxIntegrator::computeResidual(idx, C)   // no mode-specific branching
  → BoxIntegrator::computeJacobian*(idx, C)

MCMBoxModel::execute()
  → BoxIntegrator::selfDriven()
    → true:  run PETSc TS integration (PetscTSIntegrator)
    → false: no-op (MooseImplicitIntegrator — MOOSE solver handles it)
```

## Motivation

Previously, `ChemistryODEKernel` contained explicit `if (usePETScTS()) return 0`
branches in every residual/Jacobian method, creating an implicit fork between
MOOSE-driven and PETSc TS modes. The `BoxIntegrator` interface eliminates this:

- **Kernels have zero mode-specific branching** — they always delegate to the
  integrator
- **New integrator strategies** (CVODE, ARKode, custom) require only a new
  implementation of `BoxIntegrator`
- **Integration strategy is encapsulated** — `MCMBoxModel` no longer needs to
  expose `usePETScTS()` publicly

## API

| Method | Description |
|--------|-------------|
| `computeResidual(idx, C)` | dC[species_idx]/dt (molec/cm³/s) |
| `computeJacobianDiagonal(idx, C)` | ∂(dC[idx]/dt)/∂C[idx] |
| `computeJacobianOffDiagonal(i, j, C)` | ∂(dC[i]/dt)/∂C[j] |
| `reinit(time)` | Cache invalidation (MOOSE mode) or no-op (TS mode) |
| `selfDriven()` | Whether integrator handles the solve itself |
| `ppbToMolec()` | Unit conversion factor |
