# Implementation Guide

This document provides practical implementation details for the CutFEM module.

## Phase 1: Ghost Penalty Kernel Implementation

### Header Structure (CutFEMApp.h)

```cpp
#pragma once

#include "MooseApp.h"

class CutFEMApp : public MooseApp
{
public:
  static InputParameters validParams();
  CutFEMApp(const InputParameters & params);
  virtual ~CutFEMApp();
  
  static void printInfo();
};
```

### Kernel Implementation (GhostPenaltyKernel)

The kernel computes the face integral:
$$\int_F c_F h_F^{2(k-1+\gamma)} [D_n^k u] \cdot [D_n^k v] \, dS$$

**Key methods**:

1. **computeQpResidual()**
   - Calculate penalty coefficient: $c_F h_F^{2(k-1+\gamma)}$
   - Extract normal derivatives from `_grad_u[_qp]` and `_grad_u_neighbor[_qp]`
   - Compute jump: `[D_n u] = (grad_u · normal) - (grad_u_neighbor · normal)`
   - Return: penalty × jump × test function

2. **computeQpJacobian()**
   - Implement all 4 Jacobian types: ElementElement, ElementNeighbor, NeighborElement, NeighborNeighbor
   - Each type contributes to the system matrix assembly

### Parameters

| Parameter | Type | Range | Default | Meaning |
|-----------|------|-------|---------|---------|
| gamma | Real | [0, 1] | 1.0 | Stabilization strength |
| k | Integer | ≥ 1 | 1 | Derivative order |
| c_F | Real | > 0 | 1.0 | Penalty constant |

### Example Usage

```ini
[InterfaceKernels]
  [ghost_penalty]
    type = GhostPenaltyKernel
    variable = u
    neighbor_var = u
    gamma = 1.0
    k = 1
    c_F = 1.0
  []
[]
```

## Level Set Integration

### Setting up the Level Set Function

```ini
[Functions]
  [level_set_func]
    type = ParsedFunction
    expression = 'sqrt((x-0.5)^2 + (y-0.5)^2) - 0.5'
  []
[]

[AuxVariables]
  [phi]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [phi_kernel]
    type = FunctionAux
    variable = phi
    function = level_set_func
    execute_on = 'INITIAL LINEAR'
  []
[]
```

## Testing Strategy

### Convergence Study

Test across mesh sizes to verify $O(h)$ H1 convergence:

```bash
for h in 1/5 1/10 1/20 1/40; do
  # Modify mesh size and run
  mpirun -np 4 moose_opt -i test_gp.i Mesh/nx=$(5*h) Mesh/ny=$(5*h)
done
```

### Condition Number Monitoring

Enable KSP monitoring:

```ini
[Executioner]
  petsc_options_iname = '-ksp_monitor'
  petsc_options_value = 'true'
  petsc_view = true
[]
```

## Phase 2: Cut Cell Quadrature (Next Phase)

### User Object Architecture

```cpp
class CutCellQuadratureUserObject : public UserObject
{
public:
  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & uo) override;
  virtual void finalize() override;
  
private:
  bool isCutElement(const Elem* elem);
  std::vector<Point> marchingCubes(const Elem* elem);
};
```

### Marching Cubes Algorithm

1. Evaluate Level Set at all element nodes
2. Check for sign change: $\min(\phi) \times \max(\phi) < 0$
3. Find edge intersections via linear interpolation
4. Look up MC pattern from pre-computed table
5. Generate sub-element mesh and quadrature points

## Building and Testing

### Build Commands

```bash
# Compile the module
cd modules/cutfem
make -j4 METHOD=opt

# Run tests
make test

# Format code
make format

# Clean
make clobber
```

### Run Specific Test

```bash
cd test
../../../moose_test-opt -i tests/ghost_penalty/test_gp.i
```

## Debugging

### Common Issues

1. **Compilation errors**
   - Check `MOOSE_DIR` environment variable
   - Verify C++17 compiler support
   - Check header includes

2. **Test failures**
   - Compare output with gold standard
   - Check mesh size and convergence
   - Verify Level Set function evaluation

3. **Performance issues**
   - Profile with `-preload libpapi.so` (PAPI)
   - Check linear solver iteration count
   - Verify matrix assembly time

## References

For detailed mathematical derivations, see:
- IMPLEMENTATION_STRATEGY.md (root directory)
- Burman et al. (2025) Acta Numerica
- Larson & Zahedi (2020) IMA J. Numerical Analysis
