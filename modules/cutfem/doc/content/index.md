# CutFEM Module Documentation

**Cut Finite Element Methods for MOOSE Framework**

This module implements the CutFEM methodology for discretizing PDEs on domains defined by implicit functions (Level Set), without requiring the mesh to be aligned with domain boundaries.

## Quick Start

1. **Overview**: Read [overview.md](./overview.md)
2. **Implementation**: See [implementation.md](./implementation.md)  
3. **Development**: Check [workflow.md](./workflow.md)

## Key Features

- **Phase 1**: Ghost Penalty stabilization for unfitted meshes
- **Phase 2**: Non-conforming discontinuous integration via Level Set
- **Phase 3**: Dynamic interface evolution with Hamilton-Jacobi coupling

## References

- Burman et al. (2025) "Cut finite element methods" *Acta Numerica* 34:1-121
- Larson & Zahedi (2020) "Stabilization of high order cut FEM" *IMA J. Numerical Analysis* 40(3):1702-1745
- Wichrowski (2026) "Matrix-free ghost penalty evaluation" *Computers & Mathematics* 211:109-121

## Module Structure

```
include/
├── kernels/          # CutFEM kernels
├── userobjects/      # Quadrature and utility objects
├── utils/            # Helper algorithms
└── CutFEMApp.h

src/
├── kernels/
├── userobjects/
├── utils/
└── CutFEMApp.C

examples/
└── poisson_with_ghost_penalty.i

test/
└── tests/
    └── ghost_penalty/
        └── test_gp.i
```

## Contact

- MOOSE Discussions: https://github.com/idaholab/moose/discussions
- MOOSE Website: https://mooseframework.inl.gov
