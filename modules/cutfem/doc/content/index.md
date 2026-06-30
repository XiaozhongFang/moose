# CutFEM Module Documentation

**Cut Finite Element Methods for MOOSE Framework**

This module implements the CutFEM methodology for discretizing PDEs on domains defined by implicit functions (Level Set), without requiring the mesh to be aligned with domain boundaries.

## Quick Start

```bash
# Compile
make -C modules/cutfem -j2 METHOD=opt

# Run example (split-domain Poisson + LevelSet cut detection)
cd modules/cutfem
./cutfem-opt -i examples/poisson_with_ghost_penalty.i

# Run tests
cd modules/cutfem/test
./cutfem-opt -i tests/ghost_penalty/test_gp.i
./cutfem-opt -i tests/cut_cells/test_cut_detection.i
./cutfem-opt -i tests/surface_pde/test_laplace_beltrami.i
```

## Implementation Status

| Phase | Component | Files | Status |
|-------|-----------|-------|--------|
| **1** | GhostPenaltyKernel | `kernels/GhostPenaltyKernel.h/.C` | ✅ Done |
| **1** | Tests & verification | `test/tests/ghost_penalty/` (5 tests) | ✅ Done |
| **1** | Documentation | `doc/content/source/kernels/GhostPenaltyKernel.md` | ✅ Done |
| **2** | CutCellQuadratureUserObject | `userobjects/CutCellQuadratureUserObject.h/.C` | ✅ Done |
| **2** | MarchingCubes2D (16-case lookup) | `utils/MarchingCubes.h/.C` | ✅ Done |
| **2** | Triangle Gauss quadrature | `utils/MarchingCubes.h/.C` | ✅ Done |
| **2** | Interface 1D quadrature | `utils/MarchingCubes.h/.C` | ✅ Done |
| **2** | Tests | `test/tests/cut_cells/` (5 tests) | ✅ Done |
| **3** | SurfacePDEKernel (Laplace-Beltrami) | `kernels/SurfacePDEKernel.h/.C` | ✅ Done |
| **3** | SurfaceStabilizationKernel | `kernels/SurfaceStabilizationKernel.h/.C` | ✅ Done |
| **3** | H-J LevelSet advection | `test/tests/surface_pde/test_hj_evolution.i` | ✅ Done |
| **—** | CutFEMDiffusion (cut cell solve) | `kernels/CutFEMDiffusion.h/.C` | ✅ Done |
| **—** | CutFEMCombinedKernel | `kernels/CutFEMCombinedKernel.h/.C` | ✅ Done |

## Key Features

### Phase 1: Ghost Penalty Stabilization
- `GhostPenaltyKernel`: face-based gradient jump penalty (InterfaceKernel)
- Parameters: gamma ∈ [0,1], k ≥ 1 (derivative order), c_F > 0
- [Documentation](source/kernels/GhostPenaltyKernel.md)

### Phase 2: Cut Cell Integration
- `CutCellQuadratureUserObject`: LevelSet function evaluation at element nodes, sign-change cut detection
- `MarchingCubes2D`: 16-case lookup table for quadrilateral elements, sub-triangle generation
- Gauss quadrature on sub-triangles (orders 1-4, 1/3/4/6 points)
- Interface 1D Gauss quadrature (orders 1-4) for surface integration

### Phase 3: Surface PDE & Evolution
- `SurfacePDEKernel`: Laplace-Beltrami via tangential projection P = I - n⊗n
- `SurfaceStabilizationKernel`: ghost penalty on the interface (mixed stabilization)
- Hamilton-Jacobi LevelSet advection for interface evolution

### Solvers
- `CutFEMDiffusion`: single-kernel diffusion on LevelSet-defined domains using `Assembly::reinitAtPhysical`
- `CutFEMCombinedKernel`: diffusion + body force combined (avoids multi-kernel quadrature conflict)

## Complete File List

```
include/
├── base/CutFEMApp.h
├── kernels/
│   ├── GhostPenaltyKernel.h
│   ├── SurfacePDEKernel.h
│   ├── SurfaceStabilizationKernel.h
│   ├── CutFEMDiffusion.h
│   └── CutFEMCombinedKernel.h
├── userobjects/
│   └── CutCellQuadratureUserObject.h
└── utils/
    └── MarchingCubes.h

src/
├── base/CutFEMApp.C
├── main.C
├── kernels/
│   ├── GhostPenaltyKernel.C
│   ├── SurfacePDEKernel.C
│   ├── SurfaceStabilizationKernel.C
│   ├── CutFEMDiffusion.C
│   └── CutFEMCombinedKernel.C
├── userobjects/
│   └── CutCellQuadratureUserObject.C
└── utils/
    └── MarchingCubes.C

test/
├── include/base/CutFEMTestApp.h
├── src/base/CutFEMTestApp.C
└── tests/
    ├── ghost_penalty/
    │   ├── test_gp.i / test_no_gp.i
    │   ├── test_gp_fine.i / test_no_gp_fine.i
    │   ├── test_convergence.i
    │   └── tests
    ├── cut_cells/
    │   ├── test_cut_detection.i
    │   ├── test_quadrature.i
    │   ├── test_diagonal_cut.i
    │   ├── test_boundary_tangent.i
    │   ├── test_combined_kernel.i
    │   └── tests
    └── surface_pde/
        ├── test_laplace_beltrami.i
        ├── test_hj_evolution.i
        └── tests

examples/
└── poisson_with_ghost_penalty.i

doc/content/
├── index.md
├── overview.md
├── implementation.md
├── workflow.md
└── source/kernels/GhostPenaltyKernel.md
```

## Verification

| Test | Input | Result | Metric |
|------|-------|--------|--------|
| GP basic | `test_gp.i` | ✅ Converged | 2×2 split domain |
| GP fine mesh | `test_gp_fine.i` | ✅ Converged | 20×20 split domain |
| Convergence | `test_convergence.i` | ✅ L2 ~ 1e-14 | Machine precision for u=x |
| Cut detection (circle) | `test_cut_detection.i` | ✅ 28/100 | Circle r=0.3 on 10×10 |
| Cut detection (diag) | `test_diagonal_cut.i` | ✅ 10/16 | Interface y=x on 4×4 |
| Area conservation | `test_cut_detection.i` | ✅ diff=0 | Weights = sub-triangle area |
| Surface PDE | `test_laplace_beltrami.i` | ✅ Converged | Laplace-Beltrami on interface |
| H-J evolution | `test_hj_evolution.i` | ✅ 10 steps | LevelSet advection |
| Combined kernel | `test_combined_kernel.i` | ✅ Converged | CutFEMCombinedKernel |

## References

- Burman, E., Hansbo, P., Larson, M. G., & Massing, A. (2015). "CutFEM: Discretizing geometry and partial differential equations." *IJNME*, 104(7), 472-501.
- Larson, M. G., & Zahedi, S. (2020). "Stabilization of high order cut finite element methods on surfaces." *IMA J. Numerical Analysis*, 40(3), 1702-1745.
- Wichrowski, M. (2026). "Matrix-free ghost penalty evaluation via tensor product factorization." *Computers & Mathematics*, 211, 109-121.
- Burman et al. (2025). "Cut finite element methods." *Acta Numerica*, 34, 1-121.

## Contact

- MOOSE Discussions: https://github.com/idaholab/moose/discussions
- MOOSE Website: https://mooseframework.inl.gov
- GitHub: https://github.com/XiaozhongFang/moose/tree/cutfem
