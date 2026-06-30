# CutFEM Module File Structure

Current file layout for the CutFEM module under `modules/cutfem/`.

## Source Files

```
include/
├── base/
│   └── CutFEMApp.h
├── kernels/
│   ├── GhostPenaltyKernel.h              # Phase 1: gradient jump penalty
│   ├── SurfacePDEKernel.h                # Phase 3: Laplace-Beltrami
│   ├── SurfaceStabilizationKernel.h       # Phase 3: mixed stabilization
│   ├── CutFEMDiffusion.h                 # Cut cell diffusion solve
│   └── CutFEMCombinedKernel.h            # Diffusion + source combined
├── userobjects/
│   └── CutCellQuadratureUserObject.h     # Phase 2: cut detection + cache
└── utils/
    └── MarchingCubes.h                   # Phase 2: 2D MC + Gauss quadrature

src/
├── main.C                                # Entry: Moose::main<CutFEMTestApp>
├── base/CutFEMApp.C
├── kernels/
│   ├── GhostPenaltyKernel.C
│   ├── SurfacePDEKernel.C
│   ├── SurfaceStabilizationKernel.C
│   ├── CutFEMDiffusion.C
│   └── CutFEMCombinedKernel.C
├── userobjects/CutCellQuadratureUserObject.C
└── utils/MarchingCubes.C
```

## Test Files

```
test/
├── include/base/CutFEMTestApp.h
├── src/base/CutFEMTestApp.C
└── tests/
    ├── ghost_penalty/                     # Phase 1 validation
    │   ├── test_gp.i / test_no_gp.i
    │   ├── test_gp_fine.i / test_no_gp_fine.i
    │   ├── test_convergence.i
    │   └── tests
    ├── cut_cells/                         # Phase 2 validation
    │   ├── test_cut_detection.i
    │   ├── test_quadrature.i
    │   ├── test_diagonal_cut.i
    │   ├── test_boundary_tangent.i
    │   ├── test_combined_kernel.i
    │   └── tests
    └── surface_pde/                       # Phase 3 validation
        ├── test_laplace_beltrami.i
        ├── test_hj_evolution.i
        └── tests
```

## Documentation

```
doc/content/
├── index.md                              # Main module documentation
├── overview.md                           # Project overview
├── implementation.md                     # Implementation guide
├── workflow.md                           # Development workflow
└── source/kernels/
    └── GhostPenaltyKernel.md             # GP kernel reference
```

## Configuration

```
Makefile                                   # MODULE_NAME := cutfem
.clang-format                              # LLVM style, column 100
.gitignore
```

## Examples

```
examples/
└── poisson_with_ghost_penalty.i           # CutFEM demo
```

## Implementation Status

| Component | Phase | Status |
|-----------|-------|--------|
| GhostPenaltyKernel | 1 | ✅ Complete, documented |
| CutCellQuadratureUserObject | 2 | ✅ Complete (cut detection + cache) |
| MarchingCubes2D | 2 | ✅ Complete (16-case + Gauss quad) |
| SurfacePDEKernel | 3 | ✅ Complete (tangential projection) |
| SurfaceStabilizationKernel | 3 | ✅ Complete (mixed stabilization) |
| CutFEMDiffusion | — | ✅ Complete (reinitAtPhysical solve) |
| CutFEMCombinedKernel | — | ✅ Complete (combined terms) |
