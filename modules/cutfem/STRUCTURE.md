# CutFEM Module Directory Structure

This document describes the MOOSE-compliant directory structure for the CutFEM module.

## Standard MOOSE Module Layout

```
modules/cutfem/
├── doc/
│   └── content/                          # Documentation (Markdown)
│       ├── index.md                      # Main documentation page
│       ├── overview.md                   # Project overview & theory
│       ├── implementation.md             # Implementation details
│       └── workflow.md                   # Development workflow
│
├── include/                              # Header files
│   ├── base/                             # Base/application headers
│   │   └── CutFEMApp.h                   # Main application header
│   │
│   ├── kernels/                          # Kernel headers
│   │   └── GhostPenaltyKernel.h
│   │
│   ├── userobjects/                      # UserObject headers
│   │   └── CutCellQuadratureUserObject.h
│   │
│   └── utils/                            # Utility class headers
│       ├── MarchingCubes.h
│       └── CutGeometry.h
│
├── src/base/                             # Base/application implementations
│   │   └── CutFEMApp.C                   # Main application
│   │
│   ├── kernels/                          # Kernel implementations
│   │   ├── GhostPenaltyKernel.C
│   │   └── (CutFEMInterfaceKernel.C)
│   │
│   ├── userobjects/                      # UserObject implementations
│   │   └── CutCellQuadratureUserObject.C
│   │
│   └── utils/                            # Utility implementations
│       ├── MarchingCubes.C
│       └── CutGeometry.C
│   └── CutFEMApp.C                       # Main application
│
├── test/
│   └── tests/                            # Test specification directory
│       ├── ghost_penalty/                # Phase 1 tests
│       │   ├── test_gp.i                 # Main test input file
│       │   ├── test_gp_gold/             # Gold standard output
│       │   │   └── test_gp_out.e         # Exodus gold file
│       │   └── tests.txt                 # Test specification
│       │
│       ├── cut_cells/                    # Phase 2 tests (future)
│       │   ├── test_cc.i
│       │   └── tests.txt
│       │
│       └── surface_pde/                  # Phase 3 tests (future)
│           ├── test_sp.i
│           └── tests.txt
│
├── examples/
│   ├── poisson_with_ghost_penalty.i      # Phase 1 example
│   ├── circular_interface.i              # Phase 2 example (future)
│   └── evolving_interface.i              # Phase 3 example (future)
│
├── README.md                             # Project readme
├── Makefile                              # Build configuration
├── .gitignore                            # Git ignore rules
├── .clang-format                         # Code formatting rules
│
└── IMPLEMENTATION_STRATEGY.md            # Detailed technical strategy
    (root level - reference document)
```

## File Naming Conventions

### Header Files (.h)
- Location: `include/`
- Pattern: `include/category/ClassName.h`
- Example: `include/kernels/GhostPenaltyKernel.h`

### Source Files (.C)
- Location: `src/`
- Pattern: `src/category/ClassName.C`
- Must match header file structure
- Example: `src/kernels/GhostPenaltyKernel.C`

### Test Input Files (.i)
- Location: `test/tests/category/`
- Pattern: `test_category.i` or `test_specific_feature.i`
- Example: `test/tests/ghost_penalty/test_gp.i`

### Test Specification Files
- Location: `test/tests/category/`
- Name: `tests.txt`
- Defines which tests to run and their parameters

### Documentation Files (.md)
- Location: `doc/content/`
- Standard names:
  - `index.md` - Main page
  - `overview.md` - Project overview
  - `implementation.md` - Implementation guide
  - `workflow.md` - Development workflow

### Configuration Files
- `.gitignore` - Git ignore patterns
- `.clang-format` - Code formatting rules
- `Makefile` - Build system
- `README.md` - Project readme

## Directory Organization Principles

### By Category
Files are grouped by **type/category** rather than **phase**:
- `kernels/` contains all kernel-related code
- `userobjects/` contains all user object code
- `utils/` contains utility functions

### Separation of Concerns
- **Headers** go in `include/`
- **Implementations** go in `src/`
- **Tests** go in `test/tests/`
- **Documentation** goes in `doc/content/`
- **Examples** go in `examples/`

### Nested Symmetry
The directory structure in `src/` **mirrors** that in `include/`:
```
include/kernels/GhostPenaltyKernel.h    ← Header
src/kernels/GhostPenaltyKernel.C        ← Implementation
```

## Phases vs. File Organization

**Important**: Phases (1, 2, 3) do **NOT** create separate directory structures.

Instead:
- All Phase 1 files are in `include/kernels/`, `src/kernels/`, etc.
- All Phase 2 files are in `include/userobjects/`, `src/userobjects/`, etc.
- Files are organized by **functionality** not by **phase**

Exception:
- **Tests** can have phase-specific subdirectories:
  - `test/tests/ghost_penalty/` (Phase 1)
  - `test/tests/cut_cells/` (Phase 2)
  - `test/tests/surface_pde/` (Phase 3)

## Building and Testing

### Compile
```bash
cd modules/cutfem
make -j4 METHOD=opt
```

### Test Single Category
```bash
cd test/tests/ghost_penalty
../../../moose_test-opt
```

### Test All
```bash
cd modules/cutfem
make test
```

### Format Code
```bash
make format
```

## Adding New Files

When adding new functionality:

1. **Create header** in appropriate `include/` subdirectory
2. **Create source** in corresponding `src/` subdirectory
3. **Update** `include/CutFEMApp.h` to register the class
4. **Update** `Makefile` if new subdirectory was created
5. **Add tests** in `test/tests/` subdirectory
6. **Format code** with `make format`

## File Descriptions

### Base Application

| File | Type | Status | Purpose |
|------|------|--------|---------|
| `include/base/CutFEMApp.h` | Header | ✅ Complete | Application class declaration |
| `src/base/CutFEMApp.C` | Source | ✅ Complete | Application initialization & registration |

### Phase 1: Ghost Penalty (Current)

| File | Type | Status | Purpose |
|------|------|--------|---------|
| `GhostPenaltyKernel.h` | Header | ✅ Complete | Face-based penalty kernel |
| `GhostPenaltyKernel.C` | Source | 📝 To-do | Kernel implementation |
| `test_gp.i` | Input | ✅ Complete | Ghost Penalty test |
| `poisson_with_ghost_penalty.i` | Example | ✅ Complete | Usage example |

### Phase 2: Cut Cell Quadrature (Future)

| File | Type | Status | Purpose |
|------|------|--------|---------|
| `CutCellQuadratureUserObject.h` | Header | ✅ Complete | Dynamic quadrature |
| `CutCellQuadratureUserObject.C` | Source | 📝 To-do | Implementation |
| `MarchingCubes.h` | Header | 🔳 Stub | Marching Cubes algorithm |
| `MarchingCubes.C` | Source | 🔳 Stub | Algorithm implementation |

### Phase 3: Surface PDE (Future)

| File | Type | Status | Purpose |
|------|------|--------|---------|
| `SurfacePDEKernel.h` | Header | 🔳 Stub | Surface PDE kernel |
| `SurfacePDEKernel.C` | Source | 🔳 Stub | Implementation |

## References

- **MOOSE Module Naming**: https://mooseframework.inl.gov/modules/
- **Code Organization**: https://mooseframework.inl.gov/syntax/
- **Build System**: https://mooseframework.inl.gov/help/development/
