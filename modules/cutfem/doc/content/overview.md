# CutFEM Module Overview

## What is CutFEM?

Cut Finite Element Methods (CutFEM) discretize PDEs on domains defined by implicit functions (Level Set $\phi$) without requiring the mesh to be aligned with domain boundaries. This avoids the expensive remeshing required by traditional fitted FEM.

### Problem Statement

Given:
- Domain $\Omega$ defined by Level Set: $\Omega = \{\mathbf{x} : \phi(\mathbf{x}) < 0\}$
- PDE: $-\nabla \cdot (\alpha \nabla u) = f$ in $\Omega$, $u = 0$ on $\partial\Omega$
- Background mesh $\tilde{\mathcal{T}}_h$ **independent of** $\Omega$

Find: Approximate solution $u_h$ with:
- ✅ **Optimal convergence**: $\|u - u_h\|_{H^1} = O(h)$, $\|u - u_h\|_{L^2} = O(h^2)$
- ✅ **Stable conditioning**: $\kappa(\mathcal{A}_h) = O(h^{-2})$ (independent of mesh-interface position)
- ✅ **Robustness**: Works for arbitrary interface positions

## Three-Phase Implementation Plan

### Phase 1: Ghost Penalty Stabilization (2-3 months)

**Goal**: Implement Nitsche's method + Ghost Penalty to achieve optimal conditioning

**Key Formula**:
$$a_h(u,v) + s_h(u,v) = l_h(v)$$

where:
- **Nitsche part**: Enforces boundary conditions weakly
- **Ghost Penalty**: $s_{h,F}(w,v) = \sum_j c_j h^{2(j-1+\gamma)} [D_n^j w]_F \cdot [D_n^j v]_F$

**Deliverables**:
- ✅ `GhostPenaltyKernel` class
- ✅ Example: Poisson on circular unfitted domain
- ✅ Validation: Condition number improves 50%+

### Phase 2: Non-conforming Integration (3-4 months)

**Goal**: Accurate integration on cut elements via Marching Cubes

**Components**:
- Level Set integration
- Dynamic quadrature point generation
- Element sub-division algorithm

**Expected Results**:
- Integration error < 1e-8
- Same convergence rate as Phase 1

### Phase 3: Surface PDE & Dynamic Evolution (3 months)

**Goal**: Solve surface PDEs and evolve interfaces via Hamilton-Jacobi

**Features**:
- Laplace-Beltrami operator on implicit surfaces
- Mixed stabilization (face + surface terms)
- Time-dependent interface evolution

## Technical Highlights

### 1. Tensor Product Optimization (Wichrowski 2026)

For Cartesian grids, Ghost Penalty evaluation complexity:
$$O(k^{2d}) \rightarrow O(k^{d+1})$$

via factorization: $\phi_{ijk}(\mathbf{x}) = \phi_i(x_1)\phi_j(x_2)\phi_k(x_3)$

### 2. Mixed Stabilization (Larson & Zahedi 2020)

**Face stabilization**: $s_{h,F} = \sum_j c_j h^{2(j-1+\gamma)} [D_n^j w]_F [D_n^j v]_F$

**Surface stabilization**: $s_{h,\Gamma} = \sum_j c'_j h^{2(j-1+\gamma)} (D_n^j w)(D_n^j v)$

**Result**: $\kappa(\mathcal{A}_h) = O(h^{-2})$ **independent of polynomial degree $p$**

### 3. Automatic Marching Cubes

Algorithm:
1. **Detect**: If $\min_i \phi_i \cdot \max_j \phi_j < 0$ → cut element
2. **Intersect**: Linear interpolation to find $\phi = 0$ crossing points
3. **Subdivide**: Standard MC patterns
4. **Integrate**: Gauss rule on each sub-element

## Success Criteria

| Phase | Criterion | Target |
|-------|-----------|--------|
| 1 | Condition number improvement | ≥ 50% |
| 1 | H1 convergence rate | O(h) |
| 1 | L2 convergence rate | O(h²) |
| 2 | Integration accuracy | < 1e-8 error |
| 2 | Support Level Set changes | ✅ Yes |
| 3 | Long-term stability | t > 100 steps |
| 3 | Conservation laws | Satisfied |

## Key Papers

| Paper | Year | Key Contribution |
|-------|------|------------------|
| Burman et al. | 2025 | Comprehensive CutFEM review (121 pages) |
| Larson & Zahedi | 2020 | High-order stabilization theory |
| Wichrowski | 2026 | Matrix-free tensor product optimization |
| Burman et al. | 2015 | Original CutFEM method |

## Next Steps

1. Read [implementation.md](./implementation.md) for code details
2. Check [workflow.md](./workflow.md) for development process
3. Run examples: `moose_test-opt -i examples/poisson_with_ghost_penalty.i`
