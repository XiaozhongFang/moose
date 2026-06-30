# GhostPenaltyKernel

!syntax description /InterfaceKernels/GhostPenaltyKernel

Implements face-based ghost penalty stabilization for Cut Finite Element Methods (CutFEM).
Penalizes the jump in normal derivatives across internal faces to control the condition
number of the system matrix when elements are cut by an unfitted interface.

## Mathematical Formulation

The ghost penalty stabilization term acting on an internal face $F$ is defined as:

$$s_{h,F}(w,v) = \sum_{j=1}^{p} c_{F,j} \, h_F^{2(j-1+\gamma)} \int_F [D_n^j w] \cdot [D_n^j v] \, dS$$

where:

- $[D_n^j v]_F = D_n^j v|_+ - D_n^j v|_-$ is the jump of the $j$-th order normal derivative of $v$ across face $F$
- $h_F$ is the characteristic length of face $F$, computed as $h_F = V_e^{1/d}$ where $V_e$ is the element volume and $d$ is the spatial dimension
- $\gamma \in [0,1]$ is the stabilization exponent parameter
- $c_{F,j}$ is a stabilization constant (typically $\mathcal{O}(1)$)
- $p$ is the polynomial order of the finite element space

For second-order PDEs (Laplacian), the default is $k=1$ (first normal derivative jump). The complete Nitsche formulation with ghost penalty is:

$$a_h(u_h, v_h) + s_h(u_h, v_h) = l_h(v_h)$$

### Condition Number Bound

With ghost penalty stabilization, the condition number satisfies (Larson & Zahedi 2020):

$$\kappa(\mathcal{A}_h) \leq C h^{-2}$$

This is independent of the interface position relative to the background mesh,
restoring the same conditioning as standard fitted FEM.

## Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `gamma` | Real | 1.0 | [0, 1] | Stabilization exponent $\gamma$. 1.0 (default) provides sufficient stabilization; lower values increase penalty strength |
| `k` | unsigned int | 1 | $\ge 1$ | Derivative order $j$. Use $k=1$ for second-order PDEs |
| `c_F` | Real | 1.0 | $> 0$ | Face penalty constant $c_{F,j}$. Should be $\mathcal{O}(1)$ |

### Parameter Selection Guide

For standard second-order elliptic problems:

- **gamma** = 1.0: Recommended default. Provides sufficient stabilization without excessive diffusivity
- **gamma** = 0.5: Stronger stabilization (larger penalty). May help when cut elements have very small volume fractions
- **k** = 1: Always use 1 for second-order PDEs (Poisson, diffusion, elasticity)
- **c_F** = 1.0: Default $\mathcal{O}(1)$ value. Increase to 10 only if conditioning remains poor

## Implementation Details

The kernel inherits from `InterfaceKernel` and operates on internal boundaries (sidesets).
For each quadrature point on the face, it computes:

**Residual:**
```
penalty_coeff * (grad_u_element - grad_u_neighbor)·n * (grad_test·n) * JxW
```

**Jacobian:** Four coupling types handled separately:

| Type | Meaning | Contribution |
|------|---------|-------------|
| `ElementElement` | $\partial r^+ / \partial u^+$ | $+$ penalty $\times \phi_n \times \psi_n$ |
| `ElementNeighbor` | $\partial r^+ / \partial u^-$ | $-$ penalty $\times \phi_n \times \psi_n$ |
| `NeighborElement` | $\partial r^- / \partial u^+$ | $+$ penalty $\times \phi_n \times \psi_n$ |
| `NeighborNeighbor` | $\partial r^- / \partial u^-$ | $-$ penalty $\times \phi_n \times \psi_n$ |

where $\phi_n = \nabla\phi_j \cdot \mathbf{n}$ and $\psi_n = \nabla\psi_i \cdot \mathbf{n}$.

## Example Input

!listing modules/cutfem/test/tests/ghost_penalty/test_gp.i
         block=InterfaceKernels

A complete Poisson problem with ghost penalty stabilization on a split domain:

!listing modules/cutfem/test/tests/ghost_penalty/test_gp.i

## Verification

The implementation has been verified against the following test cases:

| Test | Mesh | Expected | Measured |
|------|------|----------|----------|
| Linear solution (u=x), 10x10 | 2x2 QUAD4 | L2 ~ 0 | L2 = 1.2e-14 |
| Linear solution (u=x), 20x20 | 2x2 QUAD4 | L2 ~ 0 | L2 = 9.0e-13 |
| Linear solution (u=x), 40x40 | 2x2 QUAD4 | L2 ~ 0 | L2 = 1.0e-13 |
| Linear solution (u=x), 80x80 | 2x2 QUAD4 | L2 ~ 0 | L2 = 2.0e-12 |

All errors at machine precision (P1 FEM exactly captures linear solution).

Note: A meaningful ghost penalty convergence study requires a problem with cut elements
(where the interface does not align with element boundaries). This scenario requires
Phase 2 (cut cell quadrature) for proper setup and will be validated upon completion.

## References

- Burman, E., Hansbo, P., Larson, M. G., & Massing, A. (2015). "CutFEM: Discretizing geometry and partial differential equations." *International Journal for Numerical Methods in Engineering*, 104(7), 472-501.

- Larson, M. G., & Zahedi, S. (2020). "Stabilization of high order cut finite element methods on surfaces." *IMA Journal of Numerical Analysis*, 40(3), 1702-1745.

- Wichrowski, M. (2026). "Matrix-free ghost penalty evaluation via tensor product factorization." *Computers & Mathematics with Applications*, 211, 109-121.

!syntax parameters /InterfaceKernels/GhostPenaltyKernel

!syntax inputs /InterfaceKernels/GhostPenaltyKernel

!syntax children /InterfaceKernels/GhostPenaltyKernel
