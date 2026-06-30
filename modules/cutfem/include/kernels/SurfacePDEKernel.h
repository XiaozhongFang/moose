#pragma once

#include "InterfaceKernel.h"

/**
 * Surface PDE Kernel for Laplace-Beltrami on implicit surfaces.
 *
 * Implements the weak form of the Laplace-Beltrami operator on a surface Γ:
 *
 *   ∫_Γ (P ∇u) · (P ∇v) dS
 *
 * where P = I - n⊗n is the tangential projection operator.
 * In components: (P ∇u) = ∇u - (n·∇u) n
 *
 * This kernel acts on a sideset representing the surface Γ. The surface
 * normal is taken from the sideset normal (_normals).
 *
 * For a complete surface PDE:
 *   -Δ_Γ u + u = f   on Γ
 * add a MassMatrix kernel on the same boundary.
 *
 * References:
 * - Olshanskii, Reusken, Grande (2009) "A finite element method for
 *   surface PDEs" (Trace FEM)
 * - Larson & Zahedi (2020) "Stabilization of high order cut FEM on surfaces"
 *
 * Usage:
 * \code
 * [InterfaceKernels]
 *   [surface_diffusion]
 *     type = SurfacePDEKernel
 *     variable = u
 *     neighbor_var = u
 *     boundary = 'interface'
 *   []
 * []
 * \endcode
 */
class SurfacePDEKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  SurfacePDEKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;
};
