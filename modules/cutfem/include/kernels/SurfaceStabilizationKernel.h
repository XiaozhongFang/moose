#pragma once

#include "InterfaceKernel.h"

/**
 * Surface Stabilization Kernel for CutFEM (Phase 3).
 *
 * Implements the surface ghost penalty term acting directly on Γ_h:
 *
 *   s_{h,Γ}(w,v) = Σ c'_j · h^{2(j-1+γ)} ∫_Γ D_n^j w · D_n^j v dS
 *
 * where D_n^j v = (n·∇)^j v is the j-th normal derivative on the surface.
 * This penalizes normal derivatives on the interface itself (not across it).
 *
 * Combined with GhostPenaltyKernel (face stabilization), this gives the
 * mixed stabilization that ensures condition number κ ~ O(h^{-2})
 * independent of polynomial order p (Larson & Zahedi 2020).
 *
 * References:
 * - Larson & Zahedi (2020) "Stabilization of high order cut FEM on surfaces"
 *   IMA J. Numerical Analysis, 40(3), 1702-1745.
 *
 * Usage:
 * \code
 * [InterfaceKernels]
 *   [surface_stab]
 *     type = SurfaceStabilizationKernel
 *     variable = u
 *     neighbor_var = u
 *     boundary = 'interface'
 *     gamma = 1.0
 *     k = 1
 *     c_Gamma = 1.0
 *   []
 * []
 * \endcode
 */
class SurfaceStabilizationKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  SurfaceStabilizationKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

private:
  /// Stabilization parameter gamma ∈ [0,1]
  Real _gamma;

  /// Derivative order k
  unsigned int _k;

  /// Surface penalty constant
  Real _c_Gamma;

  /// Precomputed exponent: 2(k-1+gamma)
  Real _exponent;

  /**
   * Compute the stabilization penalty coefficient.
   * Uses element volume to estimate characteristic length h.
   */
  Real computePenaltyCoeff();
};
