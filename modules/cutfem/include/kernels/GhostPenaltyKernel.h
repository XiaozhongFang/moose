#pragma once

#include "InterfaceKernel.h"

/**
 * Ghost Penalty Stabilization Kernel for CutFEM
 *
 * Implements the face-based ghost penalty stabilization term:
 *
 *   s_{h,F}(w,v) = sum_j c_{F,j} * h_F^{2(j-1+gamma)} * int_F [D_n^j w] · [D_n^j v] dS
 *
 * where:
 *   - [D_n^j v]_F is the jump in j-th order normal derivative across face F
 *   - h_F is the characteristic length scale of face F
 *   - gamma ∈ [0,1] is a stabilization parameter
 *   - c_{F,j} is a stabilization constant (typically O(1))
 *
 * References:
 * - Burman (2010), Burman & Hansbo (2012): Original ghost penalty formulation
 * - Larson & Zahedi (2020): High-order stabilization analysis
 *
 * Usage in input file:
 * \code
 * [InterfaceKernels]
 *   [ghost_penalty]
 *     type = GhostPenaltyKernel
 *     variable = u
 *     neighbor_var = u
 *     boundary = internal_surface
 *     gamma = 1.0
 *     k = 1
 *     c_F = 1.0
 *   []
 * []
 * \endcode
 */
class GhostPenaltyKernel : public InterfaceKernel
{
public:
  static InputParameters validParams();

  GhostPenaltyKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type) override;
  virtual Real computeQpJacobian(Moose::DGJacobianType type) override;

private:
  /// Stabilization parameter gamma ∈ [0,1]
  Real _gamma;

  /// Order of normal derivatives: j in D_n^j
  unsigned int _k;

  /// Face penalty constant c_{F,j}
  Real _c_F;

  /// Precomputed exponent: 2(k-1+gamma)
  Real _exponent;

  /**
   * Compute penalty coefficient: c_F * h_F^{2(k-1+gamma)}
   *
   * Uses element volume to estimate characteristic face length:
   *   h_F ≈ (element_volume)^{1/d}
   * where d is the spatial dimension.
   */
  Real computePenaltyCoeff();
};
