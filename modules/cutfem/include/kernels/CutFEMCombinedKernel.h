#pragma once

#include "Kernel.h"

class CutCellQuadratureUserObject;

/**
 * Combined diffusion + source kernel using cut cell quadrature for CutFEM.
 *
 * On cut elements, uses Assembly::reinitAtPhysical to evaluate FE data
 * at sub-triangle Gauss quadrature points. Both diffusion and source
 * terms integrate only over the physical domain (φ < 0).
 *
 * On uncut elements, falls back to standard quadrature.
 *
 * This avoids the multi-kernel issue by combining all terms into one kernel.
 *
 * Usage:
 * \code
 * [Kernels]
 *   [cutfem]
 *     type = CutFEMCombinedKernel
 *     variable = u
 *     cut_cell_quadrature = cut_cell_quad
 *     source = 1.0
 *   []
 * []
 * \endcode
 */
class CutFEMCombinedKernel : public Kernel
{
public:
  static InputParameters validParams();

  CutFEMCombinedKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;

private:
  const CutCellQuadratureUserObject & _cut_quad;
  Real _source_value;
};
