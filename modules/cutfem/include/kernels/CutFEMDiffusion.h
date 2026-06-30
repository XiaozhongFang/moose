#pragma once

#include "Kernel.h"

class CutCellQuadratureUserObject;

/**
 * Diffusion kernel using cut cell quadrature for CutFEM.
 *
 * On elements not cut by the interface, uses standard quadrature.
 * On cut elements, uses sub-triangle Gauss quadrature from
 * CutCellQuadratureUserObject via FE reinit at custom points.
 *
 * This enables solving PDEs on LevelSet-defined domains on
 * unfitted background meshes.
 *
 * Usage:
 * \code
 * [Kernels]
 *   [cut_diff]
 *     type = CutFEMDiffusion
 *     variable = u
 *     cut_cell_quadrature = cut_cell_quad
 *   []
 * []
 * \endcode
 */
class CutFEMDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  CutFEMDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;

private:
  /// Reference to the CutCellQuadratureUserObject
  const CutCellQuadratureUserObject & _cut_quad;
};
