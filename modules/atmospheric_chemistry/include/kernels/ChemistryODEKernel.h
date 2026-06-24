//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ODEKernel.h"
#include "MCMBoxModel.h"

class MooseVariableScalar;

/**
 * ODE kernel that bridges MOOSE's ScalarVariable system to the
 * MCMBoxModel computation engine for atmospheric chemistry box models.
 *
 * One instance per chemical species. Uses MCMBoxModel's cached
 * single-species interface (getDCdt / getJacobianDiagonal /
 * getJacobianOffDiagonal) to avoid redundant full-system computations
 * across kernels.
 *
 * The residual follows the ODEKernel convention:
 *   R = du/dt - f(u, t) = 0
 * so the chemical source contribution is -dC/dt (negative RHS).
 *
 * Accesses species ScalarVariable values directly from the SubProblem
 * (bypasses the ScalarCoupleable framework since species names are
 * dynamic and not known at compile time).
 */
class ChemistryODEKernel : public ODEKernel
{
public:
  static InputParameters validParams();

  ChemistryODEKernel(const InputParameters & params);

  virtual void reinit() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;

  /// Reference to the centralized box model
  const MCMBoxModel & _box_model;

  /// Index of this species in the mechanism (0..nSpecies-1)
  const unsigned int _species_idx;

  /// Pointers to all species' ScalarVariable objects (index = species index)
  std::vector<MooseVariableScalar *> _species_vars;

  /// Build the full concentration vector from scalar variable values.
  /// Returns reference to pre-allocated _C_buffer (Per.14 — no allocation per call).
  const std::vector<Real> & _buildC() const;

  /// Pre-allocated concentration buffer; reused across all residual/Jacobian calls.
  mutable std::vector<Real> _C_buffer;
};
