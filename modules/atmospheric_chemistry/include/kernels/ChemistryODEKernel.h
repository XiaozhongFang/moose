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

#include <memory>

class MooseVariableScalar;
class BoxIntegrator;

/**
 * ODE kernel that bridges MOOSE's ScalarVariable system to the
 * BoxIntegrator strategy for atmospheric chemistry box models.
 *
 * One instance per chemical species. Delegates all residual and
 * Jacobian evaluation to the BoxIntegrator — the kernel is agnostic
 * about whether the integrator runs in MOOSE-implicit or PETSc TS mode.
 * In self-driven (PETSc TS) mode, the integrator returns 0 for all
 * evaluations; MCMBoxModel::execute() handles the full-system solve.
 *
 * The residual follows the ODEKernel convention:
 *   R = du/dt - f(u, t) = 0
 * so the chemical source contribution is -dC/dt (negative RHS).
 *
 * Accesses species ScalarVariable values directly from the SubProblem
 * (bypasses the ScalarCoupleable framework since species names are
 * dynamic and not known at compile time).  Unit conversion (ppb ↔
 * molec/cm³) is applied per the `units` parameter.
 */
class ChemistryODEKernel : public ODEKernel
{
public:
  static InputParameters validParams();
  ChemistryODEKernel(const InputParameters & params);

  virtual void reinit() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;

protected:
  const BoxIntegrator & _integrator;
  const unsigned int _species_idx;
  /// ppb → molec/cm³ conversion factor (M/1e9), valid when box_model uses ppb units
  const Real _ppb_to_molec;
  std::vector<const MooseVariableScalar *> _species_vars;
  mutable std::vector<Real> _C_buffer;
  const std::vector<Real> & _buildC() const;
};
