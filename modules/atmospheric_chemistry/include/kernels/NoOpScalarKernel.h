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

/**
 * No-op scalar kernel that contributes zero residual and zero Jacobian.
 * Used in PETSc TS box mode where MCMBoxModel::execute() handles
 * the ODE integration — MOOSE needs at least one kernel per subdomain,
 * and this kernel satisfies that check without altering the solution.
 */
class NoOpScalarKernel : public ODEKernel
{
public:
  static InputParameters validParams();
  NoOpScalarKernel(const InputParameters & params);

  virtual Real computeQpResidual() override { return 0.0; }
  virtual Real computeQpJacobian() override { return 0.0; }
  virtual Real computeQpOffDiagJacobianScalar(unsigned int /*jvar*/) override { return 0.0; }
};
