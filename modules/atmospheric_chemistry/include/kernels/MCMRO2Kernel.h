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
 * ODE kernel that pins the RO2 scalar variable to the sum of peroxy radicals.
 *
 * Algebraic constraint (no TimeDerivative):
 *   R = u - RO2_sum = 0
 *
 * This makes RO2 an algebraic variable whose value at each timestep
 * is the sum of all peroxy radical species concentrations, computed
 * by MCMBoxModel::getRO2Sum(). No time integration is applied —
 * RO2 follows the instantaneous RO2 sum.
 */
class MCMRO2Kernel : public ODEKernel
{
public:
  static InputParameters validParams();
  MCMRO2Kernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const MCMBoxModel & _box_model;
  std::vector<const MooseVariableScalar *> _species_vars;
  mutable std::vector<Real> _C_buffer;
};
