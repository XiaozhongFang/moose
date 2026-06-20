//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Kernel.h"

/** First-order deposition: dC/dt = -k * C, where k is rate (/s). */
class MCMDepositionKernel : public Kernel
{
public:
  static InputParameters validParams();
  MCMDepositionKernel(const InputParameters & params);
protected:
  Real computeQpResidual() override { return _rate * _u[_qp] * _test[_i][_qp]; }
  Real computeQpJacobian() override { return _rate * _phi[_j][_qp] * _test[_i][_qp]; }
  Real _rate;
};
