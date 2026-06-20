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
#include "Function.h"

/** Emission source: dC/dt = +E(t,x) [molec/cm^3/s] */
class MCMEmissionKernel : public Kernel
{
public:
  static InputParameters validParams();
  MCMEmissionKernel(const InputParameters & params);
protected:
  Real computeQpResidual() override { return -_func.value(_t, _q_point[_qp]) * _test[_i][_qp]; }
  Real computeQpJacobian() override { return 0.0; }
  const Function & _func;
};
