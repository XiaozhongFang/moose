//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMConstraintKernel.h"

registerMooseObject("AtmosphericChemistryApp", MCMConstraintKernel);

InputParameters
MCMConstraintKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<FunctionName>("function",
      "Function providing the time-dependent constrained concentration value");
  params.addClassDescription(
      "Constrains a chemical species concentration to a prescribed function of time. "
      "Residual: R = u - f(t,x). Used for chemically 'constrained' species "
      "(e.g., NO, NO2 fixed to observations in AtChem2 mode).");
  return params;
}

MCMConstraintKernel::MCMConstraintKernel(const InputParameters & params)
  : Kernel(params), _func(getFunction("function"))
{
}

Real
MCMConstraintKernel::computeQpResidual()
{
  return _u[_qp] - _func.value(_t, _q_point[_qp]);
}

Real
MCMConstraintKernel::computeQpJacobian()
{
  return 1.0;
}
