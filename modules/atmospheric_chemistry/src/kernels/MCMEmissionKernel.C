//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMEmissionKernel.h"
registerMooseObject("AtmosphericChemistryApp", MCMEmissionKernel);

InputParameters MCMEmissionKernel::validParams()
{
  auto p = Kernel::validParams();
  p.addRequiredParam<FunctionName>("function", "Emission rate function [molec/cm^3/s]");
  p.addClassDescription("Emission source: dC/dt = +E(t,x)");
  return p;
}
MCMEmissionKernel::MCMEmissionKernel(const InputParameters & p)
  : Kernel(p), _func(getFunction("function")) {}
