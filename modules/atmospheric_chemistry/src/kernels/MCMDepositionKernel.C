//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMDepositionKernel.h"
registerMooseObject("AtmosphericChemistryApp", MCMDepositionKernel);

InputParameters MCMDepositionKernel::validParams()
{
  auto p = Kernel::validParams();
  p.addRequiredParam<Real>("rate", "First-order deposition rate constant [/s]");
  p.addClassDescription("First-order deposition loss: dC/dt = -k*C");
  return p;
}
MCMDepositionKernel::MCMDepositionKernel(const InputParameters & p)
  : Kernel(p), _rate(getParam<Real>("rate")) {}
