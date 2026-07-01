//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NoOpScalarKernel.h"

registerMooseObject("AtmosphericChemistryApp", NoOpScalarKernel);

InputParameters
NoOpScalarKernel::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addClassDescription("No-op scalar kernel with zero residual. "
                             "Satisfies MOOSE kernel coverage check.");
  return params;
}

NoOpScalarKernel::NoOpScalarKernel(const InputParameters & params)
  : ODEKernel(params)
{
}
