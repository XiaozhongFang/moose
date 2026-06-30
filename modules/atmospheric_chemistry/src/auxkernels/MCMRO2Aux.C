//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRO2Aux.h"

registerMooseObject("AtmosphericChemistryApp", MCMRO2Aux);

InputParameters
MCMRO2Aux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredCoupledVar("ro2_species",
      "List of peroxy radical species variables to sum for RO2");
  params.addClassDescription(
      "AuxKernel that computes RO2 = sum(peroxy radical species). "
      "Used in coupled (FEM) mode for diagnostic output.");
  return params;
}

MCMRO2Aux::MCMRO2Aux(const InputParameters & params)
  : AuxKernel(params)
{
  // Couple to each RO2 species and store its value pointer
  unsigned int n_ro2 = coupledComponents("ro2_species");
  _ro2_vals.resize(n_ro2);
  for (unsigned int i = 0; i < n_ro2; ++i)
    _ro2_vals[i] = &coupledValue("ro2_species", i);
}

Real
MCMRO2Aux::computeValue()
{
  Real sum = 0.0;
  for (const auto & val : _ro2_vals)
    sum += (*val)[_qp];
  return sum;
}
