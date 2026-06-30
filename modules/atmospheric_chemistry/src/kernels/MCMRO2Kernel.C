//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRO2Kernel.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", MCMRO2Kernel);

InputParameters
MCMRO2Kernel::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addRequiredParam<std::vector<VariableName>>(
      "species_variables", "Names of all species ScalarVariables (in species order)");
  params.addClassDescription(
      "ODE kernel that sets RO2 = sum(peroxy radicals). "
      "Uses an algebraic formulation (mass-matrix style) to pin RO2 "
      "to the value computed by MCMBoxModel::getRO2Sum().");
  return params;
}

MCMRO2Kernel::MCMRO2Kernel(const InputParameters & params)
  : ODEKernel(params),
    _box_model(getUserObject<MCMBoxModel>("box_model"))
{
  const auto & var_names = getParam<std::vector<VariableName>>("species_variables");
  _species_vars.reserve(var_names.size());
  for (const auto & name : var_names)
  {
    if (_sc_fe_problem.hasScalarVariable(name))
      _species_vars.push_back(&_sc_fe_problem.getScalarVariable(_sc_tid, name));
    else
      mooseError("MCMRO2Kernel: ScalarVariable '", name, "' not found");
  }
}

Real
MCMRO2Kernel::computeQpResidual()
{
  // Build concentration vector from scalar variables
  _C_buffer.assign(_species_vars.size(), 0.0);
  for (unsigned int i = 0; i < _species_vars.size(); ++i)
    _C_buffer[i] = _species_vars[i]->sln()[0];

  // If units=ppb, convert to molec/cm³ for internal RO2 computation
  Real ro2_sum = _box_model.getRO2Sum(_C_buffer);
  if (_box_model.unitsPPB())
    ro2_sum /= _box_model.ppbToMolec();

  // Algebraic constraint: R = u - RO2_sum = 0
  return _u[_i] - ro2_sum;
}

Real
MCMRO2Kernel::computeQpJacobian()
{
  // dR/du = 1
  return 1.0;
}
