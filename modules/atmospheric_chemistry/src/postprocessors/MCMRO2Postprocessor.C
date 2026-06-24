//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRO2Postprocessor.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", MCMRO2Postprocessor);

InputParameters
MCMRO2Postprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addParam<std::vector<VariableName>>(
      "species_variables", std::vector<VariableName>(),
      "Names of all species ScalarVariables (in species order). "
      "If not provided, returns 0 (box mode test compatibility).");
  params.addClassDescription(
      "Computes total RO2 (peroxy radical) concentration from MCMBoxModel.");
  return params;
}

MCMRO2Postprocessor::MCMRO2Postprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _ro2_value(0.0)
{
  const auto & var_names = getParam<std::vector<VariableName>>("species_variables");
  if (!var_names.empty())
  {
    _species_vars.reserve(var_names.size());
    for (const auto & name : var_names)
    {
      if (_fe_problem.hasScalarVariable(name))
        _species_vars.push_back(&_fe_problem.getScalarVariable(0, name));
      else
        paramError("species_variables", "ScalarVariable '", name, "' not found");
    }
  }
}

void
MCMRO2Postprocessor::execute()
{
  if (_species_vars.empty())
  {
    _ro2_value = 0.0;
    return;
  }
  const unsigned int n = _box_model.nSpecies();
  std::vector<Real> C(std::min(n, (unsigned int)_species_vars.size()), 0.0);
  for (unsigned int i = 0; i < C.size(); ++i)
    C[i] = _species_vars[i]->sln()[0];
  _ro2_value = _box_model.getRO2Sum(C);
}

Real
MCMRO2Postprocessor::getValue() const
{
  return _ro2_value;
}
