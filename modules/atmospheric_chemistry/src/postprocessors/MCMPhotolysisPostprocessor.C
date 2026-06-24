//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMPhotolysisPostprocessor.h"

registerMooseObject("AtmosphericChemistryApp", MCMPhotolysisPostprocessor);

InputParameters
MCMPhotolysisPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addParam<bool>("output_all", false,
                        "If true, output all J values (ignores j_number)");
  params.addParam<unsigned int>("j_number", 1,
                                "1-based J number to output (e.g., 1 for J1)");
  params.addClassDescription(
      "Outputs photolysis rate J-values from MCMBoxModel or SZA calculation.");
  return params;
}

MCMPhotolysisPostprocessor::MCMPhotolysisPostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _output_all(getParam<bool>("output_all")),
    _j_number(getParam<unsigned int>("j_number")),
    _j_value(0.0),
    _cosx(0.0)
{
}

void
MCMPhotolysisPostprocessor::execute()
{
  _j_value = _box_model.getJValue(_j_number);
}

Real
MCMPhotolysisPostprocessor::getValue() const
{
  return _j_value;
}
