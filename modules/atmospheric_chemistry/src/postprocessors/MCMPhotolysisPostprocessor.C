//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMPhotolysisPostprocessor.h"
#include "MCMRatesMaterial.h"

registerMooseObject("AtmosphericChemistryApp", MCMPhotolysisPostprocessor);

InputParameters
MCMPhotolysisPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<UserObjectName>("box_model", "",
      "Name of the MCMBoxModel UserObject (box mode)");
  params.addParam<MaterialPropertyName>("material", "",
      "Name of the MCMRatesMaterial (coupled mode)");
  params.addParam<unsigned int>("j_number", 1,
                                "1-based J number to output (e.g., 1 for J1)");
  params.addClassDescription(
      "Outputs photolysis rate J-values from MCMBoxModel (box mode) or MCMRatesMaterial (coupled mode).");
  return params;
}

MCMPhotolysisPostprocessor::MCMPhotolysisPostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(nullptr),
    _j_material_prop(nullptr),
    _j_number_list_prop(nullptr),
    _j_number(getParam<unsigned int>("j_number")),
    _j_value(0.0),
    _cosx(0.0)
{
  auto box_name = getParam<UserObjectName>("box_model");
  auto mat_name = getParam<MaterialPropertyName>("material");
  if (!box_name.empty())
    _box_model = &getUserObject<MCMBoxModel>("box_model");
  if (!mat_name.empty())
  {
    _j_material_prop = &getMaterialProperty<std::vector<Real>>("photolysis_rates");
    _j_number_list_prop = &getMaterialProperty<std::vector<unsigned int>>("photolysis_j_numbers");
  }
  if (!_box_model && !_j_material_prop)
    mooseError("MCMPhotolysisPostprocessor: either 'box_model' or 'material' must be specified");
}

void
MCMPhotolysisPostprocessor::execute()
{
  if (_box_model)
    _j_value = _box_model->getJValue(_j_number);
  else if (_j_material_prop && _j_number_list_prop)
  {
    // J values in photolysis_rates are sorted by J number (non-contiguous).
    // Use photolysis_j_numbers to map j_number → array index.
    const auto & jn_list = (*_j_number_list_prop)[0];
    unsigned int jn = _j_number;
    _j_value = 0.0;
    for (size_t i = 0; i < jn_list.size(); ++i)
    {
      if (jn_list[i] == jn)
      {
        if (i < (*_j_material_prop)[0].size())
          _j_value = (*_j_material_prop)[0][i];
        break;
      }
    }
  }
}

Real
MCMPhotolysisPostprocessor::getValue() const
{
  return _j_value;
}
