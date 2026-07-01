//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMLifetimePostprocessor.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", MCMLifetimePostprocessor);

InputParameters
MCMLifetimePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("box_model", "MCMBoxModel UserObject");
  params.addRequiredParam<unsigned int>("species_index", "0-based species index");
  params.addRequiredParam<std::vector<VariableName>>(
      "species_variables", "All species ScalarVariable names (in order)");
  params.addClassDescription(
      "Chemical lifetime of a species: tau = C / loss_rate [s]. "
      "Returns 1e12 s if loss rate is near zero (F0AM lifetime.m equivalent).");
  return params;
}

MCMLifetimePostprocessor::MCMLifetimePostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _species_index(getParam<unsigned int>("species_index")),
    _species_vars(getParam<std::vector<VariableName>>("species_variables"))
{
}

const std::vector<Real> &
MCMLifetimePostprocessor::_buildC() const
{
  _C_buffer.assign(_species_vars.size(), 0.0);
  for (unsigned int i = 0; i < _species_vars.size(); ++i)
  {
    if (_sc_fe_problem.hasScalarVariable(_species_vars[i]))
      _C_buffer[i] = _sc_fe_problem.getScalarVariable(_sc_tid, _species_vars[i]).sln()[0];
  }
  if (_box_model.unitsPPB())
  {
    Real conv = _box_model.ppbToMolec();
    for (auto & v : _C_buffer) v *= conv;
  }
  return _C_buffer;
}

PostprocessorValue
MCMLifetimePostprocessor::getValue() const
{
  const auto & C = _buildC();
  Real conc = (_species_index < C.size()) ? C[_species_index] : 0.0;
  Real loss = _box_model.speciesLossRate(_species_index, C);
  if (loss > 1.0e-30)
    return conc / loss;
  return 1.0e12; // sentinel: effectively infinite lifetime
}
