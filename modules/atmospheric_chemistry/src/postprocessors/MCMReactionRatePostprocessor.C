//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMReactionRatePostprocessor.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", MCMReactionRatePostprocessor);

InputParameters
MCMReactionRatePostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("box_model", "MCMBoxModel UserObject");
  params.addRequiredParam<unsigned int>("reaction_index", "0-based reaction index");
  params.addRequiredParam<std::vector<VariableName>>(
      "species_variables", "All species ScalarVariable names (in order)");
  params.addClassDescription(
      "Outputs the rate of a single chemical reaction: k[r] * prod(C[reactants]) "
      "in molec/cm³/s or ppb/s.");
  return params;
}

MCMReactionRatePostprocessor::MCMReactionRatePostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _reaction_index(getParam<unsigned int>("reaction_index")),
    _species_vars(getParam<std::vector<VariableName>>("species_variables"))
{
}

const std::vector<Real> &
MCMReactionRatePostprocessor::_buildC() const
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
MCMReactionRatePostprocessor::getValue() const
{
  return _box_model.reactionRate(_reaction_index, _buildC());
}
