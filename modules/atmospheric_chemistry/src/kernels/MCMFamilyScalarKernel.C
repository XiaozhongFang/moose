//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMFamilyScalarKernel.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", MCMFamilyScalarKernel);

InputParameters
MCMFamilyScalarKernel::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addRequiredParam<UserObjectName>("family_uo",
                                          "Name of the MCMFamilyConstraint UserObject");
  params.addRequiredParam<unsigned int>("species_index",
                                        "Index of this species in the mechanism");
  params.addParam<std::string>("family_name", "",
                               "Name of the family this kernel enforces (empty = no family). "
                               "If non-empty and species_index is the slack variable, "
                               "DAE constraint is applied.");
  params.addRequiredParam<std::vector<VariableName>>(
      "species_variables",
      "Names of all species ScalarVariables (in species order)");
  params.addClassDescription(
      "ScalarKernel for F0AM-style family conservation via DAE. "
      "Slack variable residual is corrected to enforce d(F_total)/dt = 0. "
      "Non-slack members pass through unchanged (normal ODE evolution).");
  return params;
}

MCMFamilyScalarKernel::MCMFamilyScalarKernel(const InputParameters & params)
  : ODEKernel(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _family_uo(getUserObject<MCMFamilyConstraint>("family_uo")),
    _species_idx(getParam<unsigned int>("species_index")),
    _family_name(getParam<std::string>("family_name")),
    _cached(false)
{
  const auto & var_names = getParam<std::vector<VariableName>>("species_variables");
  _species_vars.reserve(var_names.size());
  for (const auto & name : var_names)
  {
    if (_sc_fe_problem.hasScalarVariable(name))
      _species_vars.push_back(&_sc_fe_problem.getScalarVariable(_sc_tid, name));
    else
      mooseError("MCMFamilyScalarKernel: ScalarVariable '", name, "' not found");
  }
}

void
MCMFamilyScalarKernel::reinit()
{
  _cached = false;
  _box_model.markDirty();
  _box_model.setCurrentTime(_t);
}

const std::vector<Real> &
MCMFamilyScalarKernel::_buildC() const
{
  _C_buffer.assign(_species_vars.size(), 0.0);
  if (_box_model.unitsPPB())
  {
    Real conv = _box_model.ppbToMolec();
    for (unsigned int i = 0; i < _species_vars.size(); ++i)
      _C_buffer[i] = _species_vars[i]->sln()[0] * conv;
  }
  else
  {
    for (unsigned int i = 0; i < _species_vars.size(); ++i)
      _C_buffer[i] = _species_vars[i]->sln()[0];
  }
  return _C_buffer;
}

Real
MCMFamilyScalarKernel::computeQpResidual()
{
  const std::vector<Real> & C = _buildC();
  unsigned int nSp = _species_vars.size();

  // Compute full dC/dt for ALL species (cached for efficiency)
  if (!_cached)
  {
    _cached_dCdt_all.assign(nSp, 0.0);
    _box_model.computeDCdt(C, _cached_dCdt_all);

    // If this is a family slack variable, correct the source term
    // to enforce d(F_total)/dt = 0.
    // Formula: dC_slack/dt_corrected = sum(scaling_i * dC_i/dt)
    //   which equals 0 for the weighted sum of the family.
    // Implementation: R = du/dt - sum(scaling * dCdt_all[members])
    //   where the slack's own original dC/dt is excluded and replaced
    //   by the negative sum of other members' contributions.
    if (!_family_name.empty() && _family_uo.isSlack(_species_idx))
    {
      const auto & members = _family_uo.memberIndices(_family_name);
      const auto & scaling = _family_uo.scalingFactors(_family_name);
      Real total_family_source = 0.0;
      for (unsigned int k = 0; k < members.size(); ++k)
        total_family_source += scaling[k] * _cached_dCdt_all[members[k]];

      // Slack variable absorbs the algebraic constraint:
      // its effective dC/dt = 0 (total family source is zeroed)
      // Actually: the correction replaces dC_slack/dt so the
      // family total is conserved.
      // dC_slack/dt_corrected = - sum_{i != slack} scaling_i * dC_i/dt / scaling_slack
      Real slack_scaling = (!scaling.empty() && _species_idx < members.size())
                               ? scaling[_species_idx]
                               : 1.0;
      _cached_dCdt_all[_species_idx] =
          -(total_family_source - scaling[_species_idx] * _cached_dCdt_all[_species_idx]) /
          slack_scaling;
    }

    if (_box_model.unitsPPB())
    {
      Real conv = _box_model.ppbToMolec();
      for (auto & v : _cached_dCdt_all)
        v /= conv;
    }
    _cached = true;
  }

  return -_cached_dCdt_all[_species_idx];
}

Real
MCMFamilyScalarKernel::computeQpJacobian()
{
  return -_box_model.getJacobianDiagonal(_species_idx, _buildC());
}

Real
MCMFamilyScalarKernel::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  return -_box_model.getJacobianOffDiagonal(_species_idx, jvar, _buildC());
}
