//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemistryODEKernel.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"

registerMooseObject("AtmosphericChemistryApp", ChemistryODEKernel);

InputParameters
ChemistryODEKernel::validParams()
{
  InputParameters params = ODEKernel::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addRequiredParam<unsigned int>("species_index",
                                        "Index of this species in the mechanism (0-indexed)");
  params.addRequiredParam<std::vector<VariableName>>(
      "species_variables", "Names of all species ScalarVariables (in species order)");
  params.addClassDescription(
      "ODE kernel that computes chemical source terms from MCMBoxModel. "
      "One instance per species in box mode (mode=box).");
  return params;
}

ChemistryODEKernel::ChemistryODEKernel(const InputParameters & params)
  : ODEKernel(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _species_idx(getParam<unsigned int>("species_index"))
{
  // Access scalar variables directly from the FEProblemBase, bypassing the
  // ScalarCoupleable framework since species names are dynamic (not known
  // at compile time and thus cannot be declared in validParams()).
  // _sc_fe_problem is a protected member of ScalarCoupleable (FEProblemBase &).
  const auto & var_names = getParam<std::vector<VariableName>>("species_variables");
  _species_vars.reserve(var_names.size());
  for (const auto & name : var_names)
  {
    if (_sc_fe_problem.hasScalarVariable(name))
      _species_vars.push_back(&_sc_fe_problem.getScalarVariable(_sc_tid, name));
    else
      mooseError("ChemistryODEKernel: ScalarVariable '", name, "' not found");
  }
}

void
ChemistryODEKernel::reinit()
{
  // Invalidate the BoxModel cache so the next getDCdt / getJacobian*
  // call triggers a fresh full-system computation.
  _box_model.markDirty();
  _box_model.setCurrentTime(_t);
}

std::vector<Real>
ChemistryODEKernel::_buildC() const
{
  std::vector<Real> C(_species_vars.size());
  for (unsigned int i = 0; i < _species_vars.size(); ++i)
    C[i] = _species_vars[i]->sln()[0];
  return C;
}

Real
ChemistryODEKernel::computeQpResidual()
{
  // R = -dC/dt  (ODEKernel convention: R = du/dt - f, chemical source is f)
  return -_box_model.getDCdt(_species_idx, _buildC());
}

Real
ChemistryODEKernel::computeQpJacobian()
{
  return -_box_model.getJacobianDiagonal(_species_idx, _buildC());
}

Real
ChemistryODEKernel::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  return -_box_model.getJacobianOffDiagonal(_species_idx, jvar, _buildC());
}
