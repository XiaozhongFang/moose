//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemistryODEKernel.h"
#include "MCMBoxModel.h"
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
      "ODE kernel that computes chemical source terms via BoxIntegrator. "
      "One instance per species in box mode (mode=box). "
      "Delegates all residual/Jacobian evaluation to the integrator strategy.");
  return params;
}

ChemistryODEKernel::ChemistryODEKernel(const InputParameters & params)
  : ODEKernel(params),
    _integrator(getUserObject<MCMBoxModel>("box_model").getIntegrator()),
    _species_idx(getParam<unsigned int>("species_index")),
    _ppb_to_molec(getUserObject<MCMBoxModel>("box_model").ppbToMolec())
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
  // Delegate to the integrator strategy.
  //   MooseImplicitIntegrator: invalidates caches, updates time for photolysis.
  //   PetscTSIntegrator: no-op (self-driven mode, integration handled separately).
  _integrator.reinit(_t);
}

const std::vector<Real> &
ChemistryODEKernel::_buildC() const
{
  _C_buffer.assign(_species_vars.size(), 0.0);
  if (_ppb_to_molec != 1.0)
  {
    // User variable stores ppb → convert to molec/cm³ for internal computation
    for (unsigned int i = 0; i < _species_vars.size(); ++i)
      _C_buffer[i] = _species_vars[i]->sln()[0] * _ppb_to_molec;
  }
  else
  {
    // User variable is already in molec/cm³ — no conversion
    for (unsigned int i = 0; i < _species_vars.size(); ++i)
      _C_buffer[i] = _species_vars[i]->sln()[0];
  }
  return _C_buffer;
}

Real
ChemistryODEKernel::computeQpResidual()
{
  // R = -dC/dt  (ODEKernel convention: R = du/dt - f, chemical source is f)
  // In self-driven (PETSc TS) mode, the integrator returns 0.0.
  Real dCdt = _integrator.computeResidual(_species_idx, _buildC());
  // If user variable stores ppb, convert dC/dt from molec/cm³/s back to ppb/s
  if (_ppb_to_molec != 1.0)
    dCdt /= _ppb_to_molec;
  return -dCdt;
}

Real
ChemistryODEKernel::computeQpJacobian()
{
  // In self-driven (PETSc TS) mode, the integrator returns 0.0.
  return -_integrator.computeJacobianDiagonal(_species_idx, _buildC());
}

Real
ChemistryODEKernel::computeQpOffDiagJacobianScalar(unsigned int jvar)
{
  // In self-driven (PETSc TS) mode, the integrator returns 0.0.
  return -_integrator.computeJacobianOffDiagonal(_species_idx, jvar, _buildC());
}
