//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KPPChemicalSourceKernel.h"

registerMooseObject("AtmosphericChemistryApp", KPPChemicalSourceKernel);

InputParameters
KPPChemicalSourceKernel::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<std::string>("kpp_library",
      "Path to the KPP shared library (.so) containing the mechanism.");

  params.addRequiredCoupledVar("all_species", "All chemical species variables");

  params.addRequiredParam<unsigned int>("species_index",
      "Index of this species in the KPP mechanism's species ordering.");

  params.addClassDescription(
      "Chemical source kernel backed by a KPP-generated mechanism. "
      "Uses KPP Fun() for RHS and Jac_SP() for analytical Jacobian.");
  return params;
}

KPPChemicalSourceKernel::KPPChemicalSourceKernel(const InputParameters & params)
  : Kernel(params),
    _species_idx(getParam<unsigned int>("species_index")),
    _n_species(coupledComponents("all_species")),
    _all_coupled_vars(coupledIndices("all_species")),
    _all_coupled_vals(coupledValues("all_species")),
    _last_qp(std::numeric_limits<unsigned int>::max())
{
  // Load KPP mechanism shared library
  std::string lib_path = getParam<std::string>("kpp_library");
  _mechanism = std::make_shared<KPPGeneratedMechanism>(lib_path);

  if (_n_species != _mechanism->nSpecies())
    paramError("all_species",
               "Expected ", _mechanism->nSpecies(), " species from KPP library '",
               lib_path, "', but got ", _n_species, " coupled variables.");

  // Build variable number -> species index map for off-diagonal Jacobian
  _var_to_species.reserve(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _var_to_species[_all_coupled_vars[i]] = i;

  // Allocate cache vectors
  _C_cache.resize(_n_species);
  _dC_dt_cache.resize(_n_species);
}

Real
KPPChemicalSourceKernel::computeQpResidual()
{
  // Build concentration vector at this qp
  Real u_val = _u[_qp];
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    if (_all_coupled_vals[i])
      _C_cache[i] = (*_all_coupled_vals[i])[_qp];
    else
      _C_cache[i] = u_val;
  }

  // Compute dC/dt via KPP Fun
  PhysParams dummy;
  _mechanism->computeRHS(0.0, _C_cache, dummy, _dC_dt_cache);

  _last_qp = _qp;

  // Residual: -dC[this_species]/dt
  return -_test[_i][_qp] * _dC_dt_cache[_species_idx];
}

Real
KPPChemicalSourceKernel::computeQpJacobian()
{
  // Build concentration vector
  Real u_val = _u[_qp];
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    if (_all_coupled_vals[i])
      _C_cache[i] = (*_all_coupled_vals[i])[_qp];
    else
      _C_cache[i] = u_val;
  }

  // Compute sparse Jacobian via KPP Jac_SP
  PhysParams dummy;
  _J_cache.clear();
  _mechanism->computeJacobian(0.0, _C_cache, dummy, _J_cache);

  // Extract diagonal element: d(dC_i/dt)/dC_i
  for (const auto & entry : _J_cache)
  {
    if (std::get<0>(entry) == _species_idx && std::get<1>(entry) == _species_idx)
      return -_test[_i][_qp] * _phi[_j][_qp] * std::get<2>(entry);
  }

  return 0.0;
}

Real
KPPChemicalSourceKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  auto it = _var_to_species.find(jvar);
  if (it == _var_to_species.end())
    return 0.0;

  unsigned int j_species = it->second;
  if (j_species == _species_idx)
    return 0.0;  // handled by computeQpJacobian

  // Compute Jacobian if not already cached
  if (_J_cache.empty())
  {
    Real u_val = _u[_qp];
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      if (_all_coupled_vals[i])
        _C_cache[i] = (*_all_coupled_vals[i])[_qp];
      else
        _C_cache[i] = u_val;
    }

    PhysParams dummy;
    _mechanism->computeJacobian(0.0, _C_cache, dummy, _J_cache);
  }

  // Extract off-diagonal element: d(dC_i/dt)/dC_j
  for (const auto & entry : _J_cache)
  {
    if (std::get<0>(entry) == _species_idx && std::get<1>(entry) == j_species)
      return -_test[_i][_qp] * _phi[_j][_qp] * std::get<2>(entry);
  }

  return 0.0;
}
