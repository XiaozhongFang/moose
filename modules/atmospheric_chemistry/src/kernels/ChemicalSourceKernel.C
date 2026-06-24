//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalSourceKernel.h"

registerMooseObject("AtmosphericChemistryApp", ChemicalSourceKernel);

InputParameters
ChemicalSourceKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<std::vector<Real>>(
      "stoichiometric_row",
      "Stoichiometric coefficients for this species across all reactions");
  params.addRequiredCoupledVar("all_species", "All chemical species");
  params.addRequiredParam<std::vector<std::vector<Real>>>(
      "species_reactants",
      "species_reactants[k] = [rxn_0, coeff_0, rxn_1, coeff_1, ...]");
  params.addClassDescription("Chemical source term with analytical Jacobian");
  return params;
}

ChemicalSourceKernel::ChemicalSourceKernel(const InputParameters & params)
  : Kernel(params),
    _stoichiometric_row(getParam<std::vector<Real>>("stoichiometric_row")),
    _reaction_rates(getMaterialProperty<std::vector<Real>>("reaction_rates")),
    _n_species(coupledComponents("all_species")),
    _coupled_vars(coupledIndices("all_species")),
    _coupled_vals(coupledValues("all_species"))
{
  // Build CSR reactant matrix from input [rxn_0,coeff_0, rxn_1,coeff_1, ...] per species.
  // Per.16: single contiguous allocation, zero pointer-chasing vs vector<vector<Real>>.
  const auto & sr_in = getParam<std::vector<std::vector<Real>>>("species_reactants");
  _sr_row_ptr.resize(_n_species + 1);
  _sr_row_ptr[0] = 0;
  for (unsigned int k = 0; k < _n_species; ++k)
  {
    const auto & row = sr_in[k];
    for (size_t j = 0; j + 1 < row.size(); j += 2)
    {
      // Validation: row[j] is a reaction index, validate against nReactions
      if (row[j] < 0.0 || row[j] >= (Real)_stoichiometric_row.size())
        mooseError("ChemicalSourceKernel: invalid reaction index ", row[j],
                   " for species ", k, " (nReactions=", _stoichiometric_row.size(), ")");
      _sr_cols.push_back((unsigned int)row[j]);
      _sr_vals.push_back(row[j + 1]);
    }
    _sr_row_ptr[k + 1] = _sr_cols.size();
  }

  // Build O(1) reverse index: jvar → species index (Per.14)
  for (unsigned int k = 0; k < _n_species; ++k)
    _coupled_var_to_idx[_coupled_vars[k]] = k;
}

Real
ChemicalSourceKernel::computeQpResidual()
{
  Real sum = 0.0;
  const auto & rates = _reaction_rates[_qp];
  for (size_t i = 0; i < _stoichiometric_row.size(); ++i)
    sum += _stoichiometric_row[i] * rates[i];
  return -_test[_i][_qp] * sum;
}

Real
ChemicalSourceKernel::computeQpJacobian()
{
  // d(residual)/d(C_self) = -_test * phi_j * Σ_i S_i * dR_i/dC_self
  // For reactions where this species is a reactant:
  //   dR_i/dC_self = ν_self * R_i / C_self
  //   S_i is the stoichiometric coeff of THIS species
  const auto & rates = _reaction_rates[_qp];
  Real C_self = _u[_qp];
  if (C_self <= 0.0)
    return 0.0;

  Real sum = 0.0;
  for (size_t i = 0; i < _stoichiometric_row.size(); ++i)
  {
    Real S = _stoichiometric_row[i];
    if (S >= 0.0) continue; // skip products (dR/dC=0 for non-reactants)
    Real nu = -S;           // ν = -S for a pure reactant
    sum += S * nu * rates[i] / C_self;
  }

  return -_test[_i][_qp] * _phi[_j][_qp] * sum;
}

Real
ChemicalSourceKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  // O(1) lookup via pre-built reverse index (Per.14)
  auto it = _coupled_var_to_idx.find(jvar);
  if (it == _coupled_var_to_idx.end())
    return 0.0;
  unsigned int k_idx = it->second;

  // d(residual_j)/d(C_k) = -_test * phi_j * Σ_i S_{j,i} * dR_i/dC_k
  // dR_i/dC_k = ν_{i,k} * R_i / C_k  (if k is a reactant in reaction i)
  Real C_k = (*_coupled_vals[k_idx])[_qp];
  if (C_k <= 0.0)
    return 0.0;

  const auto & rates = _reaction_rates[_qp];

  Real sum = 0.0;
  for (size_t j = _sr_row_ptr[k_idx]; j < _sr_row_ptr[k_idx + 1]; ++j)
  {
    unsigned int rxn_idx = _sr_cols[j];
    Real nu = _sr_vals[j];
    sum += _stoichiometric_row[rxn_idx] * nu * rates[rxn_idx] / C_k;
  }

  // Following CoupledBEEquilibriumSub pattern: _test * _phi * derivative
  return -_test[_i][_qp] * _phi[_j][_qp] * sum;
}
