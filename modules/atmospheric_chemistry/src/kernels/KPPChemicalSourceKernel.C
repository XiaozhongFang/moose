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

  params.addRequiredParam<unsigned int>("species_index",
                                        "Row index of this species in KPP VAR order.");
  params.addRequiredCoupledVar("all_species", "All chemical species in KPP VAR order.");
  params.addParam<Real>("unit_conversion", 1.0,
                        "Variable-unit to molec/cm^3 conversion factor. Use M/1e9 for ppb.");

  params.addClassDescription(
      "KPP chemical source term with KPP-generated analytical Jacobian.");
  return params;
}

KPPChemicalSourceKernel::KPPChemicalSourceKernel(const InputParameters & params)
  : Kernel(params),
    _species_index(getParam<unsigned int>("species_index")),
    _unit_conversion(getParam<Real>("unit_conversion")),
    _kpp_rhs(getMaterialProperty<std::vector<Real>>("kpp_rhs")),
    _kpp_jacobian_dense(getMaterialProperty<std::vector<Real>>("kpp_jacobian_dense")),
    _n_species(coupledComponents("all_species")),
    _coupled_vars(coupledIndices("all_species"))
{
  if (_species_index >= _n_species)
    mooseError("KPPChemicalSourceKernel: species_index ", _species_index,
               " is outside all_species size ", _n_species);

  for (const auto i : make_range(_n_species))
    _coupled_var_to_idx[_coupled_vars[i]] = i;
}

Real
KPPChemicalSourceKernel::computeQpResidual()
{
  const auto & rhs = _kpp_rhs[_qp];
  if (rhs.size() <= _species_index)
    mooseError("KPPChemicalSourceKernel: kpp_rhs has size ", rhs.size(),
               " but species_index is ", _species_index);

  return -_test[_i][_qp] * rhs[_species_index] / _unit_conversion;
}

Real
KPPChemicalSourceKernel::computeQpJacobian()
{
  const auto & jac = _kpp_jacobian_dense[_qp];
  const auto required = static_cast<std::size_t>(_n_species) * _n_species;
  if (jac.size() < required)
    mooseError("KPPChemicalSourceKernel: kpp_jacobian_dense has size ", jac.size(),
               " but expected at least ", required);

  return -_test[_i][_qp] * _phi[_j][_qp] *
         jac[_species_index * _n_species + _species_index];
}

Real
KPPChemicalSourceKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  const auto it = _coupled_var_to_idx.find(jvar);
  if (it == _coupled_var_to_idx.end())
    return 0.0;

  const auto & jac = _kpp_jacobian_dense[_qp];
  const auto required = static_cast<std::size_t>(_n_species) * _n_species;
  if (jac.size() < required)
    mooseError("KPPChemicalSourceKernel: kpp_jacobian_dense has size ", jac.size(),
               " but expected at least ", required);

  const auto col = it->second;
  return -_test[_i][_qp] * _phi[_j][_qp] *
         jac[_species_index * _n_species + col];
}
