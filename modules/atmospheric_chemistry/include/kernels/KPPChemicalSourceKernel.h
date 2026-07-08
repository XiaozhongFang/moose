//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

#include <unordered_map>

/**
 * Chemical source kernel backed by KPP-generated RHS and analytical Jacobian
 * material properties.
 */
class KPPChemicalSourceKernel : public Kernel
{
public:
  static InputParameters validParams();

  KPPChemicalSourceKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  const unsigned int _species_index;
  const Real _unit_conversion;

  const MaterialProperty<std::vector<Real>> & _kpp_rhs;
  const MaterialProperty<std::vector<Real>> & _kpp_jacobian_dense;

  const unsigned int _n_species;
  std::vector<unsigned int> _coupled_vars;
  std::unordered_map<unsigned int, unsigned int> _coupled_var_to_idx;
};
