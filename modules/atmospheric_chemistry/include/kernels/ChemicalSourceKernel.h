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

/**
 * Kernel that implements the chemical source term for species j:
 *
 *   residual -= _test * Σ_i S_{j,i} × R_i
 *
 * with analytical Jacobian following the pattern in CoupledBEEquilibriumSub.
 */
class ChemicalSourceKernel : public Kernel
{
public:
  static InputParameters validParams();

  ChemicalSourceKernel(const InputParameters & params);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
  /// Stoichiometric row S_j for this species (one entry per reaction)
  const std::vector<Real> _stoichiometric_row;

  /// Reaction rates R_i
  const MaterialProperty<std::vector<Real>> & _reaction_rates;

  /// Number of coupled species
  const unsigned int _n_species;

  /// Variable indices of all species
  std::vector<unsigned int> _coupled_vars;

  /// Concentration values of all species
  std::vector<const VariableValue *> _coupled_vals;

  /// CSR-format reactant matrix: for species k, reactions span
  /// _sr_row_ptr[k] .. _sr_row_ptr[k+1]-1 with (reaction, coeff) pairs.
  /// Built from the input parameter in the constructor — eliminates
  /// vector<vector<Real>> overhead and static_cast<unsigned int> on indices.
  std::vector<unsigned int> _sr_cols;
  std::vector<Real>        _sr_vals;
  std::vector<size_t>      _sr_row_ptr;

  /// O(1) reverse index: jvar → species index (Per.14 — avoids O(nSpecies) scan)
  std::unordered_map<unsigned int, unsigned int> _coupled_var_to_idx;
};
