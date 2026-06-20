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

  /// _species_reactants[k] = [rxn_0, coeff_0, rxn_1, coeff_1, ...]
  /// Reactions where species k is a reactant
  const std::vector<std::vector<Real>> _species_reactants;
};
