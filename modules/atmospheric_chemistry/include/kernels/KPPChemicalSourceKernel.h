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
#include "KPPGeneratedMechanism.h"

#include <memory>
#include <vector>

/**
 * FE kernel for KPP-backed chemical source terms in coupled mode.
 *
 * For each integration point:
 *   residual -= _test * dC[this_species]/dt
 *
 * where dC/dt = Fun(C) from the KPP-generated mechanism.
 *
 * Jacobian contributions use KPP's analytical sparse Jacobian (Jac_SP).
 *
 * Usage (one kernel per species):
 *   [Kernels]
 *     [chem_O3]
 *       type = KPPChemicalSourceKernel
 *       variable = O3
 *       kpp_library = 'kpp_build_mech/libkpp_mech.so'
 *       all_species = 'O3 NO2 NO OH CO CH4 ...'
 *       species_index = 0
 *     []
 *   []
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
  /// Shared KPP mechanism instance (loaded once, shared across species via App)
  std::shared_ptr<KPPGeneratedMechanism> _mechanism;

  /// Index of this species in the KPP mechanism ordering
  unsigned int _species_idx;

  /// Number of chemical species
  unsigned int _n_species;

  /// All species variable values at the current qp
  std::vector<const VariableValue *> _all_coupled_vals;

  /// Variable numbers for all species (for off-diagonal Jacobian lookup)
  std::vector<unsigned int> _all_coupled_vars;

  /// Map from MOOSE variable number -> species index
  std::unordered_map<unsigned int, unsigned int> _var_to_species;

  /// Cached concentration vector (rebuilt each qp)
  mutable std::vector<Real> _C_cache;

  /// Cached dC/dt vector (rebuilt each qp)
  mutable std::vector<Real> _dC_dt_cache;

  /// Cached Jacobian triplets (rebuilt each qp)
  mutable std::vector<std::tuple<unsigned int, unsigned int, Real>> _J_cache;

  /// Last qp for cache validation
  mutable unsigned int _last_qp;
};
