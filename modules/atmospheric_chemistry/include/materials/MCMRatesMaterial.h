//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "FunctionParserUtils.h"

/**
 * Material that evaluates MCM reaction rate expressions at each quadrature
 * point using sequential fparser evaluation.
 *
 * Rate coefficients are evaluated in topological order using a shared
 * parameter array that holds base variables (TEMP, M, O2, N2, H2O) + all
 * coefficient values. After each coefficient is evaluated, its value is
 * stored in the array for use by later coefficients.
 *
 * Reaction rates R_i = k_i * Π[C]^ν are then computed using the evaluated
 * coefficient values and coupled species concentrations.
 */
class MCMRatesMaterial : public Material, public FunctionParserUtils<false>
{
public:
  static InputParameters validParams();

  MCMRatesMaterial(const InputParameters & params);

protected:
  virtual void computeQpProperties() override;

private:
  /// Number of chemical species
  unsigned int _n_species;

  /// Coupled species concentration values
  std::vector<const VariableValue *> _species_vals;

  /// Background atmospheric parameters
  const Real _TEMP;
  const Real _M;
  const Real _H2O_val;

  /// Coefficient names and expressions in evaluation order
  std::vector<std::string> _coeff_names;
  unsigned int _n_coefficients;

  /// Number of chemical species (for fparser variable array)
  unsigned int _n_species_material;
  /// Species concentration values at each qp (for fparser reference)
  std::vector<const VariableValue *> _species_vals_material;
  /// Species names -> index in fparser param array (constructor only)
  std::map<std::string, unsigned int> _species_name_to_index;

  /// fparser objects for each coefficient expression
  std::vector<SymFunctionPtr> _coeff_parsers;

  /// Name-to-index mapping for the big parameter array (constructor only)
  std::map<std::string, unsigned int> _name_to_index;

  /// Reaction rate expressions
  unsigned int _n_reactions;
  std::vector<SymFunctionPtr> _reaction_parsers;

  /// Reactant matrix: each row encodes [sp_idx_0, coeff_0, sp_idx_1, coeff_1, ...]
  std::vector<std::vector<Real>> _reactant_matrix;

  /// Computed reaction rates R_i
  MaterialProperty<std::vector<Real>> & _reaction_rates;

  /// J<N> variable support (constructor-computed offsets)
  unsigned int _n_j_variables;
  std::vector<std::string> _j_names;
  unsigned int _j_index_start;
  /// Pre-allocated rate constant buffer (Per.14 — reused across QP calls)
  std::vector<Real> _k_values;

  /// MCM photolysis parameters (CL, CMM, CNN) for SZA-based J calculation
  std::vector<Real> _j_CL;
  std::vector<Real> _j_CMM;
  std::vector<Real> _j_CNN;

  /// Solar zenith angle calculation parameters
  const Real _latitude;
  const Real _longitude;
  const unsigned int _day;
  const unsigned int _month;
  const unsigned int _year;
  const Real _jfac;
  /// Roof (chamber cover) open. false = CLOSED (all photolysis J=0).
  const bool _roof_open;

  /// Compute day of year from day/month/year
  unsigned int computeDayOfYear() const;
  /// Calculate cosine of solar zenith angle at simulation time t
  Real calculateCosSZA(Real t) const;

  /// fparser parameter buffer (TEMP, M, O2, N2, H2O, coeff0, coeff1, ...)
  using FunctionParserUtils<false>::_func_params;
};
