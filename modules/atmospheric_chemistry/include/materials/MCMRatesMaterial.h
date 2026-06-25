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
#include <memory>

class HybridJTableReader;

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

  /// Get a single J value by 1-based J number (e.g. 1→J1).  Returns 0 if
  /// the J number isn't in the mechanism.  For Postprocessor use.
  Real getJValue(unsigned int j_number) const;

  // ── Solar parameter getters (for MCMSolarPostprocessor coupled-mode) ──
  Real getSolarCosX()   const { return _cached_cosx; }
  Real getSolarSecX()   const { return _cached_secx; }
  Real getSolarLHA()    const { return _cached_lha; }
  Real getSolarSinLD()  const { return _cached_sinld; }
  Real getSolarCosLD()  const { return _cached_cosld; }
  Real getSolarEQT()    const { return _cached_eqt; }
  Real getLatitude()    const { return _latitude; }
  Real getLongitude()   const { return _longitude; }
  Real getDeclination() const { return _cached_dec; }

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
  /// Computed photolysis J values (exposed for Postprocessor access)
  MaterialProperty<std::vector<Real>> & _j_values;
  /// J numbers corresponding to _j_values entries (exposed for index mapping)
  MaterialProperty<std::vector<unsigned int>> & _j_number_list;
  /// Solar parameters (exposed for MCMSolarPostprocessor coupled-mode access)
  MaterialProperty<Real> & _solar_cosx;
  MaterialProperty<Real> & _solar_secx;
  MaterialProperty<Real> & _solar_lha;
  MaterialProperty<Real> & _solar_sinld;
  MaterialProperty<Real> & _solar_cosld;
  MaterialProperty<Real> & _solar_eqt;
  MaterialProperty<Real> & _solar_dec;
  MaterialProperty<Real> & _solar_lat;
  MaterialProperty<Real> & _solar_lon;

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

  /// Photolysis scheme: MCM_SZA or HYBRID
  const MooseEnum _photolysis_scheme;
  /// Optional Hybrid J-table reader (created if scheme=HYBRID)
  std::unique_ptr<HybridJTableReader> _hybrid_reader;
  /// Hybrid scheme: surface albedo [0,1]
  const Real _albedo;
  /// Hybrid scheme: O3 column [DU]
  const Real _o3column;
  /// Hybrid scheme: altitude [m]
  const Real _altitude;

  /// Cached solar parameters (computed in computeQpProperties for Postprocessor access)
  Real _cached_cosx, _cached_secx, _cached_lha;
  Real _cached_sinld, _cached_cosld, _cached_eqt, _cached_dec;

  /// Indices of species that constitute RO2 (peroxy radical sum)
  std::vector<unsigned int> _ro2_indices;

  // ── Per-parser variable indirection (ParseAndDeduceVariables) ──
  /// _func_params indices each coefficient parser actually references
  std::vector<std::vector<unsigned int>> _coeff_var_indices;
  /// Pre-allocated local param buffer per coefficient parser
  std::vector<std::vector<Real>> _coeff_local_params;
  /// _func_params indices each reaction parser actually references
  std::vector<std::vector<unsigned int>> _reaction_var_indices;
  /// Pre-allocated local param buffer per reaction parser
  std::vector<std::vector<Real>> _reaction_local_params;
  /// True if RO2 was added as an EXTRA fparser variable (not in species list)
  bool _has_ro2;

  /// Compute day of year from day/month/year
  unsigned int computeDayOfYear() const;
  /// Calculate cosine of solar zenith angle at simulation time t
  Real calculateCosSZA(Real t);

  /// fparser parameter buffer (TEMP, M, O2, N2, H2O, coeff0, coeff1, ...)
  using FunctionParserUtils<false>::_func_params;
};
