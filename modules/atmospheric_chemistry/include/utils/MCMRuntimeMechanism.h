//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IMechanism.h"
#include "FunctionParserUtils.h"
#include "StoichMatrix.h"
#include "JCalibrator.h"
#include "HybridJTableReader.h"
#include "BottomUpJIntegrator.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <array>
#include <memory>
#include <functional>

// Forward declarations
struct ParsedMechanism;
class HybridJTableReader;
class BottomUpJIntegrator;
class JCalibrator;

/**
 * Runtime chemical mechanism using MCM Facsimile (.fac) parsing.
 *
 * Parses .fac files at construction time and evaluates RHS/Jacobian
 * via runtime StoichMatrix × fparser rate coefficients.
 * This is the default mechanism backend (used when mechanism_format=MCM_FACSIMILE).
 *
 * Implements the IMechanism interface for MOOSE-agnostic chemical computation.
 * The core dC/dt evaluation is:
 *   dC/dt = StoichMatrix^T * (k .* prod(C[iG], 2))
 *
 * This class is MOOSE-agnostic except for its use of FunctionParserUtils
 * (which is a MOOSE utility header but only depends on Real/string parsing).
 * It does NOT inherit from GeneralUserObject — it is a standalone utility
 * usable from MCMBoxModel, KPP adapters, or standalone test code.
 */
class MCMRuntimeMechanism : public IMechanism, public FunctionParserUtils<false>
{
public:
  /// Type alias for fast pre-compiled handler
  using FastHandler = std::function<Real(const std::vector<Real> &)>;

  /**
   * Construct from a fully parsed mechanism.
   *
   * @param mech                Parsed mechanism (from MCMFacsimileParser)
   * @param use_limiting_reagent Enable F0AM-style limiting-reagent formulation
   *                              for RO2+RO2 termination reactions (default false)
   * @param stoich_format       Stoichiometric matrix format (default CSR)
   */
  MCMRuntimeMechanism(const ParsedMechanism & mech,
                       bool use_limiting_reagent = false,
                       StoichMatrix::Format stoich_format = StoichMatrix::CSR,
                       bool disable_fpoptimizer = false,
                       bool enable_jit = false,
                       bool fail_on_bad_deps = true,
                       unsigned int max_function_recurse = 3);

  ~MCMRuntimeMechanism() override = default;

  // ===== IMechanism interface =====

  unsigned int nSpecies() const override { return _n_species; }
  unsigned int nReactions() const override { return _n_reactions; }
  const std::vector<std::string> & speciesNames() const override { return _species_names; }

  void updateParams(const PhysParams & params) override;

  void computeRHS(Real t,
                   const std::vector<Real> & C,
                   const PhysParams & params,
                   std::vector<Real> & dC_dt) const override;

  void computeJacobian(Real t,
                        const std::vector<Real> & C,
                        const PhysParams & params,
                        std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const override;

  SpeciesRates computeSpeciesRates(Real t,
                                    const std::vector<Real> & C,
                                    const PhysParams & params) const override;

  Real ppbToMolec() const override { return _ppb_to_molec; }
  bool unitsPPB() const override { return _units_ppb; }

  // ===== Kernel-facing single-species interface (same API as MCMBoxModel) =====

  /** Mark the internal cache as dirty. Called before each residual/Jacobian eval. */
  void markDirty() const { _dirty = true; }

  /** Get dC/dt for a single species, caching the full computation. */
  Real getDCdt(unsigned int idx, const std::vector<Real> & C) const;

  /** Get diagonal Jacobian element d(dC_i)/dC_i, caching the full Jacobian. */
  Real getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const;

  /** Get off-diagonal Jacobian element d(dC_i)/dC_j, caching the full Jacobian. */
  Real getJacobianOffDiagonal(unsigned int i, unsigned int j,
                               const std::vector<Real> & C) const;

  // ===== Detailed rate diagnostics =====

  /** Compute rate (k * prod(reactants)) for a single reaction. */
  Real reactionRate(unsigned int r, const std::vector<Real> & C) const;

  /** Compute net rate contribution of reaction r to species s. */
  Real speciesReactionRate(unsigned int s, unsigned int r,
                            const std::vector<Real> & C) const;

  /** Compute total loss rate for species s (sum of consuming reactions). */
  Real speciesLossRate(unsigned int s, const std::vector<Real> & C) const;

  /** Compute total production rate for species s (sum of producing reactions). */
  Real speciesProductionRate(unsigned int s, const std::vector<Real> & C) const;

  /** Compute all reaction rates as a vector. */
  void allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const;

  /** Compute RO2 sum (peroxy radical total) from concentration vector. */
  Real getRO2Sum(const std::vector<Real> & C) const;

  // ===== Accessors =====

  const StoichMatrix & stoichiometry() const { return _stoich; }
  const std::vector<unsigned int> & ro2Indices() const { return _ro2_indices; }
  const std::vector<std::string> & ro2Species() const { return _ro2_species_names; }

  // ===== Solar / photolysis helpers (needed by evaluateCoefficients) =====

  /** Set current simulation time (seconds since midnight) for photolysis calc. */
  void setCurrentTime(Real t) const { _t = t; }

  /** Get cached solar zenith angle cosine. */
  Real getSolarCosX() const { return _solar_cosx; }
  /** Get cached solar secant (1/cosx). */
  Real getSolarSecX() const { return _solar_secx; }
  /** Get cached local hour angle (radians). */
  Real getSolarLHA() const { return _solar_lha; }
  /** Get sin(lat)*sin(dec). */
  Real getSolarSinLD() const { return _solar_sinld; }
  /** Get cos(lat)*cos(dec). */
  Real getSolarCosLD() const { return _solar_cosld; }
  /** Get equation of time. */
  Real getSolarEQT() const { return _solar_eqt; }
  /** Get solar declination (radians). */
  Real getDeclination() const { return _declination; }

  /** Set solar parameters (replaces MCMBoxModel::setSolarCycle). */
  void setSolarParams(Real lat, Real lon, int day, int month, int year);

  /** Compute cos(SZA) at given seconds since midnight (Madronich 1993). */
  Real cosSZA(Real seconds) const;

  /** Enable Hybrid J (TUV lookup) photolysis. */
  void enableHybridPhotolysis(const std::string & table_dir);

  /** Load BottomUp lamp flux and reaction map. */
  void loadBottomUpData(const std::string & data_dir, const std::string & flux_file);

  /** Set JFAC scaling factor. */
  void setJFac(Real jfac) { _jfac = jfac; }

  /** ROOF switch: false = CLOSED (all J=0), true = OPEN (normal). */
  void setRoofOpen(bool open) { _roof_open = open; }
  bool isRoofOpen() const { return _roof_open; }

  /** Invalidate BottomUp J cache. */
  void invalidatePhotolysisCache() { _bottomup_j_valid = false; }

  /** Enable/disable J-calibrator. */
  void setJCalibrator(std::unique_ptr<JCalibrator> calibrator);

  /** Get a single photolysis J value by 1-based J number. */
  Real getJValue(unsigned int j_number) const;

  /** Get the number of photolysis J variables. */
  unsigned int nJValues() const { return _n_j_vars; }

  /** Compute day of year (0-based, AtChem2 convention). */
  unsigned int computeDayOfYear() const;

  /** Calculate cos(SZA) with caching (internal). */
  Real calculateCosSZA(Real t) const;

  // ===== Public computation methods (needed by MCMBoxModel) =====

  /** Evaluate all rate coefficients (fparser evaluation + photolysis J).
   *  Called automatically by computeRHS/getDCdt when dirty, but also
   *  callable directly when coefficients need to be pre-computed before
   *  a series of evaluations (e.g., at timestep midpoint in PETSc TS mode). */
  void evaluateCoefficients();

  /** Compute dC/dt using the current (already-evaluated) rate coefficients.
   *  Unlike computeRHS(), this does NOT re-evaluate rate coefficients —
   *  it uses the cached _k from the most recent evaluateCoefficients() call.
   *  This is the hot-path for PETSc TS internal steps where coefficients
   *  are constant across the timestep. */
  void computeDCdt(const std::vector<Real> & C, std::vector<Real> & dC) const;

  /** Compute Jacobian as (row, col, val) triplets using current coefficients. */
  void computeJacobianTriplets(
      const std::vector<Real> & C,
      std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const;

  /** Build the sparse Jacobian cache from _cached_C.
   *  Called automatically by getJacobianDiagonal/getJacobianOffDiagonal. */
  void _buildJacobianCache() const;

private:
  // ===== Core computation methods =====

  /** Set up fparser for coefficient and reaction expressions. */
  void setupFparser(const ParsedMechanism & mech);

  /** Compile a fast pre-compiled handler (bypass fparser). */
  FastHandler compileFastHandler(const std::string & expr,
                                   const std::vector<unsigned int> & var_indices) const;

  /** Load mechanism data from ParsedMechanism. */
  void loadMechanism(const ParsedMechanism & mech,
                      bool use_limiting_reagent,
                      StoichMatrix::Format stoich_format);

  // ===== Member variables =====

  // -- Mechanism data --
  unsigned int _n_species;
  unsigned int _n_reactions;
  std::vector<std::string> _species_names;
  std::vector<std::string> _reaction_names;
  StoichMatrix _stoich;
  std::vector<std::array<int, 3>> _iG;
  std::vector<Real> _k;

  // -- Units --
  bool _units_ppb;
  Real _ppb_to_molec;

  // -- RO2 --
  std::vector<unsigned int> _ro2_indices;
  std::vector<std::string> _ro2_species_names;

  // -- Limiting reagent --
  bool _use_limiting_reagent;
  std::vector<bool> _limiting_reagent;
  std::vector<int> _limiting_reactant;

  // -- Physical conditions --
  Real _temperature;
  Real _air_density;
  Real _water_vapor;
  Real _press;
  Real _rh;
  Real _blheight;
  Real _jfac;

  // -- Photolysis method --
  enum PhotolysisMethod { MCM_SZA, HYBRID, BOTTOMUP } _photolysis_method;
  std::unique_ptr<HybridJTableReader> _hybrid_reader;
  std::unique_ptr<BottomUpJIntegrator> _bottomup_integrator;
  std::vector<Real> _j_CL_vals, _j_CMM_vals, _j_CNN_vals;
  std::vector<unsigned int> _j_numbers;
  std::vector<unsigned int> _j_photo_indices;
  unsigned int _j_index_start;
  unsigned int _n_j_vars;
  bool _roof_open;
  std::unique_ptr<JCalibrator> _jcalibrator;

  // -- Solar cycle --
  Real _lat, _lon;
  mutable Real _declination, _eot;
  int _day, _month, _year;

  // -- Cached solar parameters --
  mutable Real _solar_cosx, _solar_secx, _solar_lha;
  mutable Real _solar_sinld, _solar_cosld, _solar_eqt;

  // -- Current time (seconds since midnight) --
  mutable Real _t;

  // -- fparser --
  using FunctionParserUtils<false>::_func_params;
  using FunctionParserUtils<false>::evaluate;
  std::unordered_map<std::string, unsigned int> _name_to_index;
  std::vector<SymFunctionPtr> _coeff_parsers;
  std::vector<SymFunctionPtr> _reaction_parsers;
  std::vector<std::vector<unsigned int>> _coeff_var_indices;
  std::vector<std::vector<Real>> _coeff_local_params;
  std::vector<std::vector<unsigned int>> _reaction_var_indices;
  std::vector<std::vector<Real>> _reaction_local_params;
  std::vector<FastHandler> _coeff_fast;
  std::vector<FastHandler> _reaction_fast;

  // -- Cache for single-species interface --
  mutable bool _dirty;
  mutable std::vector<Real> _cached_dC;
  mutable std::vector<Real> _cached_C;
  mutable std::vector<Real> _cached_diag_J;
  mutable std::vector<unsigned int> _cached_od_cols;
  mutable std::vector<Real> _cached_od_vals;
  mutable std::vector<size_t> _cached_od_row_ptr;

  // -- Scratch buffers --
  mutable std::vector<Real> _scratch_G;
  mutable std::vector<Real> _scratch_rates;

  // -- BottomUp cache --
  mutable std::map<std::string, Real> _cached_bottomup_j;
  mutable Real _cached_bottomup_T;
  mutable Real _cached_bottomup_P;
  mutable bool _bottomup_j_valid;
};
