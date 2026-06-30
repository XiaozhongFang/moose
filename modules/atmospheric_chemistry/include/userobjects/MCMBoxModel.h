//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "MCMFacsimileParser.h"
#include "HybridJTableReader.h"
#include "BottomUpJIntegrator.h"
#include "FunctionParserUtils.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <map>

/**
 * Lightweight sparse stoichiometric matrix with pluggable storage format.
 *
 * Follows PETSc's parameter-driven format selection pattern (cf. -mat_type aij|dense):
 * a MooseEnum parameter "stoich_format" selects the storage backend at construction
 * time.  All compute methods use a single template iteration interface
 * (forEachInRow) that the compiler fully inlines — zero virtual-call overhead
 * regardless of the format chosen.
 *
 * Supported formats:
 *   CSR   — Compressed Sparse Row, net stoichiometric coefficients.
 *           PETSc AIJ-compatible, compact, optimal for HPC with large mechanisms.
 *   COO   — AtChem2-style split reactant / product vectors (clhs/crhs).
 *           Enables separate loss-rate / production-rate diagnostics.
 *   DENSE — Dense 2D array.  Simple, cache-friendly, zero indirect-addressing
 *           overhead.  Best for tiny mechanisms (< ~50 species) where memory is
 *           irrelevant and instruction-level efficiency matters.
 *   CSC   — Compressed Sparse Column: species-major storage.  Each column
 *           lists (reaction, coeff) for all reactions involving that species.
 *           Enables column queries ("which reactions involve species X?").
 *           Row iteration builds a lightweight forward index at build() time.
 *
 * Adding a new format:
 *   1. Add enum value to Format.
 *   2. Add data members.
 *   3. Add branches in forEachInRow(), get(), build().
 */
struct StoichMatrix
{
  enum Format { CSR, COO, DENSE, CSC };

  Format format;

  unsigned int nSpecies = 0, nReactions = 0;

  // ---- CSR data ----
  std::vector<int>    csr_cols;
  std::vector<Real>   csr_vals;
  std::vector<size_t> csr_row_ptr;

  // ---- COO data (AtChem2-style clhs / crhs) ----
  std::vector<int>    lhs_species, rhs_species;
  std::vector<Real>   lhs_coeff,   rhs_coeff;
  std::vector<size_t> lhs_row_ptr, rhs_row_ptr;

  // ---- DENSE data ----
  /// dense[r][s] = net stoichiometric coefficient (0.0 for non-participating)
  std::vector<std::vector<Real>> dense;

  // ---- CSC data (species-major, column queries) ----
  /// Column s spans csc_col_ptr[s] .. csc_col_ptr[s+1]-1.
  /// Element k: reaction = csc_rows[k], coefficient = csc_c_vals[k].
  /// Row iteration uses the CSR fields (populated alongside CSC as a forward
  /// index), so forEachInRow simply delegates to the CSR branch.
  std::vector<int>    csc_rows;
  std::vector<Real>   csc_c_vals;
  std::vector<size_t> csc_col_ptr;

  /// Build from parsed mechanism.
  void build(const ParsedMechanism & mech, Format fmt);

  /**
   * Iterate non-zero stoichiometric entries for reaction r.
   * Calls fn(species_index, net_coefficient) for each participating species.
   *
   * CSR: net_coefficient = product - reactant (stored directly).
   * COO: reactants (-coeff) then products (+coeff), matching AtChem2 resid().
   *       NOTE: if a species appears on BOTH sides of a reaction, fn() is
   *       called twice — once per occurrence.  The net effect sums correctly.
   *
   * The fn lambda is fully inlined — zero dispatch overhead.
   */
  template <typename F>
  void forEachInRow(unsigned int r, F && fn) const
  {
    switch (format)
    {
      case CSR:
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          fn(csr_cols[k], csr_vals[k]);
        break;
      case COO:
        for (size_t k = lhs_row_ptr[r]; k < lhs_row_ptr[r + 1]; ++k)
          fn(lhs_species[k], -lhs_coeff[k]);
        for (size_t k = rhs_row_ptr[r]; k < rhs_row_ptr[r + 1]; ++k)
          fn(rhs_species[k], rhs_coeff[k]);
        break;
      case DENSE:
        for (unsigned int s = 0; s < nSpecies; ++s)
          if (std::abs(dense[r][s]) > 1e-30)
            fn((int)s, dense[r][s]);
        break;
      case CSC:
        // CSC stores species-major; row iteration uses the CSR forward index
        // populated during build().  Same performance as pure CSR.
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          fn(csr_cols[k], csr_vals[k]);
        break;
    }
  }

  /// O(k) lookup — k = entries per row (typically 2-10).  Diagnostic use only.
  /// Returns 0.0 for out-of-bounds indices or if species s does not participate
  /// in reaction r.
  Real get(unsigned int r, unsigned int s) const
  {
    if (r >= nReactions || s >= nSpecies)
      return 0.0;

    switch (format)
    {
      case CSR:
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          if ((unsigned int)csr_cols[k] == s) return csr_vals[k];
        return 0.0;
      case COO:
        for (size_t k = lhs_row_ptr[r]; k < lhs_row_ptr[r + 1]; ++k)
          if ((unsigned int)lhs_species[k] == s) return -lhs_coeff[k];
        for (size_t k = rhs_row_ptr[r]; k < rhs_row_ptr[r + 1]; ++k)
          if ((unsigned int)rhs_species[k] == s) return rhs_coeff[k];
        return 0.0;
      case DENSE:
        return dense[r][s];
      case CSC:
        // Species-major: scan column s for reaction r
        for (size_t k = csc_col_ptr[s]; k < csc_col_ptr[s + 1]; ++k)
          if ((unsigned int)csc_rows[k] == r) return csc_c_vals[k];
        return 0.0;
    }
    return 0.0;
  }
};

/**
 * Centralized box model for atmospheric chemistry ODE systems.
 *
 * Manages the chemical system matrices (stoichiometric coefficients,
 * reactant indices, rate constants) and provides the core dC/dt
 * computation.  Designed to support arbitrary numbers of species and
 * reactions (tested up to full MCM v3.3.1: ~5832 species, ~17224 reactions).
 *
 * The stoichiometric matrix is stored in a swappable StoichMatrix backend
 * (default: CSR for HPC; optional: COO for AtChem2-style diagnostics).
 * The "stoich_format" parameter selects the backend at construction time —
 * analogous to PETSc's -mat_type parameter for runtime format selection.
 *
 * Core computation (equivalent to F0AM dydt_eval.m and AtChem2 resid()):
 *   dC/dt = f^T * (k .* prod(C[iG], 2))
 *
 * Usage:
 *   This object is typically created by AtmosphericChemistryAction during
 *   the "add_user_object" task (box mode).
 */
class MCMBoxModel : public GeneralUserObject, public FunctionParserUtils<false>
{
public:
  static InputParameters validParams();
  MCMBoxModel(const InputParameters & params);

  using FunctionParserUtils<false>::evaluate;
  using FunctionParserUtils<false>::_func_params;

  // -- GeneralUserObject interface --
  void initialize() override;
  void execute() override {}
  void finalize() override {}

  // -- Core chemistry computation --
  /**
   * Compute dC/dt for given species concentrations.
   * @param C  Species concentrations [molec/cm^3] (length = nSpecies)
   * @param dC Output: time derivatives [molec/cm^3/s] (length = nSpecies)
   */
  void computeDCdt(const std::vector<Real> & C, std::vector<Real> & dC) const;

  /**
   * Compute the analytical Jacobian matrix d(dC/dt)/dC.
   * @param C  Species concentrations at which to evaluate
   * @param J  Output: Jacobian as triplet list (row, col, value)
   */
  void computeJacobianTriplets(const std::vector<Real> & C,
                               std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const;

  // -- Cached single-species interface (for ScalarKernel / ODEKernel) --
  /**
   * Mark the internal cache as dirty. Called by ChemistryODEKernel::reinit()
   * before each residual/Jacobian evaluation. Idempotent — multiple calls
   * between evaluations are safe.
   */
  void markDirty() const { _dirty = true; }

  /** Set current simulation time (seconds since midnight) for photolysis calculation. */
  void setCurrentTime(Real t) const { _t = t; }

  /** Compute RO2 sum (peroxy radical total) from concentration vector. */
  Real getRO2Sum(const std::vector<Real> & C) const;

  /** Get a single photolysis J value by 1-based J number (e.g., 1 for J1). */
  Real getJValue(unsigned int j_number) const;

  /** Get the number of photolysis J variables. */
  unsigned int nJValues() const { return _n_j_vars; }

  /** Get cached solar zenith angle cosine (computed during evaluateCoefficients). */
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
  /** Get latitude (degrees). */
  Real getLatitude() const { return _lat; }
  /** Get longitude (degrees). */
  Real getLongitude() const { return _lon; }
  /** Get solar declination (radians). */
  Real getDeclination() const { return _declination; }

  /**
   * Get dC/dt for a single species, caching the full computation.
   * On first call after markDirty(), computes the complete dC/dt vector
   * and caches it. Subsequent calls return cached values.
   * @param idx Species index (0..nSpecies-1)
   * @param C   Current concentration vector [molec/cm^3]
   */
  Real getDCdt(unsigned int idx, const std::vector<Real> & C) const;

  /**
   * Get diagonal Jacobian element d(dC_i)/dC_i, caching the full Jacobian.
   * On first call after markDirty(), builds the sparse Jacobian cache.
   */
  Real getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const;

  /**
   * Get off-diagonal Jacobian element d(dC_i)/dC_j, caching the full Jacobian.
   * Returns 0.0 for entries that are structurally zero in the Jacobian.
   */
  Real getJacobianOffDiagonal(unsigned int i, unsigned int j, const std::vector<Real> & C) const;

  // -- Query interface --
  /// Whether the model uses ppb units for input/output (vs molec/cm³)
  bool unitsPPB() const { return _units_ppb; }
  /// ppb → molec/cm³ conversion factor
  Real ppbToMolec() const { return _ppb_to_molec; }

  unsigned int nSpecies() const { return _n_species; }
  unsigned int nReactions() const { return _n_reactions; }
  const std::vector<std::string> & speciesNames() const { return _species_names; }
  const std::vector<std::string> & reactionNames() const { return _reaction_names; }

  /** Populate internal matrices from a ParsedMechanism. */
  void loadMechanism(const ParsedMechanism & mech);

  // -- Constrained species (AtChem2 mode) --
  /** Mark species as constrained (fixed to observed values, not solved). */
  void setConstrainedSpecies(const std::vector<std::string> & names);

  /** Update the fixed values of constrained species at current time. */
  void updateConstrainedValues(const std::vector<Real> & values);

  /** Compute dC/dt for the FULL concentration vector (constrained included). */
  void computeDCdtFull(const std::vector<Real> & C_full, std::vector<Real> & dC_full) const;

  /** Number of constrained (non-solved) species. */
  unsigned int nConstrained() const { return _constrained_set.size(); }

  // -- Detailed reaction rate output (AtChem2 lossRates/productionRates) --
  /** Compute rate (k * prod(reactants)) for a single reaction. */
  Real reactionRate(unsigned int r, const std::vector<Real> & C) const;

  /** Compute net rate contribution of reaction r to species s. */
  Real speciesReactionRate(unsigned int s, unsigned int r, const std::vector<Real> & C) const;

  /** Compute all reaction rates as a vector. */
  void allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const;

  // -- Photolysis --
  /** Enable Hybrid J (TUV lookup) photolysis. Call before loadMechanism. */
  void enableHybridPhotolysis(const std::string & table_dir);

  /** Map J<N> names to reaction indices and set up photolysis rate pointer. */
  void mapPhotolysisReactions();

  /** Update photolysis rates for current SZA/alt/albedo/O3. */
  void updatePhotolysis(Real sza, Real albedo, Real o3col, Real altitude);

  /** Load BottomUp lamp flux and reaction map. */
  void loadBottomUpData(const std::string & data_dir, const std::string & flux_file);

  /** Update photolysis via BottomUp integration (constant T,P for chamber). */
  void updatePhotolysisBottomUp();

  /** Update photolysis using SZA formula (default MCM method). */
  void updatePhotolysisSZA(Real sza, Real jfac = 1.0);

  /** ROOF (chamber cover) switch. CLOSED = all photolysis rates forced to zero. */
  void setRoofOpen(bool open) { _roof_open = open; }
  bool isRoofOpen() const { return _roof_open; }

  /** Invalidate BottomUp J-value cache. Call when photolysis-relevant
   *  parameters (T, P, lamp flux, reaction map) change during a simulation.
   *  Next evaluateCoefficients() will recompute all J-values. */
  void invalidatePhotolysisCache() { _bottomup_j_valid = false; }

  /**
   * Auto-calculate JFAC from a reference species (e.g., J4=NO2).
   * JFAC = constrained_value / parameterized_value.
   * @param ref_j_name Reference J name (e.g. "J4"), must have constrained data
   * @param constrained_val Reference species' measured (constrained) J value
   */
  void calcJFAC(const std::string & ref_j_name, Real constrained_val);

  // -- Solar cycle + convergence (F0AM nDays/Converge) --
  /** Set solar cycle params for multi-day simulation. */
  void setSolarCycle(Real lat, Real lon, int day, int month, int year);

  /** Compute cos(SZA) at given seconds since midnight (Madronich 1993). */
  Real cosSZA(Real seconds) const;

  /** Run one solar cycle day: update photolysis at each time step. */
  void advanceSolarCycle(Real seconds_from_midnight);

  // -- Dilution --
  /** Set dilution rate (kdil, /s) and background concentrations. */
  void setDilution(Real kdil, const std::vector<Real> & conc_bkgd);

  /** Compute dC/dt with dilution: dC/dt = chemistry - kdil*(C - C_bkgd). */
  void computeDCdtWithDilution(const std::vector<Real> & C, std::vector<Real> & dC) const;

  /** Check solar cycle convergence: max relative change in species. */
  Real checkConvergence(const std::vector<Real> & C_prev, const std::vector<Real> & C_curr,
                        const std::vector<std::string> & conv_species = {}) const;

protected:
  unsigned int _n_species;
  unsigned int _n_reactions;
  std::vector<std::string> _species_names;
  std::vector<std::string> _reaction_names;

  /// Stoichiometric matrix with swappable storage backend.
  /// Format selected via "stoich_format" parameter (CSR or COO).
  StoichMatrix _stoich;

  /// Reactant indices: _iG[reaction][0..2], padded with ONE index (0).
  /// Fixed-size inner arrays (3 ints) — single contiguous allocation, cache-friendly.
  std::vector<std::array<int, 3>> _iG;

  /// Pre-computed (or template) rate constants
  std::vector<Real> _k;

  // -- Constrained species tracking --
  std::set<unsigned int> _constrained_set;
  std::vector<Real> _constrained_values;

  // -- Photolysis --
  enum PhotolysisMethod { MCM_SZA, HYBRID, BOTTOMUP } _photolysis_method;
  std::unique_ptr<HybridJTableReader> _hybrid_reader;
  std::unique_ptr<BottomUpJIntegrator> _bottomup_integrator;
  std::vector<unsigned int> _j_reaction_indices;
  std::vector<Real> _j_CL, _j_CMM, _j_CNN;

  /// Units: true = ppb, false = molec/cm³
  bool _units_ppb;
  /// ppb → molec/cm³ conversion factor, set during evaluateCoefficients
  Real _ppb_to_molec;

  // -- Solar cycle --
  Real _lat, _lon;
  mutable Real _declination, _eot;
  int _day, _month, _year;

  /// Cached solar parameters (computed once per evaluateCoefficients call)
  mutable Real _solar_cosx, _solar_secx, _solar_lha;
  mutable Real _solar_sinld, _solar_cosld, _solar_eqt;

  // -- Dilution --
  Real _kdil;
  std::vector<Real> _conc_bkgd;

  /// ROOF chamber cover: false = CLOSED (all J=0), true = OPEN (normal)
  bool _roof_open;

  // -- fparser for complex rate coefficients --
  std::vector<SymFunctionPtr> _coeff_parsers;
  std::vector<SymFunctionPtr> _reaction_parsers;
  /// O(1) average lookup — evaluated every species/J in evaluateCoefficients().
  std::unordered_map<std::string, unsigned int> _name_to_index;
  /// Pre-computed J photo index offsets into _func_params (avoid string+map in hot path)
  std::vector<unsigned int> _j_photo_indices;
  unsigned int _j_index_start;
  Real _temperature;
  Real _air_density;
  Real _water_vapor;
  Real _press;       // pressure (mbar), 0 = unused
  Real _rh;          // relative humidity (%), -1 = unused
  Real _blheight;    // boundary layer height (m)
  Real _jfac;
  // Photolysis J parameters (CL/CMM/CNN per J number)
  std::vector<unsigned int> _j_numbers;
  std::vector<Real> _j_CL_vals, _j_CMM_vals, _j_CNN_vals;
  unsigned int _n_j_vars;
  void evaluateCoefficients();
  void setupFparser(const ParsedMechanism & mech);
  unsigned int computeDayOfYear() const;
  Real calculateCosSZA(Real t) const;

  /// Current simulation time (seconds since midnight), set before evaluateCoefficients
  mutable Real _t;

  // -- Cache for single-species ODEKernel interface --
  /// Dirty flag: true when cache needs recomputation
  mutable bool _dirty;
  /// Cached dC/dt vector (length = nSpecies)
  mutable std::vector<Real> _cached_dC;
  /// Cached last concentration vector for Jacobian cache keying
  mutable std::vector<Real> _cached_C;
  /// Cached diagonal Jacobian elements
  mutable std::vector<Real> _cached_diag_J;
  /// Cached off-diagonal Jacobian elements: key = (row << 32) | col
  mutable std::unordered_map<uint64_t, Real> _cached_offdiag_J;

  /// Scratch buffers for computeDCdt (Per.14 — reused across calls, no per-call allocation)
  mutable std::vector<Real> _scratch_G;
  mutable std::vector<Real> _scratch_rates;

  /// RO2 species indices in the species vector (built from parser's ro2_species)
  std::vector<unsigned int> _ro2_indices;

  /// Cached BottomUp J-values (recomputed only on T/P change)
  mutable std::map<std::string, Real> _cached_bottomup_j;
  mutable Real _cached_bottomup_T;
  mutable Real _cached_bottomup_P;
  mutable bool _bottomup_j_valid;

  // ── Per-parser variable indirection (ParseAndDeduceVariables) ──
  /// _func_params indices each coefficient parser actually references
  std::vector<std::vector<unsigned int>> _coeff_var_indices;
  /// Pre-allocated local param buffer per coefficient parser
  std::vector<std::vector<Real>> _coeff_local_params;
  /// _func_params indices each reaction parser actually references
  std::vector<std::vector<unsigned int>> _reaction_var_indices;
  /// Pre-allocated local param buffer per reaction parser
  std::vector<std::vector<Real>> _reaction_local_params;

  /// Build the sparse Jacobian cache from the current _cached_C
  void _buildJacobianCache() const;
};
