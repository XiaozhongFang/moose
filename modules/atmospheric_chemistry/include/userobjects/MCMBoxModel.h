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
#include "FunctionParserUtils.h"
#include "IMechanism.h"
#include <petscts.h>

#include "BoxIntegrator.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <map>
#include <memory>

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
  ~MCMBoxModel();

  using FunctionParserUtils<false>::evaluate;
  using FunctionParserUtils<false>::_func_params;

  // -- GeneralUserObject interface --
  void initialize() override;
  void execute() override;
  void finalize() override {}

  // -- Box integrator strategy (MOOSE implicit or PETSc TS) --
  /** Access the box integrator for per-species residual/Jacobian evaluation */
  const BoxIntegrator & getIntegrator() const { return *_integrator; }

  // -- PETSc TS standalone integrator (box mode only) --
  /** Initialize PETSc TS context */
  void setupPETScTS();
  /** Run one PETSc TS step from t0 to t1, update _ts_X in-place */
  void runPETScStep(PetscReal t0, PetscReal t1);

  // -- PETSc TS callback wrappers (static, used as PETSc function pointers) --
protected:
  static PetscErrorCode tsRHSFunction(TS ts, PetscReal t, Vec C, Vec F, void *ctx);
  static PetscErrorCode tsRHSJacobian(TS ts, PetscReal t, Vec C, Mat Amat, Mat Pmat, void *ctx);

public:
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
  void markDirty() const { if (_mechanism) _mechanism->markDirty(); }

  /** Set current simulation time (seconds since midnight) for photolysis calculation. */
  void setCurrentTime(Real t) const { if (_mechanism) _mechanism->setCurrentTime(t); }

  /** Compute RO2 sum (peroxy radical total) from concentration vector. */
  Real getRO2Sum(const std::vector<Real> & C) const;

  /** Get a single photolysis J value by 1-based J number (e.g., 1 for J1). */
  Real getJValue(unsigned int j_number) const;

  /** Get the number of photolysis J variables. */
  unsigned int nJValues() const { return _mechanism ? _mechanism->nJValues() : 0; }

  /** Get cached solar zenith angle cosine (computed during evaluateCoefficients). */
  Real getSolarCosX() const { return _mechanism ? _mechanism->getSolarCosX() : _solar_cosx; }
  /** Get cached solar secant (1/cosx). */
  Real getSolarSecX() const { return _mechanism ? _mechanism->getSolarSecX() : _solar_secx; }
  /** Get cached local hour angle (radians). */
  Real getSolarLHA() const { return _mechanism ? _mechanism->getSolarLHA() : _solar_lha; }
  /** Get sin(lat)*sin(dec). */
  Real getSolarSinLD() const { return _mechanism ? _mechanism->getSolarSinLD() : _solar_sinld; }
  /** Get cos(lat)*cos(dec). */
  Real getSolarCosLD() const { return _mechanism ? _mechanism->getSolarCosLD() : _solar_cosld; }
  /** Get equation of time. */
  Real getSolarEQT() const { return _mechanism ? _mechanism->getSolarEQT() : _solar_eqt; }
  /** Get latitude (degrees). */
  Real getLatitude() const { return _lat; }
  /** Get longitude (degrees). */
  Real getLongitude() const { return _lon; }
  /** Get solar declination (radians). */
  Real getDeclination() const { return _mechanism ? _mechanism->getDeclination() : _declination; }

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
  bool unitsPPB() const { return _mechanism ? _mechanism->unitsPPB() : _units_ppb; }
  /// ppb → molec/cm³ conversion factor
  Real ppbToMolec() const { return _mechanism ? _mechanism->ppbToMolec() : _ppb_to_molec; }

  unsigned int nSpecies() const { return _n_species; }
  unsigned int nReactions() const { return _n_reactions; }
  const std::vector<std::string> & speciesNames() const { return _species_names; }
  const std::vector<std::string> & reactionNames() const { return _reaction_names; }

  /** Get RO2 (peroxy radical) species name list. */
  const std::vector<std::string> & getRO2Species() const;
  /** Get RO2 species indices in the species vector. */
  const std::vector<unsigned int> & getRO2Indices() const;



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

  /** Compute total loss rate for species s (sum of consuming reactions). */
  Real speciesLossRate(unsigned int s, const std::vector<Real> & C) const;

  /** Compute total production rate for species s (sum of producing reactions). */
  Real speciesProductionRate(unsigned int s, const std::vector<Real> & C) const;

  /** Compute all reaction rates as a vector. */
  void allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const;

  // -- Photolysis --
  /** Enable Hybrid J (TUV lookup) photolysis. Call before loadMechanism. */
  void enableHybridPhotolysis(const std::string & table_dir);

  /** Load BottomUp lamp flux and reaction map. */
  void loadBottomUpData(const std::string & data_dir, const std::string & flux_file);

  /** ROOF (chamber cover) switch. CLOSED = all photolysis rates forced to zero. */
  void setRoofOpen(bool open) { _roof_open = open; if (_mechanism) _mechanism->setRoofOpen(open); }
  bool isRoofOpen() const { return _mechanism ? _mechanism->isRoofOpen() : _roof_open; }

  /** Set JFAC scaling factor (light intensity multiplier). */
  void setJFac(Real jfac) { _jfac = jfac; if (_mechanism) _mechanism->setJFac(jfac); }
  Real getJFac() const { return _jfac; }

  /** Invalidate BottomUp J-value cache. Call when photolysis-relevant
   *  parameters (T, P, lamp flux, reaction map) change during a simulation.
   *  Next evaluateCoefficients() will recompute all J-values. */
  void invalidatePhotolysisCache() { _bottomup_j_valid = false; if (_mechanism) _mechanism->invalidatePhotolysisCache(); }

  // -- Solar cycle + convergence (F0AM nDays/Converge) --
  /** Set solar cycle params for multi-day simulation. */
  void setSolarCycle(Real lat, Real lon, int day, int month, int year);

  /** Compute cos(SZA) at given seconds since midnight (Madronich 1993). */
  Real cosSZA(Real seconds) const;

  // -- Dilution --
  /** Set dilution rate (kdil, /s) and background concentrations. */
  void setDilution(Real kdil, const std::vector<Real> & conc_bkgd);

  /** Enable Gaussian dispersion dilution (F0AM tgauss model).
   *  dilrate = -1/(tgauss + 2*(t+t_start)) * (conc - conc_bkgd) */
  void setGaussianDispersion(Real tgauss, const std::vector<Real> & conc_bkgd,
                             Real t_start = 0.0);

  /** Compute dC/dt with dilution: dC/dt = chemistry - kdil*(C - C_bkgd). */
  void computeDCdtWithDilution(const std::vector<Real> & C, std::vector<Real> & dC) const;

  /** Check solar cycle convergence: max relative change in species. */
  Real checkConvergence(const std::vector<Real> & C_prev, const std::vector<Real> & C_curr,
                        const std::vector<std::string> & conv_species = {}) const;

protected:
  // -- Chemical mechanism delegate --
  std::unique_ptr<IMechanism> _mechanism;

  // -- Mechanism metadata (cached for fast query) --
  unsigned int _n_species;
  unsigned int _n_reactions;
  std::vector<std::string> _species_names;
  std::vector<std::string> _reaction_names;

  /// Saved initial concentrations (reloaded for self-driven integrators
  /// before the TS solve, because the MOOSE Newton solver may clear ICs)
  std::vector<Real> _initial_conc;

  /// Units: true = ppb, false = molec/cm³ (cached fallback)
  bool _units_ppb;
  /// ppb → molec/cm³ conversion factor (cached fallback)
  Real _ppb_to_molec;

  // -- Constrained species tracking --
  std::set<unsigned int> _constrained_set;
  std::vector<Real> _constrained_values;

  // -- Solar cycle (cached for postprocessor access) --
  Real _lat, _lon;
  mutable Real _declination, _eot;
  int _day, _month, _year;

  /// Cached solar parameters (fallback when mechanism not yet created)
  mutable Real _solar_cosx, _solar_secx, _solar_lha;
  mutable Real _solar_sinld, _solar_cosld, _solar_eqt;

  // -- Dilution --
  Real _kdil;
  std::vector<Real> _conc_bkgd;
  Real _tgauss;          // Gaussian time constant (s)
  Real _t_start_dil;     // Start time for Gaussian dispersion
  bool _use_gaussian;    // true = use Gaussian, false = simple first-order

  /// ROOF chamber cover: false = CLOSED (all J=0), true = OPEN (normal)
  bool _roof_open;
  Real _jfac;

  // -- Photolysis config (cached, applied to mechanism in initialize()) --
  std::string _photolysis_scheme;
  std::string _hybrid_table_dir;
  std::string _bottomup_data_dir;
  std::string _lamp_flux_file;

  /// Current simulation time (seconds since midnight)
  mutable Real _t;

  /// BottomUp cache validity flag (fallback)
  mutable bool _bottomup_j_valid;

  /// RO2 species (fallback when mechanism not yet created)
  std::vector<unsigned int> _ro2_indices;
  std::vector<std::string> _ro2_species_names;

  // -- Box integrator strategy --
  std::unique_ptr<BoxIntegrator> _integrator;

  // -- PETSc TS members --
  bool _use_box_solver = false;
  TS _ts = nullptr;
  Vec _ts_X = nullptr;
  Mat _ts_J = nullptr;
  PetscReal _solver_rtol = 1e-6;
  PetscReal _solver_atol = 1e-10;
  std::string _solver_type = "bdf";
  /// Chemical solver backend name (from chem_solver param)
  std::string _chem_solver;

  // -- SUNDIALS direct solver members --
  /// True when user selected solver_type=sundials; triggers solveSundialsCVODE()
  /// instead of PETSc TS steps inside execute().
  bool _use_sundials = false;
  /// True when chem_solver=kpp_* triggers KPP integration path in execute()
  bool _use_kpp = false;
  /// True when chem_solver=petsc_ts uses PETSc TS as self-driven integrator
  bool _use_petsc_ts = false;

public:
  // SUNDIALS direct integration entry point.  Called from execute() when
  // _use_sundials is true and _integrator->selfDriven() is true.
  void solveSundialsCVODEWrapper(Real t0, Real t1, std::vector<Real> & C);
};

