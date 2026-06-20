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

#include <string>
#include <vector>

/**
 * Centralized box model for atmospheric chemistry ODE systems.
 *
 * Manages the chemical system matrices (stoichiometric coefficients,
 * reactant indices, rate constants) and provides the core dC/dt
 * computation.  Designed to support arbitrary numbers of species and
 * reactions (tested up to full MCM v3.3.1: ~5832 species, ~17224 reactions).
 *
 * The matrix layout follows F0AM's convention:
 *   - _f[i][j] : stoichiometric coefficient for reaction i, species j
 *   - _iG[i][0..1] : reactant indices for reaction i (-1 for pseudo-first-order)
 *   - _k[i] : pre-computed rate constant for reaction i
 *
 * Core computation (equivalent to F0AM dydt_eval.m):
 *   dC/dt = f^T * (k .* prod(C[iG], 2))
 *
 * Usage:
 *   This object is typically created by MCMFacsimileAction during
 *   the "add_user_object" task.
 */
class MCMBoxModel : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MCMBoxModel(const InputParameters & params);

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

  // -- Query interface --
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

  /** Update photolysis using SZA formula (default MCM method). */
  void updatePhotolysisSZA(Real sza, Real jfac = 1.0);

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

  /// Stoichiometric coefficients: _f[reaction][species]
  std::vector<std::vector<Real>> _f;

  /// Reactant indices: _iG[reaction][0..1], -1 for pseudo-first-order
  std::vector<std::vector<int>> _iG;

  /// Pre-computed (or template) rate constants
  std::vector<Real> _k;

  // -- Constrained species tracking --
  std::set<unsigned int> _constrained_set;
  std::vector<Real> _constrained_values;

  // -- Photolysis --
  enum PhotolysisMethod { MCM_SZA, HYBRID } _photolysis_method;
  std::unique_ptr<HybridJTableReader> _hybrid_reader;
  std::vector<unsigned int> _j_reaction_indices;
  std::vector<Real> _j_CL, _j_CMM, _j_CNN;

  // -- Solar cycle --
  Real _lat, _lon, _declination, _eot;
  int _day, _month, _year;

  // -- Dilution --
  Real _kdil;
  std::vector<Real> _conc_bkgd;
};
