//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"

#include <vector>
#include <string>
#include <tuple>

/**
 * Physical parameters for chemical rate coefficient evaluation.
 *
 * Bundles 20+ scattered member variables from MCMBoxModel into a single
 * struct that can be passed through the IMechanism interface.
 * This struct is MOOSE-agnostic — only Real and std::vector types.
 */
struct PhysParams
{
  Real temperature = 298.15;     // K
  Real air_density = 2.46e19;   // molec/cm³ (used when press == 0)
  Real water_vapor = 2.46e17;   // molec/cm³ (used when rh < 0)
  Real pressure = 0.0;          // mbar (>0 → air_density computed dynamically)
  Real rh = -1.0;               // % (>=0 → water_vapor computed dynamically)
  Real jfac = 1.0;              // Photolysis scaling factor
  Real latitude = 51.51;        // deg N
  Real longitude = 0.13;        // deg E
  Real albedo = 0.1;
  Real o3column = 350.0;        // Dobson Units (for HYBRID photolysis)
  Real altitude = 0.0;          // m (for HYBRID photolysis)
  Real cos_sza = 1.0;           // Solar zenith angle cosine
  Real blheight = 0.0;          // Boundary layer height (m)

  // Photolysis rate values J[N] (computed externally by PhotolysisEngine)
  std::vector<Real> j_vals;
};

/**
 * Per-species production and loss rates for diagnostic output.
 */
struct SpeciesRates
{
  std::vector<Real> production;
  std::vector<Real> loss;
};

/**
 * Abstract mechanism evaluation interface.
 *
 * Pure virtual interface for computing chemical RHS, Jacobian, and
 * diagnostic rates.  Two concrete implementations:
 *
 *   1. MCMRuntimeMechanism — parses .fac at runtime, evaluates via
 *      StoichMatrix + fparser.  Flexible, no pre-generation step.
 *
 *   2. KPPGeneratedMechanism — wraps KPP-generated C code loaded
 *      via dlopen.  Maximum performance, mechanism-specific.
 *
 * This interface is deliberately MOOSE-agnostic: no dependence on
 * FEProblem, MooseVariable, or other MOOSE framework objects.
 * Only std::vector<Real>, plain Real, and the MOOSE Real typedef.
 */
class IMechanism
{
public:
  virtual ~IMechanism() = default;

  ///@{
  /// Mechanism dimensions and metadata.
  virtual unsigned int nSpecies() const = 0;
  virtual unsigned int nReactions() const = 0;
  virtual const std::vector<std::string> & speciesNames() const = 0;
  ///@}

  /// Update mechanism-internal state when physical parameters change.
  /// Triggers rate coefficient re-evaluation (fparser re-eval or KPP Update_RCONST).
  virtual void updateParams(const PhysParams & params) = 0;

  /// Compute dC/dt = f(C, params) for all species.
  /// @param[in]  t     Current time (seconds since midnight, for photolysis)
  /// @param[in]  C     Species concentrations [molec/cm³] (length nSpecies)
  /// @param[in]  params Physical parameters (temperature, pressure, J-values, etc.)
  /// @param[out] dC_dt Time derivatives [molec/cm³/s] (length nSpecies)
  virtual void computeRHS(Real t,
                           const std::vector<Real> & C,
                           const PhysParams & params,
                           std::vector<Real> & dC_dt) const = 0;

  /// Compute analytical Jacobian as (row, col, value) triplets.
  /// Default implementation uses finite differences — override with
  /// analytical version for better performance.
  virtual void computeJacobian(
      Real t,
      const std::vector<Real> & C,
      const PhysParams & params,
      std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const;

  /// Compute per-species production and loss rates.
  virtual SpeciesRates computeSpeciesRates(
      Real t,
      const std::vector<Real> & C,
      const PhysParams & params) const = 0;

  /// ppb → molec/cm³ conversion factor.
  virtual Real ppbToMolec() const { return 1.0; }
  virtual bool unitsPPB() const { return false; }
};
