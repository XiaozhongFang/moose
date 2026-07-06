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

#include <string>
#include <vector>

/**
 * IMechanism implementation wrapping KPP-generated C code via dlopen.
 *
 * Loads a KPP-generated shared library (e.g. libkpp_small_test.so) at
 * construction time, resolves the KPP C API symbols (Fun, Jac_SP, etc.),
 * and implements the IMechanism interface by calling those KPP functions
 * directly.
 *
 * This provides a high-performance alternative to MCMRuntimeMechanism
 * when a KPP-generated mechanism is available.  The mechanism-specific
 * .so is compiled from KPP-generated .c files together with the
 * kpp_adapter.c bridge.
 *
 * Reuses KppBoxIntegrator for the self-driven integration path
 * (which calls KPP's Rosenbrock integrator via kpp_integrate).
 * This class is the IMechanism front-end for RHS/Jacobian evaluations
 * used by MOOSE's implicit time integration (PETSc TS).
 */
class KPPGeneratedMechanism : public IMechanism
{
public:
  /**
   * Load a KPP-generated shared library.
   *
   * @param lib_path  Path to the KPP shared library (e.g. "libkpp_small_test.so").
   *                  If not an absolute path, searched via system dlopen rules
   *                  (LD_LIBRARY_PATH, RPATH, etc.).
   */
  KPPGeneratedMechanism(const std::string & lib_path);

  ~KPPGeneratedMechanism() override;

  // ===== IMechanism interface =====

  unsigned int nSpecies() const override { return _n_species; }
  unsigned int nReactions() const override { return _n_reactions; }
  const std::vector<std::string> & speciesNames() const override { return _species_names; }

  void updateParams(const PhysParams & params) override;
  void computeRHS(Real t,
                   const std::vector<Real> & C,
                   const PhysParams & params,
                   std::vector<Real> & dC_dt) const override;
  void computeJacobian(
      Real t,
      const std::vector<Real> & C,
      const PhysParams & params,
      std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const override;
  SpeciesRates computeSpeciesRates(
      Real t,
      const std::vector<Real> & C,
      const PhysParams & params) const override;

  // ===== Accessors / configuration =====

  void setRoofOpen(bool open) override { _roof_open = open; }
  bool isRoofOpen() const override { return _roof_open; }
  void setJFac(Real jfac) override { _jfac = jfac; }
  void invalidatePhotolysisCache() override {}
  void markDirty() const override {}
  void setCurrentTime(Real t) const override { _t = t; }

  Real getJValue(unsigned int j_number) const override;
  unsigned int nJValues() const override { return _n_j_vals; }

  Real getRO2Sum(const std::vector<Real> & C) const override;

  Real reactionRate(unsigned int r, const std::vector<Real> & C) const override;
  Real speciesReactionRate(unsigned int s,
                            unsigned int r,
                            const std::vector<Real> & C) const override;
  Real speciesLossRate(unsigned int s, const std::vector<Real> & C) const override;
  Real speciesProductionRate(unsigned int s, const std::vector<Real> & C) const override;
  void allReactionRates(const std::vector<Real> & C,
                         std::vector<Real> & rates) const override;

  Real getDCdt(unsigned int idx, const std::vector<Real> & C) const override;
  Real getJacobianDiagonal(unsigned int idx,
                            const std::vector<Real> & C) const override;
  Real getJacobianOffDiagonal(unsigned int i,
                               unsigned int j,
                               const std::vector<Real> & C) const override;

private:
  // ===== KPP C API function pointer types =====

  /// Fun(Y, FIX, RCONST, Ydot) — compute RHS
  using KppFunFn = void (*)(double[], double[], double[], double[]);
  /// Jac_SP(Y, FIX, RCONST, JVS) — compute sparse Jacobian
  using KppJacFn = void (*)(double[], double[], double[], double[]);
  /// Initialize() — set up KPP global state
  using KppInitFn = void (*)(void);

  // ===== Loaded function pointers =====

  void * _lib_handle;
  KppFunFn _fun;
  KppJacFn _jac;
  KppInitFn _init;

  // ===== Pointers to KPP global variables (resolved via dlsym) =====

  int * _kpp_nspec;
  int * _kpp_nvar;
  int * _kpp_nreact;
  double * _kpp_C;       // C[NSPEC] — concentration array
  double * _kpp_VAR;     // VAR[NVAR] — variable species array
  double * _kpp_FIX;     // FIX[] — fixed species
  double * _kpp_RCONST;  // RCONST[] — rate constants

  // ===== Mechanism metadata =====

  unsigned int _n_species;
  unsigned int _n_reactions;
  std::vector<std::string> _species_names;

  // ===== Runtime state =====

  bool _roof_open;
  Real _jfac;
  mutable Real _t;
  unsigned int _n_j_vals;
};
