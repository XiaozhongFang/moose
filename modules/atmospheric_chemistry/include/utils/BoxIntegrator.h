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
#include "ConsoleStreamInterface.h"

#include <vector>
#include <memory>

class MooseApp;
class MCMBoxModel;

/**
 * Abstract interface for box-model chemical integration.
 *
 * Provides two interaction modes:
 *   - Self-driven mode: the integrator handles the full ODE solve internally
 *     (PetscTSIntegrator). Residual/Jacobian evaluation returns 0 — MOOSE's
 *     nonlinear solver is bypassed.
 *   - MOOSE-driven mode: the integrator provides per-species residual and
 *     Jacobian evaluations for MOOSE's Newton solver (MooseImplicitIntegrator).
 *
 * ChemistryODEKernel holds a reference to BoxIntegrator and calls
 * computeResidual/Jacobian*() unconditionally — the kernel has no mode-specific
 * branching.  The integrator's selfDriven() flag controls whether MCMBoxModel's
 * execute() runs a full-system integration step.
 */
class BoxIntegrator
{
public:
  virtual ~BoxIntegrator() = default;

  ///@{
  /// Chemical source term for species idx: dC[species_idx]/dt (molec/cm³/s).
  virtual Real computeResidual(unsigned int species_idx,
                                const std::vector<Real> & C) const = 0;

  /// Diagonal Jacobian: ∂(dC[species_idx]/dt) / ∂C[species_idx].
  virtual Real computeJacobianDiagonal(unsigned int species_idx,
                                        const std::vector<Real> & C) const = 0;

  /// Off-diagonal Jacobian: ∂(dC[species_idx]/dt) / ∂C[jvar].
  virtual Real computeJacobianOffDiagonal(unsigned int species_idx,
                                           unsigned int jvar,
                                           const std::vector<Real> & C) const = 0;
  ///@}

  /**
   * Reinitialize integrator state at the start of a new timestep.
   * In MOOSE-driven mode, this invalidates caches and updates the time for
   * photolysis computation.  In self-driven mode, this is a no-op.
   */
  virtual void reinit(Real time) const = 0;

  /**
   * Whether this integrator drives the ODE solve itself.
   * If true, MCMBoxModel::execute() runs a full-system integration step;
   * ChemistryODEKernel residuals return 0.  If false, MOOSE's nonlinear
   * solver owns the solve and ChemistryODEKernel provides real residuals.
   */
  virtual bool selfDriven() const = 0;

  /**
   * ppb → molec/cm³ conversion factor.
   * Used by ChemistryODEKernel for unit conversion of the concentration vector.
   */
  virtual Real ppbToMolec() const = 0;
};

/**
 * MOOSE-implicit integrator: wraps MCMBoxModel for per-species evaluation.
 *
 * computeResidual/Jacobian*() forward to MCMBoxModel's getDCdt() /
 * getJacobian*() methods.  solve() is a no-op — MOOSE's Newton solver
 * drives the timestep via ChemistryODEKernel.
 */
class MooseImplicitIntegrator : public BoxIntegrator
{
public:
  MooseImplicitIntegrator(const MCMBoxModel & box_model);
  ~MooseImplicitIntegrator() override = default;

  Real computeResidual(unsigned int species_idx,
                        const std::vector<Real> & C) const override;

  Real computeJacobianDiagonal(unsigned int species_idx,
                                const std::vector<Real> & C) const override;

  Real computeJacobianOffDiagonal(unsigned int species_idx,
                                   unsigned int jvar,
                                   const std::vector<Real> & C) const override;

  void reinit(Real time) const override;
  bool selfDriven() const override { return false; }
  Real ppbToMolec() const override;

private:
  const MCMBoxModel & _box;
};

/**
 * PETSc TS integrator: wraps MCMBoxModel for full-system integration.
 *
 * computeResidual/Jacobian*() return a tiny damping term (1e-20) to keep
 * MOOSE's Jacobian non-singular.  The actual chemistry integration is
 * handled by MCMBoxModel::runPETScStep() in execute().
 * Used in box mode when petsc_ts=true.
 */
class PetscTSIntegrator : public BoxIntegrator
{
public:
  PetscTSIntegrator(const MCMBoxModel & box_model) : _box(box_model) {}
  ~PetscTSIntegrator() override = default;

  Real computeResidual(unsigned int /*species_idx*/,
                        const std::vector<Real> & /*C*/) const override
  { return 0.0; }  // chemistry handled by execute(); no-op for MOOSE Newton

  Real computeJacobianDiagonal(unsigned int /*species_idx*/,
                                const std::vector<Real> & /*C*/) const override
  { return 0.0; }  // zero Jacobian — ChemistryODEKernel is a true no-op

  Real computeJacobianOffDiagonal(unsigned int /*species_idx*/,
                                   unsigned int /*jvar*/,
                                   const std::vector<Real> & /*C*/) const override { return 0.0; }

  void reinit(Real /*time*/) const override {}
  bool selfDriven() const override { return true; }
  Real ppbToMolec() const override { return 1.0; }

private:
  const MCMBoxModel & _box;
};

// The integrator remains declared without SUNDIALS so its existing runtime
// stub can report that the selected solver is unavailable.  When SUNDIALS is
// available, pull in the types needed by the callback signatures.
#if defined(HAVE_SUNDIALS)
#include <sundials/sundials_types.h>    // sunrealtype, sunindextype, SUNContext
#include <sundials/sundials_matrix.h>  // SUNMatrix
#include <nvector/nvector_serial.h>    // N_Vector (defines N_Vector)
#endif

/**
 * SUNDIALS CVODE/ARKODE integrator: wraps MCMBoxModel for full-system integration
 * using the standalone SUNDIALS library (bypassing PETSc TS).
 *
 * Without SUNDIALS support, solveSundialsCVODE() reports a runtime build error.
 * computeResidual/Jacobian*() return 0, matching the PetscTSIntegrator no-op contract.
 * selfDriven() === true so MCMBoxModel::execute() dispatches to
 * solveSundialsCVODE() instead of runPETScStep().
 *
 * The actual SUNDIALS solve is implemented in BoxIntegrator.C.
 */
class SundialsBoxIntegrator : public BoxIntegrator, public ConsoleStreamInterface
{
public:
  SundialsBoxIntegrator(const MCMBoxModel & box_model,
                         MooseApp & app,
                         Real rtol = 1e-6,
                         Real atol = 1e-10)
    : BoxIntegrator(),
      ConsoleStreamInterface(app),
      _box(box_model),
      _rtol(rtol),
      _atol(atol)
  {}
  ~SundialsBoxIntegrator() override = default;

  Real computeResidual(unsigned int /*species_idx*/,
                        const std::vector<Real> & /*C*/) const override { return 0.0; }

  Real computeJacobianDiagonal(unsigned int /*species_idx*/,
                                const std::vector<Real> & /*C*/) const override { return 0.0; }

  Real computeJacobianOffDiagonal(unsigned int /*species_idx*/,
                                   unsigned int /*jvar*/,
                                   const std::vector<Real> & /*C*/) const override { return 0.0; }

  void reinit(Real /*time*/) const override {}
  bool selfDriven() const override { return true; }
  Real ppbToMolec() const override { return 1.0; }

  /**
   * Advance the chemistry by one MOOSE timestep using SUNDIALS CVODE.
   *
   * @param t0  Start time of the current step (s).
   * @param t1  End time of the current step (s).
   * @param C   On entry: current species concentrations (molec/cm³).
   *            On exit: integrated species concentrations at t1.
   */
  void solveSundialsCVODE(Real t0, Real t1, std::vector<Real> & C) const;

private:
  const MCMBoxModel & _box;
  const Real _rtol;
  const Real _atol;

#if defined(HAVE_SUNDIALS)
  /** SUNDIALS RHS callback: reads N_Vector y, writes dy/dt into dy. */
  static int sundialsRHSF(sunrealtype t, N_Vector y, N_Vector dy, void *user_data);

  /**
   * Analytical Jacobian callback for SUNDIALS (CVodeSetJacFn).
   *
   * Matches the KPP/F0AM/AtChem2 pattern: the chemistry-aware Jacobian
   * ∂(dC_i/dt)/∂C_j is computed by chain rule over the reaction network
   * (_iG stoichiometry + _k rate coefficients) and written directly into
   * the dense SUNMatrix J.
   *
   * Signature required by SUNDIALS:
   *   int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
   *            void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
   *
   * This replaces SUNDIALS' built-in finite-difference Jacobian, which is
   * too fragile on the extremely stiff MCM mechanism (600+ species, 30+
   * orders-of-magnitude spread in lifetimes) and caused Newton divergence.
   */
  static int sundialsJacFn(sunrealtype t, N_Vector y, N_Vector fy,
                            SUNMatrix J, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
};
