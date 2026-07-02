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
#include <memory>

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
 * computeResidual/Jacobian*() return 0 — ChemistryODEKernel is a no-op.
 * solve() runs the PETSc TS integration via MCMBoxModel::runPETScStep().
 * Used in box mode when petsc_ts=true.
 */
class PetscTSIntegrator : public BoxIntegrator
{
public:
  PetscTSIntegrator(const MCMBoxModel & box_model) : _box(box_model) {}
  ~PetscTSIntegrator() override = default;

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

private:
  const MCMBoxModel & _box;
};
