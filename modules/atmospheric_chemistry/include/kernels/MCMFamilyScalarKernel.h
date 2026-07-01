//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ODEKernel.h"
#include "MCMBoxModel.h"
#include "MCMFamilyConstraint.h"

/**
 * ScalarKernel implementing F0AM-style family conservation via DAE.
 *
 * For the DAE slack variable (first member of each family):
 *   The chemical source is corrected so that d(F_total)/dt = 0.
 *   R = du/dt - (sum_{members} scaling_i * dC_i/dt - dC_slack/dt_original)
 *   where dC_i/dt_original is the chemical source, and the correction
 *   ensures F_total is conserved.
 *
 * For non-slack family members: the residual is unchanged (normal ODE).
 * The slack variable absorbs the algebraic constraint.
 *
 * Reference: F0AM dydt_eval.m family section (lines ~85-87):
 *   dydt(:,j(m)) = Fc - Ft    where Fc = sum(scaling * conc(members))
 *                                    Ft = family_target
 */
class MCMFamilyScalarKernel : public ODEKernel
{
public:
  static InputParameters validParams();
  MCMFamilyScalarKernel(const InputParameters & params);

  virtual void reinit() override;
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobianScalar(unsigned int jvar) override;

protected:
  const MCMBoxModel & _box_model;
  const MCMFamilyConstraint & _family_uo;
  /// Species index this kernel is attached to
  const unsigned int _species_idx;
  /// Family name this kernel belongs to (empty = no family)
  const std::string _family_name;
  /// Species variables for building concentration vector
  std::vector<const MooseVariableScalar *> _species_vars;
  mutable std::vector<Real> _C_buffer;
  mutable bool _cached;
  mutable std::vector<Real> _cached_dCdt_all;
  const std::vector<Real> & _buildC() const;
};
