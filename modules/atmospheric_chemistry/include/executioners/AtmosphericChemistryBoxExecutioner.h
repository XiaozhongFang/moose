//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Executioner.h"
#include "BoxIntegrator.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Executioner for 0-D ODE box-model atmospheric chemistry simulations.
 *
 * For self-driven solvers (petsc_ts, sundials, kpp_*), this Executioner
 * takes over the full time integration loop — reading concentrations from
 * ScalarVariables, advancing with the chemistry integrator, and writing
 * results back — bypassing MOOSE's normal nonlinear solve.
 *
 * For moose_implicit mode, the Executioner delegates to the standard
 * Transient solve path (MOOSE's Newton solver owns the integration).
 *
 * Replaces the standard [Executioner] type=Transient for box-mode
 * atmospheric chemistry problems.
 */
class AtmosphericChemistryBoxExecutioner : public Executioner
{
public:
  static InputParameters validParams();

  AtmosphericChemistryBoxExecutioner(const InputParameters & params);

  virtual void execute() override;
  virtual bool lastSolveConverged() const override { return _last_converged; }

protected:
  /// Read species concentrations from ScalarVariables into a vector
  std::vector<Real> readConcentrations() const;

  /// Write integrated concentrations back to ScalarVariables and solution vectors
  void writeConcentrations(const std::vector<Real> & C);

  /// Perform one time step with the self-driven integrator
  void stepSelfDriven(Real t0, Real t1);

  /// Perform one time step using MOOSE's standard Transient solve
  void stepMooseImplicit(Real dt);

  /// Output at current time/step
  void outputStep();

  /// Print timestep summary
  void printStep(Real t, Real dt, bool converged);

  /// Time management
  Real _start_time;
  Real _end_time;
  Real _dt;
  Real _dtmin;
  Real _dtmax;
  int _num_steps;
  int _t_step;
  Real _time;

  /// UserObject name for the MCMBoxModel
  std::string _box_model_name;

  /// Species variable names (from the box model)
  std::vector<std::string> _species_names;

  /// Whether the last solve converged
  bool _last_converged;

  /// Solver mode
  std::string _chem_solver;
  bool _self_driven;
};
