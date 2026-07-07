//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryBoxExecutioner.h"
#include "FEProblem.h"
#include "MooseEnum.h"
#include "MCMBoxModel.h"

#include <iomanip>

registerMooseObject("AtmosphericChemistryApp", AtmosphericChemistryBoxExecutioner);

InputParameters
AtmosphericChemistryBoxExecutioner::validParams()
{
  InputParameters params = Executioner::validParams();

  params.addRequiredParam<std::string>("box_model",
      "Name of the MCMBoxModel UserObject that manages the chemistry integration.");

  params.addParam<Real>("start_time", 0.0, "Start time of the simulation");
  params.addParam<Real>("end_time", 1.0, "End time of the simulation");
  params.addParam<Real>("dt", 1.0, "Initial time step size");
  params.addParam<Real>("dtmin", 1e-10, "Minimum time step size");
  params.addParam<Real>("dtmax", 1e10, "Maximum time step size");
  params.addParam<int>("num_steps", 100000, "Maximum number of time steps");

  params.addParam<std::string>("chem_solver", "petsc_ts",
      "Chemical solver backend: petsc_ts, sundials, kpp_rosenbrock, "
      "kpp_sdirk, kpp_runge_kutta, or moose_implicit.");

  params.addClassDescription(
      "Executioner for 0-D ODE box-model atmospheric chemistry. "
      "For self-driven solvers (petsc_ts, sundials, kpp_*), manages the "
      "time-stepping loop directly. For moose_implicit, the standard "
      "MOOSE nonlinear solve is used.");
  return params;
}

AtmosphericChemistryBoxExecutioner::AtmosphericChemistryBoxExecutioner(
    const InputParameters & params)
  : Executioner(params),
    _start_time(getParam<Real>("start_time")),
    _end_time(getParam<Real>("end_time")),
    _dt(getParam<Real>("dt")),
    _dtmin(getParam<Real>("dtmin")),
    _dtmax(getParam<Real>("dtmax")),
    _num_steps(getParam<int>("num_steps")),
    _t_step(0),
    _time(_start_time),
    _box_model_name(getParam<std::string>("box_model")),
    _last_converged(true),
    _chem_solver(getParam<std::string>("chem_solver")),
    _self_driven(_chem_solver != "moose_implicit")
{
  if (_start_time > _end_time)
    mooseError("AtmosphericChemistryBoxExecutioner: start_time > end_time");
  if (_dt <= 0.0)
    mooseError("AtmosphericChemistryBoxExecutioner: dt must be positive");
  if (_num_steps <= 0)
    mooseError("AtmosphericChemistryBoxExecutioner: num_steps must be positive");
}

void
AtmosphericChemistryBoxExecutioner::execute()
{
  // For moose_implicit, we currently rely on the standard Transient path.
  // The Box Executioner with moose_implicit would need FEProblemSolve.
  if (!_self_driven)
    mooseError("AtmosphericChemistryBoxExecutioner: moose_implicit solver "
               "is not yet supported. Use [Executioner] type=Transient instead.");

  _console << "\nAtmosphericChemistryBoxExecutioner: Starting box chemistry"
           << "\n  solver: " << _chem_solver
           << "\n  time: [" << _start_time << ", " << _end_time << "]"
           << "\n  dt: " << _dt
           << "\n  max steps: " << _num_steps
           << std::endl;

  // Get reference to the MCMBoxModel UserObject
  auto & box_model = const_cast<MCMBoxModel &>(
      _fe_problem.getUserObject<MCMBoxModel>(_box_model_name));

  // Initialize time at start
  _fe_problem.time() = _start_time;
  _fe_problem.dt() = _dt;
  _fe_problem.timeStep() = 0;

  // Output initial state
  outputStep();

  // Main time loop
  _t_step = 0;
  _time = _start_time;
  _last_converged = true;

  while (_t_step < _num_steps && _time < _end_time && _last_converged)
  {
    _t_step++;
    Real t_start = _time;
    Real dt_eff = std::min(_dt, _end_time - _time);
    Real t_end = t_start + dt_eff;

    // Update FEProblem time state for this step
    _fe_problem.time() = t_end;
    _fe_problem.timeOld() = t_start;
    _fe_problem.dt() = dt_eff;
    _fe_problem.timeStep() = _t_step;

    // Perform chemistry integration
    try
    {
      box_model.stepChemistry(t_start, t_end);
      _last_converged = true;
    }
    catch (std::exception &)
    {
      _last_converged = false;
      _console << "AtmosphericChemistryBoxExecutioner: Step " << _t_step
               << " FAILED at t=[" << t_start << "," << t_end << "]"
               << std::endl;
      break;
    }

    // Output
    printStep(t_end, dt_eff, _last_converged);
    outputStep();

    // Advance time
    _time = t_end;
  }

  // Final output at end_time
  if (_time >= _end_time)
  {
    _console << "\nAtmosphericChemistryBoxExecutioner: Finished at t=" << _time
             << " (" << _t_step << " steps)" << std::endl;
  }
  else if (!_last_converged)
  {
    _console << "\nAtmosphericChemistryBoxExecutioner: ABORTED at t=" << _time
             << " (step " << _t_step << " failed)" << std::endl;
  }
  else
  {
    _console << "\nAtmosphericChemistryBoxExecutioner: STOPPED at t=" << _time
             << " (max steps reached)" << std::endl;
  }
}

void
AtmosphericChemistryBoxExecutioner::outputStep()
{
  // Trigger CSV/Exodus output at specified intervals
  _fe_problem.outputStep(EXEC_TIMESTEP_END);
}

void
AtmosphericChemistryBoxExecutioner::printStep(Real t, Real dt, bool converged)
{
  if (converged)
    _console << "  Step " << std::setw(5) << _t_step
             << ": t=" << std::setw(12) << std::setprecision(6) << std::scientific << t
             << " dt=" << std::setw(12) << dt
             << std::endl;
  else
    _console << "  Step " << std::setw(5) << _t_step
             << ": t=" << std::setw(12) << t
             << " ** FAILED **" << std::endl;
}
