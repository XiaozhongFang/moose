//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "MechanismLoader.h"

#include <string>
#include <vector>

/**
 * Action for 0-D ODE box-model atmospheric chemistry simulations.
 *
 * Creates ScalarVariable (family=SCALAR) for each species,
 * ODETimeDerivative + ChemistryODEKernel, MCMBoxModel UserObject,
 * and optional family-conservation DAE constraints.
 *
 * Suitable for large mechanisms (up to full MCM ~5832 species).
 *
 * Replaces [AtmosphericChemistry] mode=box.
 */
class AtmosphericChemistryBoxAction : public Action
{
public:
  static InputParameters validParams();

  AtmosphericChemistryBoxAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// Create scalar variables for each chemical species
  void actAddVariable();
  /// Create the MCMBoxModel UserObject
  void actAddUserObject();
  /// Create ChemistryODEKernel + ODETimeDerivative for each species
  void actAddScalarKernel();
  /// Create family-conservation UserObject (if families configured)
  void actAddFamilyUserObject();
  /// Create family-conservation scalar kernels (if families configured)
  void actAddFamilyScalarKernel();

  /// Parsed species names in mechanism order
  std::vector<std::string> _species;
  /// Full mechanism data (for sparse coupling matrix in moose_implicit mode)
  MechanismData _mech_data;
  /// Chemical solver backend (from chem_solver param)
  std::string _chem_solver;
  /// Whether to use a self-driven integrator (PETSc TS, SUNDIALS, KPP)
  bool _use_box_solver;
  /// Whether the selected mechanism is handled by KPP-generated code.
  bool _use_kpp;
  /// Whether RO2 diagnostic variable should be created
  bool _ro2_diagnostic_enabled;
  /// Family conservation (F0AM DAE method)
  std::vector<std::string> _family_names;
  std::vector<std::vector<std::string>> _family_members;
  std::vector<std::vector<Real>> _family_scaling;
};
