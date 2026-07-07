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
#include <map>
#include <set>

/**
 * Unified Action for atmospheric chemistry simulations.
 *
 * Supports two modes via the 'mode' parameter:
 *
 *   mode = box (0-D ODE box model)
 *     Creates ScalarVariable (family=SCALAR) for each species,
 *     ODETimeDerivative + ChemistryODEKernel, and an MCMBoxModel
 *     UserObject for the computation engine.
 *     Suitable for large mechanisms (up to full MCM ~5832 species).
 *
 *   mode = coupled (FEM transport + chemistry)
 *     Creates MooseVariableFE (family=LAGRANGE) for each species,
 *     MCMRatesMaterial for rate evaluation, TimeDerivative +
 *     ChemicalSourceKernel for each species, and optional Diffusion
 *     kernels via include_transport.
 *     Suitable for spatially-resolved simulations (5--50 species).
 *
 * Parses MCM Facsimile (.fac) mechanism files via MechanismLoader.
 * Replaces the deprecated MCMFacsimileAction.
 */
class AtmosphericChemistryAction : public Action
{
public:
  static InputParameters validParams();

  AtmosphericChemistryAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// Build the reactant index matrix for the Material (coupled mode only)
  std::vector<std::vector<Real>> buildReactantMatrix() const;

  /// Box mode tasks
  void actBoxAddVariable();
  void actBoxAddUserObject();
  void actBoxAddScalarKernel();
  void actBoxAddFamilyUO();
  void actBoxAddFamilyScalarKernel();

  /// Coupled mode tasks (equivalent to old MCMFacsimileAction)
  void actCoupledAddVariable();
  void actCoupledAddMaterial();
  void actCoupledAddKernel();

  /// System configuration parsed from .fac
  const std::string _mechanism_file;
  std::vector<std::string> _species;
  std::vector<std::string> _ro2_species;

  /// Parse KPP .spc files to extract species names (for mechanism_format=KPP)
  std::vector<std::string> parseKPPSpecies(const std::string & kpp_file) const;

  /// Full mechanism data (from MechanismLoader), used by coupled mode tasks
  MechanismData _mech_data;

  /// Mode: "box" or "coupled"
  const MooseEnum _mode;

  /// Include transport in coupled mode
  const bool _include_transport;
  /// Chemical solver backend (derived from chem_solver + box_solver fallback)
  MooseEnum _chem_solver;
  /// PETSc TS standalone integrator for box mode (derived from _chem_solver)
  bool _use_box_solver;

  /// Whether RO2 diagnostic variable should be created (computed once in constructor)
  bool _ro2_diagnostic_enabled;
  /// Whether the RO2 conflict warning has been printed
  bool _ro2_warning_printed;

  /// Family conservation (F0AM DAE method)
  std::vector<std::string> _family_names;
  std::vector<std::vector<std::string>> _family_members;
  std::vector<std::vector<Real>> _family_scaling;
};
