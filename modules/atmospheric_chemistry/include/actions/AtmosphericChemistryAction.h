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
 * Parses MCM Facsimile (.fac) mechanism files via MCMFacsimileParser.
 * Replaces the deprecated MCMFacsimileAction.
 */
class AtmosphericChemistryAction : public Action
{
public:
  static InputParameters validParams();

  AtmosphericChemistryAction(const InputParameters & params);

  virtual void act() override;

  const std::string & getRateCoefficient(const std::string & name) const;

protected:
  struct Reaction
  {
    std::string rate_expression;
    std::vector<std::pair<Real, std::string>> reactants;
    std::vector<std::pair<Real, std::string>> products;
  };

  /// Build the reactant index matrix for the Material (coupled mode only)
  std::vector<std::vector<Real>> buildReactantMatrix() const;

  /// Box mode tasks
  void actBoxAddVariable();
  void actBoxAddUserObject();
  void actBoxAddScalarKernel();

  /// Coupled mode tasks (equivalent to old MCMFacsimileAction)
  void actCoupledAddVariable();
  void actCoupledAddMaterial();
  void actCoupledAddKernel();

  /// System configuration parsed from .fac
  const std::string _mechanism_file;
  std::vector<std::string> _species;
  std::map<std::string, std::string> _rate_coefficients;
  std::map<std::string, Real> _photolysis_rates;
  std::map<std::string, Real> _j_CL, _j_CMM, _j_CNN;
  std::vector<std::string> _ro2_species;
  std::vector<Reaction> _reactions;
  std::vector<std::vector<Real>> _stoichiometric_matrix;
  std::map<std::string, std::string> _converted_coefficients;
  std::vector<std::string> _eval_order;
  std::vector<std::string> _reaction_rate_expressions;
  std::set<std::string> _coefficient_names;
  std::set<std::string> _base_variables;

  /// Mode: "box" or "coupled"
  const MooseEnum _mode;

  /// Include transport in coupled mode
  const bool _include_transport;

  /// Whether RO2 diagnostic variable should be created (computed once in constructor)
  bool _ro2_diagnostic_enabled;
  /// Whether the RO2 conflict warning has been printed
  bool _ro2_warning_printed;
};
