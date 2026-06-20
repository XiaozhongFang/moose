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
 * Action that parses an MCM (Master Chemical Mechanism) Facsimile-format
 * mechanism file and sets up the corresponding ODE system.
 *
 * Delegates all .fac parsing to MCMFacsimileParser (shared with MCMBoxModel).
 *
 * Creates:
 *   - Nonlinear variables for each chemical species
 *   - MCMRatesMaterial for evaluating reaction rates
 *   - TimeDerivative + ChemicalSourceKernel for each species
 *   - Optional transport kernels (include_transport = true)
 */
class MCMFacsimileAction : public Action
{
public:
  static InputParameters validParams();

  MCMFacsimileAction(const InputParameters & params);

  virtual void act() override;

  const std::string & getRateCoefficient(const std::string & name) const;

protected:
  struct Reaction
  {
    std::string rate_expression;
    std::vector<std::pair<Real, std::string>> reactants;
    std::vector<std::pair<Real, std::string>> products;
  };

  /// Build the reactant index matrix for the Material
  std::vector<std::vector<Real>> buildReactantMatrix() const;

  const std::string _mechanism_file;
  std::vector<std::string> _species;
  std::map<std::string, std::string> _rate_coefficients;
  std::map<std::string, Real> _photolysis_rates;
  std::map<std::string, Real> _j_CL, _j_CMM, _j_CNN;
  std::vector<Reaction> _reactions;
  std::vector<std::vector<Real>> _stoichiometric_matrix;
  std::map<std::string, std::string> _converted_coefficients;
  std::vector<std::string> _eval_order;
  std::vector<std::string> _reaction_rate_expressions;
  std::set<std::string> _coefficient_names;
  std::set<std::string> _base_variables;
  bool _include_transport;
};
