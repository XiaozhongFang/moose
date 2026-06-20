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

#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <queue>

/**
 * Parsed reaction from an MCM Facsimile (.fac) file.
 */
struct ParsedReaction
{
  std::string rate_expression;
  std::vector<std::pair<Real, std::string>> reactants;
  std::vector<std::pair<Real, std::string>> products;
};

/**
 * Complete parsed mechanism — ready for building an ODE system.
 */
struct ParsedMechanism
{
  std::vector<std::string> species;
  std::vector<std::string> reaction_names;
  std::vector<ParsedReaction> reactions;

  /// Full stoichiometric matrix: stoichiometry[species][reaction]
  std::vector<std::vector<Real>> stoichiometry;

  /// Reactant indices: reactant_indices[reaction] = {idx0, idx1} (-1 for pseudo-first-order)
  std::vector<std::vector<int>> reactant_indices;

  /// Rate coefficient expressions (in topological evaluation order)
  std::vector<std::string> coefficient_names;
  std::vector<std::string> coefficient_expressions;

  /// Photolysis parameter vectors (sorted by J number)
  std::vector<unsigned int> j_numbers;
  std::vector<Real> j_CL, j_CMM, j_CNN;

  /// RO2 (peroxy radical) species list (AtChem2-compatible)
  std::vector<std::string> ro2_species;
};

/**
 * Standalone parser for MCM Facsimile-format (.fac) mechanism files.
 *
 * Extracts species, reactions, rate coefficients, and stoichiometric
 * relationships independently of any MOOSE Action.  Both
 * MCMFacsimileAction (for variable/kernel setup) and MCMBoxModel (for
 * direct computation) use this parser.
 */
class MCMFacsimileParser
{
public:
  /**
   * Parse a .fac file and return the complete mechanism.
   * @param filename Path to the .fac file
   * @param photolysis_file Optional path to photolysis parameter file
   */
  ParsedMechanism parse(const std::string & filename,
                        const std::string & photolysis_file = "");

private:
  void parseFile(const std::string & filename);
  void processStatement(const std::string & statement);
  void parseReactionLine(const std::string & line);
  void extractSpecies(const std::string & side,
                      std::vector<std::pair<Real, std::string>> & out);
  std::string convertRateExpression(const std::string & expr) const;
  void scanPhotolysisReferences();
  void loadPhotolysisParameters(const std::string & filename);
  void buildStoichiometricMatrix();
  void buildReactantIndices();
  void expandRateCoefficients();

  std::vector<std::string> _species;
  std::vector<ParsedReaction> _reactions;
  std::map<std::string, std::string> _rate_coefficients;
  std::map<std::string, Real> _j_CL, _j_CMM, _j_CNN;
  std::set<std::string> _photolysis_rates;
  std::set<std::string> _coefficient_names;
  std::set<std::string> _base_variables;
  std::vector<std::string> _eval_order;
  std::map<std::string, std::string> _converted_coefficients;

  /// Transposed stoichiometric matrix: stoichiometry[species][reaction]
  std::vector<std::vector<Real>> _stoichiometric;

  /// Reactant indices per reaction: _reactant_indices[reaction][0..1]
  std::vector<std::vector<int>> _reactant_indices;
};
