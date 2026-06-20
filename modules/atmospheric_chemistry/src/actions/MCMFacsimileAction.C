//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMFacsimileAction.h"
#include "MCMFacsimileParser.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "pcrecpp.h"

#include <fstream>
#include <sstream>
#include <algorithm>

registerMooseAction("AtmosphericChemistryApp", MCMFacsimileAction, "add_variable");
registerMooseAction("AtmosphericChemistryApp", MCMFacsimileAction, "add_material");
registerMooseAction("AtmosphericChemistryApp", MCMFacsimileAction, "add_kernel");

InputParameters
MCMFacsimileAction::validParams()
{
  InputParameters params = Action::validParams();
  params.addRequiredParam<std::string>(
      "mechanism_file", "Path to the MCM Facsimile-format mechanism file (.fac)");
  params.addParam<bool>("include_transport", false,
      "Add Diffusion kernels for non-placeholder species");
  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19, "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17, "Background water vapor concentration (molecules/cm^3)");
  params.addParam<std::string>(
      "mcm_photolysis_file",
      "doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat",
      "Path to the MCM photolysis-rates parameter file for SZA-based J calculation.");
  params.addParam<Real>("latitude", 51.51, "Latitude in degrees (North positive)");
  params.addParam<Real>("longitude", 0.13, "Longitude in degrees (East positive)");
  params.addParam<unsigned int>("day", 21, "Day of month for solar zenith angle calculation");
  params.addParam<unsigned int>("month", 6, "Month for solar zenith angle calculation");
  params.addParam<unsigned int>("year", 2010, "Year for solar zenith angle calculation");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addClassDescription(
      "Parses an MCM Facsimile-format mechanism file and sets up the box model ODE system");
  return params;
}

MCMFacsimileAction::MCMFacsimileAction(const InputParameters & params)
  : Action(params),
    _mechanism_file(getParam<std::string>("mechanism_file")),
    _include_transport(getParam<bool>("include_transport"))
{
  // Delegate parsing to MCMFacsimileParser (shared with MCMBoxModel)
  MCMFacsimileParser parser;

  // Resolve photolysis file path (relative to input file directory)
  std::string photo_path = getParam<std::string>("mcm_photolysis_file");
  {
    std::ifstream test_file(photo_path);
    if (!test_file.good())
    {
      auto input_files = _app.getInputFileNames();
      for (auto & input_file : input_files)
      {
        auto pos = input_file.find_last_of("/\\");
        if (pos != std::string::npos)
        {
          std::string resolved = input_file.substr(0, pos) + "/" + photo_path;
          test_file.open(resolved);
          if (test_file.good()) { photo_path = resolved; break; }
        }
      }
      // Fallback: look in mechanism file's directory
      if (!test_file.good())
      {
        auto pos = _mechanism_file.find_last_of("/\\");
        if (pos != std::string::npos)
        {
          auto bname_pos = photo_path.find_last_of("/\\");
          std::string base = (bname_pos != std::string::npos) ? photo_path.substr(bname_pos + 1) : photo_path;
          std::string resolved = _mechanism_file.substr(0, pos) + "/" + base;
          test_file.open(resolved);
          if (test_file.good()) photo_path = resolved;
        }
      }
    }
  }

  ParsedMechanism mech = parser.parse(_mechanism_file, photo_path);

  _species = mech.species;
  for (auto & r : mech.reactions)
  {
    Reaction rx;
    rx.rate_expression = r.rate_expression;
    rx.reactants = r.reactants;
    rx.products = r.products;
    _reactions.push_back(rx);
  }

  _stoichiometric_matrix = mech.stoichiometry;

  for (unsigned int i = 0; i < mech.coefficient_names.size(); ++i)
  {
    _rate_coefficients[mech.coefficient_names[i]] = mech.coefficient_expressions[i];
    _converted_coefficients[mech.coefficient_names[i]] = mech.coefficient_expressions[i];
    _coefficient_names.insert(mech.coefficient_names[i]);
  }
  _eval_order = mech.coefficient_names;

  _reaction_rate_expressions.resize(mech.reactions.size());
  for (unsigned int i = 0; i < mech.reactions.size(); ++i)
    _reaction_rate_expressions[i] = mech.reactions[i].rate_expression;

  for (unsigned int i = 0; i < mech.j_numbers.size(); ++i)
  {
    std::string jkey = "J<" + std::to_string(mech.j_numbers[i]) + ">";
    _photolysis_rates[jkey] = 0.0;
    _j_CL[jkey] = mech.j_CL[i];
    _j_CMM[jkey] = mech.j_CMM[i];
    _j_CNN[jkey] = mech.j_CNN[i];
  }

  _base_variables = {"TEMP", "M", "O2", "N2", "H2O"};
  for (auto & [jname, _] : _photolysis_rates)
    _base_variables.insert(jname);

  _console << "MCMFacsimileAction: Parsed " << _species.size() << " species, "
           << _rate_coefficients.size() << " rate coefficients, " << _reactions.size()
           << " reactions, " << _photolysis_rates.size() << " photolysis J<N> references"
           << " from " << _mechanism_file << std::endl;
}

const std::string &
MCMFacsimileAction::getRateCoefficient(const std::string & name) const
{
  static const std::string empty;
  auto it = _rate_coefficients.find(name);
  return (it != _rate_coefficients.end()) ? it->second : empty;
}

std::vector<std::vector<Real>>
MCMFacsimileAction::buildReactantMatrix() const
{
  std::vector<std::vector<Real>> matrix;
  for (size_t i = 0; i < _reactions.size(); ++i)
  {
    std::vector<Real> row;
    for (auto & [coeff, name] : _reactions[i].reactants)
    {
      auto it = std::find(_species.begin(), _species.end(), name);
      if (it != _species.end())
      {
        row.push_back(static_cast<Real>(it - _species.begin()));
        row.push_back(coeff);
      }
      else
        mooseWarning("MCMFacsimileAction: Species '", name, "' not found in species list");
    }
    matrix.push_back(row);
  }
  return matrix;
}

void
MCMFacsimileAction::act()
{
  if (_current_task == "add_variable")
  {
    auto type = AddVariableAction::variableType(FEType(0, LAGRANGE));
    auto var_params = _factory.getValidParams(type);
    for (const auto & sp : _species)
      _problem->addVariable(type, sp, var_params);
    _console << "MCMFacsimileAction: Created " << _species.size() << " nonlinear variables"
             << std::endl;
  }
  else if (_current_task == "add_material")
  {
    auto params = _factory.getValidParams("MCMRatesMaterial");
    params.set<Real>("temperature") = getParam<Real>("temperature");
    params.set<Real>("air_density") = getParam<Real>("air_density");
    params.set<Real>("water_vapor") = getParam<Real>("water_vapor");
    params.set<std::vector<std::string>>("species_list") = _species;
    params.set<std::vector<VariableName>>("species_variables") =
        std::vector<VariableName>(_species.begin(), _species.end());
    params.set<std::vector<std::string>>("reaction_rate_expressions") = _reaction_rate_expressions;
    params.set<std::vector<std::vector<Real>>>("reactant_matrix") = buildReactantMatrix();

    std::vector<std::string> coeff_names, coeff_exprs;
    for (auto & name : _eval_order)
    {
      coeff_names.push_back(name);
      coeff_exprs.push_back(_converted_coefficients[name]);
    }
    params.set<std::vector<std::string>>("coefficient_names") = coeff_names;
    params.set<std::vector<std::string>>("coefficient_expressions") = coeff_exprs;

    std::vector<Real> j_cl_vals, j_cmm_vals, j_cnn_vals;
    {
      std::set<unsigned int> j_nums;
      for (auto & [jname, _] : _photolysis_rates)
      {
        unsigned int num;
        pcrecpp::RE("J<([0-9]+)>").FullMatch(jname, &num);
        j_nums.insert(num);
      }
      for (auto & jn : j_nums)
      {
        std::string jkey = "J<" + std::to_string(jn) + ">";
        j_cl_vals.push_back(_j_CL.count(jkey) ? _j_CL[jkey] : 0.0);
        j_cmm_vals.push_back(_j_CMM.count(jkey) ? _j_CMM[jkey] : 0.0);
        j_cnn_vals.push_back(_j_CNN.count(jkey) ? _j_CNN[jkey] : 0.0);
      }
    }
    params.set<std::vector<Real>>("j_cl_values") = j_cl_vals;
    params.set<std::vector<Real>>("j_cmm_values") = j_cmm_vals;
    params.set<std::vector<Real>>("j_cnn_values") = j_cnn_vals;

    params.set<Real>("latitude") = getParam<Real>("latitude");
    params.set<Real>("longitude") = getParam<Real>("longitude");
    params.set<unsigned int>("day") = getParam<unsigned int>("day");
    params.set<unsigned int>("month") = getParam<unsigned int>("month");
    params.set<unsigned int>("year") = getParam<unsigned int>("year");
    params.set<Real>("jfac") = getParam<Real>("jfac");

    _problem->addMaterial("MCMRatesMaterial", "mcm_rates_material", params);
    _console << "MCMFacsimileAction: Created MCMRatesMaterial with " << _reactions.size()
             << " reactions" << std::endl;
  }
  else if (_current_task == "add_kernel")
  {
    for (unsigned int j = 0; j < _species.size(); ++j)
    {
      auto td_params = _factory.getValidParams("TimeDerivative");
      td_params.set<NonlinearVariableName>("variable") = _species[j];
      _problem->addKernel("TimeDerivative", "td_" + _species[j], td_params);

      std::vector<std::vector<Real>> species_reactants(_species.size());
      for (unsigned int r = 0; r < _reactions.size(); ++r)
        for (auto & [coeff, name] : _reactions[r].reactants)
        {
          auto it = std::find(_species.begin(), _species.end(), name);
          if (it != _species.end())
          {
            unsigned int sidx = it - _species.begin();
            species_reactants[sidx].push_back(static_cast<Real>(r));
            species_reactants[sidx].push_back(coeff);
          }
        }

      auto src_params = _factory.getValidParams("ChemicalSourceKernel");
      src_params.set<NonlinearVariableName>("variable") = _species[j];
      src_params.set<std::vector<Real>>("stoichiometric_row") = _stoichiometric_matrix[j];
      src_params.set<std::vector<VariableName>>("all_species") =
          std::vector<VariableName>(_species.begin(), _species.end());
      src_params.set<std::vector<std::vector<Real>>>("species_reactants") = species_reactants;
      _problem->addKernel("ChemicalSourceKernel", "src_" + _species[j], src_params);
    }
    _console << "MCMFacsimileAction: Created TimeDerivative + ChemicalSourceKernel for "
             << _species.size() << " species" << std::endl;

    if (_include_transport)
    {
      for (unsigned int j = 0; j < _species.size(); ++j)
      {
        bool is_placeholder = true;
        for (unsigned int r = 0; r < _reactions.size(); ++r)
          if (std::abs(_stoichiometric_matrix[j][r]) > 1e-30)
          { is_placeholder = false; break; }
        if (is_placeholder) continue;
        auto diff = _factory.getValidParams("Diffusion");
        diff.set<NonlinearVariableName>("variable") = _species[j];
        _problem->addKernel("Diffusion", "diff_" + _species[j], diff);
      }
    }
  }
}
