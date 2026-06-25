//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryAction.h"
#include "MCMFacsimileParser.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "pcrecpp.h"

#include <fstream>
#include <sstream>
#include <algorithm>

registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryAction, "add_variable");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryAction, "add_user_object");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryAction, "add_material");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryAction, "add_kernel");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryAction, "add_scalar_kernel");

InputParameters
AtmosphericChemistryAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::string>(
      "mechanism_file", "Path to the MCM Facsimile-format mechanism file (.fac)");

  MooseEnum mode_enum("box coupled", "box");
  params.addRequiredParam<MooseEnum>("mode", mode_enum,
      "Simulation mode: 'box' for 0-D ODE ScalarVariable, 'coupled' for FEM transport+chemistry");

  MooseEnum format_enum("MCM_FACSIMILE", "MCM_FACSIMILE");
  params.addParam<MooseEnum>("mechanism_format", format_enum,
      "Chemical mechanism format");

  params.addParam<bool>("include_transport", false,
      "Add Diffusion kernels for non-placeholder species (coupled mode only)");

  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19, "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17,
      "Background water vapor concentration (molecules/cm^3)");
  params.addParam<Real>("press", 0.0,
      "Pressure (mbar).  If >0, M computed dynamically via ideal gas law.");
  params.addParam<std::string>(
      "mcm_photolysis_file",
      "doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat",
      "Path to the MCM photolysis-rates parameter file for SZA-based J calculation.");

  MooseEnum mcm_version("v3.1 v3.2 v3.3.1", "v3.3.1");
  params.addParam<MooseEnum>("mcm_version", mcm_version,
      "Version of the Master Chemical Mechanism (v3.1, v3.2, v3.3.1)");

  params.addParam<Real>("latitude", 51.51, "Latitude in degrees (North positive)");
  params.addParam<Real>("longitude", 0.13, "Longitude in degrees (East positive)");
  params.addParam<unsigned int>("day", 21, "Day of month for solar zenith angle calculation");
  params.addParam<unsigned int>("month", 6, "Month for solar zenith angle calculation");
  params.addParam<unsigned int>("year", 2010, "Year for solar zenith angle calculation");
  MooseEnum photo_scheme("MCM_SZA HYBRID", "MCM_SZA");
  params.addParam<MooseEnum>("photolysis_scheme", photo_scheme,
      "Photolysis scheme: MCM_SZA (empirical CL*cos^CMM*exp(-CNN*sec)), "
      "HYBRID (4D TUV lookup table interpolation)");

  params.addParam<std::string>("hybrid_table_dir", "",
      "Directory containing F0AM Hybrid J-value table files "
      "(table_J<N>.dat, axis_*.dat, index.txt). Required if photolysis_scheme=HYBRID.");

  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addParam<bool>("roof_open", true, "Roof (chamber cover) open. false = CLOSED (all J=0)");

  params.addClassDescription(
      "Unified atmospheric chemistry Action. Supports box mode (0-D ODE via "
      "ScalarVariable + ChemistryODEKernel + MCMBoxModel) and coupled mode "
      "(FEM transport + chemistry via MooseVariableFE + ChemicalSourceKernel + "
      "MCMRatesMaterial). Replaces the deprecated MCMFacsimileAction.");
  return params;
}

AtmosphericChemistryAction::AtmosphericChemistryAction(const InputParameters & params)
  : Action(params),
    _mechanism_file(getParam<std::string>("mechanism_file")),
    _mode(getParam<MooseEnum>("mode")),
    _include_transport(getParam<bool>("include_transport"))
{
  // Reject absolute paths — prevents reading arbitrary system files via
  // malicious mechanism_file / photolysis_file parameters.
  // Relative paths (including "..") are allowed — MOOSE test harnesses
  // routinely use "../../../doc/..." to reference database files.
  if (!_mechanism_file.empty() && _mechanism_file[0] == '/')
    mooseError("AtmosphericChemistry: mechanism_file must be relative, got absolute: ", _mechanism_file);

  // Parse the .fac mechanism file via MCMFacsimileParser (shared with MCMBoxModel)
  MCMFacsimileParser parser;

  std::string mcm_ver = getParam<MooseEnum>("mcm_version");
  parser.setMCMVersion(mcm_ver);

  std::string photo_path = getParam<std::string>("mcm_photolysis_file");
  if (!photo_path.empty() && photo_path[0] == '/')
    mooseError("AtmosphericChemistry: mcm_photolysis_file must be relative, got absolute: ", photo_path);
  std::string peroxy_path =
      "doc/content/modules/atmospheric_chemistry/database/mcm_peroxy_radicals_" + mcm_ver + ".dat";
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
      if (!test_file.good())
      {
        auto pos = _mechanism_file.find_last_of("/\\");
        if (pos != std::string::npos)
        {
          auto bname_pos = photo_path.find_last_of("/\\");
          std::string base = (bname_pos != std::string::npos)
                               ? photo_path.substr(bname_pos + 1)
                               : photo_path;
          std::string resolved = _mechanism_file.substr(0, pos) + "/" + base;
          test_file.open(resolved);
          if (test_file.good()) photo_path = resolved;
        }
      }
    }
  }

  ParsedMechanism mech = parser.parse(_mechanism_file, photo_path, peroxy_path);

  _species = mech.species;
  _ro2_species = mech.ro2_species;
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

  _console << "AtmosphericChemistry: Parsed " << _species.size() << " species, "
           << _rate_coefficients.size() << " rate coefficients, " << _reactions.size()
           << " reactions, " << _photolysis_rates.size() << " photolysis J<N> references"
           << " from " << _mechanism_file << " (mode=" << _mode << ")" << std::endl;
}

const std::string &
AtmosphericChemistryAction::getRateCoefficient(const std::string & name) const
{
  static const std::string empty;
  auto it = _rate_coefficients.find(name);
  return (it != _rate_coefficients.end()) ? it->second : empty;
}

std::vector<std::vector<Real>>
AtmosphericChemistryAction::buildReactantMatrix() const
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
        mooseWarning("AtmosphericChemistry: Species '", name, "' not found in species list");
    }
    matrix.push_back(row);
  }
  return matrix;
}

void
AtmosphericChemistryAction::act()
{
  if (_current_task == "add_variable")
  {
    if (_mode == "box")
      actBoxAddVariable();
    else
      actCoupledAddVariable();
  }
  else if (_current_task == "add_user_object")
  {
    if (_mode == "box")
      actBoxAddUserObject();
  }
  else if (_current_task == "add_scalar_kernel")
  {
    if (_mode == "box")
      actBoxAddScalarKernel();
  }
  else if (_current_task == "add_material")
  {
    if (_mode == "coupled")
      actCoupledAddMaterial();
  }
  else if (_current_task == "add_kernel")
  {
    if (_mode == "coupled")
      actCoupledAddKernel();
  }
}

// ===== Box mode tasks =====

void
AtmosphericChemistryAction::actBoxAddVariable()
{
  auto var_params = _factory.getValidParams("MooseVariableScalar");
  for (const auto & sp : _species)
    _problem->addVariable("MooseVariableScalar", sp, var_params);
  _console << "AtmosphericChemistry (box): Created " << _species.size()
           << " ScalarVariable(s)" << std::endl;
}

void
AtmosphericChemistryAction::actBoxAddUserObject()
{
  auto params = _factory.getValidParams("MCMBoxModel");
  params.set<std::string>("mechanism_file") = _mechanism_file;
  params.set<Real>("temperature") = getParam<Real>("temperature");
  params.set<Real>("air_density") = getParam<Real>("air_density");
  params.set<Real>("water_vapor") = getParam<Real>("water_vapor");
  params.set<Real>("press") = getParam<Real>("press");
  params.set<Real>("latitude") = getParam<Real>("latitude");
  params.set<Real>("longitude") = getParam<Real>("longitude");
  params.set<unsigned int>("day") = getParam<unsigned int>("day");
  params.set<unsigned int>("month") = getParam<unsigned int>("month");
  params.set<unsigned int>("year") = getParam<unsigned int>("year");
  params.set<Real>("jfac") = getParam<Real>("jfac");
  params.set<std::string>("photolysis_file") = getParam<std::string>("mcm_photolysis_file");
  _problem->addUserObject("MCMBoxModel", "box_model", params);
  _console << "AtmosphericChemistry (box): Created MCMBoxModel UserObject" << std::endl;
}

void
AtmosphericChemistryAction::actBoxAddScalarKernel()
{
  // Build the species_variables list for ChemistryODEKernel
  std::vector<VariableName> species_vars(_species.begin(), _species.end());

  for (unsigned int j = 0; j < _species.size(); ++j)
  {
    // ODETimeDerivative: contributes du/dt to the residual
    auto td_params = _factory.getValidParams("ODETimeDerivative");
    td_params.set<NonlinearVariableName>("variable") = _species[j];
    _problem->addScalarKernel("ODETimeDerivative", "td_" + _species[j], td_params);

    // ChemistryODEKernel: contributes -dC/dt (chemical source)
    auto chem_params = _factory.getValidParams("ChemistryODEKernel");
    chem_params.set<NonlinearVariableName>("variable") = _species[j];
    chem_params.set<UserObjectName>("box_model") = "box_model";
    chem_params.set<unsigned int>("species_index") = j;
    chem_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("ChemistryODEKernel", "chem_" + _species[j], chem_params);
  }
  _console << "AtmosphericChemistry (box): Created ODETimeDerivative + ChemistryODEKernel for "
           << _species.size() << " species" << std::endl;
}

// ===== Coupled mode tasks (equivalent to old MCMFacsimileAction) =====

void
AtmosphericChemistryAction::actCoupledAddVariable()
{
  auto type = AddVariableAction::variableType(FEType(0, LAGRANGE));
  auto var_params = _factory.getValidParams(type);
  for (const auto & sp : _species)
    _problem->addVariable(type, sp, var_params);
  _console << "AtmosphericChemistry (coupled): Created " << _species.size()
           << " nonlinear variable(s)" << std::endl;
}

void
AtmosphericChemistryAction::actCoupledAddMaterial()
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
  std::vector<unsigned int> j_numbers_all;
  {
    // Load ALL photolysis parameters from the MCM photolysis-rates file.
    // The parser only transfers mechanism-referenced J<N>, but we need the
    // full set for MCMPhotolysisPostprocessor output (e.g. J11-J61).
    // Re-read the file here to get every J<N> entry.
    std::string photo_path = getParam<std::string>("mcm_photolysis_file");
    // Resolve relative paths (same logic as in the constructor)
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
      }
    }
    std::ifstream pfile(photo_path);
    if (pfile.good())
    {
      std::string line;
      std::getline(pfile, line); // skip header
      while (std::getline(pfile, line))
      {
        if (line.empty() || line[0] == '#') continue;
        // Convert Fortran D-notation (6.073D-05) → E-notation (6.073E-05)
        std::replace(line.begin(), line.end(), 'D', 'E');
        std::replace(line.begin(), line.end(), 'd', 'e');
        // Parse "j l m n name tau" columns (we only need first 4)
        std::istringstream iss(line);
        unsigned int jn;
        double cl, cmm, cnn;
        std::string unused1, unused2;
        if (iss >> jn >> cl >> cmm >> cnn)
        {
          j_numbers_all.push_back(jn);
          j_cl_vals.push_back(cl);
          j_cmm_vals.push_back(cmm);
          j_cnn_vals.push_back(cnn);
        }
      }
    }
  }
  params.set<std::vector<unsigned int>>("j_numbers") = j_numbers_all;
  params.set<std::vector<Real>>("j_cl_values") = j_cl_vals;
  params.set<std::vector<Real>>("j_cmm_values") = j_cmm_vals;
  params.set<std::vector<Real>>("j_cnn_values") = j_cnn_vals;

  // Pass RO2 lumped-species list for fparser variable setup (box mode
  // equivalent: MCMBoxModel::setupFparser adds RO2 as a derived variable)
  if (!_ro2_species.empty())
    params.set<std::vector<std::string>>("ro2_species") = _ro2_species;

  params.set<Real>("latitude") = getParam<Real>("latitude");
  params.set<Real>("longitude") = getParam<Real>("longitude");
  params.set<unsigned int>("day") = getParam<unsigned int>("day");
  params.set<unsigned int>("month") = getParam<unsigned int>("month");
  params.set<unsigned int>("year") = getParam<unsigned int>("year");
  params.set<Real>("jfac") = getParam<Real>("jfac");
  params.set<bool>("roof_open") = getParam<bool>("roof_open");
  params.set<MooseEnum>("photolysis_scheme") = getParam<MooseEnum>("photolysis_scheme");
  params.set<std::string>("hybrid_table_dir") = getParam<std::string>("hybrid_table_dir");

  _problem->addMaterial("MCMRatesMaterial", "mcm_rates_material", params);
  _console << "AtmosphericChemistry (coupled): Created MCMRatesMaterial with "
           << _reactions.size() << " reactions" << std::endl;
}

void
AtmosphericChemistryAction::actCoupledAddKernel()
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
  _console << "AtmosphericChemistry (coupled): Created TimeDerivative + ChemicalSourceKernel for "
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
