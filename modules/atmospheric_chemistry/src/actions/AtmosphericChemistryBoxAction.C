//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryBoxAction.h"
#include "ChemistryMechanismSpec.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "MooseEnum.h"
#include "libmesh/coupling_matrix.h"

#include <algorithm>
#include <map>

#include <algorithm>

registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryBoxAction, "add_variable");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryBoxAction, "add_user_object");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryBoxAction, "add_scalar_kernel");

InputParameters
AtmosphericChemistryBoxAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::string>(
      "mechanism_file", "Path to the mechanism file (.fac for MCM, .kpp for KPP)");

  // Mechanism format auto-detected from chem_solver and file extension.

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
  MooseEnum photo_scheme("MCM_SZA HYBRID BOTTOMUP", "MCM_SZA");
  params.addParam<MooseEnum>("photolysis_scheme", photo_scheme,
      "Photolysis scheme: MCM_SZA (empirical), HYBRID (4D TUV lookup), or "
      "BOTTOMUP (cross-section x QY x lamp-flux)");

  params.addParam<std::string>("hybrid_table_dir", "",
      "Directory containing F0AM Hybrid J-value table files");
  params.addParam<std::string>("lamp_flux_file", "",
      "Path to lamp/actinic flux file (required if photolysis_scheme=BOTTOMUP)");
  params.addParam<std::string>("bottomup_data_dir",
      "../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup",
      "Directory containing BottomUp photolysis data files");
  params.addParam<Real>("albedo", 0.1, "Surface albedo (0-1), used by HYBRID scheme");
  params.addParam<Real>("o3column", 350.0, "O3 column in Dobson Units");
  params.addParam<Real>("altitude", 0.0, "Altitude in meters, used by HYBRID scheme");

  MooseEnum units_enum("molec_cm3 ppb", "molec_cm3");
  params.addParam<MooseEnum>("units", units_enum,
      "Concentration units for input/output: 'molec_cm3' (default) or 'ppb'.");

  params.addParam<bool>("output_ro2_sum", false,
      "Create a diagnostic variable 'RO2' = sum of peroxy radical concentrations.");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addParam<Real>("default_ic", 0.0,
      "Default initial concentration (molec/cm3) for species without explicit ICs.");
  params.addParam<bool>("roof_open", true, "Roof (chamber cover) open.");
  params.addParam<bool>("use_limiting_reagent", false,
      "Enable F0AM-style limiting-reagent formulation for RO2+RO2 termination.");

  // --- Chemical solver selection ---
  MooseEnum solver_enum(
      "moose_implicit petsc_ts sundials kpp_rosenbrock kpp_sdirk kpp_runge_kutta",
      "moose_implicit");
  params.addParam<MooseEnum>("chem_solver", solver_enum,
      "Chemical ODE solver backend:\n"
      "  moose_implicit  — MOOSE Newton solver\n"
      "  petsc_ts        — PETSc TS (BDF/ARKIMEX)\n"
      "  sundials        — SUNDIALS CVODE\n"
      "  kpp_rosenbrock  — KPP Rosenbrock shared-library solver\n"
      "  kpp_sdirk       — KPP SDIRK shared-library solver\n"
      "  kpp_runge_kutta — KPP implicit Runge-Kutta shared-library solver");

  params.addParam<Real>("chem_solver_rtol", 1e-6,
      "Relative tolerance for the chemical ODE solver.");
  params.addParam<Real>("chem_solver_atol", 1e-10,
      "Absolute tolerance for the chemical ODE solver.");
  MooseEnum chem_type_enum(
      "bdf arkimex eimex rosw mimex beuler cn rk theta ssp", "bdf");
  params.addParam<MooseEnum>("chem_solver_type", chem_type_enum,
      "ODE solver type (petsc_ts only): 'bdf' (default), 'arkimex', etc.");

  // --- Family conservation (F0AM DAE method) ---
  params.addParam<std::vector<std::string>>("family_names", {},
      "Names of chemical families for DAE conservation (e.g. 'NOx', 'Ox').");
  params.addParam<std::vector<std::vector<std::string>>>(
      "family_members", {},
      "Member species for each family. First member is the DAE slack variable.");
  params.addParam<std::vector<std::vector<Real>>>("family_scaling", {},
      "Scaling/weighting factors for each family member.");

  params.addClassDescription(
      "Action for 0-D ODE box-model atmospheric chemistry. Creates scalar "
      "variables, ODE kernels, and an MCMBoxModel UserObject.");
  return params;
}

AtmosphericChemistryBoxAction::AtmosphericChemistryBoxAction(const InputParameters & params)
  : Action(params),
    _chem_solver(getParam<MooseEnum>("chem_solver")),
    _use_box_solver(false),
    _ro2_diagnostic_enabled(false)
{
  // ---- Validate mechanism file path ----
  std::string mech_file = getParam<std::string>("mechanism_file");
  if (!mech_file.empty() && mech_file[0] == '/')
    mooseError("AtmosphericChemistryBox: mechanism_file must be relative, got absolute: ", mech_file);

  // ---- Parse mechanism ----
  std::string mcm_ver = getParam<MooseEnum>("mcm_version");
  std::string photo_path = getParam<std::string>("mcm_photolysis_file");
  if (!photo_path.empty() && photo_path[0] == '/')
    mooseError("AtmosphericChemistryBox: mcm_photolysis_file must be relative, got absolute: ", photo_path);
  std::string peroxy_path =
      "doc/content/modules/atmospheric_chemistry/database/mcm_peroxy_radicals_" + mcm_ver + ".dat";

  auto input_files = _app.getInputFileNames();
  ChemistryMechanismSpec spec(mech_file, _chem_solver, mcm_ver,
                               photo_path, peroxy_path, input_files);

  _species = spec.species();
  _mech_data = spec.mechanismData();

  // ---- Derive use_box_solver ----
  _use_box_solver = (_chem_solver != "moose_implicit");

  // ---- RO2 diagnostic ----
  bool want_ro2 = getParam<bool>("output_ro2_sum");
  bool has_ro2 = std::find(_species.begin(), _species.end(), "RO2") != _species.end();
  _ro2_diagnostic_enabled = want_ro2 && !has_ro2;

  // ---- Family conservation ----
  _family_names = getParam<std::vector<std::string>>("family_names");
  if (!_family_names.empty())
  {
    _family_members = getParam<std::vector<std::vector<std::string>>>("family_members");
    _family_scaling = getParam<std::vector<std::vector<Real>>>("family_scaling");
    if (_family_members.size() != _family_names.size())
      mooseError("AtmosphericChemistryBox: family_names (", _family_names.size(),
                 ") and family_members (", _family_members.size(), ") must have same length");
    if (_family_scaling.size() != _family_names.size())
      mooseError("AtmosphericChemistryBox: family_names (", _family_names.size(),
                 ") and family_scaling (", _family_scaling.size(), ") must have same length");
    for (unsigned int i = 0; i < _family_names.size(); ++i)
      if (_family_members[i].size() != _family_scaling[i].size())
        mooseError("AtmosphericChemistryBox: Family '", _family_names[i],
                   "' has ", _family_members[i].size(),
                   " members but ", _family_scaling[i].size(), " scalings.");
  }
}

void
AtmosphericChemistryBoxAction::act()
{
  // Delegate to task-specific methods based on which MOOSE task triggered this call
  const std::string & task = _current_task;

  if (task == "add_variable")
    actAddVariable();
  else if (task == "add_user_object")
    actAddUserObject();
  else if (task == "add_scalar_kernel")
    actAddScalarKernel();
}

void
AtmosphericChemistryBoxAction::actAddVariable()
{
  // Create scalar variables for all species
  auto var_params = _factory.getValidParams("MooseVariableScalar");
  for (const auto & sp : _species)
    _problem->addVariable("MooseVariableScalar", sp, var_params);

  // RO2 diagnostic variable
  if (_ro2_diagnostic_enabled)
    _problem->addVariable("MooseVariableScalar", "RO2", var_params);

  _console << "AtmosphericChemistryBox: Created " << _species.size()
           << " scalar variable(s)" << std::endl;
}

void
AtmosphericChemistryBoxAction::actAddUserObject()
{
  // Create MCMBoxModel
  auto uo_params = _factory.getValidParams("MCMBoxModel");
  uo_params.set<std::string>("mechanism_file") = getParam<std::string>("mechanism_file");
  uo_params.set<Real>("temperature") = getParam<Real>("temperature");
  uo_params.set<Real>("air_density") = getParam<Real>("air_density");
  uo_params.set<Real>("water_vapor") = getParam<Real>("water_vapor");
  uo_params.set<Real>("press") = getParam<Real>("press");
  uo_params.set<Real>("latitude") = getParam<Real>("latitude");
  uo_params.set<Real>("longitude") = getParam<Real>("longitude");
  uo_params.set<unsigned int>("day") = getParam<unsigned int>("day");
  uo_params.set<unsigned int>("month") = getParam<unsigned int>("month");
  uo_params.set<unsigned int>("year") = getParam<unsigned int>("year");
  uo_params.set<MooseEnum>("units") = getParam<MooseEnum>("units");
  uo_params.set<Real>("jfac") = getParam<Real>("jfac");
  uo_params.set<Real>("default_ic") = getParam<Real>("default_ic");
  uo_params.set<bool>("use_limiting_reagent") = getParam<bool>("use_limiting_reagent");
  // Photolysis file: empty for BOTTOMUP (parsed from data files), 
  // otherwise the MCM photolysis rates file (mcm_photolysis_file).
  {
    auto scheme = getParam<MooseEnum>("photolysis_scheme");
    uo_params.set<std::string>("photolysis_file") =
        (scheme == "BOTTOMUP") ? "" : getParam<std::string>("mcm_photolysis_file");
  }
  uo_params.set<MooseEnum>("photolysis_scheme") = getParam<MooseEnum>("photolysis_scheme");
  uo_params.set<std::string>("hybrid_table_dir") = getParam<std::string>("hybrid_table_dir");
  uo_params.set<std::string>("lamp_flux_file") = getParam<std::string>("lamp_flux_file");
  uo_params.set<std::string>("bottomup_data_dir") = getParam<std::string>("bottomup_data_dir");
  uo_params.set<Real>("albedo") = getParam<Real>("albedo");
  uo_params.set<Real>("o3column") = getParam<Real>("o3column");
  uo_params.set<Real>("altitude") = getParam<Real>("altitude");
  // Solver parameters
  uo_params.set<MooseEnum>("chem_solver") = MooseEnum(
      "moose_implicit petsc_ts sundials kpp_rosenbrock kpp_sdirk kpp_runge_kutta", _chem_solver);
  uo_params.set<Real>("solver_rtol") = getParam<Real>("chem_solver_rtol");
  uo_params.set<Real>("solver_atol") = getParam<Real>("chem_solver_atol");
  {
    MooseEnum ts_type("bdf arkimex eimex rosw mimex beuler cn rk theta ssp sundials", "bdf");
    ts_type = getParam<MooseEnum>("chem_solver_type");
    uo_params.set<MooseEnum>("solver_type") = ts_type;
  }

  _problem->addUserObject("MCMBoxModel", "box_model", uo_params);

  _console << "AtmosphericChemistryBox: Created MCMBoxModel 'box_model'" << std::endl;
}

void
AtmosphericChemistryBoxAction::actAddScalarKernel()
{
  // Build the species_variables list
  std::vector<VariableName> species_vars(_species.begin(), _species.end());

  // If families are active, identify slack variable indices to skip
  std::set<std::string> slack_species;
  for (unsigned int f = 0; f < _family_names.size(); ++f)
    if (!_family_members[f].empty())
      slack_species.insert(_family_members[f][0]);

  for (unsigned int j = 0; j < _species.size(); ++j)
  {
    // ODETimeDerivative: only in moose_implicit mode
    if (_chem_solver == "moose_implicit")
    {
      auto td_params = _factory.getValidParams("ODETimeDerivative");
      td_params.set<NonlinearVariableName>("variable") = _species[j];
      _problem->addScalarKernel("ODETimeDerivative", "td_" + _species[j], td_params);
    }

    // ChemistryODEKernel (skip family slack variables)
    if (slack_species.count(_species[j]))
      continue;

    auto chem_params = _factory.getValidParams("ChemistryODEKernel");
    chem_params.set<NonlinearVariableName>("variable") = _species[j];
    chem_params.set<UserObjectName>("box_model") = "box_model";
    chem_params.set<unsigned int>("species_index") = j;
    chem_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("ChemistryODEKernel", "chem_" + _species[j], chem_params);
  }

  // RO2 algebraic kernel
  if (_ro2_diagnostic_enabled)
  {
    auto ro2_params = _factory.getValidParams("MCMRO2Kernel");
    ro2_params.set<NonlinearVariableName>("variable") = "RO2";
    ro2_params.set<UserObjectName>("box_model") = "box_model";
    ro2_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("MCMRO2Kernel", "ro2_kernel", ro2_params);
  }

  // Logging
  {
    std::string mode_label =
        (_chem_solver == "moose_implicit") ? "MOOSE implicit"
      : (_chem_solver == "sundials")      ? "SUNDIALS CVODE"
      : (_chem_solver.find("kpp_") == 0)  ? "KPP " + _chem_solver
      : "PETSc TS";
    _console << "AtmosphericChemistryBox: " << mode_label << " — "
             << "created ChemistryODEKernel + ODETimeDerivative for "
             << _species.size() << " species" << std::endl;
  }

  // Set sparse Jacobian coupling pattern (MOOSE implicit mode only)
  if (_chem_solver == "moose_implicit")
  {
    std::map<std::string, unsigned int> sp_idx;
    for (unsigned int i = 0; i < _species.size(); ++i)
      sp_idx[_species[i]] = i;

    auto cm = std::make_unique<libMesh::CouplingMatrix>(_species.size());
    for (unsigned int i = 0; i < _species.size(); ++i)
    {
      (*cm)(i, i) = 1;
      for (unsigned int r = 0; r < _mech_data.reactions.size(); ++r)
      {
        if (std::abs(_mech_data.stoichiometric_matrix[i][r]) < 1.0e-30)
          continue;
        for (const auto & reactant_pair : _mech_data.reactions[r].reactants)
        {
          auto it = sp_idx.find(reactant_pair.second);
          if (it != sp_idx.end())
            (*cm)(i, it->second) = 1;
        }
      }
    }
    _problem->setCouplingMatrix(std::move(cm), 0);
    _console << "AtmosphericChemistryBox: Set sparse Jacobian coupling pattern"
             << " (" << _species.size() << "x" << _species.size() << ")" << std::endl;
  }
}
