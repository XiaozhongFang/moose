//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryAction.h"
#include "MechanismLoader.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "libmesh/coupling_matrix.h"

#include <algorithm>
#include <set>
#include <unordered_map>

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
  MooseEnum photo_scheme("MCM_SZA HYBRID BOTTOMUP", "MCM_SZA");
  params.addParam<MooseEnum>("photolysis_scheme", photo_scheme,
      "Photolysis scheme: MCM_SZA (empirical CL*cos^CMM*exp(-CNN*sec)), "
      "HYBRID (4D TUV lookup table interpolation), or "
      "BOTTOMUP (lab-chamber cross-section × quantum-yield × lamp-flux integration)");

  params.addParam<std::string>("hybrid_table_dir", "",
      "Directory containing F0AM Hybrid J-value table files "
      "(table_J<N>.dat, axis_*.dat, index.txt). Required if photolysis_scheme=HYBRID.");

  params.addParam<std::string>("lamp_flux_file", "",
      "Path to lamp/actinic flux file (relative to bottomup_data_dir). "
      "Required if photolysis_scheme=BOTTOMUP.");

  params.addParam<std::string>("bottomup_data_dir",
      "../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup",
      "Directory containing BottomUp photolysis data files "
      "(CrossSections/, QuantumYields/, bottomup_jmap.dat).");

  params.addParam<Real>("albedo", 0.1, "Surface albedo (0-1), used by HYBRID scheme");
  params.addParam<Real>("o3column", 350.0, "O3 column in Dobson Units, used by HYBRID scheme");
  params.addParam<Real>("altitude", 0.0, "Altitude in meters, used by HYBRID scheme");

  MooseEnum units_enum("molec_cm3 ppb", "molec_cm3");
  params.addParam<MooseEnum>("units", units_enum,
      "Concentration units for input/output: 'molec_cm3' (default) or 'ppb'. "
      "When 'ppb', ICs are in ppb and ChemistryODEKernel converts internally.");

  params.addParam<bool>("output_ro2_sum", false,
      "Create a diagnostic variable 'RO2' = sum of peroxy radical concentrations. "
      "Only effective when no species named 'RO2' exists in the mechanism. "
      "Recommended default: false (use Postprocessors for precise RO2 output).");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addParam<Real>("default_ic", 0.0,
      "Default initial concentration (molec/cm³) for species without explicit ICs.\n"
      "Set to 1e6 for large mechanisms to prevent Jacobian singularities.");
  params.addParam<bool>("roof_open", true, "Roof (chamber cover) open. false = CLOSED (all J=0)");
  params.addParam<bool>("use_limiting_reagent", false,
      "Enable F0AM-style limiting-reagent formulation for RO2+RO2 termination "
      "reactions: rate = k * min([A],[B])² instead of k * [A] * [B]. "
      "Default false (standard MCM chemistry). Set true for F0AM-compatible "
      "RO2 termination or when comparing against F0AM reference outputs.");

  // --- 化学求解器选择（替代 box_solver*，v2.0 后移除旧参数） ---
  MooseEnum solver_enum(
      "moose_implicit petsc_ts sundials kpp_rosenbrock kpp_sdirk kpp_runge_kutta",
      "petsc_ts");
  params.addParam<MooseEnum>("chem_solver", solver_enum,
      "Chemical ODE solver backend for box mode:\n"
      "  moose_implicit  — MOOSE Newton solver owns the integration\n"
      "  petsc_ts        — PETSc TS (BDF/ARKIMEX, default)\n"
      "  sundials        — SUNDIALS CVODE (if compiled with SUNDIALS)\n"
      "  kpp_rosenbrock  — KPP Rosenbrock (if compiled with KPP)\n"
      "  kpp_sdirk       — KPP SDIRK (if compiled with KPP)\n"
      "  kpp_runge_kutta — KPP Runge-Kutta (if compiled with KPP)");

  params.addParam<Real>("chem_solver_rtol", 1e-6,
      "Relative tolerance for the chemical ODE solver's adaptive integrator.\n"
      "Formerly box_solver_rtol.");
  params.addParam<Real>("chem_solver_atol", 1e-10,
      "Absolute tolerance for the chemical ODE solver's adaptive integrator.\n"
      "Formerly box_solver_atol.");
  MooseEnum chem_type_enum(
      "bdf arkimex eimex rosw mimex beuler cn rk theta ssp", "bdf");
  params.addParam<MooseEnum>("chem_solver_type", chem_type_enum,
      "ODE solver type (petsc_ts only): 'bdf' (default), 'arkimex', etc.\n"
      "Ignored for sundials and kpp_* solvers.\n"
      "Formerly box_solver_type.");

  // --- (已废弃) 旧式 box ODE 求解器参数 ---
  MooseEnum ts_type_enum("bdf arkimex eimex rosw mimex beuler cn rk theta ssp sundials", "bdf");
  params.addDeprecatedParam<bool>("box_solver", false,
      "Use 'chem_solver' instead. This parameter is deprecated.",
      "box_solver is deprecated. Use 'chem_solver' instead.");
  params.addDeprecatedParam<MooseEnum>("box_solver_type", ts_type_enum,
      "Use 'chem_solver_type' instead.",
      "box_solver_type is deprecated. Use 'chem_solver_type' instead.");
  params.addDeprecatedParam<Real>("box_solver_rtol", 1e-6,
      "Use 'chem_solver_rtol' instead.",
      "box_solver_rtol is deprecated. Use 'chem_solver_rtol' instead.");
  params.addDeprecatedParam<Real>("box_solver_atol", 1e-10,
      "Use 'chem_solver_atol' instead.",
      "box_solver_atol is deprecated. Use 'chem_solver_atol' instead.");

  // --- Family conservation (F0AM DAE method) ---
  params.addParam<std::vector<std::string>>("family_names", {},
      "Names of chemical families for DAE conservation (e.g. 'NOx', 'Ox'). "
      "If non-empty, family_{members,scaling} must also be provided.");
  params.addParam<std::vector<std::vector<std::string>>>(
      "family_members", {},
      "Member species for each family. "
      "First member is the DAE slack variable (algebraically constrained).");
  params.addParam<std::vector<std::vector<Real>>>("family_scaling", {},
      "Scaling/weighting factors for each family member. "
      "E.g. Ox = O3 + NO2 + 2*NO3 → scaling = [1, 1, 2]");

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
    _include_transport(getParam<bool>("include_transport")),
    _use_box_solver(getParam<bool>("box_solver"))
{
  // Reject absolute paths — prevents reading arbitrary system files via
  // malicious mechanism_file / photolysis_file parameters.
  // Relative paths (including "..") are allowed — MOOSE test harnesses
  // routinely use "../../../doc/..." to reference database files.
  if (!_mechanism_file.empty() && _mechanism_file[0] == '/')
    mooseError("AtmosphericChemistry: mechanism_file must be relative, got absolute: ", _mechanism_file);

  // Validate: box_solver requires box mode
  if (_use_box_solver && _mode != "box")
    mooseError("AtmosphericChemistry: box_solver=true requires mode=box. "
               "Coupled mode cannot use standalone ODE solver integration.");

  // Parse the .fac mechanism file via MechanismLoader (shared with MCMBoxModel)
  std::string mcm_ver = getParam<MooseEnum>("mcm_version");
  std::string photo_path = getParam<std::string>("mcm_photolysis_file");
  if (!photo_path.empty() && photo_path[0] == '/')
    mooseError("AtmosphericChemistry: mcm_photolysis_file must be relative, got absolute: ", photo_path);
  std::string peroxy_path =
      "doc/content/modules/atmospheric_chemistry/database/mcm_peroxy_radicals_" + mcm_ver + ".dat";

  auto input_files = _app.getInputFileNames();
  MechanismData data = MechanismLoader::load(
      _mechanism_file, photo_path, mcm_ver, peroxy_path, input_files);

  // ── Copy from MechanismData to Action members ────────────────────────────
  _species = data.species;
  _ro2_species = data.ro2_species;

  // Convert MechanismData::Reaction → Action::Reaction (same fields)
  for (auto & r : data.reactions)
  {
    Reaction rx;
    rx.rate_expression = r.rate_expression;
    rx.reactants = r.reactants;
    rx.products = r.products;
    _reactions.push_back(rx);
  }
  _stoichiometric_matrix = data.stoichiometric_matrix;

  for (unsigned int i = 0; i < data.eval_order.size(); ++i)
  {
    _rate_coefficients[data.eval_order[i]] = data.rate_coefficients.at(data.eval_order[i]);
    _converted_coefficients[data.eval_order[i]] = data.converted_coefficients.at(data.eval_order[i]);
    _coefficient_names.insert(data.eval_order[i]);
  }
  _eval_order = data.eval_order;
  _reaction_rate_expressions = data.reaction_rate_expressions;

  // Photolysis (mechanism-referenced J<N> only)
  for (unsigned int i = 0; i < data.j_numbers.size(); ++i)
  {
    std::string jkey = "J<" + std::to_string(data.j_numbers[i]) + ">";
    _j_CL[jkey] = data.j_CL[i];
    _j_CMM[jkey] = data.j_CMM[i];
    _j_CNN[jkey] = data.j_CNN[i];
  }

  // Base variables (TEMP, M, O2, N2, H2O + photolysis J<N> references)
  _base_variables = {"TEMP", "M", "O2", "N2", "H2O"};
  for (auto & [jname, _] : _j_CL)
    _base_variables.insert(jname);

  // Resolved photo path + full photolysis set
  _resolved_photo_path = data.resolved_photo_path;
  _j_numbers_all = data.j_numbers_all;
  _j_cl_all = data.j_cl_values;
  _j_cmm_all = data.j_cmm_values;
  _j_cnn_all = data.j_cnn_values;

  _console << "AtmosphericChemistry: Parsed " << _species.size() << " species, "
           << _rate_coefficients.size() << " rate coefficients, " << _reactions.size()
           << " reactions, " << _j_CL.size() << " photolysis J<N> references"
           << " from " << _mechanism_file << " (mode=" << _mode << ")" << std::endl;

  // Compute RO2 diagnostic availability at construction time
  _ro2_warning_printed = false;
  bool want_ro2 = getParam<bool>("output_ro2_sum");
  bool has_ro2 = std::find(_species.begin(), _species.end(), "RO2") != _species.end();
  _ro2_diagnostic_enabled = want_ro2 && !has_ro2;
  if (want_ro2 && has_ro2)
    _console << "AtmosphericChemistry Warning: 'RO2' is already a species "
             << "in the mechanism. Skipping automatic RO2 diagnostic variable.\n"
             << "To output RO2 concentration manually:\n"
             << (_mode == "box"
                 ? "  [Postprocessors]\n"
                   "    [ro2_sum]\n"
                   "      type = MCMRO2Postprocessor\n"
                   "      box_model = box_model\n"
                   "      species_variables = '<all species in order>'\n"
                   "    []\n"
                   "  []\n"
                 : "  [AuxKernels]\n"
                   "    [ro2_aux]\n"
                   "      type = MCMRO2Aux\n"
                   "      variable = RO2_out\n"
                   "      ro2_species = '<auto-detected peroxy radicals>'\n"
                   "    []\n"
                   "  []\n")
             << std::endl;

  // Load family conservation configuration
  _family_names = getParam<std::vector<std::string>>("family_names");
  if (!_family_names.empty())
  {
    _family_members = getParam<std::vector<std::vector<std::string>>>("family_members");
    _family_scaling = getParam<std::vector<std::vector<Real>>>("family_scaling");
    if (_family_members.size() != _family_names.size())
      mooseError("AtmosphericChemistry: family_names (", _family_names.size(),
                 ") and family_members (", _family_members.size(), ") must have same length");
    if (_family_scaling.size() != _family_names.size())
      mooseError("AtmosphericChemistry: family_names (", _family_names.size(),
                 ") and family_scaling (", _family_scaling.size(), ") must have same length");
    for (unsigned int i = 0; i < _family_names.size(); ++i)
    {
      if (_family_members[i].size() != _family_scaling[i].size())
        mooseError("AtmosphericChemistry: Family '", _family_names[i],
                   "' has ", _family_members[i].size(), " members but ",
                   _family_scaling[i].size(), " scaling factors");
    }
    _console << "AtmosphericChemistry: Loaded " << _family_names.size()
             << " family conservation groups: ";
    for (const auto & fn : _family_names)
      _console << fn << " ";
    _console << std::endl;
  }
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
    {
      actBoxAddUserObject();
      if (!_family_names.empty())
        actBoxAddFamilyUO();
    }
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
  // RO2 diagnostic variable (sum of peroxy radicals)
  if (_ro2_diagnostic_enabled)
  {
    _problem->addVariable("MooseVariableScalar", "RO2", var_params);
    _console << "AtmosphericChemistry (box): Created " << _species.size()
             << " ScalarVariable(s) + RO2" << std::endl;
  }
  else
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
  params.set<MooseEnum>("units") = getParam<MooseEnum>("units");
  params.set<Real>("jfac") = getParam<Real>("jfac");
  params.set<Real>("default_ic") = getParam<Real>("default_ic");
  params.set<bool>("use_limiting_reagent") = getParam<bool>("use_limiting_reagent");
  {
    auto scheme = getParam<MooseEnum>("photolysis_scheme");
    params.set<std::string>("photolysis_file") =
        (scheme == "BOTTOMUP") ? "" : getParam<std::string>("mcm_photolysis_file");
  }
  params.set<MooseEnum>("photolysis_scheme") = getParam<MooseEnum>("photolysis_scheme");
  params.set<std::string>("hybrid_table_dir") = getParam<std::string>("hybrid_table_dir");
  params.set<std::string>("lamp_flux_file") = getParam<std::string>("lamp_flux_file");
  params.set<std::string>("bottomup_data_dir") = getParam<std::string>("bottomup_data_dir");
  params.set<Real>("albedo") = getParam<Real>("albedo");
  params.set<Real>("o3column") = getParam<Real>("o3column");
  params.set<Real>("altitude") = getParam<Real>("altitude");
  // Box ODE solver parameters (forwarded to MCMBoxModel)
  params.set<MooseEnum>("box_solver_mode") = _use_box_solver ? MooseEnum("moose_implicit petsc_ts", "petsc_ts") : MooseEnum("moose_implicit petsc_ts", "moose_implicit");
  if (_use_box_solver)
  {
    params.set<MooseEnum>("solver_type") = getParam<MooseEnum>("box_solver_type");
    params.set<Real>("solver_rtol") = getParam<Real>("box_solver_rtol");
    params.set<Real>("solver_atol") = getParam<Real>("box_solver_atol");
  }
  _problem->addUserObject("MCMBoxModel", "box_model", params);
  _console << "AtmosphericChemistry (box): Created MCMBoxModel UserObject" << std::endl;
  if (_use_box_solver)
    _console << "AtmosphericChemistry (box): Box ODE solver enabled "
             << "(type=" << getParam<MooseEnum>("box_solver_type")
             << ", rtol=" << getParam<Real>("box_solver_rtol")
             << ", atol=" << getParam<Real>("box_solver_atol") << ")" << std::endl;
}

void
AtmosphericChemistryAction::actBoxAddFamilyUO()
{
  auto fam_params = _factory.getValidParams("MCMFamilyConstraint");
  fam_params.set<std::vector<std::string>>("family_names") = _family_names;
  fam_params.set<std::vector<std::vector<std::string>>>("family_members") = _family_members;
  fam_params.set<std::vector<std::vector<Real>>>("family_scaling") = _family_scaling;
  fam_params.set<std::vector<std::string>>("species_list") = _species;
  _problem->addUserObject("MCMFamilyConstraint", "family_uo", fam_params);

  // Build species variables list for family kernel
  std::vector<VariableName> species_vars(_species.begin(), _species.end());

  unsigned int fam_count = _family_names.size();
  for (unsigned int f = 0; f < fam_count; ++f)
  {
    const auto & members = _family_members[f];
    if (members.empty()) continue;

    // Find the first member (slack variable) index in _species
    unsigned int slack_idx = 0;
    bool found = false;
    for (unsigned int j = 0; j < _species.size(); ++j)
    {
      if (_species[j] == members[0])
      {
        slack_idx = j;
        found = true;
        break;
      }
    }
    if (!found)
    {
      _console << "AtmosphericChemistry Warning: Family '"
               << _family_names[f] << "' slack species '"
               << members[0] << "' not found in mechanism, skipping"
               << std::endl;
      continue;
    }

    // Create MCMFamilyScalarKernel for the slack variable
    // (replaces ChemistryODEKernel to enforce DAE constraint)
    auto fam_k_params = _factory.getValidParams("MCMFamilyScalarKernel");
    fam_k_params.set<NonlinearVariableName>("variable") = _species[slack_idx];
    fam_k_params.set<UserObjectName>("box_model") = "box_model";
    fam_k_params.set<UserObjectName>("family_uo") = "family_uo";
    fam_k_params.set<unsigned int>("species_index") = slack_idx;
    fam_k_params.set<std::string>("family_name") = _family_names[f];
    fam_k_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("MCMFamilyScalarKernel",
                              "fam_" + _family_names[f] + "_" + _species[slack_idx],
                              fam_k_params);

    _console << "AtmosphericChemistry (box): Family '"
             << _family_names[f] << "' slack = "
             << members[0] << " (idx=" << slack_idx << ")"
             << std::endl;
  }
}

void
AtmosphericChemistryAction::actBoxAddScalarKernel()
{
  // Build the species_variables list for ChemistryODEKernel (includes RO2 for coupling)
  std::vector<VariableName> species_vars(_species.begin(), _species.end());

  // If families are active, identify slack variable indices to skip
  // (MCMFamilyScalarKernel replaces ChemistryODEKernel for slack species)
  std::set<std::string> slack_species;
  for (unsigned int f = 0; f < _family_names.size(); ++f)
    if (!_family_members[f].empty())
      slack_species.insert(_family_members[f][0]); // first member = slack

  for (unsigned int j = 0; j < _species.size(); ++j)
  {
    // ODETimeDerivative: contributes du/dt to the residual.
    // Skipped in PETSc TS mode — the TS handles time integration directly.
    if (!_use_box_solver)
    {
      auto td_params = _factory.getValidParams("ODETimeDerivative");
      td_params.set<NonlinearVariableName>("variable") = _species[j];
      _problem->addScalarKernel("ODETimeDerivative", "td_" + _species[j], td_params);
    }

    // ChemistryODEKernel: contributes -dC/dt (chemical source).
    // Created unconditionally — the BoxIntegrator strategy handles
    // mode-specific behavior (real evaluation in MOOSE implicit mode,
    // zero return in PETSc TS mode).
    // For family slack variables, MCMFamilyScalarKernel is created
    // separately.  Skip ChemistryODEKernel here to avoid double-counting.
    if (slack_species.count(_species[j]))
      continue;

    auto chem_params = _factory.getValidParams("ChemistryODEKernel");
    chem_params.set<NonlinearVariableName>("variable") = _species[j];
    chem_params.set<UserObjectName>("box_model") = "box_model";
    chem_params.set<unsigned int>("species_index") = j;
    chem_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("ChemistryODEKernel", "chem_" + _species[j], chem_params);
  }

  // RO2 algebraic kernel: RO2 = sum(peroxy radicals), no time derivative
  if (_ro2_diagnostic_enabled)
  {
    auto ro2_params = _factory.getValidParams("MCMRO2Kernel");
    ro2_params.set<NonlinearVariableName>("variable") = "RO2";
    ro2_params.set<UserObjectName>("box_model") = "box_model";
    ro2_params.set<std::vector<VariableName>>("species_variables") = species_vars;
    _problem->addScalarKernel("MCMRO2Kernel", "ro2_kernel", ro2_params);
  }

  // ── Logging ──
  {
    std::string mode_label = _use_box_solver ? "PETSc TS" : "MOOSE implicit";
    _console << "AtmosphericChemistry (box): " << mode_label << " — "
             << "created ChemistryODEKernel for " << _species.size() << " species";
    if (_ro2_diagnostic_enabled)
      _console << " + RO2";
    if (!_use_box_solver)
      _console << " + ODETimeDerivative";
    _console << std::endl;
  }

  // ── Set sparse Jacobian coupling pattern (MOOSE implicit mode only) ──
  // PETSc TS mode bypasses MOOSE's nonlinear solver — no coupling matrix needed.
  if (!_use_box_solver)
  {
    std::map<std::string, unsigned int> sp_idx;
    for (unsigned int i = 0; i < _species.size(); ++i)
      sp_idx[_species[i]] = i;

    auto cm = std::make_unique<libMesh::CouplingMatrix>(_species.size());
    for (unsigned int i = 0; i < _species.size(); ++i)
    {
      (*cm)(i, i) = 1;
      for (unsigned int r = 0; r < _reactions.size(); ++r)
      {
        if (std::abs(_stoichiometric_matrix[i][r]) < 1.0e-30)
          continue;
        for (const auto & reactant_pair : _reactions[r].reactants)
        {
          auto it = sp_idx.find(reactant_pair.second);
          if (it != sp_idx.end())
            (*cm)(i, it->second) = 1;
        }
      }
    }
    _problem->setCouplingMatrix(std::move(cm), 0);
    _console << "AtmosphericChemistry (box): Set sparse Jacobian coupling pattern"
             << " (" << _species.size() << "×" << _species.size() << ")" << std::endl;
  }
}

// ===== Coupled mode tasks (equivalent to old MCMFacsimileAction) =====

void
AtmosphericChemistryAction::actCoupledAddVariable()
{
  auto type = AddVariableAction::variableType(FEType(0, LAGRANGE));
  auto var_params = _factory.getValidParams(type);
  for (const auto & sp : _species)
    _problem->addVariable(type, sp, var_params);
  // RO2 AuxVariable: constant monomial (only when no species conflict)
  if (_ro2_diagnostic_enabled)
  {
    auto aux_params = _factory.getValidParams("MooseVariable");
    aux_params.set<MooseEnum>("family") = MooseEnum("LAGRANGE MONOMIAL", "MONOMIAL");
    aux_params.set<MooseEnum>("order") = MooseEnum("FIRST SECOND THIRD FOURTH CONSTANT", "CONSTANT");
    _problem->addAuxVariable("MooseVariable", "RO2", aux_params);
  }
  _console << "AtmosphericChemistry (coupled): Created " << _species.size()
           << " nonlinear variable(s) + RO2" << std::endl;
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

  // Use the pre-resolved photolysis file path and pre-loaded full J<N>
  // parameter set from MechanismLoader.  No need to re-resolve or re-read
  // the photolysis file — the constructor already did this.
  std::vector<unsigned int> j_numbers_all;
  std::vector<Real> j_cl_vals, j_cmm_vals, j_cnn_vals;
  {
    auto scheme = getParam<MooseEnum>("photolysis_scheme");
    if (scheme != "BOTTOMUP")
    {
      j_numbers_all = _j_numbers_all;
      j_cl_vals = _j_cl_all;
      j_cmm_vals = _j_cmm_all;
      j_cnn_vals = _j_cnn_all;
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
  params.set<MooseEnum>("units") = getParam<MooseEnum>("units");
  params.set<bool>("roof_open") = getParam<bool>("roof_open");
  params.set<MooseEnum>("photolysis_scheme") = getParam<MooseEnum>("photolysis_scheme");
  params.set<std::string>("hybrid_table_dir") = getParam<std::string>("hybrid_table_dir");
  params.set<std::string>("lamp_flux_file") = getParam<std::string>("lamp_flux_file");
  params.set<std::string>("bottomup_data_dir") = getParam<std::string>("bottomup_data_dir");

  _problem->addMaterial("MCMRatesMaterial", "mcm_rates_material", params);
  _console << "AtmosphericChemistry (coupled): Created MCMRatesMaterial with "
           << _reactions.size() << " reactions" << std::endl;
}

void
AtmosphericChemistryAction::actCoupledAddKernel()
{
  // ── Precompute species name → index mapping ──
  // O(1) lookup per reactant instead of O(N_species) std::find.
  std::unordered_map<std::string, unsigned int> species_name_to_idx;
  species_name_to_idx.reserve(_species.size());
  for (unsigned int i = 0; i < _species.size(); ++i)
    species_name_to_idx[_species[i]] = i;

  // ── Precompute species_reactants matrix (once, not per species) ──
  // species_reactants[k] = [rxn_0, coeff_0, rxn_1, coeff_1, ...]
  // Lists which reactions have species k as a reactant, used by
  // ChemicalSourceKernel's off-diagonal Jacobian.
  std::vector<std::vector<Real>> species_reactants(_species.size());
  for (unsigned int r = 0; r < _reactions.size(); ++r)
    for (auto & [coeff, name] : _reactions[r].reactants)
    {
      auto it = species_name_to_idx.find(name);
      if (it != species_name_to_idx.end())
      {
        unsigned int sidx = it->second;
        species_reactants[sidx].push_back(static_cast<Real>(r));
        species_reactants[sidx].push_back(coeff);
      }
    }

  // ── Precompute unit conversion and all_species list ──
  auto u = getParam<MooseEnum>("units");
  Real M = getParam<Real>("air_density");
  Real unit_conversion = (u == "ppb") ? M / 1.0e9 : 1.0;
  std::vector<VariableName> all_species(_species.begin(), _species.end());

  for (unsigned int j = 0; j < _species.size(); ++j)
  {
    auto td_params = _factory.getValidParams("TimeDerivative");
    td_params.set<NonlinearVariableName>("variable") = _species[j];
    _problem->addKernel("TimeDerivative", "td_" + _species[j], td_params);

    auto src_params = _factory.getValidParams("ChemicalSourceKernel");
    src_params.set<NonlinearVariableName>("variable") = _species[j];
    src_params.set<std::vector<Real>>("stoichiometric_row") = _stoichiometric_matrix[j];
    src_params.set<std::vector<VariableName>>("all_species") = all_species;
    src_params.set<std::vector<std::vector<Real>>>("species_reactants") = species_reactants;
    src_params.set<Real>("unit_conversion") = unit_conversion;
    _problem->addKernel("ChemicalSourceKernel", "src_" + _species[j], src_params);
  }
  _console << "AtmosphericChemistry (coupled): Created TimeDerivative + ChemicalSourceKernel for "
           << _species.size() << " species" << std::endl;

  // RO2 AuxKernel: sum of peroxy radical species (coupled mode diagnostic)
  if (_ro2_diagnostic_enabled && !_ro2_species.empty())
  {
    std::vector<VariableName> ro2_coupled(_ro2_species.begin(), _ro2_species.end());
    auto ro2_params = _factory.getValidParams("MCMRO2Aux");
    ro2_params.set<AuxVariableName>("variable") = "RO2";
    ro2_params.set<std::vector<VariableName>>("ro2_species") = ro2_coupled;
    _problem->addAuxKernel("MCMRO2Aux", "ro2_aux", ro2_params);
    _console << "  Added RO2 AuxKernel (" << ro2_coupled.size() << " RO2 species)" << std::endl;
  }
}
