//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryCoupledAction.h"
#include "ChemistryMechanismSpec.h"
#include "KPPGeneratedMechanism.h"
#include "AddVariableAction.h"
#include "FEProblem.h"
#include "MooseEnum.h"

#include <algorithm>
#include <cstdlib>
#include <unordered_map>

registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryCoupledAction, "add_variable");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryCoupledAction, "add_material");
registerMooseAction("AtmosphericChemistryApp", AtmosphericChemistryCoupledAction, "add_kernel");

namespace
{
const char * const transport_velocity_property = "atmospheric_chemistry_transport_velocity";
}

InputParameters
AtmosphericChemistryCoupledAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::string>(
      "mechanism_file", "Path to the mechanism file (.fac for MCM, .kpp for KPP)");

  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19, "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17,
      "Background water vapor concentration (molecules/cm^3)");
  params.addParam<Real>("press", 0.0,
      "Pressure (mbar). If >0, M computed dynamically via ideal gas law.");
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
      "Photolysis scheme: MCM_SZA, HYBRID (4D TUV lookup), or BOTTOMUP");

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

  params.addParam<std::vector<VariableName>>(
      "transport_velocity",
      {},
      "Scalar velocity component variables used to add a ConservativeAdvection kernel for every "
      "species.");
  params.addParam<MaterialPropertyName>(
      "transport_velocity_material",
      "RealVectorValue material property used to add a ConservativeAdvection kernel for every "
      "species.");
  MooseEnum upwinding_type("none full", "none");
  params.addParam<MooseEnum>(
      "advection_upwinding_type", upwinding_type, "Upwinding type for species advection kernels.");
  params.addParam<FunctionName>(
      "diffusivity",
      "Diffusion coefficient/function used to add a FunctionDiffusion kernel for every species. "
      "Numeric constants are accepted.");
  params.addParam<FunctionName>(
      "density_weighted_diffusivity",
      "Turbulent diffusivity K used to add an AtmosphericDensityWeightedDiffusion kernel for every "
      "species.");
  params.addParam<FunctionName>(
      "density_weighted_air_density",
      "Air number-density function rho used with density_weighted_diffusivity.");
  MooseEnum density_weighted_component("x y z", "z");
  params.addParam<MooseEnum>("density_weighted_component",
                             density_weighted_component,
                             "Coordinate direction for density-weighted diffusion.");
  params.addParam<Real>(
      "density_weighted_coordinate_scale",
      1.0,
      "Scale from the density-weighted diffusion mesh coordinate unit to meters. Use 1000 when "
      "height is in km and diffusivity is in m^2/s.");

  MooseEnum solver_enum("moose_implicit kpp_rosenbrock kpp_sdirk kpp_runge_kutta",
                        "moose_implicit");
  params.addParam<MooseEnum>(
      "chem_solver", solver_enum,
      "Chemical mechanism backend. Use moose_implicit for MCM .fac files; KPP values "
      "select a generated KPP shared library for .kpp mechanisms.");

  params.addParam<bool>("output_ro2_sum", false,
      "Create a diagnostic variable 'RO2' = sum of peroxy radical concentrations.");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addParam<bool>("roof_open", true, "Roof (chamber cover) open.");
  params.addParam<bool>("use_limiting_reagent", false,
      "Enable F0AM-style limiting-reagent formulation for RO2+RO2 termination.");

  params.addClassDescription(
      "Action for FEM transport + chemistry (coupled mode). Creates FE variables, "
      "chemistry material properties, source kernels, and optional transport kernels for each "
      "species.");
  return params;
}

AtmosphericChemistryCoupledAction::AtmosphericChemistryCoupledAction(const InputParameters & params)
  : Action(params),
    _chem_solver(getParam<MooseEnum>("chem_solver")),
    _use_kpp(false),
    _ro2_diagnostic_enabled(false)
{
  if (!getParam<std::vector<VariableName>>("transport_velocity").empty() &&
      isParamValid("transport_velocity_material"))
    paramError("transport_velocity",
               "Specify only one of transport_velocity or transport_velocity_material.");
  if (getParam<std::vector<VariableName>>("transport_velocity").size() > LIBMESH_DIM)
    paramError("transport_velocity",
               "The number of velocity component variables cannot exceed LIBMESH_DIM.");
  if (isParamValid("diffusivity") && isParamValid("density_weighted_diffusivity"))
    paramError("diffusivity",
               "Specify only one of diffusivity or density_weighted_diffusivity.");
  if (isParamValid("density_weighted_diffusivity") !=
      isParamValid("density_weighted_air_density"))
    paramError("density_weighted_diffusivity",
               "density_weighted_diffusivity and density_weighted_air_density must be specified "
               "together.");
  if (getParam<Real>("density_weighted_coordinate_scale") <= 0.0)
    paramError("density_weighted_coordinate_scale",
               "The density-weighted coordinate scale must be positive.");

  std::string mech_file = getParam<std::string>("mechanism_file");
  if (!mech_file.empty() && mech_file[0] == '/')
    mooseError("AtmosphericChemistryCoupled: mechanism_file must be relative, got absolute: ", mech_file);

  std::string mcm_ver = getParam<MooseEnum>("mcm_version");
  std::string photo_path = getParam<std::string>("mcm_photolysis_file");
  if (!photo_path.empty() && photo_path[0] == '/')
    mooseError("AtmosphericChemistryCoupled: mcm_photolysis_file must be relative, got absolute: ", photo_path);
  std::string peroxy_path =
      "doc/content/modules/atmospheric_chemistry/database/mcm_peroxy_radicals_" + mcm_ver + ".dat";

  _use_kpp = (_chem_solver == "kpp_rosenbrock" || _chem_solver == "kpp_sdirk" ||
              _chem_solver == "kpp_runge_kutta");
  if (!_use_kpp && mech_file.size() > 4 &&
      mech_file.substr(mech_file.size() - 4) == ".kpp")
    _use_kpp = true;

  if (_use_kpp)
  {
    KPPGeneratedMechanism mechanism(kppLibraryPath(mech_file));
    _species = mechanism.speciesNames();
    _ro2_species.clear();
  }
  else
  {
    auto input_files = _app.getInputFileNames();
    ChemistryMechanismSpec spec(mech_file, _chem_solver, mcm_ver,
                                 photo_path, peroxy_path, input_files);
    _species = spec.species();
    _ro2_species = spec.ro2Species();
    _mech_data = spec.mechanismData();
  }

  // RO2 diagnostic
  bool want_ro2 = getParam<bool>("output_ro2_sum");
  bool has_ro2 = std::find(_species.begin(), _species.end(), "RO2") != _species.end();
  _ro2_diagnostic_enabled = want_ro2 && !has_ro2;
  if (_use_kpp && _ro2_diagnostic_enabled)
  {
    mooseWarning("AtmosphericChemistryCoupled: output_ro2_sum is not available for KPP "
                 "mechanisms because RO2 metadata is not exported by the KPP shared library.");
    _ro2_diagnostic_enabled = false;
  }
}

void
AtmosphericChemistryCoupledAction::act()
{
  const std::string & task = _current_task;

  if (task == "add_variable")
    actAddVariable();
  else if (task == "add_material")
    actAddMaterial();
  else if (task == "add_kernel")
    actAddKernel();
}

std::string
AtmosphericChemistryCoupledAction::kppLibraryPath(const std::string & mechanism_file) const
{
  const char * kpp_lib_env = std::getenv("KPP_LIB");
  if (kpp_lib_env)
    return kpp_lib_env;

  if (mechanism_file.empty())
    return "libkpp_generated.so";

  auto slash = mechanism_file.find_last_of('/');
  std::string dir = (slash != std::string::npos) ? mechanism_file.substr(0, slash + 1) : "";
  std::string base = mechanism_file.substr(slash + 1);
  auto dot = base.find_last_of('.');
  if (dot != std::string::npos)
    base = base.substr(0, dot);

  return dir + "kpp_build_" + base + "/libkpp_" + base + ".so";
}

bool
AtmosphericChemistryCoupledAction::hasAdvection() const
{
  return !getParam<std::vector<VariableName>>("transport_velocity").empty() ||
         isParamValid("transport_velocity_material");
}

bool
AtmosphericChemistryCoupledAction::hasDiffusion() const
{
  return isParamValid("diffusivity");
}

bool
AtmosphericChemistryCoupledAction::hasDensityWeightedDiffusion() const
{
  return isParamValid("density_weighted_diffusivity");
}

void
AtmosphericChemistryCoupledAction::addTransportKernels(const std::string & species_name)
{
  if (hasAdvection())
  {
    auto adv_params = _factory.getValidParams("ConservativeAdvection");
    adv_params.set<NonlinearVariableName>("variable") = species_name;
    adv_params.set<MooseEnum>("upwinding_type") = getParam<MooseEnum>("advection_upwinding_type");

    const auto & velocity = getParam<std::vector<VariableName>>("transport_velocity");
    if (!velocity.empty())
      adv_params.set<MaterialPropertyName>("velocity_material") = transport_velocity_property;
    else
      adv_params.set<MaterialPropertyName>("velocity_material") =
          getParam<MaterialPropertyName>("transport_velocity_material");

    _problem->addKernel("ConservativeAdvection", "advect_" + species_name, adv_params);
  }

  if (hasDiffusion())
  {
    auto diff_params = _factory.getValidParams("FunctionDiffusion");
    diff_params.set<NonlinearVariableName>("variable") = species_name;
    diff_params.set<FunctionName>("function") = getParam<FunctionName>("diffusivity");
    _problem->addKernel("FunctionDiffusion", "diff_" + species_name, diff_params);
  }

  if (hasDensityWeightedDiffusion())
  {
    auto diff_params = _factory.getValidParams("AtmosphericDensityWeightedDiffusion");
    diff_params.set<NonlinearVariableName>("variable") = species_name;
    diff_params.set<FunctionName>("diffusivity") =
        getParam<FunctionName>("density_weighted_diffusivity");
    diff_params.set<FunctionName>("density") =
        getParam<FunctionName>("density_weighted_air_density");
    diff_params.set<MooseEnum>("component") = getParam<MooseEnum>("density_weighted_component");
    diff_params.set<Real>("coordinate_scale") =
        getParam<Real>("density_weighted_coordinate_scale");
    _problem->addKernel(
        "AtmosphericDensityWeightedDiffusion", "rho_diff_" + species_name, diff_params);
  }
}

void
AtmosphericChemistryCoupledAction::actAddVariable()
{
  // Create FE (LAGRANGE) variables for all species
  auto type = AddVariableAction::variableType(FEType(0, LAGRANGE));
  auto var_params = _factory.getValidParams(type);
  for (const auto & sp : _species)
    _problem->addVariable(type, sp, var_params);

  // RO2 AuxVariable (constant monomial)
  if (_ro2_diagnostic_enabled)
  {
    auto aux_params = _factory.getValidParams("MooseVariable");
    aux_params.set<MooseEnum>("family") = MooseEnum("LAGRANGE MONOMIAL", "MONOMIAL");
    aux_params.set<MooseEnum>("order") = MooseEnum("CONSTANT");
    _problem->addAuxVariable("MooseVariable", "RO2", aux_params);
  }

  _console << "AtmosphericChemistryCoupled: Created " << _species.size()
           << " FE variable(s)" << std::endl;
}

std::vector<std::vector<Real>>
AtmosphericChemistryCoupledAction::buildReactantMatrix() const
{
  std::map<std::string, unsigned int> sp_idx;
  for (unsigned int i = 0; i < _species.size(); ++i)
    sp_idx[_species[i]] = i;

  std::vector<std::vector<Real>> mat(_mech_data.reactions.size());
  for (unsigned int r = 0; r < _mech_data.reactions.size(); ++r)
  {
    for (auto & [coeff, name] : _mech_data.reactions[r].reactants)
    {
      auto it = sp_idx.find(name);
      if (it != sp_idx.end())
      {
        mat[r].push_back(static_cast<Real>(it->second));
        mat[r].push_back(coeff);
      }
    }
  }
  return mat;
}

void
AtmosphericChemistryCoupledAction::actAddMaterial()
{
  const auto & velocity = getParam<std::vector<VariableName>>("transport_velocity");
  if (!velocity.empty())
  {
    auto velocity_params = _factory.getValidParams("VectorFromComponentVariablesMaterial");
    velocity_params.set<MaterialPropertyName>("vector_prop_name") = transport_velocity_property;
    velocity_params.set<std::vector<VariableName>>("u") = {velocity[0]};
    if (velocity.size() > 1)
      velocity_params.set<std::vector<VariableName>>("v") = {velocity[1]};
    if (velocity.size() > 2)
      velocity_params.set<std::vector<VariableName>>("w") = {velocity[2]};
    _problem->addMaterial("VectorFromComponentVariablesMaterial",
                          "atmospheric_chemistry_transport_velocity",
                          velocity_params);
  }

  if (_use_kpp)
  {
    auto params = _factory.getValidParams("KPPMechanismMaterial");
    params.set<std::string>("lib_path") = kppLibraryPath(getParam<std::string>("mechanism_file"));
    params.set<std::vector<VariableName>>("species_variables") =
        std::vector<VariableName>(_species.begin(), _species.end());
    params.set<Real>("temperature") = getParam<Real>("temperature");
    params.set<Real>("air_density") = getParam<Real>("air_density");
    params.set<Real>("water_vapor") = getParam<Real>("water_vapor");
    params.set<Real>("press") = getParam<Real>("press");
    params.set<Real>("jfac") = getParam<Real>("jfac");
    params.set<bool>("roof_open") = getParam<bool>("roof_open");
    params.set<MooseEnum>("units") = getParam<MooseEnum>("units");
    params.set<MooseEnum>("photolysis_scheme") = getParam<MooseEnum>("photolysis_scheme");
    params.set<std::string>("lamp_flux_file") = getParam<std::string>("lamp_flux_file");
    params.set<std::string>("bottomup_data_dir") = getParam<std::string>("bottomup_data_dir");

    _problem->addMaterial("KPPMechanismMaterial", "kpp_mechanism_material", params);
    _console << "AtmosphericChemistryCoupled: Created KPPMechanismMaterial with "
             << _species.size() << " variable species" << std::endl;
    return;
  }

  auto params = _factory.getValidParams("MCMRatesMaterial");
  params.set<Real>("temperature") = getParam<Real>("temperature");
  params.set<Real>("air_density") = getParam<Real>("air_density");
  params.set<Real>("water_vapor") = getParam<Real>("water_vapor");
  params.set<std::vector<std::string>>("species_list") = _species;
  params.set<std::vector<VariableName>>("species_variables") =
      std::vector<VariableName>(_species.begin(), _species.end());
  params.set<std::vector<std::string>>("reaction_rate_expressions") =
      _mech_data.reaction_rate_expressions;
  params.set<std::vector<std::vector<Real>>>("reactant_matrix") = buildReactantMatrix();

  std::vector<std::string> coeff_names, coeff_exprs;
  for (auto & name : _mech_data.eval_order)
  {
    coeff_names.push_back(name);
    coeff_exprs.push_back(_mech_data.converted_coefficients.at(name));
  }
  params.set<std::vector<std::string>>("coefficient_names") = coeff_names;
  params.set<std::vector<std::string>>("coefficient_expressions") = coeff_exprs;

  // J-value parameters
  std::vector<unsigned int> j_numbers_all;
  std::vector<Real> j_cl_vals, j_cmm_vals, j_cnn_vals;
  {
    auto scheme = getParam<MooseEnum>("photolysis_scheme");
    if (scheme != "BOTTOMUP")
    {
      j_numbers_all = _mech_data.j_numbers_all;
      j_cl_vals = _mech_data.j_cl_values;
      j_cmm_vals = _mech_data.j_cmm_values;
      j_cnn_vals = _mech_data.j_cnn_values;
    }
  }
  params.set<std::vector<unsigned int>>("j_numbers") = j_numbers_all;
  params.set<std::vector<Real>>("j_cl_values") = j_cl_vals;
  params.set<std::vector<Real>>("j_cmm_values") = j_cmm_vals;
  params.set<std::vector<Real>>("j_cnn_values") = j_cnn_vals;

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
  _console << "AtmosphericChemistryCoupled: Created MCMRatesMaterial with "
           << _mech_data.reactions.size() << " reactions" << std::endl;
}

void
AtmosphericChemistryCoupledAction::actAddKernel()
{
  if (_use_kpp)
  {
    auto u = getParam<MooseEnum>("units");
    Real M = getParam<Real>("air_density");
    const Real pressure = getParam<Real>("press");
    if (pressure > 0.0)
      M = (pressure * 100.0) / (1.380649e-23 * getParam<Real>("temperature")) * 1.0e-6;
    const Real unit_conversion = (u == "ppb") ? M / 1.0e9 : 1.0;
    std::vector<VariableName> all_species(_species.begin(), _species.end());

    for (const auto j : make_range(_species.size()))
    {
      auto td_params = _factory.getValidParams("TimeDerivative");
      td_params.set<NonlinearVariableName>("variable") = _species[j];
      _problem->addKernel("TimeDerivative", "td_" + _species[j], td_params);

      auto src_params = _factory.getValidParams("KPPChemicalSourceKernel");
      src_params.set<NonlinearVariableName>("variable") = _species[j];
      src_params.set<unsigned int>("species_index") = j;
      src_params.set<std::vector<VariableName>>("all_species") = all_species;
      src_params.set<Real>("unit_conversion") = unit_conversion;
      _problem->addKernel("KPPChemicalSourceKernel", "src_" + _species[j], src_params);

      addTransportKernels(_species[j]);
    }

    _console << "AtmosphericChemistryCoupled: Created TimeDerivative + KPPChemicalSourceKernel for "
             << _species.size() << " species" << std::endl;
    return;
  }

  // Build species index mapping
  std::unordered_map<std::string, unsigned int> species_name_to_idx;
  species_name_to_idx.reserve(_species.size());
  for (const auto i : make_range(_species.size()))
    species_name_to_idx[_species[i]] = i;

  // Build species_reactants matrix
  std::vector<std::vector<Real>> species_reactants(_species.size());
  for (const auto r : make_range(_mech_data.reactions.size()))
    for (auto & [coeff, name] : _mech_data.reactions[r].reactants)
    {
      auto it = species_name_to_idx.find(name);
      if (it != species_name_to_idx.end())
      {
        unsigned int sidx = it->second;
        species_reactants[sidx].push_back(static_cast<Real>(r));
        species_reactants[sidx].push_back(coeff);
      }
    }

  // Unit conversion
  auto u = getParam<MooseEnum>("units");
  Real M = getParam<Real>("air_density");
  Real unit_conversion = (u == "ppb") ? M / 1.0e9 : 1.0;
  std::vector<VariableName> all_species(_species.begin(), _species.end());

  for (const auto j : make_range(_species.size()))
  {
    auto td_params = _factory.getValidParams("TimeDerivative");
    td_params.set<NonlinearVariableName>("variable") = _species[j];
    _problem->addKernel("TimeDerivative", "td_" + _species[j], td_params);

    auto src_params = _factory.getValidParams("ChemicalSourceKernel");
    src_params.set<NonlinearVariableName>("variable") = _species[j];
    src_params.set<std::vector<Real>>("stoichiometric_row") = _mech_data.stoichiometric_matrix[j];
    src_params.set<std::vector<VariableName>>("all_species") = all_species;
    src_params.set<std::vector<std::vector<Real>>>("species_reactants") = species_reactants;
    src_params.set<Real>("unit_conversion") = unit_conversion;
    _problem->addKernel("ChemicalSourceKernel", "src_" + _species[j], src_params);

    addTransportKernels(_species[j]);
  }

  _console << "AtmosphericChemistryCoupled: Created TimeDerivative + ChemicalSourceKernel for "
           << _species.size() << " species" << std::endl;

  // RO2 AuxKernel
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
