//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KPPMechanismMaterial.h"

#include <algorithm>
#include <cmath>
#include <limits>

registerMooseObject("AtmosphericChemistryApp", KPPMechanismMaterial);

InputParameters
KPPMechanismMaterial::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<std::string>("lib_path", "Path to the KPP-generated shared library.");
  params.addRequiredCoupledVar("species_variables",
                               "Chemical species variables in KPP-generated VAR order.");
  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19, "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17,
                        "Background water vapor concentration (molecules/cm^3)");
  params.addParam<Real>("press", 0.0,
                        "Pressure (mbar). If >0, M is computed from the ideal gas law.");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for KPP photolysis rates.");
  params.addParam<bool>("roof_open", true,
                        "Roof (chamber cover) open. false forces KPP SUN to zero.");
  MooseEnum photo_scheme("MCM_SZA HYBRID BOTTOMUP", "MCM_SZA");
  params.addParam<MooseEnum>(
      "photolysis_scheme",
      photo_scheme,
      "Photolysis scheme. KPP coupled mode currently uses BOTTOMUP to set exported J globals.");
  params.addParam<std::string>(
      "lamp_flux_file",
      "",
      "Lamp/actinic flux file relative to bottomup_data_dir. Required for BOTTOMUP.");
  params.addParam<std::string>(
      "bottomup_data_dir",
      "../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup",
      "Directory containing BottomUp photolysis data files.");

  MooseEnum units_enum("molec_cm3 ppb", "molec_cm3");
  params.addParam<MooseEnum>("units", units_enum,
                             "Concentration units for coupled variables.");

  params.addClassDescription(
      "Evaluates a KPP-generated mechanism RHS and analytical Jacobian for coupled mode.");
  return params;
}

KPPMechanismMaterial::KPPMechanismMaterial(const InputParameters & params)
  : Material(params),
    _mechanism(getParam<std::string>("lib_path")),
    _temperature(getParam<Real>("temperature")),
    _air_density(getParam<Real>("air_density")),
    _water_vapor(getParam<Real>("water_vapor")),
    _pressure(getParam<Real>("press")),
    _jfac(getParam<Real>("jfac")),
    _roof_open(getParam<bool>("roof_open")),
    _units(getParam<MooseEnum>("units")),
    _photolysis_scheme(getParam<MooseEnum>("photolysis_scheme")),
    _lamp_flux_file(getParam<std::string>("lamp_flux_file")),
    _bottomup_data_dir(getParam<std::string>("bottomup_data_dir")),
    _cached_bottomup_temperature(std::numeric_limits<Real>::quiet_NaN()),
    _cached_bottomup_pressure(std::numeric_limits<Real>::quiet_NaN()),
    _bottomup_j_valid(false),
    _kpp_rhs(declareProperty<std::vector<Real>>("kpp_rhs")),
    _kpp_jacobian_dense(declareProperty<std::vector<Real>>("kpp_jacobian_dense"))
{
  _mechanism.setRoofOpen(_roof_open);
  _mechanism.setJFac(_jfac);

  if (_photolysis_scheme != "MCM_SZA" && _photolysis_scheme != "BOTTOMUP")
    mooseError("KPPMechanismMaterial: photolysis_scheme=", _photolysis_scheme,
               " is not supported in KPP coupled mode. Use BOTTOMUP for chamber mechanisms.");

  if (_photolysis_scheme == "BOTTOMUP")
  {
    if (_lamp_flux_file.empty())
      mooseError("KPPMechanismMaterial: lamp_flux_file is required when "
                 "photolysis_scheme=BOTTOMUP.");
    if (!_bottomup_data_dir.empty() && _bottomup_data_dir[0] == '/')
      mooseError("KPPMechanismMaterial: bottomup_data_dir must be relative, got absolute: ",
                 _bottomup_data_dir);

    _bottomup_integrator = std::make_unique<BottomUpJIntegrator>(_bottomup_data_dir);
    _bottomup_integrator->loadLampFlux(_lamp_flux_file);
    _bottomup_integrator->loadReactionMap("bottomup_jmap.dat");
  }

  _n_species = coupledComponents("species_variables");
  if (_n_species != _mechanism.nSpecies())
    mooseError("KPPMechanismMaterial: species_variables has ", _n_species,
               " entries but the KPP library has ", _mechanism.nSpecies(),
               " variable species.");

  _species_vals.resize(_n_species);
  for (const auto i : make_range(_n_species))
    _species_vals[i] = &coupledValue("species_variables", i);
}

Real
KPPMechanismMaterial::airDensity() const
{
  if (_pressure <= 0.0)
    return _air_density;

  return (_pressure * 100.0) / (1.380649e-23 * _temperature) * 1.0e-6;
}

void
KPPMechanismMaterial::applyPhotolysisGlobals()
{
  if (_photolysis_scheme != "BOTTOMUP" || !_bottomup_integrator)
    return;

  const Real current_pressure = _pressure > 0.0 ? _pressure : 1013.25;
  if (!_bottomup_j_valid ||
      std::abs(_temperature - _cached_bottomup_temperature) > 1.0e-12 ||
      std::abs(current_pressure - _cached_bottomup_pressure) > 1.0e-12)
  {
    _cached_bottomup_j = _bottomup_integrator->computeAllJ(_temperature, current_pressure);
    _cached_bottomup_temperature = _temperature;
    _cached_bottomup_pressure = current_pressure;
    _bottomup_j_valid = true;
  }

  const Real chamber_factor = _roof_open ? _jfac : 0.0;
  for (const auto & [jname, value] : _cached_bottomup_j)
    _mechanism.setGlobal(jname, value * chamber_factor);
}

void
KPPMechanismMaterial::computeQpProperties()
{
  const Real current_air_density = airDensity();
  const Real unit_conversion = (_units == "ppb") ? current_air_density / 1.0e9 : 1.0;

  std::vector<Real> concentrations(_n_species);
  for (const auto i : make_range(_n_species))
    concentrations[i] = std::max((*_species_vals[i])[_qp], 0.0) * unit_conversion;

  PhysParams phys;
  phys.temperature = _temperature;
  phys.air_density = _air_density;
  phys.water_vapor = _water_vapor;
  phys.pressure = _pressure;
  phys.jfac = 1.0;

  _kpp_rhs[_qp].assign(_n_species, 0.0);
  applyPhotolysisGlobals();
  _mechanism.computeRHS(_t, concentrations, phys, _kpp_rhs[_qp]);

  std::vector<std::tuple<unsigned int, unsigned int, Real>> jacobian_triplets;
  applyPhotolysisGlobals();
  _mechanism.computeJacobian(_t, concentrations, phys, jacobian_triplets);

  _kpp_jacobian_dense[_qp].assign(_n_species * _n_species, 0.0);
  for (const auto & entry : jacobian_triplets)
  {
    const auto row = std::get<0>(entry);
    const auto col = std::get<1>(entry);
    if (row < _n_species && col < _n_species)
      _kpp_jacobian_dense[_qp][row * _n_species + col] = std::get<2>(entry);
  }
}
