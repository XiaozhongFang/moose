//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRatesMaterial.h"
#include "pcrecpp.h"

registerMooseObject("AtmosphericChemistryApp", MCMRatesMaterial);

InputParameters
MCMRatesMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params += FunctionParserUtils<false>::validParams();

  params.addRequiredParam<Real>("temperature", "Ambient temperature (K)");
  params.addRequiredParam<Real>("air_density", "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17, "Background water vapor concentration (molecules/cm^3)");

  params.addRequiredParam<std::vector<std::string>>(
      "species_list", "List of chemical species names");
  params.addRequiredCoupledVar("species_variables",
                               "All chemical species concentrations as coupled variables");

  params.addRequiredParam<std::vector<std::string>>(
      "reaction_rate_expressions",
      "Rate expressions for each reaction (may reference coefficient names)");

  params.addRequiredParam<std::vector<std::vector<Real>>>(
      "reactant_matrix",
      "Reactant encoding per reaction: flattened [sp_idx, coeff, sp_idx, coeff, ...]");

  params.addRequiredParam<std::vector<std::string>>(
      "coefficient_names",
      "Rate coefficient names in evaluation order");
  params.addRequiredParam<std::vector<std::string>>(
      "coefficient_expressions",
      "Rate coefficient expressions in evaluation order (converted)");

  params.addRequiredParam<std::vector<Real>>(
      "j_cl_values", "MCM photolysis CL parameter per J<N>, sorted by J number.");
  params.addRequiredParam<std::vector<Real>>(
      "j_cmm_values", "MCM photolysis CMM parameter per J<N>, sorted by J number.");
  params.addRequiredParam<std::vector<Real>>(
      "j_cnn_values", "MCM photolysis CNN parameter per J<N>, sorted by J number.");

  params.addParam<Real>("latitude", 51.51, "Latitude in degrees (North positive)");
  params.addParam<Real>("longitude", 0.13, "Longitude in degrees (East positive)");
  params.addParam<unsigned int>("day", 21, "Day of month for solar zenith angle calculation");
  params.addParam<unsigned int>("month", 6, "Month for solar zenith angle calculation");
  params.addParam<unsigned int>("year", 2010, "Year for solar zenith angle calculation");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor for photolysis rates");
  params.addParam<bool>("roof_open", true, "Roof (chamber cover) open. false = CLOSED (all J=0)");

  params.addClassDescription("Material that evaluates MCM reaction rates via sequential fparser evaluation");
  return params;
}

MCMRatesMaterial::MCMRatesMaterial(const InputParameters & params)
  : Material(params),
    FunctionParserUtils<false>(params),
    _TEMP(getParam<Real>("temperature")),
    _M(getParam<Real>("air_density")),
    _H2O_val(getParam<Real>("water_vapor")),
    _coeff_names(getParam<std::vector<std::string>>("coefficient_names")),
    _n_coefficients(_coeff_names.size()),
    _n_reactions(getParam<std::vector<std::string>>("reaction_rate_expressions").size()),
    _reactant_matrix(getParam<std::vector<std::vector<Real>>>("reactant_matrix")),
    _reaction_rates(declareProperty<std::vector<Real>>("reaction_rates")),
    _j_CL(getParam<std::vector<Real>>("j_cl_values")),
    _j_CMM(getParam<std::vector<Real>>("j_cmm_values")),
    _j_CNN(getParam<std::vector<Real>>("j_cnn_values")),
    _latitude(getParam<Real>("latitude")),
    _longitude(getParam<Real>("longitude")),
    _day(getParam<unsigned int>("day")),
    _month(getParam<unsigned int>("month")),
    _year(getParam<unsigned int>("year")),
    _jfac(getParam<Real>("jfac")),
    _roof_open(getParam<bool>("roof_open"))
{
  // Couple to species variables
  _n_species = coupledComponents("species_variables");
  _species_vals.resize(_n_species);
  for (unsigned int i = 0; i < _n_species; ++i)
    _species_vals[i] = &coupledValue("species_variables", i);

  // Read all expression data first (for J<N> detection and conversion)
  auto coeff_exprs = getParam<std::vector<std::string>>("coefficient_expressions");
  auto rxn_exprs = getParam<std::vector<std::string>>("reaction_rate_expressions");

  // Convert J<N> to PHOTOJ_N_ in all expressions (fparser doesn't allow < and >)
  auto replace_j = [](std::string & s) {
    std::string result;
    std::string::size_type i = 0;
    while (i < s.size())
    {
      if (i + 1 < s.size() && s[i] == 'J' && s[i+1] == '<')
      {
        std::string num;
        i += 2; // skip J<
        while (i < s.size() && s[i] >= '0' && s[i] <= '9')
        { num += s[i]; i++; }
        result += "PHOTOJ" + num;
        if (i < s.size() && s[i] == '>') i++; // skip >
      }
      else
      { result += s[i]; i += 1; }
    }
    s = result;
  };
  for (auto & e : coeff_exprs) replace_j(e);
  for (auto & e : rxn_exprs) replace_j(e);

  // Initialize species count (needed before J index calculation)
  _n_species_material = coupledComponents("species_variables");

  // Detect PHOTOJ numbers from expressions
  _j_index_start = 5 + _n_coefficients + _n_species_material;
  pcrecpp::RE re_j("PHOTOJ([0-9]+)");
  std::set<int> j_numbers;
  for (auto & e : coeff_exprs)
  {
    int jnum; pcrecpp::StringPiece sp(e);
    while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum);
  }
  for (auto & e : rxn_exprs)
  {
    int jnum; pcrecpp::StringPiece sp(e);
    while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum);
  }
  _n_j_variables = j_numbers.size();
  _j_names.clear();

  // Build master variable list with all variables
  std::string vars = "TEMP,M,O2,N2,H2O";
  _name_to_index["TEMP"] = 0; _name_to_index["M"] = 1;
  _name_to_index["O2"] = 2; _name_to_index["N2"] = 3;
  _name_to_index["H2O"] = 4;

  for (unsigned int i = 0; i < _n_coefficients; ++i)
  {
    vars += "," + _coeff_names[i];
    _name_to_index[_coeff_names[i]] = 5 + i;
  }

  // Add species names (for RO2 = CH3O2 etc.)
  _species_vals_material.resize(_n_species_material);
  auto species_list = getParam<std::vector<std::string>>("species_list");
  for (unsigned int i = 0; i < _n_species_material; ++i)
  {
    vars += "," + species_list[i];
    _species_name_to_index[species_list[i]] = 5 + _n_coefficients + i;
    _species_vals_material[i] = &coupledValue("species_variables", i);
  }

  // Add PHOTOJ variables
  for (auto & n : j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(n);
    vars += "," + jname;
    _j_names.push_back(jname);
  }

  // Resize parameter buffer
  _func_params.resize(5 + _n_coefficients + _n_species_material + _n_j_variables, 0.0);

  // The photolysis-rates file contains ALL MCM J values (~35); the mechanism
  // may only reference a subset.  Require at least as many entries as the
  // mechanism uses — extra entries in the file are harmless.
  if (_j_CL.size() < _n_j_variables || _j_CMM.size() < _n_j_variables || _j_CNN.size() < _n_j_variables)
    mooseError("MCMRatesMaterial: MCM photolysis parameter vectors (j_cl/cmm/cnn_values) have size ",
               _j_CL.size(), "/", _j_CMM.size(), "/", _j_CNN.size(),
               " but the mechanism references ", _n_j_variables,
               " J<N> variables. The photolysis-rates file must contain at least that many entries.");

  // Note: J values are not initialized here; they will be calculated from
  // solar zenith angle in computeQpProperties().

  // Parse coefficient expressions
  _coeff_parsers.resize(_n_coefficients);
  for (unsigned int i = 0; i < _n_coefficients; ++i)
  {
    _coeff_parsers[i] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_coeff_parsers[i]);
    if (_coeff_parsers[i]->Parse(coeff_exprs[i], vars) >= 0)
      mooseError("MCMRatesMaterial: Invalid coefficient expression for '",
                 _coeff_names[i], "':\n", coeff_exprs[i], "\n",
                 _coeff_parsers[i]->ErrorMsg());
    if (!_disable_fpoptimizer)
      _coeff_parsers[i]->Optimize();
  }

  // Parse reaction rate expressions
  _reaction_parsers.resize(_n_reactions);
  for (unsigned int i = 0; i < _n_reactions; ++i)
  {
    _reaction_parsers[i] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_reaction_parsers[i]);
    const std::string & rxn_expr = rxn_exprs[i];
    if (_reaction_parsers[i]->Parse(rxn_expr, vars) >= 0)
      mooseError("MCMRatesMaterial: Invalid reaction expression for reaction ", i, ":\n",
                 rxn_exprs[i], "\n",
                 _reaction_parsers[i]->ErrorMsg());
    if (!_disable_fpoptimizer)
      _reaction_parsers[i]->Optimize();
  }
}

void
MCMRatesMaterial::computeQpProperties()
{
  // Step 1: Set base variable values
  _func_params[0] = _TEMP;
  _func_params[1] = _M;
  _func_params[2] = 0.21 * _M; // O2 volume fraction
  _func_params[3] = 0.78 * _M; // N2 volume fraction
  _func_params[4] = _H2O_val;

  // Step 2: Calculate photolysis rates from solar zenith angle (MCM formula)
  // J = CL * cos(theta)^CMM * exp(-CNN / cos(theta)) * JFAC
  {
    Real cosx = calculateCosSZA(_t);
    Real secx = (cosx > 1.0e-10) ? (1.0 / cosx) : 1.0e2;
    for (unsigned int i = 0; i < _n_j_variables; ++i)
    {
      if (cosx > 1.0e-10)
        _func_params[_j_index_start + i] =
            _j_CL[i] * std::pow(cosx, _j_CMM[i]) * std::exp(-_j_CNN[i] * secx) * _jfac;
      else
        _func_params[_j_index_start + i] = 0.0;
    }
    // Roof (chamber cover) CLOSED: all J values forced to zero
    if (!_roof_open)
      for (unsigned int i = 0; i < _n_j_variables; ++i)
        _func_params[_j_index_start + i] = 0.0;
  }

  // Step 3: Update species concentration variables (for RO2 = CH3O2 etc.)
  for (unsigned int i = 0; i < _n_species_material; ++i)
    _func_params[5 + _n_coefficients + i] = (*_species_vals_material[i])[_qp];

  // Step 4: Evaluate rate coefficients in topological order.
  for (unsigned int i = 0; i < _n_coefficients; ++i)
  {
    Real val = evaluate(_coeff_parsers[i]);
    if (std::isnan(val))
      val = 0.0;
    _func_params[5 + i] = val;
  }

  // Step 5: Compute reaction rate k_i for each reaction.
  // _k_values is pre-allocated member buffer (Per.14 — no per-QP allocation).
  _k_values.assign(_n_reactions, 0.0);
  for (unsigned int i = 0; i < _n_reactions; ++i)
    _k_values[i] = evaluate(_reaction_parsers[i]);

  // Step 6: Compute R_i = k_i * Π [C_reactant]^ν
  _reaction_rates[_qp].resize(_n_reactions);
  for (unsigned int i = 0; i < _n_reactions; ++i)
  {
    Real rate = _k_values[i];
    if (std::isnan(rate))
      rate = 0.0;
    const auto & row = _reactant_matrix[i];
    for (size_t j = 0; j + 1 < row.size(); j += 2)
    {
      // Defense against negative/corrupt indices in reactant_matrix
      if (row[j] < 0.0 || row[j] >= (Real)_n_species)
        mooseError("MCMRatesMaterial: invalid species index ", row[j],
                   " in reaction ", i, " (nSpecies=", _n_species, ")");
      unsigned int sp_idx = static_cast<unsigned int>(row[j]);
      Real coeff = row[j + 1];
      rate *= std::pow(std::max((*_species_vals[sp_idx])[_qp], 0.0), coeff);
    }
    _reaction_rates[_qp][i] = rate;
  }
}

unsigned int
MCMRatesMaterial::computeDayOfYear() const
{
  unsigned int days_in_months[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0)
    days_in_months[1] = 29;

  unsigned int doy = 0;
  for (unsigned int m = 0; m < _month - 1; ++m)
    doy += days_in_months[m];
  doy += _day;
  return doy;
}

Real
MCMRatesMaterial::calculateCosSZA(Real t) const
{
  constexpr Real pi = 3.14159265358979323846;

  // Day angle (radians)
  unsigned int doy = computeDayOfYear();
  unsigned int days_in_year =
      ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0) ? 366 : 365;
  Real theta = 2.0 * pi * static_cast<Real>(doy) / static_cast<Real>(days_in_year);

  // Solar declination (radians) — Madronich (1993)
  Real dec = 0.006918 - 0.399912 * cos(theta) + 0.070257 * sin(theta) -
             0.006758 * cos(2.0 * theta) + 0.000907 * sin(2.0 * theta) -
             0.002697 * cos(3.0 * theta) + 0.001480 * sin(3.0 * theta);

  // Equation of time
  Real eqt = 0.000075 + 0.001868 * cos(theta) - 0.032077 * sin(theta) -
             0.014615 * cos(2.0 * theta) - 0.040849 * sin(2.0 * theta);

  // Local hour angle (radians)
  Real current_frac_hour = std::fmod(t / 3600.0, 24.0);
  Real lha = pi * ((current_frac_hour / 12.0) - (1.0 + _longitude / 180.0)) + eqt;

  // Cosine of solar zenith angle
  Real lat_rad = _latitude * pi / 180.0;
  Real cosx = cos(lha) * cos(lat_rad) * cos(dec) + sin(lat_rad) * sin(dec);

  if (cosx <= 0.0)
    cosx = 0.0;

  return cosx;
}
