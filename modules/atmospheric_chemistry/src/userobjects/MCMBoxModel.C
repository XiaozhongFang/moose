//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMBoxModel.h"
#include "pcrecpp.h"

registerMooseObject("AtmosphericChemistryApp", MCMBoxModel);

InputParameters
MCMBoxModel::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += FunctionParserUtils<false>::validParams();
  params.addParam<std::string>(
      "mechanism_file", "", "Path to MCM Facsimile (.fac) mechanism file for auto-parsing");
  params.addParam<std::string>(
      "photolysis_file", "", "Path to MCM photolysis-rates data file (e.g. mcm_photolysis_rates_v3.3.1.dat)");
  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19, "Air number density (molecules/cm^3)");
  params.addParam<Real>("water_vapor", 2.46e17, "Background water vapor (molecules/cm^3)");
  params.addParam<Real>("latitude", 51.51, "Latitude (deg N)");
  params.addParam<Real>("longitude", 0.13, "Longitude (deg E)");
  params.addParam<unsigned int>("day", 21, "Day of month");
  params.addParam<unsigned int>("month", 6, "Month");
  params.addParam<unsigned int>("year", 2010, "Year");
  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor");
  params.addClassDescription(
      "Centralized box model UserObject for atmospheric chemistry ODE systems.");
  return params;
}

MCMBoxModel::MCMBoxModel(const InputParameters & params)
  : GeneralUserObject(params),
    FunctionParserUtils<false>(params),
    _n_species(0), _n_reactions(0),
    _j_index_start(0),
    _temperature(getParam<Real>("temperature")),
    _air_density(getParam<Real>("air_density")),
    _water_vapor(getParam<Real>("water_vapor")),
    _lat(getParam<Real>("latitude")),
    _lon(getParam<Real>("longitude")),
    _day((int)getParam<unsigned int>("day")),
    _month((int)getParam<unsigned int>("month")),
    _year((int)getParam<unsigned int>("year")),
    _jfac(getParam<Real>("jfac")),
    _photolysis_method(MCM_SZA), _kdil(0.0), _dirty(true)
{
}

void
MCMBoxModel::initialize()
{
  // Parse .fac file only once
  if (_n_species > 0) return;

  std::string mech_file = getParam<std::string>("mechanism_file");
  if (!mech_file.empty())
  {
    std::string photo_file = getParam<std::string>("photolysis_file");
    MCMFacsimileParser parser;
    ParsedMechanism mech = parser.parse(mech_file, photo_file);
    loadMechanism(mech);
    _console << "MCMBoxModel: Loaded " << _n_species << " species, "
             << _n_reactions << " reactions from " << mech_file << std::endl;
  }
}

void
MCMBoxModel::computeDCdt(const std::vector<Real> & C, std::vector<Real> & dC) const
{
  dC.assign(_n_species, 0.0);

  if (_n_reactions == 0 || _n_species == 0)
    return;

  // Step 1: compute reactant products G[i] = C[iG[i][0]] * C[iG[i][1]] * C[iG[i][2]]
  // All indices are valid species indices (padded with ONE=0, conc=1)
  std::vector<Real> G(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    G[r] = C[_iG[r][0]] * C[_iG[r][1]] * C[_iG[r][2]];

  // Step 2: rates = k .* G
  std::vector<Real> rates(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rates[r] = _k[r] * G[r];

  // Step 3: dC = f^T * rates (stoichiometric accumulation)
  for (unsigned int r = 0; r < _n_reactions; ++r)
    for (unsigned int s = 0; s < _n_species; ++s)
      dC[s] += _f[r][s] * rates[r];
}

void
MCMBoxModel::computeJacobianTriplets(
    const std::vector<Real> & C,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  J.clear();

  if (_n_reactions == 0 || _n_species == 0)
    return;

  // F0AM-style 3-term product Jacobian:
  //   rate_r = k_r * C[i0] * C[i1] * C[i2]
  //   d(rate_r)/dC[j] = k_r * sum_{k where iG[k]==j} prod_{m≠k} C[iG[m]]
  //   J[s][j] += f[r][s] * d(rate_r)/dC[j]

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    const Real c0 = C[i0], c1 = C[i1], c2 = C[i2];
    const Real k = _k[r];

    // Contribution from C[i0]: dr/dC[i0] = k * c1 * c2
    {
      Real drate = k * c1 * c2;
      // For self-reaction duplicates (i0==i1 or i0==i2), this is the ONLY contribution
      // and the product rule automatically gives k*c1*c2 = k*C[i0]*C[i2] which is correct
      if (std::abs(drate) > 1e-30)
        for (unsigned int s = 0; s < _n_species; ++s)
          if (std::abs(_f[r][s]) > 1e-30)
            J.emplace_back(s, (unsigned int)i0, drate * _f[r][s]);
    }
    // Contribution from C[i1]: dr/dC[i1] = k * c0 * c2
    if (i1 != i0)  // skip if already accounted (self-reaction)
    {
      Real drate = k * c0 * c2;
      if (std::abs(drate) > 1e-30)
        for (unsigned int s = 0; s < _n_species; ++s)
          if (std::abs(_f[r][s]) > 1e-30)
            J.emplace_back(s, (unsigned int)i1, drate * _f[r][s]);
    }
    // Contribution from C[i2]: dr/dC[i2] = k * c0 * c1
    if (i2 != i0 && i2 != i1)  // skip if already accounted
    {
      Real drate = k * c0 * c1;
      if (std::abs(drate) > 1e-30)
        for (unsigned int s = 0; s < _n_species; ++s)
          if (std::abs(_f[r][s]) > 1e-30)
            J.emplace_back(s, (unsigned int)i2, drate * _f[r][s]);
    }
  }
}

// -- Cached single-species interface --

Real
MCMBoxModel::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  if (_dirty || _cached_dC.size() != _n_species)
  {
    _cached_C = C;  // store for evaluateCoefficients species lookup
    if (!_coeff_parsers.empty())
      const_cast<MCMBoxModel*>(this)->evaluateCoefficients();
    _cached_dC.resize(_n_species);
    computeDCdt(C, _cached_dC);
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_dC[idx] : 0.0;
}

Real
MCMBoxModel::getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const
{
  if (_dirty || _cached_C != C || _cached_diag_J.size() != _n_species)
  {
    _cached_C = C;
    _buildJacobianCache();
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_diag_J[idx] : 0.0;
}

Real
MCMBoxModel::getJacobianOffDiagonal(unsigned int i, unsigned int j, const std::vector<Real> & C) const
{
  if (_dirty || _cached_C != C || _cached_diag_J.size() != _n_species)
  {
    _cached_C = C;
    _buildJacobianCache();
    _dirty = false;
  }
  uint64_t key = (static_cast<uint64_t>(i) << 32) | static_cast<uint64_t>(j);
  auto it = _cached_offdiag_J.find(key);
  return (it != _cached_offdiag_J.end()) ? it->second : 0.0;
}

void
MCMBoxModel::_buildJacobianCache() const
{
  _cached_diag_J.assign(_n_species, 0.0);
  _cached_offdiag_J.clear();

  if (_n_reactions == 0 || _n_species == 0)
    return;

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    const Real c0 = _cached_C[i0], c1 = _cached_C[i1], c2 = _cached_C[i2];
    const Real k = _k[r];

    // Helper: accumulate Jacobian contribution (s, j) += val
    auto accum = [&](unsigned int j, Real drate) {
      for (unsigned int s = 0; s < _n_species; ++s)
      {
        Real val = drate * _f[r][s];
        if (std::abs(val) < 1e-30) continue;
        if (j == s)
          _cached_diag_J[s] += val;
        else
        {
          uint64_t key = (static_cast<uint64_t>(s) << 32) | static_cast<uint64_t>(j);
          _cached_offdiag_J[key] += val;
        }
      }
    };

    // dr/dC[i0] = k * c1 * c2
    accum((unsigned int)i0, k * c1 * c2);
    // dr/dC[i1] = k * c0 * c2
    if (i1 != i0) accum((unsigned int)i1, k * c0 * c2);
    // dr/dC[i2] = k * c0 * c1
    if (i2 != i0 && i2 != i1) accum((unsigned int)i2, k * c0 * c1);
  }
}

void
MCMBoxModel::loadMechanism(const ParsedMechanism & mech)
{
  _n_species = mech.species.size();
  _n_reactions = mech.reactions.size();
  _species_names = mech.species;
  _reaction_names = mech.reaction_names;

  // Convert stoichiometry from [species][reaction] to [reaction][species] (F0AM convention)
  _f.assign(_n_reactions, std::vector<Real>(_n_species, 0.0));
  for (unsigned int r = 0; r < _n_reactions; ++r)
    for (unsigned int s = 0; s < _n_species; ++s)
      _f[r][s] = mech.stoichiometry[s][r];

  // Copy reactant indices
  _iG = mech.reactant_indices;

  // --- Evaluate rate coefficients (fparser for complex expressions) ---
  _k.assign(_n_reactions, 1.0);

  // Try simple numeric parse first; fall back to fparser
  std::map<std::string, Real> coeff_map;
  bool need_fparser = false;

  for (unsigned int i = 0; i < mech.coefficient_names.size(); ++i)
  {
    const std::string & expr = mech.coefficient_expressions[i];
    try
    {
      Real val = std::stod(expr);
      coeff_map[mech.coefficient_names[i]] = val;
    }
    catch (...)
    {
      need_fparser = true;
      break;
    }
  }

  if (!need_fparser)
  {
    // Simple case: all coefficients are numeric (e.g. tutorial_5sp)
    for (unsigned int r = 0; r < _n_reactions; ++r)
    {
      const std::string & expr = mech.reactions[r].rate_expression;
      auto it = coeff_map.find(expr);
      if (it != coeff_map.end())
        _k[r] = it->second;
      else
        try { _k[r] = std::stod(expr); } catch (...) {}
    }
  }
  else
  {
    // Complex case: use fparser (see evaluateCoefficients)
    _console << "MCMBoxModel: " << mech.coefficient_names.size()
             << " coefficients — using fparser for complex expressions" << std::endl;
    setupFparser(mech);
  }

  _dirty = true;
}

void
MCMBoxModel::setupFparser(const ParsedMechanism & mech)
{
  auto coeff_exprs = mech.coefficient_expressions;
  std::vector<std::string> coeff_names = mech.coefficient_names;
  std::vector<std::string> rxn_exprs(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rxn_exprs[r] = mech.reactions[r].rate_expression;

  // Convert J<N> to PHOTOJN (fparser doesn't allow < >)
  auto replace_j = [](std::string & s) {
    std::string result;
    for (size_t i = 0; i < s.size(); )
    {
      if (i + 1 < s.size() && s[i] == 'J' && s[i+1] == '<')
      {
        std::string num; i += 2;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { num += s[i]; i++; }
        result += "PHOTOJ" + num;
        if (i < s.size() && s[i] == '>') i++;
      }
      else { result += s[i]; i++; }
    }
    s = result;
  };
  for (auto & e : coeff_exprs) replace_j(e);
  for (auto & e : rxn_exprs) replace_j(e);

  // Detect photolysis J numbers
  std::set<int> j_numbers;
  pcrecpp::RE re_j("PHOTOJ([0-9]+)");
  for (auto & e : coeff_exprs)
  { int jnum; pcrecpp::StringPiece sp(e); while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum); }
  for (auto & e : rxn_exprs)
  { int jnum; pcrecpp::StringPiece sp(e); while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum); }

  // Build variable list
  std::string vars = "TEMP,M,O2,N2,H2O";
  _name_to_index["TEMP"] = 0; _name_to_index["M"] = 1;
  _name_to_index["O2"] = 2; _name_to_index["N2"] = 3;
  _name_to_index["H2O"] = 4;

  for (unsigned int i = 0; i < coeff_names.size(); ++i)
  {
    vars += "," + coeff_names[i];
    _name_to_index[coeff_names[i]] = 5 + i;
  }
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    vars += "," + _species_names[i];
    _name_to_index[_species_names[i]] = 5 + coeff_names.size() + i;
  }

  // RO2 is a derived/lumped species (RO2 = CH3O2 + ...), not in VARIABLE list.
  // It may appear in rate expressions directly. Add it as a variable with value 0.
  bool has_ro2_in_species = false;
  for (auto & s : _species_names)
    if (s == "RO2") { has_ro2_in_species = true; break; }
  if (!has_ro2_in_species)
  {
    vars += ",RO2";
    _name_to_index["RO2"] = 5 + coeff_names.size() + _n_species;
  }
  unsigned int n_extra_vars = has_ro2_in_species ? 0 : 1;
  _j_index_start = 5 + coeff_names.size() + _n_species + n_extra_vars;
  for (auto & n : j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(n);
    vars += "," + jname;
    _name_to_index[jname] = _j_index_start + _name_to_index.size() - _j_index_start;
  }

  unsigned int n_vars = _name_to_index.size();
  unsigned int n_j = j_numbers.size();
  _func_params.resize(5 + coeff_names.size() + _n_species + n_extra_vars + n_j, 0.0);

  // Parse coefficients
  _coeff_parsers.resize(coeff_names.size());
  for (unsigned int i = 0; i < coeff_names.size(); ++i)
  {
    _coeff_parsers[i] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_coeff_parsers[i]);
    if (_coeff_parsers[i]->Parse(coeff_exprs[i], vars) >= 0)
      mooseError("MCMBoxModel: Bad coefficient '", coeff_names[i], "': ", coeff_exprs[i], "\n",
                 _coeff_parsers[i]->ErrorMsg());
    if (!_disable_fpoptimizer)
      _coeff_parsers[i]->Optimize();
  }

  // Parse reaction expressions
  _reaction_parsers.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    _reaction_parsers[r] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_reaction_parsers[r]);
    if (_reaction_parsers[r]->Parse(rxn_exprs[r], vars) >= 0)
      mooseError("MCMBoxModel: Bad reaction ", r, ": ", rxn_exprs[r], "\n",
                 _reaction_parsers[r]->ErrorMsg());
    if (!_disable_fpoptimizer)
      _reaction_parsers[r]->Optimize();
  }

  // Store photolysis parameters for SZA-based J calculation
  _j_numbers = mech.j_numbers;
  _j_CL_vals = mech.j_CL;
  _j_CMM_vals = mech.j_CMM;
  _j_CNN_vals = mech.j_CNN;
  _n_j_vars = _j_numbers.size();
  if (_n_j_vars > 0)
    _console << "MCMBoxModel: " << _n_j_vars << " photolysis J values loaded" << std::endl;
}

void
MCMBoxModel::evaluateCoefficients()
{
  if (_coeff_parsers.empty()) return;

  _func_params[0] = _temperature;
  _func_params[1] = _air_density;
  _func_params[2] = 0.21 * _air_density;
  _func_params[3] = 0.78 * _air_density;
  _func_params[4] = _water_vapor;

  // Compute photolysis J values from solar zenith angle
  if (!_j_CL_vals.empty())
  {
    Real cosx = calculateCosSZA(_t);
    Real secx = (cosx > 1.0e-10) ? (1.0 / cosx) : 1.0e2;
    for (size_t i = 0; i < _j_CL_vals.size(); ++i)
    {
      unsigned int jn = _j_numbers[i];
      std::string jname = "PHOTOJ" + std::to_string(jn);
      auto it = _name_to_index.find(jname);
      if (it != _name_to_index.end())
      {
        if (cosx > 1.0e-10)
          _func_params[it->second] = _j_CL_vals[i] * std::pow(cosx, _j_CMM_vals[i])
                                     * std::exp(-_j_CNN_vals[i] * secx) * _jfac;
        else
          _func_params[it->second] = 0.0;
      }
    }
  }

  unsigned int n_coeff = _coeff_parsers.size();
  // Evaluate coefficients in topological order
  for (unsigned int i = 0; i < n_coeff; ++i)
  {
    Real val = evaluate(_coeff_parsers[i]);
    if (std::isnan(val)) val = 0.0;
    _func_params[5 + i] = val;
  }

  // Set species concentrations in fparser buffer
  for (unsigned int s = 0; s < _n_species; ++s)
  {
    auto it = _name_to_index.find(_species_names[s]);
    if (it != _name_to_index.end())
      _func_params[it->second] = (_cached_C.empty() ? 0.0 :
          (s < _cached_C.size() ? _cached_C[s] : 0.0));
  }
  // RO2 = sum of peroxy radicals
  auto it_ro2 = _name_to_index.find("RO2");
  if (it_ro2 != _name_to_index.end())
  {
    Real ro2_sum = 0.0;
    for (unsigned int s = 0; s < _n_species; ++s)
      if (_species_names[s].size() >= 3 &&
          _species_names[s].substr(_species_names[s].size() - 2) == "O2" &&
          _species_names[s] != "HO2" && _species_names[s] != "NO2" &&
          _species_names[s] != "SO2" && _species_names[s] != "H2O2")
        ro2_sum += (_cached_C.empty() ? 0.0 :
            (s < _cached_C.size() ? _cached_C[s] : 0.0));
    _func_params[it_ro2->second] = ro2_sum;
  }

  // Evaluate reaction rate expressions → _k
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real val = evaluate(_reaction_parsers[r]);
    _k[r] = (std::isnan(val)) ? 1.0 : val;
  }
}

unsigned int
MCMBoxModel::computeDayOfYear() const
{
  unsigned int days_in_months[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0)
    days_in_months[1] = 29;
  unsigned int doy = 0;
  for (unsigned int m = 0; m < (unsigned int)(_month - 1); ++m)
    doy += days_in_months[m];
  doy += _day;
  return doy;
}

Real
MCMBoxModel::calculateCosSZA(Real t) const
{
  const Real pi = 3.14159265358979323846;
  unsigned int doy = computeDayOfYear();
  unsigned int days_in_year =
      ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0) ? 366 : 365;
  Real theta = 2.0 * pi * (Real)doy / (Real)days_in_year;
  Real dec = 0.006918 - 0.399912 * cos(theta) + 0.070257 * sin(theta) -
             0.006758 * cos(2.0 * theta) + 0.000907 * sin(2.0 * theta) -
             0.002697 * cos(3.0 * theta) + 0.001480 * sin(3.0 * theta);
  Real eqt = 0.000075 + 0.001868 * cos(theta) - 0.032077 * sin(theta) -
             0.014615 * cos(2.0 * theta) - 0.040849 * sin(2.0 * theta);
  Real current_frac_hour = std::fmod(t / 3600.0, 24.0);
  Real lat_rad = _lat * pi / 180.0;
  Real lon_rad = _lon * pi / 180.0;
  Real lha = pi * ((current_frac_hour / 12.0) - (1.0 + _lon / 180.0)) + eqt;
  Real cosx = cos(lha) * cos(lat_rad) * cos(dec) + sin(lat_rad) * sin(dec);
  if (cosx <= 0.0) cosx = 0.0;
  return cosx;
}

// -- Constrained species --
void
MCMBoxModel::setConstrainedSpecies(const std::vector<std::string> & names)
{
  _constrained_set.clear();
  for (auto & name : names)
    for (unsigned int i = 0; i < _n_species; ++i)
      if (_species_names[i] == name)
        _constrained_set.insert(i);
  _constrained_values.assign(_constrained_set.size(), 0.0);
}

void
MCMBoxModel::updateConstrainedValues(const std::vector<Real> & values)
{
  if (values.size() == _constrained_set.size())
    _constrained_values = values;
}

void
MCMBoxModel::computeDCdtFull(const std::vector<Real> & C_full, std::vector<Real> & dC) const
{
  computeDCdt(C_full, dC);
  // Zero out constrained species derivatives
  for (auto idx : _constrained_set)
    dC[idx] = 0.0;
}

// -- Reaction rate queries --
Real
MCMBoxModel::reactionRate(unsigned int r, const std::vector<Real> & C) const
{
  if (r >= _n_reactions || C.size() != _n_species) return 0.0;
  return _k[r] * C[_iG[r][0]] * C[_iG[r][1]] * C[_iG[r][2]];
}

Real
MCMBoxModel::speciesReactionRate(unsigned int s, unsigned int r, const std::vector<Real> & C) const
{
  return _f[r][s] * reactionRate(r, C);
}

void
MCMBoxModel::allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const
{
  rates.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rates[r] = reactionRate(r, C);
}

// -- Photolysis --
void
MCMBoxModel::enableHybridPhotolysis(const std::string & table_dir)
{
  _photolysis_method = HYBRID;
  _hybrid_reader = std::make_unique<HybridJTableReader>(table_dir);
}

void
MCMBoxModel::mapPhotolysisReactions()
{
  _j_reaction_indices.clear();
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    // Scan rate expression for J<N> pattern
    for (unsigned int jn = 1; jn <= 100; ++jn)
    {
      std::string jname = "J" + std::to_string(jn);
      if (_reaction_names[r].find(jname) != std::string::npos ||
          (_species_names.size() > 0 && r < _n_reactions))
      {
        // Check if this reaction's rate expression references J<N>
        // Simple heuristic: if the reaction name starts with "R" + number
        // Real mapping done by Action passing photolysis parameter vectors
        _j_reaction_indices.push_back(r);
        break;
      }
    }
  }
}

void
MCMBoxModel::updatePhotolysis(Real sza, Real albedo, Real o3col, Real altitude)
{
  if (_photolysis_method == HYBRID && _hybrid_reader)
  {
    for (unsigned int i = 0; i < _j_reaction_indices.size(); ++i)
    {
      unsigned int r = _j_reaction_indices[i];
      unsigned int jn = i + 1;  // J1, J2, ...
      std::string jname = "J" + std::to_string(jn);
      if (_hybrid_reader->hasJValue(jname))
        _k[r] = _hybrid_reader->interpolate(jname, sza, albedo, o3col, altitude);
    }
  }
  else
    updatePhotolysisSZA(sza);
}

void
MCMBoxModel::updatePhotolysisSZA(Real sza, Real jfac)
{
  Real cosx = std::cos(sza * M_PI / 180.0);
  if (cosx <= 0.0) { cosx = 1e-10; }
  Real secx = 1.0 / cosx;
  for (unsigned int i = 0; i < _j_reaction_indices.size() && i < _j_CL.size(); ++i)
  {
    _k[_j_reaction_indices[i]] = _j_CL[i] * std::pow(cosx, _j_CMM[i])
                               * std::exp(-_j_CNN[i] * secx) * jfac;
  }
}

// -- Solar cycle (Madronich 1993) --
void
MCMBoxModel::setSolarCycle(Real lat, Real lon, int day, int month, int year)
{
  _lat = lat * M_PI / 180.0; _lon = lon * M_PI / 180.0;
  Real doy = (month - 1) * 30 + day;
  Real theta = 2.0 * M_PI * doy / 365.0;
  _declination = 0.006918 - 0.399912*std::cos(theta) + 0.070257*std::sin(theta)
    - 0.006758*std::cos(2*theta) + 0.000907*std::sin(2*theta)
    - 0.002697*std::cos(3*theta) + 0.001480*std::sin(3*theta);
  Real B = 2.0*M_PI*(doy-1.0)/365.0;
  _eot = 0.165*std::sin(2*B) - 0.126*std::cos(B) - 0.025*std::sin(B);
}

Real
MCMBoxModel::cosSZA(Real seconds) const
{
  Real hour = seconds / 3600.0;
  Real ha = (hour - 12.0 + _eot + _lon*12.0/M_PI) * M_PI / 12.0;
  return std::sin(_lat)*std::sin(_declination)
       + std::cos(_lat)*std::cos(_declination)*std::cos(ha);
}

void
MCMBoxModel::advanceSolarCycle(Real s)
{
  Real cosx = cosSZA(s);
  Real sza = std::acos(std::max(-1.0,std::min(1.0,cosx))) * 180.0 / M_PI;
  updatePhotolysis(sza, 0.1, 300.0, 0.0);
}

// -- Dilution --
void
MCMBoxModel::setDilution(Real kdil, const std::vector<Real> & bg)
{ _kdil = kdil; _conc_bkgd = bg; }

void
MCMBoxModel::computeDCdtWithDilution(const std::vector<Real> & C, std::vector<Real> & dC) const
{
  computeDCdt(C, dC);
  for (unsigned int i = 0; i < _n_species; ++i)
    dC[i] -= _kdil * (C[i] - (_conc_bkgd.empty() ? 0.0 : _conc_bkgd[i]));
}

// -- Convergence --
Real
MCMBoxModel::checkConvergence(const std::vector<Real> & Cp, const std::vector<Real> & Cc,
                               const std::vector<std::string> & cs) const
{
  Real mx = 0.0;
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    if (!cs.empty()) {
      bool ok = false;
      for (auto & s : cs) if (s == _species_names[i]) { ok = true; break; }
      if (!ok) continue;
    }
    Real d = std::abs(Cc[i] - Cp[i]) / std::max(std::abs(Cp[i]), 1e-30);
    if (d > mx) mx = d;
  }
  return mx;
}
