//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMBoxModel.h"
#include <iomanip>

registerMooseObject("AtmosphericChemistryApp", MCMBoxModel);

InputParameters
MCMBoxModel::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addParam<std::string>(
      "mechanism_file", "", "Path to MCM Facsimile (.fac) mechanism file for auto-parsing");
  params.addClassDescription(
      "Centralized box model UserObject for atmospheric chemistry ODE systems. "
      "Manages stoichiometric coefficients, reactant indices, and rate constants, "
      "and provides dC/dt computation following F0AM's dydt_eval algorithm.");
  return params;
}

MCMBoxModel::MCMBoxModel(const InputParameters & params)
  : GeneralUserObject(params), _n_species(0), _n_reactions(0),
    _photolysis_method(MCM_SZA), _kdil(0.0), _dirty(true)
{
}

void
MCMBoxModel::initialize()
{
  // Parse .fac file if provided (auto-populate from mechanism_file)
  std::string mech_file = getParam<std::string>("mechanism_file");
  if (!mech_file.empty())
  {
    MCMFacsimileParser parser;
    ParsedMechanism mech = parser.parse(mech_file);
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

  // --- Evaluate rate coefficients ---
  // Build map from coefficient name → numeric value
  std::map<std::string, Real> coeff_map;
  for (unsigned int i = 0; i < mech.coefficient_names.size(); ++i)
  {
    const std::string & expr = mech.coefficient_expressions[i];
    Real val = 1.0;
    try
    {
      val = std::stod(expr);
    }
    catch (...)
    {
      mooseWarning("MCMBoxModel: coefficient '", mech.coefficient_names[i],
                   "' has non-numeric expression '", expr, "'; using 1.0");
    }
    coeff_map[mech.coefficient_names[i]] = val;
  }

  // Resolve reaction rate expressions
  _k.assign(_n_reactions, 1.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const std::string & expr = mech.reactions[r].rate_expression;
    auto it = coeff_map.find(expr);
    if (it != coeff_map.end())
      _k[r] = it->second;
    else
    {
      // Try direct numeric parse
      try { _k[r] = std::stod(expr); }
      catch (...)
      {
        mooseWarning("MCMBoxModel: cannot evaluate rate '", expr, "' for reaction ", r,
                     "; using 1.0");
      }
    }
  }

  _dirty = true;

  // --- Debug: print matrices for comparison with F0AM ---
  _console << "\n===== MCMBoxModel matrices =====\n";
  _console << "Species (" << _n_species << "):";
  for (unsigned int s = 0; s < _n_species; ++s)
    _console << " [" << s << "]=" << _species_names[s];
  _console << "\n\nStoichiometric matrix f (nRx×nSp, F0AM convention):\n";
  _console << "         ";
  for (unsigned int s = 0; s < _n_species; ++s)
    _console << std::setw(6) << _species_names[s];
  _console << "\n";
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    _console << "Rx" << std::setw(2) << r << ":  ";
    for (unsigned int s = 0; s < _n_species; ++s)
      _console << std::setw(6) << _f[r][s];
    _console << "\n";
  }
  _console << "\niG (reactant indices per reaction):\n";
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _console << "  Rx" << r << ": [" << _iG[r][0] << ", " << _iG[r][1] << ", " << _iG[r][2] << "]\n";
  _console << "\nk (rate constants):\n";
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _console << "  Rx" << r << ": " << std::scientific << _k[r] << "\n";
  _console << "===== end matrices =====\n\n";
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
