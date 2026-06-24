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

// ---- StoichMatrix::build ----------------------------------------------------

void
StoichMatrix::build(const ParsedMechanism & mech, Format fmt)
{
  format = fmt;
  nSpecies = mech.species.size();
  nReactions = mech.reactions.size();

  // Clear all format vectors unconditionally — prevents stale data if build()
  // is ever called more than once (e.g. mechanism reload with different format).
  csr_cols.clear(); csr_vals.clear(); csr_row_ptr.clear();
  lhs_species.clear(); lhs_coeff.clear(); lhs_row_ptr.clear();
  rhs_species.clear(); rhs_coeff.clear(); rhs_row_ptr.clear();
  dense.clear();
  csc_rows.clear(); csc_c_vals.clear(); csc_col_ptr.clear();

  // Build name -> index map (needed by COO format for reactant/product lookup)
  std::unordered_map<std::string, unsigned int> name_to_idx;
  for (unsigned int i = 0; i < nSpecies; ++i)
    name_to_idx[mech.species[i]] = i;

  switch (format)
  {
    case CSR:
    {
      csr_row_ptr.resize(nReactions + 1);
      csr_row_ptr[0] = 0;
      csr_cols.clear();
      csr_vals.clear();
      for (unsigned int r = 0; r < nReactions; ++r)
      {
        for (unsigned int s = 0; s < nSpecies; ++s)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csr_cols.push_back((int)s);
            csr_vals.push_back(val);
          }
        }
        csr_row_ptr[r + 1] = csr_cols.size();
      }
      break;
    }
    case COO:
    {
      lhs_row_ptr.resize(nReactions + 1);
      rhs_row_ptr.resize(nReactions + 1);
      lhs_row_ptr[0] = rhs_row_ptr[0] = 0;
      lhs_species.clear(); lhs_coeff.clear();
      rhs_species.clear(); rhs_coeff.clear();

      for (unsigned int r = 0; r < nReactions; ++r)
      {
        const auto & rx = mech.reactions[r];
        for (const auto & [coeff, name] : rx.reactants)
        {
          auto it = name_to_idx.find(name);
          if (it != name_to_idx.end())
          {
            lhs_species.push_back((int)it->second);
            lhs_coeff.push_back(coeff);
          }
        }
        lhs_row_ptr[r + 1] = lhs_species.size();

        for (const auto & [coeff, name] : rx.products)
        {
          auto it = name_to_idx.find(name);
          if (it != name_to_idx.end())
          {
            rhs_species.push_back((int)it->second);
            rhs_coeff.push_back(coeff);
          }
        }
        rhs_row_ptr[r + 1] = rhs_species.size();
      }
      break;
    }
    case DENSE:
    {
      dense.assign(nReactions, std::vector<Real>(nSpecies, 0.0));
      for (unsigned int s = 0; s < nSpecies; ++s)
        for (unsigned int r = 0; r < nReactions; ++r)
          dense[r][s] = mech.stoichiometry[s][r];
      break;
    }
    case CSC:
    {
      // Build species-major (column) storage from dense stoichiometry.
      // Simultaneously populate CSR fields as a row-iteration forward index.
      csc_col_ptr.resize(nSpecies + 1);
      csc_col_ptr[0] = 0;
      csc_rows.clear();
      csc_c_vals.clear();

      csr_row_ptr.resize(nReactions + 1);
      csr_row_ptr[0] = 0;
      csr_cols.clear();
      csr_vals.clear();

      // Two-pass: first count nonzeros per species, then fill
      for (unsigned int s = 0; s < nSpecies; ++s)
      {
        for (unsigned int r = 0; r < nReactions; ++r)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csc_rows.push_back((int)r);
            csc_c_vals.push_back(val);
          }
        }
        csc_col_ptr[s + 1] = csc_rows.size();
      }

      // Build CSR forward index for row iteration
      for (unsigned int r = 0; r < nReactions; ++r)
      {
        for (unsigned int s = 0; s < nSpecies; ++s)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csr_cols.push_back((int)s);
            csr_vals.push_back(val);
          }
        }
        csr_row_ptr[r + 1] = csr_cols.size();
      }
      break;
    }
  }
}

// ---- MCMBoxModel ------------------------------------------------------------

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
  MooseEnum stoich_fmt("CSR COO DENSE CSC", "CSR");
  params.addParam<MooseEnum>(
      "stoich_format", stoich_fmt,
      "Stoichiometric matrix storage format.  CSR = compressed sparse row (PETSc AIJ-compatible, "
      "HPC default).  COO = AtChem2-style split reactant/product.  DENSE = dense 2D array "
      "(best for < ~50 species).  CSC = compressed sparse column (species-major; enables "
      "column queries 'which reactions involve species X?').  Analogous to PETSc -mat_type.");
  params.addClassDescription(
      "Centralized box model UserObject for atmospheric chemistry ODE systems.");
  return params;
}

MCMBoxModel::MCMBoxModel(const InputParameters & params)
  : GeneralUserObject(params),
    FunctionParserUtils<false>(params),
    _n_species(0), _n_reactions(0),
    _photolysis_method(MCM_SZA),
    _lat(getParam<Real>("latitude")),
    _lon(getParam<Real>("longitude")),
    _day((int)getParam<unsigned int>("day")),
    _month((int)getParam<unsigned int>("month")),
    _year((int)getParam<unsigned int>("year")),
    _kdil(0.0),
    _roof_open(true),
    _j_index_start(0),
    _temperature(getParam<Real>("temperature")),
    _air_density(getParam<Real>("air_density")),
    _water_vapor(getParam<Real>("water_vapor")),
    _jfac(getParam<Real>("jfac")),
    _t(0.0),
    _dirty(true)
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
  // All indices are valid species indices (padded with ONE=0, conc=1).
  // Per.14: _scratch_G / _scratch_rates are mutable members pre-allocated once,
  // avoiding heap allocation on every dC/dt call (~272 KiB per call for full MCM).
  _scratch_G.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _scratch_G[r] = C[_iG[r][0]] * C[_iG[r][1]] * C[_iG[r][2]];

  // Step 2: rates = k .* G
  _scratch_rates.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _scratch_rates[r] = _k[r] * _scratch_G[r];

  // Step 3: dC = f^T * rates  — format-agnostic iteration via StoichMatrix.
  // forEachInRow dispatches based on _stoich.format (CSR or COO); the lambda
  // is fully inlined, so there is zero virtual-call overhead.
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real rate = _scratch_rates[r];
    _stoich.forEachInRow(r, [&](int s, Real coeff) {
      dC[s] += coeff * rate;
    });
  }
}

void
MCMBoxModel::computeJacobianTriplets(
    const std::vector<Real> & C,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  J.clear();

  if (_n_reactions == 0 || _n_species == 0)
    return;

  // F0AM-style 3-term product Jacobian — format-agnostic via StoichMatrix.
  //   rate_r = k_r * C[i0] * C[i1] * C[i2]
  //   d(rate_r)/dC[j] = k_r * sum_{k where iG[k]==j} prod_{m≠k} C[iG[m]]
  //   J[s][j] += f[r][s] * d(rate_r)/dC[j]

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    const Real c0 = C[i0], c1 = C[i1], c2 = C[i2];
    const Real k = _k[r];

    // Helper: for reactant j with derivative drate, emit J(s,j) for every
    // species s that participates in reaction r.
    auto emit_contrib = [&](unsigned int j, Real drate) {
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        J.emplace_back((unsigned int)s, j, drate * coeff);
      });
    };

    // Contribution from C[i0]: dr/dC[i0] = k * c1 * c2
    emit_contrib((unsigned int)i0, k * c1 * c2);
    // Contribution from C[i1]: dr/dC[i1] = k * c0 * c2
    if (i1 != i0) emit_contrib((unsigned int)i1, k * c0 * c2);
    // Contribution from C[i2]: dr/dC[i2] = k * c0 * c1
    if (i2 != i0 && i2 != i1) emit_contrib((unsigned int)i2, k * c0 * c1);
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
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        Real val = drate * coeff;
        unsigned int us = (unsigned int)s;
        if (j == us)
          _cached_diag_J[us] += val;
        else
        {
          uint64_t key = (static_cast<uint64_t>(us) << 32) | static_cast<uint64_t>(j);
          _cached_offdiag_J[key] += val;
        }
      });
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

  // Build RO2 species index list from parser's explicit ro2_species
  _ro2_indices.clear();
  for (const auto & ro2_name : mech.ro2_species)
  {
    auto it = std::find(_species_names.begin(), _species_names.end(), ro2_name);
    if (it != _species_names.end())
      _ro2_indices.push_back((unsigned int)(it - _species_names.begin()));
  }

  // Build stoichiometric matrix in the selected format (parameter "stoich_format").
  // CSR: compact, PETSc AIJ-compatible, optimal for HPC.
  // COO: AtChem2-style split reactant/product — enables loss/production diagnostics.
  StoichMatrix::Format fmt = StoichMatrix::CSR;
  if (isParamValid("stoich_format"))
  {
    MooseEnum fmt_enum = getParam<MooseEnum>("stoich_format");
    if (fmt_enum == "COO")        fmt = StoichMatrix::COO;
    else if (fmt_enum == "DENSE") fmt = StoichMatrix::DENSE;
    else if (fmt_enum == "CSC")   fmt = StoichMatrix::CSC;
  }
  _stoich.build(mech, fmt);

  // Copy reactant indices — convert vector<vector<int>> to flat array<int,3>
  // Per.16: single contiguous allocation vs 17k independent heap allocations.
  _iG.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const auto & src = mech.reactant_indices[r];
    _iG[r] = {src[0], src[1], src[2]};
  }

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

  // Pre-compute J photo indices into _func_params (Per.14 — avoids string+map in evaluateCoefficients).
  // Only store indices for J numbers that were actually registered in _name_to_index.
  // (The parser may detect J numbers in mech.j_numbers that don't appear in the
  // coefficient/reaction expressions scanned by setupFparser, so we match the old
  // defensive find()-based approach.)
  _j_photo_indices.clear();
  _j_photo_indices.reserve(_j_numbers.size());
  for (auto & jn : _j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(jn);
    auto it = _name_to_index.find(jname);
    if (it != _name_to_index.end())
      _j_photo_indices.push_back(it->second);
    else
      _j_photo_indices.push_back((unsigned int)-1); // sentinel: skip in evaluateCoefficients
  }

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
  _func_params[2] = 0.21 * _air_density;  // O2 volume fraction
  _func_params[3] = 0.78 * _air_density;  // N2 volume fraction
  _func_params[4] = _water_vapor;

  // Compute photolysis J values from solar zenith angle (MCM formula)
  // J = l * cosx^m * exp(-n * secx) * JFAC * ROOF
  // ROOF = CLOSED (0) or OPEN (1); JFAC ∈ [0,1]
  if (!_j_CL_vals.empty())
  {
    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    Real cosx = calculateCosSZA(_t);
    Real secx = (cosx > 1.0e-10) ? (1.0 / cosx) : 1.0e2;
    // Per.14: use pre-computed _j_photo_indices instead of string+map lookup.
    // Sentinel value (unsigned int)-1 means J number was not registered; skip it.
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      if (cosx > 1.0e-10)
        _func_params[idx] = _j_CL_vals[i] * std::pow(cosx, _j_CMM_vals[i])
                          * std::exp(-_j_CNN_vals[i] * secx) * _jfac * roof_factor;
      else
        _func_params[idx] = 0.0;
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
  // RO2 = sum of peroxy radicals (using explicit species list from parser)
  auto it_ro2 = _name_to_index.find("RO2");
  if (it_ro2 != _name_to_index.end())
  {
    Real ro2_sum = 0.0;
    for (auto idx : _ro2_indices)
      ro2_sum += (_cached_C.empty() ? 0.0 :
          (idx < _cached_C.size() ? _cached_C[idx] : 0.0));
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
  constexpr Real pi = 3.14159265358979323846;
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
  return _stoich.get(r, s) * reactionRate(r, C);
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
    updatePhotolysisSZA(sza, _jfac);
}

void
MCMBoxModel::calcJFAC(const std::string & ref_j_name, Real constrained_val)
{
  if (constrained_val <= 0.0)
  {
    _jfac = 0.0;
    return;
  }

  // Parse J number from name like "J4"
  unsigned int ref_jn = 0;
  if (ref_j_name.size() > 1 && ref_j_name[0] == 'J')
    ref_jn = (unsigned int)std::stoi(ref_j_name.substr(1));

  // Find the photolysis parameters for this reference J
  Real cl_val = 0.0, cmm_val = 1.0, cnn_val = 0.0;
  for (size_t i = 0; i < _j_numbers.size(); ++i)
    if (_j_numbers[i] == ref_jn)
    {
      cl_val = _j_CL_vals[i];
      cmm_val = _j_CMM_vals[i];
      cnn_val = _j_CNN_vals[i];
      break;
    }

  if (cl_val == 0.0)
  {
    mooseWarning("MCMBoxModel: JFAC reference ", ref_j_name, " has no photolysis parameters");
    _jfac = 1.0;
    return;
  }

  // Compute parameterized value at current SZA
  Real cosx = calculateCosSZA(_t);
  if (cosx <= 1.0e-10)
  {
    _jfac = 0.0;
    return;
  }
  Real secx = 1.0 / cosx;
  Real j_calc = cl_val * std::pow(cosx, cmm_val) * std::exp(-cnn_val * secx);

  if (j_calc <= 0.0)
    _jfac = 0.0;
  else
    _jfac = constrained_val / j_calc;  // JFAC = observed / calculated
}

void
MCMBoxModel::updatePhotolysisSZA(Real sza, Real jfac)
{
  if (!_roof_open) return;  // ROOF CLOSED: J rates handled in evaluateCoefficients
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
  _lat = lat * M_PI / 180.0;
  _lon = lon * M_PI / 180.0;
  _day = day;
  _month = month;
  _year = year;

  unsigned int doy = computeDayOfYear();
  unsigned int days_in_year =
      ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) ? 366 : 365;
  Real theta = 2.0 * M_PI * (Real)doy / (Real)days_in_year;
  Real B = 2.0 * M_PI * ((Real)doy - 1.0) / (Real)days_in_year;

  _declination = 0.006918 - 0.399912 * std::cos(theta) + 0.070257 * std::sin(theta) -
                 0.006758 * std::cos(2.0 * theta) + 0.000907 * std::sin(2.0 * theta) -
                 0.002697 * std::cos(3.0 * theta) + 0.001480 * std::sin(3.0 * theta);
  _eot = 0.165 * std::sin(2.0 * B) - 0.126 * std::cos(B) - 0.025 * std::sin(B);
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

Real
MCMBoxModel::getRO2Sum(const std::vector<Real> & C) const
{
  Real sum = 0.0;
  for (auto idx : _ro2_indices)
    if (idx < C.size())
      sum += C[idx];
  return sum;
}

Real
MCMBoxModel::getJValue(unsigned int j_number) const
{
  // Look up J value from fparser parameter buffer using PHOTOJ<N> key
  std::string jname = "PHOTOJ" + std::to_string(j_number);
  auto it = _name_to_index.find(jname);
  if (it != _name_to_index.end() && it->second < _func_params.size())
    return _func_params[it->second];
  return 0.0;
}

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
