//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRuntimeMechanism.h"
#include "MCMFacsimileParser.h"
#include "HybridJTableReader.h"
#include "BottomUpJIntegrator.h"
#include "JCalibrator.h"
#include "pcrecpp.h"

#include <algorithm>
#include <regex>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

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

// ---- MCMRuntimeMechanism ----------------------------------------------------

MCMRuntimeMechanism::MCMRuntimeMechanism(const ParsedMechanism & mech,
                                           bool use_limiting_reagent,
                                           StoichMatrix::Format stoich_format,
                                           bool disable_fpoptimizer,
                                           bool enable_jit,
                                           bool fail_on_bad_deps,
                                           unsigned int max_function_recurse)
  : IMechanism(),
    FunctionParserUtils<false>([=]() {
      InputParameters p = FunctionParserUtils<false>::validParams();
      p.set<bool>("disable_fpoptimizer") = disable_fpoptimizer;
      p.set<bool>("enable_jit") = enable_jit;
      p.set<bool>("enable_ad_cache") = true;
      p.set<bool>("enable_auto_optimize") = true && !disable_fpoptimizer;
      return p;
    }()),
    _n_species(0),
    _n_reactions(0),
    _units_ppb(false),
    _ppb_to_molec(1.0),
    _use_limiting_reagent(use_limiting_reagent),
    _temperature(298.15),
    _air_density(2.46e19),
    _water_vapor(2.46e17),
    _press(0.0),
    _rh(-1.0),
    _blheight(0.0),
    _jfac(1.0),
    _photolysis_method(MCM_SZA),
    _j_index_start(0),
    _n_j_vars(0),
    _roof_open(true),
    _lat(51.51),
    _lon(0.13),
    _declination(0.0),
    _eot(0.0),
    _day(21),
    _month(6),
    _year(2010),
    _solar_cosx(0.0),
    _solar_secx(1.0e2),
    _solar_lha(0.0),
    _solar_sinld(0.0),
    _solar_cosld(0.0),
    _solar_eqt(0.0),
    _t(0.0),
    _dirty(true),
    _cached_bottomup_T(0.0),
    _cached_bottomup_P(0.0),
    _bottomup_j_valid(false)
{
  loadMechanism(mech, use_limiting_reagent, stoich_format);
}

void
MCMRuntimeMechanism::loadMechanism(const ParsedMechanism & mech,
                                    bool use_limiting_reagent,
                                    StoichMatrix::Format stoich_format)
{
  _n_species = mech.species.size();
  _n_reactions = mech.reactions.size();
  _species_names = mech.species;
  _reaction_names = mech.reaction_names;

  // Build RO2 species index list and name list from parser's explicit ro2_species
  _ro2_indices.clear();
  _ro2_species_names = mech.ro2_species;
  for (const auto & ro2_name : mech.ro2_species)
  {
    auto it = std::find(_species_names.begin(), _species_names.end(), ro2_name);
    if (it != _species_names.end())
      _ro2_indices.push_back((unsigned int)(it - _species_names.begin()));
  }

  // Build stoichiometric matrix in the selected format.
  _stoich.build(mech, stoich_format);

  // Copy reactant indices — convert vector<vector<int>> to flat array<int,3>
  _iG.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const auto & src = mech.reactant_indices[r];
    _iG[r] = {src[0], src[1], src[2]};
  }

  // Load limiting-reagent flags (F0AM RO2 termination reactions)
  _use_limiting_reagent = use_limiting_reagent;
  _limiting_reagent = mech.is_limiting_reagent;
  _limiting_reactant = mech.limiting_reactant;

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
    // Complex case: use fparser
    setupFparser(mech);
  }

  _dirty = true;
}

void
MCMRuntimeMechanism::setupFparser(const ParsedMechanism & mech)
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

  // Build name-to-index mapping for ALL possible variables.
  _name_to_index["TEMP"] = 0; _name_to_index["M"] = 1;
  _name_to_index["O2"] = 2; _name_to_index["N2"] = 3;
  _name_to_index["H2O"] = 4;

  for (unsigned int i = 0; i < coeff_names.size(); ++i)
    _name_to_index[coeff_names[i]] = 5 + i;
  for (unsigned int i = 0; i < _n_species; ++i)
    _name_to_index[_species_names[i]] = 5 + coeff_names.size() + i;

  bool has_ro2_in_species = false;
  for (auto & s : _species_names)
    if (s == "RO2") { has_ro2_in_species = true; break; }
  if (!has_ro2_in_species)
    _name_to_index["RO2"] = 5 + coeff_names.size() + _n_species;
  unsigned int n_extra_vars = has_ro2_in_species ? 0 : 1;
  _j_index_start = 5 + coeff_names.size() + _n_species + n_extra_vars;

  for (auto & n : j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(n);
    _name_to_index[jname] = _j_index_start + _name_to_index.size() - _j_index_start;
  }

  unsigned int n_j = j_numbers.size();
  _func_params.resize(5 + coeff_names.size() + _n_species + n_extra_vars + n_j, 0.0);

  // ---- Parse coefficients using ParseAndDeduceVariables ----
  _coeff_parsers.resize(coeff_names.size());
  _coeff_var_indices.resize(coeff_names.size());
  _coeff_local_params.resize(coeff_names.size());

  for (unsigned int i = 0; i < coeff_names.size(); ++i)
  {
    _coeff_parsers[i] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_coeff_parsers[i]);

    std::string discoveredVars;
    int nFound = 0;
    if (_coeff_parsers[i]->ParseAndDeduceVariables(coeff_exprs[i], discoveredVars, &nFound) >= 0)
      throw std::runtime_error("MCMRuntimeMechanism: Bad coefficient '" + coeff_names[i] +
                               "': " + coeff_exprs[i]);

    // Build per-parser index into _func_params
    auto & indices = _coeff_var_indices[i];
    indices.reserve(nFound);
    std::istringstream vstream(discoveredVars);
    std::string vname;
    while (std::getline(vstream, vname, ','))
    {
      auto it = _name_to_index.find(vname);
      if (it != _name_to_index.end())
        indices.push_back(it->second);
      else
      {
        std::cerr << "Warning: MCMRuntimeMechanism: coefficient '" << coeff_names[i]
                  << "' references undefined variable '" << vname
                  << "'. Auto-registering with value 0." << std::endl;
        unsigned int idx = _func_params.size();
        _func_params.push_back(0.0);
        _name_to_index[vname] = idx;
        indices.push_back(idx);
      }
    }
    _coeff_local_params[i].resize(indices.size());

    if (!_disable_fpoptimizer)
      _coeff_parsers[i]->Optimize();
  }

  // ---- Parse reaction expressions using ParseAndDeduceVariables ----
  _reaction_parsers.resize(_n_reactions);
  _reaction_var_indices.resize(_n_reactions);
  _reaction_local_params.resize(_n_reactions);

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    _reaction_parsers[r] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_reaction_parsers[r]);

    std::string discoveredVars;
    int nFound = 0;
    if (_reaction_parsers[r]->ParseAndDeduceVariables(rxn_exprs[r], discoveredVars, &nFound) >= 0)
      throw std::runtime_error("MCMRuntimeMechanism: Bad reaction " + std::to_string(r) +
                               ": " + rxn_exprs[r]);

    auto & indices = _reaction_var_indices[r];
    indices.reserve(nFound);
    std::istringstream vstream(discoveredVars);
    std::string vname;
    while (std::getline(vstream, vname, ','))
    {
      auto it = _name_to_index.find(vname);
      if (it != _name_to_index.end())
        indices.push_back(it->second);
      else
      {
        std::cerr << "Warning: MCMRuntimeMechanism: reaction " << r
                  << " references undefined variable '" << vname
                  << "'. Auto-registering with value 0." << std::endl;
        unsigned int idx = _func_params.size();
        _func_params.push_back(0.0);
        _name_to_index[vname] = idx;
        indices.push_back(idx);
      }
    }
    _reaction_local_params[r].resize(indices.size());

    if (!_disable_fpoptimizer)
      _reaction_parsers[r]->Optimize();
  }

  // Store photolysis parameters
  if (!mech.j_numbers.empty())
  {
    _j_numbers = mech.j_numbers;
    _j_CL_vals = mech.j_CL;
    _j_CMM_vals = mech.j_CMM;
    _j_CNN_vals = mech.j_CNN;
  }
  else
  {
    _j_numbers.assign(j_numbers.begin(), j_numbers.end());
    _j_CL_vals.assign(_j_numbers.size(), 0.0);
    _j_CMM_vals.assign(_j_numbers.size(), 1.0);
    _j_CNN_vals.assign(_j_numbers.size(), 0.0);
  }
  _n_j_vars = _j_numbers.size();

  // Pre-compute J photo indices into _func_params
  _j_photo_indices.clear();
  _j_photo_indices.reserve(_j_numbers.size());
  for (auto & jn : _j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(jn);
    auto it = _name_to_index.find(jname);
    if (it != _name_to_index.end())
      _j_photo_indices.push_back(it->second);
    else
      _j_photo_indices.push_back((unsigned int)-1);
  }

  // Build fast pre-compiled handlers for coefficient and reaction expressions.
  _coeff_fast.resize(_coeff_parsers.size());
  for (unsigned int i = 0; i < _coeff_parsers.size(); ++i)
    _coeff_fast[i] = compileFastHandler(coeff_exprs[i], _coeff_var_indices[i]);

  _reaction_fast.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _reaction_fast[r] = compileFastHandler(rxn_exprs[r], _reaction_var_indices[r]);
}

MCMRuntimeMechanism::FastHandler
MCMRuntimeMechanism::compileFastHandler(const std::string & expr,
                                         const std::vector<unsigned int> & var_indices) const
{
  // Pattern 1: Single variable reference — most reaction expressions (e.g. "KMT01").
  if (var_indices.size() == 1 && std::regex_match(expr, std::regex("^[A-Za-z_][A-Za-z0-9_]*$")))
  {
    unsigned int idx = var_indices[0];
    return [idx](const std::vector<Real> & p) { return p[idx]; };
  }

  // Pattern 2: Simple numeric constant.
  {
    char * end = nullptr;
    double val = std::strtod(expr.c_str(), &end);
    if (end && *end == '\0')
      return [val](const std::vector<Real> &) { return val; };
  }

  // Pattern 3: Simple Arrhenius: A*exp(B/TEMP)
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) { return A * std::exp(B / p[0]); };
    }
  }

  // Pattern 4: Modified Arrhenius: A*(TEMP/300)^B*exp(C/TEMP)
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*\(TEMP/300\)\^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      double C = std::stod(m[3].str());
      return [A, B, C](const std::vector<Real> & p) {
        return A * std::pow(p[0] / 300.0, B) * std::exp(C / p[0]);
      };
    }
  }

  // Pattern 5: Power temp: A*(TEMP/300)^B
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*\(TEMP/300\)\^([0-9.eE+\-]+)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) {
        return A * std::pow(p[0] / 300.0, B);
      };
    }
  }

  // Pattern 6: Arrhenius with M: A*exp(B/TEMP)*M
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)\*M$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) {
        return A * std::exp(B / p[0]) * p[1];
      };
    }
  }

  // No pattern matched — fall back to fparser
  return nullptr;
}

void
MCMRuntimeMechanism::evaluateCoefficients()
{
  if (_coeff_parsers.empty()) return;

  // Compute M (air density) dynamically from press/temp if press > 0
  Real M_val = _air_density;
  if (_press > 0.0)
  {
    constexpr Real NA_over_R = 6.02214129e23 / 8.3144621;
    M_val = 1.0e-6 * NA_over_R * (_press * 100.0 / _temperature);
  }

  // Update ppb conversion factor from the actual air density used
  _ppb_to_molec = (_units_ppb) ? M_val / 1.0e9 : 1.0;

  // Compute H2O dynamically from rh/temp/press if rh >= 0
  Real H2O_val = _water_vapor;
  if (_rh >= 0.0)
  {
    Real temp_c = _temperature - 273.15;
    Real wvp = (_rh / 100.0) * 6.116441 * std::pow(10.0, (7.591386 * temp_c) / (temp_c + 240.7263));
    Real press_mbar = (_press > 0.0) ? _press : 1013.25;
    Real h2o_ppu = wvp / (press_mbar - wvp);
    H2O_val = h2o_ppu * M_val;
  }

  _func_params[0] = _temperature;
  _func_params[1] = M_val;
  _func_params[2] = 0.21 * M_val;  // O2 volume fraction
  _func_params[3] = 0.78 * M_val;  // N2 volume fraction
  _func_params[4] = H2O_val;

  // Compute photolysis J values into _func_params.
  if (_photolysis_method == BOTTOMUP && _bottomup_integrator)
  {
    const Real T_cur = _temperature;
    const Real P_cur = _press > 0 ? _press : 1013.25;
    if (!_bottomup_j_valid ||
        std::abs(T_cur - _cached_bottomup_T) > 1.0e-6 ||
        std::abs(P_cur - _cached_bottomup_P) > 1.0e-6)
    {
      _cached_bottomup_j = _bottomup_integrator->computeAllJ(T_cur, P_cur);
      _cached_bottomup_T = T_cur;
      _cached_bottomup_P = P_cur;
      _bottomup_j_valid = true;
    }
    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    for (size_t i = 0; i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      std::string jname = "J" + std::to_string(jn);
      auto it = _cached_bottomup_j.find(jname);
      _func_params[idx] = (it != _cached_bottomup_j.end()) ? it->second * _jfac * roof_factor : 0.0;
    }
  }
  else if (!_j_CL_vals.empty())
  {
    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    Real cosx = calculateCosSZA(_t);
    Real secx = _solar_secx;

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

  // Auto-calibrate J values from observed data
  if (_jcalibrator && _jcalibrator->hasObservedData())
  {
    std::map<unsigned int, Real> param_J;
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      param_J[jn] = _func_params[idx];
    }
    _jcalibrator->calibrate(param_J);
    _jcalibrator->applyTo(param_J);
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      auto it = param_J.find(jn);
      if (it != param_J.end())
        _func_params[idx] = it->second;
    }
  }

  // Populate ALL _func_params slots with species and RO2 values
  for (unsigned int s = 0; s < _n_species; ++s)
  {
    auto it = _name_to_index.find(_species_names[s]);
    if (it != _name_to_index.end())
      _func_params[it->second] = (_cached_C.empty() ? 0.0 :
          (s < _cached_C.size() ? _cached_C[s] : 0.0));
  }
  auto it_ro2 = _name_to_index.find("RO2");
  if (it_ro2 != _name_to_index.end())
  {
    Real ro2_sum = 0.0;
    for (auto idx : _ro2_indices)
      ro2_sum += (_cached_C.empty() ? 0.0 :
          (idx < _cached_C.size() ? _cached_C[idx] : 0.0));
    _func_params[it_ro2->second] = ro2_sum;
  }

  // Evaluate coefficients in topological order
  unsigned int n_coeff = _coeff_parsers.size();
  for (unsigned int i = 0; i < n_coeff; ++i)
  {
    Real val;
    if (_coeff_fast[i])
      val = _coeff_fast[i](_func_params);
    else
    {
      const auto & indices = _coeff_var_indices[i];
      auto & local = _coeff_local_params[i];
      for (size_t j = 0; j < indices.size(); ++j)
        local[j] = _func_params[indices[j]];
      val = evaluate(_coeff_parsers[i], local);
    }
    if (std::isnan(val) || std::isinf(val)) val = 0.0;
    _func_params[5 + i] = val;
  }

  // Evaluate reaction rate expressions → _k
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real val;
    if (_reaction_fast[r])
      val = _reaction_fast[r](_func_params);
    else
    {
      const auto & indices = _reaction_var_indices[r];
      auto & local = _reaction_local_params[r];
      for (size_t j = 0; j < indices.size(); ++j)
        local[j] = _func_params[indices[j]];
      val = evaluate(_reaction_parsers[r], local);
    }
    _k[r] = (std::isnan(val) || std::isinf(val)) ? 0.0 : val;
  }
}

void
MCMRuntimeMechanism::evaluateCoefficients(const std::vector<Real> & C)
{
  _cached_C = C;
  evaluateCoefficients();
}

// ---- IMechanism interface implementations ----

void
MCMRuntimeMechanism::updateParams(const PhysParams & params)
{
  _temperature = params.temperature;
  _air_density = params.air_density;
  _water_vapor = params.water_vapor;
  _press = params.pressure;
  _rh = params.rh;
  _jfac = params.jfac;
  _blheight = params.blheight;
  _lat = params.latitude;
  _lon = params.longitude;
}

void
MCMRuntimeMechanism::computeRHS(Real t,
                                 const std::vector<Real> & C,
                                 const PhysParams & params,
                                 std::vector<Real> & dC_dt) const
{
  // Update physical parameters
  const_cast<MCMRuntimeMechanism *>(this)->updateParams(params);
  const_cast<MCMRuntimeMechanism *>(this)->_t = t;

  // Compute dC/dt
  computeDCdt(C, dC_dt);
}

void
MCMRuntimeMechanism::computeJacobian(
    Real t,
    const std::vector<Real> & C,
    const PhysParams & params,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  // Update physical parameters
  const_cast<MCMRuntimeMechanism *>(this)->updateParams(params);
  const_cast<MCMRuntimeMechanism *>(this)->_t = t;

  // Compute Jacobian
  computeJacobianTriplets(C, J);
}

SpeciesRates
MCMRuntimeMechanism::computeSpeciesRates(Real t,
                                          const std::vector<Real> & C,
                                          const PhysParams & params) const
{
  // Update physical parameters
  const_cast<MCMRuntimeMechanism *>(this)->updateParams(params);
  const_cast<MCMRuntimeMechanism *>(this)->_t = t;

  // Evaluate rate coefficients
  const_cast<MCMRuntimeMechanism *>(this)->evaluateCoefficients(C);

  SpeciesRates rates;
  rates.production.resize(_n_species, 0.0);
  rates.loss.resize(_n_species, 0.0);

  for (unsigned int s = 0; s < _n_species; ++s)
  {
    rates.production[s] = speciesProductionRate(s, C);
    rates.loss[s] = speciesLossRate(s, C);
  }

  return rates;
}

// ---- Core computation methods ----

void
MCMRuntimeMechanism::computeDCdt(const std::vector<Real> & C, std::vector<Real> & dC) const
{
  dC.assign(_n_species, 0.0);

  if (_n_reactions == 0 || _n_species == 0)
    return;

  const_cast<MCMRuntimeMechanism *>(this)->evaluateCoefficients(C);

  // Step 1: compute reactant products G[i] = C[iG[i][0]] * C[iG[i][1]] * C[iG[i][2]]
  _scratch_G.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real c0 = (_iG[r][0] >= 0) ? C[_iG[r][0]] : 1.0;
    Real c1 = (_iG[r][1] >= 0) ? C[_iG[r][1]] : 1.0;
    Real c2 = (_iG[r][2] >= 0) ? C[_iG[r][2]] : 1.0;
    _scratch_G[r] = c0 * c1 * c2;
  }

  // Limiting-reagent override (F0AM RO2 termination reactions)
  if (_use_limiting_reagent && !_limiting_reagent.empty())
  {
    for (unsigned int r = 0; r < _n_reactions; ++r)
    {
      if (!_limiting_reagent[r])
        continue;
      const int i0 = _iG[r][0], i1 = _iG[r][1];
      if (i0 >= 0 && i1 >= 0)
      {
        Real c0 = C[i0], c1 = C[i1];
        Real minc = std::min(c0, c1);
        _scratch_G[r] = minc;
      }
    }
  }

  // Step 2: rates = k .* G
  _scratch_rates.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _scratch_rates[r] = _k[r] * _scratch_G[r];

  // Step 3: dC = f^T * rates
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real rate = _scratch_rates[r];
    _stoich.forEachInRow(r, [&](int s, Real coeff) {
      dC[s] += coeff * rate;
    });
  }
}

void
MCMRuntimeMechanism::computeJacobianTriplets(
    const std::vector<Real> & C,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  J.clear();

  if (_n_reactions == 0 || _n_species == 0)
    return;

  const_cast<MCMRuntimeMechanism *>(this)->evaluateCoefficients(C);

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    const Real c0 = (i0 >= 0) ? C[i0] : 1.0;
    const Real c1 = (i1 >= 0) ? C[i1] : 1.0;
    const Real c2 = (i2 >= 0) ? C[i2] : 1.0;
    const Real k = _k[r];

    auto emit_contrib = [&](unsigned int j, Real drate) {
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        J.emplace_back((unsigned int)s, j, drate * coeff);
      });
    };

    // Limiting-reagent (LR) Jacobian
    if (_use_limiting_reagent && !_limiting_reagent.empty() && r < _limiting_reagent.size() && _limiting_reagent[r])
    {
      if (i0 >= 0 && i1 >= 0)
      {
        unsigned int min_idx = (c0 <= c1) ? (unsigned int)i0 : (unsigned int)i1;
        Real drate = k;
        emit_contrib(min_idx, drate);
      }
      continue;
    }

    // Standard Jacobian: 3-term product derivative
    if (i0 >= 0) emit_contrib((unsigned int)i0, k * c1 * c2);
    if (i1 >= 0 && i1 != i0) emit_contrib((unsigned int)i1, k * c0 * c2);
    if (i2 >= 0 && i2 != i0 && i2 != i1) emit_contrib((unsigned int)i2, k * c0 * c1);
  }
}

// ---- Cached single-species interface ----

Real
MCMRuntimeMechanism::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  if ((_dirty) || _cached_dC.size() != _n_species)
  {
    _cached_C = C;
    _cached_dC.resize(_n_species);
    computeDCdt(C, _cached_dC);
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_dC[idx] : 0.0;
}

Real
MCMRuntimeMechanism::getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const
{
  if (_dirty || _cached_diag_J.size() != _n_species)
  {
    _cached_C = C;
    const_cast<MCMRuntimeMechanism *>(this)->evaluateCoefficients(C);
    _buildJacobianCache();
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_diag_J[idx] : 0.0;
}

Real
MCMRuntimeMechanism::getJacobianOffDiagonal(unsigned int i, unsigned int j, const std::vector<Real> & C) const
{
  if (_dirty || _cached_od_row_ptr.size() != _n_species + 1)
  {
    _cached_C = C;
    const_cast<MCMRuntimeMechanism *>(this)->evaluateCoefficients(C);
    _buildJacobianCache();
    _dirty = false;
  }
  size_t lo = _cached_od_row_ptr[i];
  size_t hi = _cached_od_row_ptr[i + 1];
  while (lo < hi)
  {
    size_t mid = lo + (hi - lo) / 2;
    unsigned int col = _cached_od_cols[mid];
    if (col == j)
      return _cached_od_vals[mid];
    else if (col < j)
      lo = mid + 1;
    else
      hi = mid;
  }
  return 0.0;
}

void
MCMRuntimeMechanism::_buildJacobianCache() const
{
  _cached_diag_J.assign(_n_species, 0.0);

  if (_n_reactions == 0 || _n_species == 0)
  {
    _cached_od_row_ptr.assign(_n_species + 1, 0);
    _cached_od_cols.clear();
    _cached_od_vals.clear();
    return;
  }

  std::vector<std::vector<std::pair<unsigned int, Real>>> temp(_n_species);

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    const Real c0 = (i0 >= 0) ? _cached_C[i0] : 1.0;
    const Real c1 = (i1 >= 0) ? _cached_C[i1] : 1.0;
    const Real c2 = (i2 >= 0) ? _cached_C[i2] : 1.0;
    const Real k = _k[r];

    auto accum = [&](unsigned int j, Real drate) {
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        Real val = drate * coeff;
        unsigned int us = (unsigned int)s;
        if (j == us)
          _cached_diag_J[us] += val;
        else
          temp[us].emplace_back(j, val);
      });
    };

    bool is_lr = (_use_limiting_reagent && !_limiting_reagent.empty() && r < _limiting_reagent.size() && _limiting_reagent[r]);
    if (is_lr && i0 >= 0 && i1 >= 0)
    {
      unsigned int min_idx = (c0 <= c1) ? (unsigned int)i0 : (unsigned int)i1;
      Real drate = k;
      accum(min_idx, drate);
      continue;
    }

    accum((unsigned int)i0, k * c1 * c2);
    if (i1 != i0) accum((unsigned int)i1, k * c0 * c2);
    if (i2 != i0 && i2 != i1) accum((unsigned int)i2, k * c0 * c1);
  }

  // Flatten temp to CSR
  _cached_od_row_ptr.resize(_n_species + 1);
  _cached_od_row_ptr[0] = 0;
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    auto & row = temp[i];
    if (row.size() <= 1)
    {
      if (row.size() == 1)
        _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i] + 1;
      else
        _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i];
      continue;
    }
    std::sort(row.begin(), row.end());
    size_t w = 0;
    for (size_t k = 1; k < row.size(); ++k)
    {
      if (row[w].first == row[k].first)
        row[w].second += row[k].second;
      else
        row[++w] = row[k];
    }
    row.resize(w + 1);
    _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i] + row.size();
  }

  size_t nnz = _cached_od_row_ptr[_n_species];
  _cached_od_cols.resize(nnz);
  _cached_od_vals.resize(nnz);
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    size_t base = _cached_od_row_ptr[i];
    for (size_t k = 0; k < temp[i].size(); ++k)
    {
      _cached_od_cols[base + k] = temp[i][k].first;
      _cached_od_vals[base + k] = temp[i][k].second;
    }
  }
}

// ---- Reaction rate queries ----

Real
MCMRuntimeMechanism::reactionRate(unsigned int r, const std::vector<Real> & C) const
{
  if (r >= _n_reactions || C.size() != _n_species) return 0.0;
  Real c0 = (_iG[r][0] >= 0) ? C[_iG[r][0]] : 1.0;
  Real c1 = (_iG[r][1] >= 0) ? C[_iG[r][1]] : 1.0;
  Real c2 = (_iG[r][2] >= 0) ? C[_iG[r][2]] : 1.0;
  // Limiting-reagent branch: rate = k * min(c0, c1) (F0AM-compatible)
  if (_use_limiting_reagent && !_limiting_reagent.empty() &&
      r < _limiting_reagent.size() && _limiting_reagent[r])
    return _k[r] * std::min(c0, c1);
  return _k[r] * c0 * c1 * c2;
}

Real
MCMRuntimeMechanism::speciesReactionRate(unsigned int s, unsigned int r, const std::vector<Real> & C) const
{
  return _stoich.get(r, s) * reactionRate(r, C);
}

void
MCMRuntimeMechanism::allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const
{
  rates.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rates[r] = reactionRate(r, C);
}

Real
MCMRuntimeMechanism::speciesLossRate(unsigned int s, const std::vector<Real> & C) const
{
  Real total = 0.0;
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real coeff = _stoich.get(r, s);
    if (coeff < 0.0)
      total += (-coeff) * reactionRate(r, C);
  }
  return total;
}

Real
MCMRuntimeMechanism::speciesProductionRate(unsigned int s, const std::vector<Real> & C) const
{
  Real total = 0.0;
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real coeff = _stoich.get(r, s);
    if (coeff > 0.0)
      total += coeff * reactionRate(r, C);
  }
  return total;
}

Real
MCMRuntimeMechanism::getRO2Sum(const std::vector<Real> & C) const
{
  Real sum = 0.0;
  for (auto idx : _ro2_indices)
    if (idx < C.size())
      sum += C[idx];
  return sum;
}

// ---- Solar / time helpers ----

unsigned int
MCMRuntimeMechanism::computeDayOfYear() const
{
  unsigned int days_in_months[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0)
    days_in_months[1] = 29;
  unsigned int doy = 0;
  for (unsigned int m = 0; m < (unsigned int)(_month - 1); ++m)
    doy += days_in_months[m];
  doy += _day - 1;
  return doy;
}

Real
MCMRuntimeMechanism::calculateCosSZA(Real t) const
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

  _declination = dec;
  _eot = eqt;
  _solar_lha = lha;
  _solar_sinld = sin(lat_rad) * sin(dec);
  _solar_cosld = cos(lat_rad) * cos(dec);
  _solar_eqt = eqt;

  if (cosx <= 0.0) cosx = 0.0;
  Real secx = (cosx > 1.0e-10) ? (1.0 / cosx) : 1.0e2;
  _solar_cosx = cosx;
  _solar_secx = secx;
  return cosx;
}

void
MCMRuntimeMechanism::setSolarParams(Real lat, Real lon, int day, int month, int year)
{
  _lat = lat;
  _lon = lon;
  _day = day;
  _month = month;
  _year = year;
}

Real
MCMRuntimeMechanism::cosSZA(Real seconds) const
{
  Real hour = seconds / 3600.0;
  Real ha = (hour - 12.0 + _eot + _lon * 12.0 / M_PI) * M_PI / 12.0;
  return std::sin(_lat * M_PI / 180.0) * std::sin(_declination)
       + std::cos(_lat * M_PI / 180.0) * std::cos(_declination) * std::cos(ha);
}

// ---- Photolysis setup ----

void
MCMRuntimeMechanism::enableHybridPhotolysis(const std::string & table_dir)
{
  _photolysis_method = HYBRID;
  _hybrid_reader = std::make_unique<HybridJTableReader>(table_dir);
}

void
MCMRuntimeMechanism::loadBottomUpData(const std::string & data_dir, const std::string & flux_file)
{
  _photolysis_method = BOTTOMUP;
  _bottomup_integrator = std::make_unique<BottomUpJIntegrator>(data_dir);
  _bottomup_integrator->loadLampFlux(flux_file);
  _bottomup_integrator->loadReactionMap("bottomup_jmap.dat");
}

Real
MCMRuntimeMechanism::getJValue(unsigned int j_number) const
{
  std::string jname = "PHOTOJ" + std::to_string(j_number);
  auto it = _name_to_index.find(jname);
  if (it != _name_to_index.end() && it->second < _func_params.size())
    return _func_params[it->second];

  for (size_t i = 0; i < _j_numbers.size(); ++i)
  {
    if (_j_numbers[i] == j_number)
    {
      const Real roof_factor = _roof_open ? 1.0 : 0.0;
      Real cosx = _solar_cosx;
      if (cosx > 1.0e-10)
        return _j_CL_vals[i] * std::pow(cosx, _j_CMM_vals[i])
               * std::exp(-_j_CNN_vals[i] / cosx) * _jfac * roof_factor;
      return 0.0;
    }
  }
  return 0.0;
}

void
MCMRuntimeMechanism::setJCalibrator(std::unique_ptr<JCalibrator> calibrator)
{
  _jcalibrator = std::move(calibrator);
}
