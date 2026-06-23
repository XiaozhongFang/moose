//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMFacsimileParser.h"
#include "MooseUtils.h"
#include "pcrecpp.h"

#include <fstream>
#include <sstream>
#include <algorithm>

ParsedMechanism
MCMFacsimileParser::parse(const std::string & filename, const std::string & photolysis_file)
{
  _species.clear();
  _reactions.clear();
  _rate_coefficients.clear();
  _j_CL.clear();
  _j_CMM.clear();
  _j_CNN.clear();
  _photolysis_rates.clear();
  _coefficient_names.clear();
  _eval_order.clear();
  _converted_coefficients.clear();

  // 1. Parse the .fac file
  parseFile(filename);

  // 2. Auto-detect photolysis J<N> references
  scanPhotolysisReferences();

  // 3. Load photolysis parameters
  if (!photolysis_file.empty())
    loadPhotolysisParameters(photolysis_file);

  // 4. Build coefficient name set
  for (auto & [name, _] : _rate_coefficients)
    _coefficient_names.insert(name);

  // 5. Define base variables
  _base_variables = {"TEMP", "M", "O2", "N2", "H2O", "O", "N", "H", "OH", "HO2"};
  for (auto & jname : _photolysis_rates)
  {
    _base_variables.insert(jname);
    std::string alt = jname;
    std::replace(alt.begin(), alt.end(), '<', '(');
    std::replace(alt.begin(), alt.end(), '>', ')');
    _base_variables.insert(alt);
  }

  // 6. Build stoichiometric matrix
  buildStoichiometricMatrix();

  // 7. Build reactant indices (for F0AM-style dC/dt)
  buildReactantIndices();

  // 8. Expand rate coefficients
  expandRateCoefficients();

  // Build output
  ParsedMechanism mech;
  mech.species = _species;

  const unsigned int n_rxn = _reactions.size();
  const unsigned int n_sp = _species.size();
  for (unsigned int i = 0; i < n_rxn; ++i)
  {
    std::string rxn_name = "R" + std::to_string(i + 1);
    mech.reaction_names.push_back(rxn_name);
  }
  mech.reactions = _reactions;
  mech.stoichiometry.assign(n_sp, std::vector<Real>(n_rxn, 0.0));
  for (unsigned int r = 0; r < n_rxn; ++r)
    for (unsigned int s = 0; s < n_sp; ++s)
      mech.stoichiometry[s][r] = _stoichiometric[s][r];

  mech.reactant_indices.assign(n_rxn, std::vector<int>(2, -1));
  for (unsigned int r = 0; r < n_rxn; ++r)
    for (int k = 0; k < 2; ++k)
      mech.reactant_indices[r][k] = _reactant_indices[r][k];

  mech.coefficient_names = _eval_order;
  for (auto & name : _eval_order)
    mech.coefficient_expressions.push_back(
        _converted_coefficients.count(name) ? _converted_coefficients[name] : "");

  // Detect RO2 species (O2-suffix or RO2 pattern, AtChem2-compatible)
  for (auto & sp : _species)
    if ((sp.size() >= 3 && sp.substr(sp.size()-2) == "O2") ||
        sp.find("RO2") != std::string::npos)
      mech.ro2_species.push_back(sp);

  std::cout << "MCMFacsimileParser: Detected " << mech.ro2_species.size()
            << " RO2 species" << std::endl;

  // Transfer photolysis data
  for (auto & [jkey, cl_val] : _j_CL)
  {
    unsigned int jn;
    pcrecpp::RE("J<([0-9]+)>").FullMatch(jkey, &jn);
    mech.j_numbers.push_back(jn);
    mech.j_CL.push_back(cl_val);
    mech.j_CMM.push_back(_j_CMM[jkey]);
    mech.j_CNN.push_back(_j_CNN[jkey]);
  }

  return mech;
}

// ---------------------------------------------------------------------------
// Private parsing methods
// ---------------------------------------------------------------------------

void
MCMFacsimileParser::parseFile(const std::string & filename)
{
  std::ifstream file(filename);
  if (!file.good())
    mooseError("MCMFacsimileParser: Cannot open mechanism file: ", filename);

  std::string line, buffer;
  bool has_buffer = false;

  while (std::getline(file, line))
  {
    line.erase(line.find_last_not_of(" \t\r") + 1);
    if (line.empty())
      continue;

    std::string trimmed = line;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    if (trimmed.empty() || trimmed[0] == '*')
      continue;

    if (!has_buffer)
    {
      buffer = line;
      if (line.find(';') == std::string::npos)
        has_buffer = true;
      else
      {
        processStatement(buffer);
        buffer.clear();
      }
    }
    else
    {
      buffer += " " + line;
      if (line.find(';') != std::string::npos)
      {
        processStatement(buffer);
        buffer.clear();
        has_buffer = false;
      }
    }
  }

  file.close();
}

void
MCMFacsimileParser::processStatement(const std::string & statement)
{
  std::string s = statement;
  s.erase(0, s.find_first_not_of(" \t"));
  s.erase(s.find_last_not_of(" \t") + 1);
  if (s.empty())
    return;

  // Strip trailing inline * ... ; comment:
  //   "KMT01 = 2.0D-12 ; * text ;" → "KMT01 = 2.0D-12"
  auto star_pos = s.rfind(" * ");
  if (star_pos != std::string::npos)
  {
    auto trail_sc = s.find(';', star_pos);
    if (trail_sc != std::string::npos)
      s.erase(star_pos);
  }

  auto sc_pos = s.rfind(';');
  if (sc_pos != std::string::npos)
    s.erase(sc_pos);
  s.erase(s.find_last_not_of(" \t") + 1);
  if (s.empty())
    return;

  // Strip CONSTANT prefix: "CONSTANT TEMP 298.15" → "TEMP 298.15"
  if (s.find("CONSTANT ") == 0)
    s = s.substr(9);

  if (s.find("VARIABLE") == 0)
  {
    std::string var_part = s.substr(8);
    std::istringstream iss(var_part);
    std::string sp;
    while (iss >> sp)
      if (!sp.empty())
        _species.push_back(sp);
  }
  else if (s[0] == '%')
  {
    parseReactionLine(s);
  }
  else if (s.find('=') != std::string::npos)
  {
    // Skip RO2 sum declaration: "RO2 = CH3O2 + C2H5O2" is a section marker, not a rate coefficient
    auto eq_pos = s.find('=');
    std::string name = s.substr(0, eq_pos);
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);
    if (name == "RO2")
      return;

    std::string expr = s.substr(eq_pos + 1);
    expr.erase(0, expr.find_first_not_of(" \t"));
    expr.erase(expr.find_last_not_of(" \t") + 1);
    _rate_coefficients[name] = expr;
  }
}

void
MCMFacsimileParser::parseReactionLine(const std::string & line)
{
  std::string work = line;
  work.erase(0, 1);
  work.erase(0, work.find_first_not_of(" \t"));

  auto colon_pos = work.find(':');
  if (colon_pos == std::string::npos)
    return;

  std::string rate = work.substr(0, colon_pos);
  std::string equation = work.substr(colon_pos + 1);
  rate.erase(0, rate.find_first_not_of(" \t"));
  rate.erase(rate.find_last_not_of(" \t") + 1);
  equation.erase(0, equation.find_first_not_of(" \t"));
  equation.erase(equation.find_last_not_of(" \t") + 1);

  std::string rate_expr = convertRateExpression(rate);

  auto eq_pos = equation.find('=');
  if (eq_pos == std::string::npos)
    return;

  ParsedReaction rxn;
  rxn.rate_expression = rate_expr;
  extractSpecies(equation.substr(0, eq_pos), rxn.reactants);
  extractSpecies(equation.substr(eq_pos + 1), rxn.products);
  _reactions.push_back(rxn);
}

void
MCMFacsimileParser::extractSpecies(const std::string & side,
                                   std::vector<std::pair<Real, std::string>> & out)
{
  std::string trimmed = side;
  trimmed.erase(0, trimmed.find_first_not_of(" \t"));
  trimmed.erase(trimmed.find_last_not_of(" \t") + 1);
  if (trimmed.empty())
    return;

  std::istringstream stream(trimmed);
  std::string term;
  // Match coefficient + species: "2 NO2", "0.5O3", "A", "2.0D-3*C5H8"
  // Space between coefficient and name is allowed (AtChem2-compatible)
  // Handle "A", "2H2O", "2 H2O" — coefficient is optional
  pcrecpp::RE re_species("^\\s*(?:(\\d*\\.?\\d*(?:[dDeE][+-]?\\d+)?)\\s+)?([a-zA-Z_]\\w*)\\s*$");

  while (std::getline(stream, term, '+'))
  {
    term.erase(0, term.find_first_not_of(" \t"));
    term.erase(term.find_last_not_of(" \t") + 1);
    if (term.empty())
      continue;

    // Try space-separated format first (AtChem2: "2 H2O"), then compact
    std::string coeff_str, name;
    if (re_species.FullMatch(term, &coeff_str, &name))
    {
      Real coeff = 1.0;
      if (!coeff_str.empty())
      {
        std::string num = coeff_str;
        std::replace(num.begin(), num.end(), 'D', 'E');
        std::replace(num.begin(), num.end(), 'd', 'e');
        coeff = std::stod(num);
      }
      auto it = std::find_if(out.begin(), out.end(),
          [&name](const auto & p) { return p.second == name; });
      if (it != out.end())
        it->first += coeff;
      else
        out.emplace_back(coeff, name);
    }
  }
}

std::string
MCMFacsimileParser::convertRateExpression(const std::string & expr) const
{
  std::string s = expr;
  {
    pcrecpp::RE re_fortran_D("(\\d)D([+-]?)");
    re_fortran_D.GlobalReplace("\\1E\\2", &s);
    pcrecpp::RE re_fortran_d("(\\d)d([+-]?)");
    re_fortran_d.GlobalReplace("\\1e\\2", &s);
  }
  std::replace(s.begin(), s.end(), '@', '^');
  {
    auto pos = s.find("**");
    while (pos != std::string::npos) { s.replace(pos, 2, "^"); pos = s.find("**", pos + 1); }
  }
  {
    auto pos = s.find("EXP(");
    while (pos != std::string::npos) { s.replace(pos, 4, "exp("); pos = s.find("EXP(", pos + 4); }
  }
  {
    auto pos = s.find("LOG10(");
    while (pos != std::string::npos)
    {
      s.replace(pos, 6, "log10(");
      pos = s.find("LOG10(", pos + 6);
    }
  }
  return s;
}

void
MCMFacsimileParser::scanPhotolysisReferences()
{
  pcrecpp::RE re_j("J<([0-9]+)>");
  for (auto & [name, expr] : _rate_coefficients)
  {
    unsigned int jn;
    pcrecpp::StringPiece sp(expr);
    while (re_j.FindAndConsume(&sp, &jn))
      _photolysis_rates.insert("J<" + std::to_string(jn) + ">");
  }
  for (auto & rxn : _reactions)
  {
    unsigned int jn;
    pcrecpp::StringPiece sp(rxn.rate_expression);
    while (re_j.FindAndConsume(&sp, &jn))
      _photolysis_rates.insert("J<" + std::to_string(jn) + ">");
  }
}

void
MCMFacsimileParser::loadPhotolysisParameters(const std::string & filename)
{
  std::ifstream file(filename);
  if (!file.good())
  {
    mooseWarning("MCMFacsimileParser: Photolysis file not found: ", filename);
    return;
  }

  std::string line;
  while (std::getline(file, line))
  {
    auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos) continue;
    char fc = line[first];
    if (fc == '#' || fc == '!' || fc == 'j') continue;
    break;
  }

  auto parse_line = [&](const std::string & l) {
    std::string trimmed = l;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == '!') return;
    int jnum;
    Real cl_val, cmm_val, cnn_val;
    std::string name;
    int tau;
    {
      pcrecpp::RE re_f("(\\d)D([+-]?)");
      re_f.GlobalReplace("\\1E\\2", &trimmed);
    }
    std::istringstream iss(trimmed);
    if (iss >> jnum >> cl_val >> cmm_val >> cnn_val >> name >> tau)
    {
      std::string jkey = "J<" + std::to_string(jnum) + ">";
      if (_photolysis_rates.count(jkey))
      {
        _j_CL[jkey] = cl_val;
        _j_CMM[jkey] = cmm_val;
        _j_CNN[jkey] = cnn_val;
      }
    }
  };

  parse_line(line);
  while (std::getline(file, line))
    parse_line(line);
  file.close();
}

void
MCMFacsimileParser::buildStoichiometricMatrix()
{
  unsigned int n_sp = _species.size(), n_rx = _reactions.size();
  _stoichiometric.assign(n_sp, std::vector<Real>(n_rx, 0.0));

  for (unsigned int r = 0; r < n_rx; ++r)
  {
    for (auto & [coeff, name] : _reactions[r].reactants)
    {
      auto it = std::find(_species.begin(), _species.end(), name);
      if (it != _species.end())
        _stoichiometric[it - _species.begin()][r] -= coeff;
    }
    for (auto & [coeff, name] : _reactions[r].products)
    {
      auto it = std::find(_species.begin(), _species.end(), name);
      if (it != _species.end())
        _stoichiometric[it - _species.begin()][r] += coeff;
    }
  }
}

void
MCMFacsimileParser::buildReactantIndices()
{
  unsigned int n_rx = _reactions.size();
  _reactant_indices.assign(n_rx, std::vector<int>(2, -1));

  for (unsigned int r = 0; r < n_rx; ++r)
  {
    unsigned int k = 0;
    for (auto & [coeff, name] : _reactions[r].reactants)
    {
      auto it = std::find(_species.begin(), _species.end(), name);
      if (it == _species.end())
        continue;

      int idx = (int)(it - _species.begin());
      // Unroll merged coefficients (e.g. B+B → coeff=2 → push B twice)
      for (int c = 0; c < (int)coeff && k < 2; ++c)
      {
        _reactant_indices[r][k] = idx;
        ++k;
      }
    }
  }
}

void
MCMFacsimileParser::expandRateCoefficients()
{
  std::map<std::string, std::set<std::string>> deps;
  pcrecpp::RE re_token("[A-Z][A-Za-z0-9_]*");

  for (auto & [name, expr] : _rate_coefficients)
  {
    std::set<std::string> expr_deps;
    std::string token;
    pcrecpp::StringPiece input(expr);
    while (re_token.FindAndConsume(&input, &token))
      if (_coefficient_names.count(token) && !_base_variables.count(token) && token != name)
        expr_deps.insert(token);
    deps[name] = expr_deps;
  }

  std::map<std::string, unsigned> in_degree;
  for (auto & [name, _] : deps)
    in_degree[name] = deps[name].size();

  std::map<std::string, std::set<std::string>> dependents;
  for (auto & [name, expr_deps] : deps)
    for (auto & dep : expr_deps)
      dependents[dep].insert(name);

  std::queue<std::string> q;
  for (auto & [name, deg] : in_degree)
    if (deg == 0)
      q.push(name);

  _eval_order.clear();
  while (!q.empty())
  {
    auto name = q.front();
    q.pop();
    _eval_order.push_back(name);
    for (auto & dep : dependents[name])
      if (--in_degree[dep] == 0)
        q.push(dep);
  }

  _converted_coefficients = _rate_coefficients;
  for (auto & [name, expr] : _converted_coefficients)
    expr = convertRateExpression(expr);
}


