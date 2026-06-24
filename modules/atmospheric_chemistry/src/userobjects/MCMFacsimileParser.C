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
MCMFacsimileParser::parse(const std::string & filename, const std::string & photolysis_file,
                          const std::string & peroxy_file)
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
  _parsed_ro2_names.clear();
  if (_mcm_version.empty()) _mcm_version = "v3.3.1";

  // 1. Parse the mechanism file (.fac or .kpp)
  if (filename.size() >= 4 && filename.substr(filename.size() - 4) == ".kpp")
    parseKPP(filename);
  else
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

  mech.reactant_indices.assign(n_rxn, std::vector<int>(3, 0));
  for (unsigned int r = 0; r < n_rxn; ++r)
    for (int k = 0; k < 3; ++k)
      mech.reactant_indices[r][k] = _reactant_indices[r][k];

  mech.coefficient_names = _eval_order;
  for (auto & name : _eval_order)
    mech.coefficient_expressions.push_back(
        _converted_coefficients.count(name) ? _converted_coefficients[name] : "");

  // Load RO2 reference list for validation (if peroxy_file provided)
  std::set<std::string> ro2_reference;
  if (!peroxy_file.empty())
  {
    std::ifstream ref_file(peroxy_file);
    if (ref_file.good())
    {
      std::string ref_line;
      while (std::getline(ref_file, ref_line))
      {
        ref_line.erase(0, ref_line.find_first_not_of(" \t\r"));
        ref_line.erase(ref_line.find_last_not_of(" \t\r") + 1);
        if (!ref_line.empty())
          ro2_reference.insert(ref_line);
      }
      ref_file.close();
    }
  }

  // Build RO2 species list from explicit "Peroxy radicals" section.
  // Fall back to name-based heuristic if the section was absent.
  mech.ro2_species.clear();
  if (!_parsed_ro2_names.empty())
  {
    // Use explicit list from the "Peroxy radicals" section
    for (const auto & ro2_name : _parsed_ro2_names)
    {
      auto it = std::find(_species.begin(), _species.end(), ro2_name);
      if (it != _species.end())
        mech.ro2_species.push_back(ro2_name);
      else
        mooseWarning("MCMFacsimileParser: RO2 species '", ro2_name,
                     "' from Peroxy radicals section not found in species list");
    }
  }
  else
  {
    // Fallback: name-based heuristic (for .fac files without explicit section)
    static const std::set<std::string> non_ro2 = {
      "HO2", "NO2", "SO2", "H2O2", "O2", "N2O2",
      "NO3", "HNO3", "CO2", "CLO2", "CL2O2", "BRO2"
    };
    for (auto & sp : _species)
    {
      if (non_ro2.count(sp)) continue;
      if ((sp.size() >= 3 && sp.substr(sp.size()-2) == "O2") ||
          sp.find("RO2") != std::string::npos)
        mech.ro2_species.push_back(sp);
    }
  }

  // Validate RO2 species against reference list (AtChem2-compatible warning)
  if (!ro2_reference.empty())
  {
    unsigned int missing = 0, unknown = 0;
    for (const auto & ro2_name : mech.ro2_species)
      if (!ro2_reference.count(ro2_name))
      {
        if (missing < 10)
          mooseWarning("MCMFacsimileParser: RO2 species '", ro2_name,
                       "' not in MCM ", _mcm_version, " reference list");
        missing++;
      }
    for (const auto & ref_name : ro2_reference)
      if (std::find(mech.ro2_species.begin(), mech.ro2_species.end(), ref_name) ==
          mech.ro2_species.end())
        unknown++;
    if (missing > 10)
      mooseWarning("MCMFacsimileParser: ... and ", missing - 10,
                   " more RO2 species not in reference list");
    if (unknown > 0)
      mooseWarning("MCMFacsimileParser: ", unknown,
                   " species in MCM ", _mcm_version,
                   " reference list are NOT in the mechanism — reactions may be missing");
  }

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

  // --- Read all lines ---
  std::vector<std::string> raw_lines;
  std::string line;
  while (std::getline(file, line))
  {
    line.erase(line.find_last_not_of(" \t\r") + 1);
    if (!line.empty())
      raw_lines.push_back(line);
  }
  file.close();

  // --- Pre-processing: fix broken continuation lines (AtChem2 fix_mechanism_fac.py logic) ---
  // In the "Reaction definitions" section, lines not starting with '%' or '*'
  // are continuation fragments and should be appended to the previous line.
  bool in_rxn_section = false;
  std::vector<std::string> fixed_lines;
  for (const auto & l : raw_lines)
  {
    std::string trimmed = l;
    trimmed.erase(0, trimmed.find_first_not_of(" \t"));
    if (trimmed.empty()) continue;

    if (!in_rxn_section && trimmed.find("Reaction definitions") != std::string::npos)
      in_rxn_section = true;

    if (in_rxn_section && trimmed[0] != '%' && trimmed[0] != '*' && !fixed_lines.empty())
    {
      // Continuation line: append to previous
      fixed_lines.back() += " " + trimmed;
    }
    else
    {
      fixed_lines.push_back(l);
    }
  }

  // --- Split multi-statement lines (lines with multiple ';') into separate lines ---
  // Find end of header (before "Generic Rate Coefficients") — keep header intact
  size_t header_end = 0;
  for (size_t i = 0; i < fixed_lines.size(); ++i)
    if (fixed_lines[i].find("Generic Rate Coefficients") != std::string::npos)
    { header_end = i; break; }

  std::vector<std::string> split_lines;
  for (size_t i = 0; i < fixed_lines.size(); ++i)
  {
    if (i < header_end)
    {
      split_lines.push_back(fixed_lines[i]);
      continue;
    }
    // Count semicolons outside of the header
    std::string l = fixed_lines[i];
    size_t sc_count = std::count(l.begin(), l.end(), ';');
    if (sc_count <= 1)
    {
      split_lines.push_back(l);
    }
    else
    {
      // Split by ';', keeping each segment as a separate line
      std::istringstream sc_stream(l);
      std::string segment;
      while (std::getline(sc_stream, segment, ';'))
      {
        segment.erase(0, segment.find_first_not_of(" \t"));
        segment.erase(segment.find_last_not_of(" \t") + 1);
        if (!segment.empty())
        {
          // Restore the ';' as it's the statement terminator
          split_lines.push_back(segment + " ;");
        }
      }
    }
  }

  // --- Process lines ---
  std::string buffer;
  bool has_buffer = false;

  for (const auto & raw_line : split_lines)
  {
    std::string line = raw_line;
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
    {
      // Parse "Peroxy radicals" section: RO2 = A + B + C + ...
      std::string ro2_rhs = s.substr(eq_pos + 1);
      ro2_rhs.erase(0, ro2_rhs.find_first_not_of(" \t"));
      ro2_rhs.erase(ro2_rhs.find_last_not_of(" \t") + 1);
      std::istringstream ro2_ss(ro2_rhs);
      std::string ro2_term;
      while (std::getline(ro2_ss, ro2_term, '+'))
      {
        ro2_term.erase(0, ro2_term.find_first_not_of(" \t"));
        ro2_term.erase(ro2_term.find_last_not_of(" \t") + 1);
        if (!ro2_term.empty())
          _parsed_ro2_names.push_back(ro2_term);
      }
      return;
    }

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

  // Find ONE species index (FACSIMILE placeholder, must have net coeff = 0)
  int one_idx = -1;
  for (unsigned int s = 0; s < n_sp; ++s)
    if (_species[s] == "ONE")
    { one_idx = (int)s; break; }

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
    // ONE is a placeholder (conc=1) — must have net stoichiometry 0
    if (one_idx >= 0)
      _stoichiometric[one_idx][r] = 0.0;
  }
}

void
MCMFacsimileParser::buildReactantIndices()
{
  unsigned int n_rx = _reactions.size();
  // 3-column iG matching F0AM convention. Default to ONE index (0).
  _reactant_indices.assign(n_rx, std::vector<int>(3, 0));

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
      for (int c = 0; c < (int)coeff && k < 3; ++c)
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

// ---------------------------------------------------------------------------
// KPP format parser (AtChem2 kpp_conversion.py equivalent)
// ---------------------------------------------------------------------------
void
MCMFacsimileParser::parseKPP(const std::string & filename)
{
  std::ifstream file(filename);
  if (!file.good())
    mooseError("MCMFacsimileParser: Cannot open KPP file: ", filename);

  std::vector<std::string> lines;
  std::string line;
  while (std::getline(file, line))
  {
    line.erase(0, line.find_first_not_of(" \t"));
    line.erase(line.find_last_not_of(" \t\r") + 1);
    if (!line.empty())
      lines.push_back(line);
  }
  file.close();

  // --- Helper: extract section between markers ---
  auto section = [&](const std::string & start, const std::string & end)
      -> std::vector<std::string>
  {
    std::vector<std::string> result;
    bool in = false;
    for (const auto & l : lines)
    {
      if (!in && l.find(start) != std::string::npos) { in = true; continue; }
      if (in)
      {
        if (!end.empty() && l.find(end) != std::string::npos) break;
        result.push_back(l);
      }
    }
    return result;
  };

  // --- Section 1: Peroxy radicals (RO2 sum) ---
  // Format: RO2 = &
  //           C(ind_CH3O2) + C(ind_C2H5O2) + ...  &
  //           )
  auto ro2_lines = section("RO2 =", ")");
  for (auto & l : ro2_lines)
  {
    // Remove continuation markers
    l.erase(std::remove(l.begin(), l.end(), '&'), l.end());
    // Replace C(ind_NAME) with NAME
    pcrecpp::RE re_c("C\\(ind_([A-Za-z0-9_]+)\\)");
    re_c.GlobalReplace("\\1", &l);
    // Split by '+' and add to parsed RO2 list
    std::istringstream ro2_ss(l);
    std::string ro2_term;
    while (std::getline(ro2_ss, ro2_term, '+'))
    {
      ro2_term.erase(0, ro2_term.find_first_not_of(" \t"));
      ro2_term.erase(ro2_term.find_last_not_of(" \t") + 1);
      if (!ro2_term.empty())
        _parsed_ro2_names.push_back(ro2_term);
    }
  }

  // --- Section 2 & 3: Generic Rate Coefficients + Complex reactions ---
  auto rates_lines = section(")", "#ENDINLINE");
  // Remove `CALL mcm_constants()` line if present
  if (!rates_lines.empty())
  {
    auto & last = rates_lines.back();
    if (last.find("mcm_constants") != std::string::npos || last.find("CALL") != std::string::npos)
      rates_lines.pop_back();
  }

  // Known simple rate coefficients (MCM convention)
  static const std::set<std::string> simple_list = {
    "KRO2NO", "KRO2HO2", "KAPHO2", "KAPNO", "KRO2NO3",
    "KNO3AL", "KDEC", "KROPRIM", "KROSEC", "KCH3O2",
    "K298CH3O2", "K14ISOM1"
  };

  for (const auto & l : rates_lines)
  {
    auto eq_pos = l.find('=');
    if (eq_pos == std::string::npos) continue;
    std::string name = l.substr(0, eq_pos);
    std::string expr = l.substr(eq_pos + 1);
    name.erase(0, name.find_first_not_of(" \t"));
    name.erase(name.find_last_not_of(" \t") + 1);
    expr.erase(0, expr.find_first_not_of(" \t"));
    expr.erase(expr.find_last_not_of(" \t") + 1);
    // KPP uses ** for exponentiation; convert to @ (FACSIMILE convention)
    std::string::size_type pos = 0;
    while ((pos = expr.find("**", pos)) != std::string::npos)
    { expr.replace(pos, 2, "^"); pos++; }
    _rate_coefficients[name] = expr;
  }

  // --- Section 4: Reaction definitions ---
  auto rxn_lines = section("#EQUATIONS", "");
  pcrecpp::RE re_kpp_rxn("\\{\\s*(\\d+\\.?\\d*(?:[dDeE][+-]?\\d+)?)\\s*\\}\\s*");
  pcrecpp::RE re_kpp_j("J\\((\\d+)\\)");

  for (const auto & l : rxn_lines)
  {
    std::string rate_str, eqn_str;

    // Match KPP format: {rate} A + B = C + D : ... ;
    std::string work = l;
    // Remove leading {rate}
    std::string rate_num;
    if (re_kpp_rxn.PartialMatch(work, &rate_num))
    {
      // Reaction uses KPP rate number; we'll handle differently
    }

    // Split by ':' (rate : equation)
    auto colon_pos = work.find(':');
    if (colon_pos != std::string::npos)
    {
      rate_str = work.substr(0, colon_pos);
      eqn_str = work.substr(colon_pos + 1);
    }
    else
    {
      continue; // not a reaction line
    }

    // Clean up
    rate_str.erase(0, rate_str.find_first_not_of(" \t{}"));
    rate_str.erase(rate_str.find_last_not_of(" \t") + 1);
    eqn_str.erase(0, eqn_str.find_first_not_of(" \t"));
    auto sc_pos = eqn_str.rfind(';');
    if (sc_pos != std::string::npos) eqn_str.erase(sc_pos);
    eqn_str.erase(eqn_str.find_last_not_of(" \t") + 1);

    // Convert J(N) to J<N>
    re_kpp_j.GlobalReplace("J<\\1>", &rate_str);

    // Convert ** to ^
    std::string::size_type pos2 = 0;
    while ((pos2 = rate_str.find("**", pos2)) != std::string::npos)
    { rate_str.replace(pos2, 2, "^"); pos2++; }

    // Parse equation
    auto eq_pos2 = eqn_str.find('=');
    if (eq_pos2 == std::string::npos) continue;

    std::string reactants_str = eqn_str.substr(0, eq_pos2);
    std::string products_str = eqn_str.substr(eq_pos2 + 1);

    ParsedReaction rxn;
    rxn.rate_expression = convertRateExpression(rate_str);
    extractSpecies(reactants_str, rxn.reactants);
    extractSpecies(products_str, rxn.products);
    _reactions.push_back(rxn);

    // Collect species
    for (auto & [_, name] : rxn.reactants)
      if (std::find(_species.begin(), _species.end(), name) == _species.end())
        _species.push_back(name);
    for (auto & [_, name] : rxn.products)
      if (std::find(_species.begin(), _species.end(), name) == _species.end())
        _species.push_back(name);
  }
}


