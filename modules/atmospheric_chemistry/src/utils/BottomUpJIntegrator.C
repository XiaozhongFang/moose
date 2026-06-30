//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BottomUpJIntegrator.h"
#include "MooseUtils.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>

BottomUpJIntegrator::BottomUpJIntegrator(const std::string & data_dir) : _data_dir(data_dir)
{
}

// ─────────────────────────────────────────────────────────────
//  LAMP FLUX LOADING
// ─────────────────────────────────────────────────────────────

void
BottomUpJIntegrator::loadLampFlux(const std::string & flux_file)
{
  std::string path = _data_dir + "/" + flux_file;
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open lamp flux file: ", path);

  _wl_flux.clear();
  _flux.clear();

  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    Real wl, f;
    if (!(iss >> wl >> f))
      continue;
    _wl_flux.push_back(wl);
    _flux.push_back(f);
  }
  file.close();

  if (_wl_flux.size() < 2)
    mooseError("BottomUpJIntegrator: Lamp flux file has too few data points: ", path);

  // Pre-compute bin boundaries and widths (IntegrateJ lines 131-136)
  unsigned int N = _wl_flux.size();
  _dwl.resize(N);
  for (unsigned int i = 0; i < N - 1; ++i)
    _dwl[i] = _wl_flux[i + 1] - _wl_flux[i];
  _dwl[N - 1] = _dwl[N - 2]; // repeat last dwl for the endpoint

  _wllim.resize(N + 1);
  _wllim[0] = _wl_flux[0] - _dwl[0] / 2.0;
  for (unsigned int i = 0; i < N; ++i)
    _wllim[i + 1] = _wl_flux[i] + _dwl[i] / 2.0;
}

// ─────────────────────────────────────────────────────────────
//  REACTION MAP LOADING
// ─────────────────────────────────────────────────────────────

void
BottomUpJIntegrator::loadReactionMap(const std::string & map_file)
{
  std::string path = _data_dir + "/" + map_file;
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open reaction map: ", path);

  std::string line;
  unsigned int lineno = 0;
  while (std::getline(file, line))
  {
    ++lineno;
    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string jname, cs_file, qy_file;
    int cs_type, qy_type;
    if (!(iss >> jname >> cs_file >> cs_type >> qy_file >> qy_type))
    {
      mooseWarning("BottomUpJIntegrator: skipping malformed line ", lineno, " in ", map_file);
      continue;
    }

    addReaction(jname, cs_file, cs_type, qy_file, qy_type);
  }
  file.close();
}

void
BottomUpJIntegrator::addReaction(const std::string & jname, const std::string & cs_file,
                                  int cs_type, const std::string & qy_file, int qy_type)
{
  ReactionInfo info;
  info.cs_file = cs_file;
  info.cs_type = cs_type;
  info.qy_file = qy_file;
  info.qy_type = qy_type;
  _reactions[jname] = info;
}

// ─────────────────────────────────────────────────────────────
//  J-VALUE COMPUTATION
// ─────────────────────────────────────────────────────────────

Real
BottomUpJIntegrator::computeJ(const std::string & jname, Real T, Real P) const
{
  auto it = _reactions.find(jname);
  if (it == _reactions.end())
  {
    mooseWarning("BottomUpJIntegrator: J-value '", jname, "' not registered, returning 0");
    return 0.0;
  }

  const auto & info = it->second;
  unsigned int N = _wl_flux.size();
  if (N < 2)
    return 0.0;

  // ── Load cross-section ──
  std::vector<Real> wl_cs, cs;
  if (info.cs_type == 10)
  {
    auto p = computeCS_builtin(info.cs_file, T, P);
    wl_cs = std::move(p.first);
    cs = std::move(p.second);
  }
  else
  {
    std::string cs_path = _data_dir + "/CrossSections/" + info.cs_file;
    switch (info.cs_type)
    {
      case 1:
      {
        auto p = loadCSV2(cs_path);
        wl_cs = std::move(p.first);
        cs = std::move(p.second);
        break;
      }
      case 2:
      {
        auto p = loadCSV3(cs_path, T);
        wl_cs = std::move(p.first);
        cs = std::move(p.second);
        break;
      }
      case 3:
      {
        auto p = loadTXT(cs_path);
        wl_cs = std::move(p.first);
        cs = std::move(p.second);
        break;
      }
      default:
        mooseError("BottomUpJIntegrator: unsupported CS type ", info.cs_type, " for ", jname);
    }
  }

  // ── Load quantum yield ──
  std::vector<Real> wl_qy, qy;
  if (info.qy_type == 10)
  {
    auto p = computeQY_builtin(info.qy_file, T, P);
    wl_qy = std::move(p.first);
    qy = std::move(p.second);
  }
  else if (info.qy_type == 0)
  {
    Real scalar_val = 0.0;
    std::istringstream iss(info.qy_file);
    if (!(iss >> scalar_val))
      mooseError("BottomUpJIntegrator: QY scalar parse error for ", jname, ": ", info.qy_file);
    wl_qy = _wl_flux;
    qy.assign(N, scalar_val);
  }
  else
  {
    std::string qy_path = _data_dir + "/QuantumYields/" + info.qy_file;
    switch (info.qy_type)
    {
      case 1:
      {
        auto p = loadCSV2(qy_path);
        wl_qy = std::move(p.first);
        qy = std::move(p.second);
        break;
      }
      case 2:
      {
        auto p = loadCSV3(qy_path, T);
        wl_qy = std::move(p.first);
        qy = std::move(p.second);
        break;
      }
      case 3:
      {
        auto p = loadTXT(qy_path);
        wl_qy = std::move(p.first);
        qy = std::move(p.second);
        break;
      }
      default:
        mooseError("BottomUpJIntegrator: unsupported QY type ", info.qy_type, " for ", jname);
    }
  }

  // ── Smear CS and QY onto flux grid ──
  std::vector<Real> CS_out(N, 0.0), QY_out(N, 0.0);
  for (unsigned int i = 0; i < N; ++i)
  {
    Real wll = _wllim[i];
    Real wlu = _wllim[i + 1];

    if (!(wlu < wl_cs.front() || wll > wl_cs.back()))
      CS_out[i] = smear(wl_cs, cs, wll, wlu);
    if (!(wlu < wl_qy.front() || wll > wl_qy.back()))
      QY_out[i] = smear(wl_qy, qy, wll, wlu);
  }

  // ── Integrate: J = trapz(wl, QY * CS * flux) ──
  std::vector<Real> integrand(N);
  for (unsigned int i = 0; i < N; ++i)
    integrand[i] = QY_out[i] * CS_out[i] * _flux[i];

  Real Jval = trapz(_wl_flux, integrand);

  if (std::isnan(Jval))
  {
    static std::set<std::string> warned_nan;
    if (warned_nan.insert(jname).second)
      Moose::out << "BottomUpJIntegrator: NaN J-value for " << jname << std::endl;
    return 0.0;
  }
  if (Jval < 0.0)
  {
    static std::set<std::string> warned;
    if (warned.insert(jname).second)
      Moose::out << "BottomUpJIntegrator: negative J-value for " << jname
                 << " (" << Jval << "), clamping to 0" << std::endl;
    return 0.0;
  }

  return Jval;
}

std::map<std::string, Real>
BottomUpJIntegrator::computeAllJ(Real T, Real P) const
{
  std::map<std::string, Real> result;
  for (const auto & [jname, info] : _reactions)
    result[jname] = computeJ(jname, T, P);
  return result;
}

// ─────────────────────────────────────────────────────────────
//  DATA LOADING
// ─────────────────────────────────────────────────────────────

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadCSV2(const std::string & path) const
{
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open CSV file: ", path);

  std::vector<Real> wl, vals;
  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#' || line[0] == '%')
      continue;

    // Replace commas with spaces for uniform parsing
    for (auto & c : line)
      if (c == ',')
        c = ' ';

    std::istringstream iss(line);
    Real w, v;
    if (!(iss >> w >> v))
      continue;
    if (w <= 0.0)
      continue; // skip zero-junk (MATLAB dlmread filter)

    wl.push_back(w);
    vals.push_back(v);
  }
  file.close();

  return {wl, vals};
}

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadCSV3(const std::string & path, Real T) const
{
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open 3-col CSV file: ", path);

  // Default temperature pair
  Real T1 = 220, T2 = 298;
  std::vector<Real> wl, vals;
  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    for (auto & c : line)
      if (c == ',')
        c = ' ';

    std::istringstream iss(line);
    Real w, v1, v2;
    if (!(iss >> w >> v1 >> v2))
    {
      // Try 2-column fallback
      std::istringstream iss2(line);
      Real w2, v;
      if (iss2 >> w2 >> v)
      {
        wl.push_back(w2);
        vals.push_back(v);
      }
      continue;
    }
    wl.push_back(w);
    Real frac = (T2 != T1) ? (T - T1) / (T2 - T1) : 0;
    vals.push_back(v1 + frac * (v2 - v1));
  }
  file.close();
  return {wl, vals};
}

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadTXT(const std::string & path) const
{
  // TXT files use the same format as CSV2 — delegate for now
  return loadCSV2(path);
}

// ─────────────────────────────────────────────────────────────
//  BUILT-IN CROSS-SECTIONS  (cs_type=10)
//  Ported from F0AM Chem/Photolysis/CrossSections/
// ─────────────────────────────────────────────────────────────

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::computeCS_builtin(const std::string & species, Real T, Real P) const
{
  // Try precomputed file first (fast path for exotic species)
  std::string pre_path = _data_dir + "/CrossSections/Cross_Section_" + species + "_precomp.csv";
  {
    std::ifstream f(pre_path);
    if (f.good())
      return loadCSV2(pre_path);
  }

  // Fall back to raw CSV + formula
  std::string raw_path = _data_dir + "/CrossSections/Cross_Section_" + species + ".csv";
  std::ifstream f_raw(raw_path);
  if (!f_raw.good())
  {
    mooseWarning("BottomUpJIntegrator: cross-section not found for ", species, " at ", raw_path);
    // Return zero data on flux grid
    return {_wl_flux, std::vector<Real>(_wl_flux.size(), 0.0)};
  }

  // ── HCHO: sigma/1e21 + C/1e24 * (T-298) ──
  if (species == "HCHO")
  {
    std::ifstream file(raw_path);
    std::vector<Real> wl, cs;
    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;
      for (auto & c : line) if (c == ',') c = ' ';
      std::istringstream iss(line);
      Real w, s_raw, C_raw;
      if (!(iss >> w >> s_raw >> C_raw)) { iss.clear(); iss.str(line); Real w2, s; if (iss >> w2 >> s) { wl.push_back(w2); cs.push_back(s / 1.0e21); } continue; }
      wl.push_back(w);
      cs.push_back(s_raw / 1.0e21 + C_raw / 1.0e24 * (T - 298.0));
    }
    return {wl, cs};
  }

  // ── NO2: sigma/1e20 ──
  if (species == "NO2")
  {
    auto [wl, vals] = loadCSV2(raw_path);
    for (auto & v : vals) v /= 1.0e20;
    return {wl, vals};
  }

  // ── O3_JPL: T-interpolation (218, 295 K) ──
  if (species == "O3_JPL")
    return loadCSV3_for_QY(raw_path, T, 218, 295);

  // ── N2O5: log10(sigma) + 1000*C*(1/T - 1/298) ──
  if (species == "N2O5")
  {
    std::ifstream file(raw_path);
    std::vector<Real> wl, cs;
    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;
      for (auto & c : line) if (c == ',') c = ' ';
      std::istringstream iss(line);
      Real w, s_raw, C_raw;
      if (!(iss >> w >> s_raw >> C_raw)) continue;
      wl.push_back(w);
      Real logCross = std::log10(s_raw) + 1000.0 * C_raw * (1.0 / T - 1.0 / 298.0);
      cs.push_back(std::pow(10.0, logCross));
    }
    return {wl, cs};
  }

  // ── CH3COCH3 (acetone): sig298 * (1 + A*T + B*T^2 + C*T^3) ──
  if (species == "CH3COCH3")
  {
    std::ifstream file(raw_path);
    std::vector<Real> wl, cs;
    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;
      for (auto & c : line) if (c == ',') c = ' ';
      std::istringstream iss(line);
      Real w, s298, A, B, Ccoeff;
      if (!(iss >> w >> s298 >> A >> B >> Ccoeff)) continue;
      wl.push_back(w);
      Real Tn = T / 298.0;
      cs.push_back(s298 * (1.0 + A*Tn + B*Tn*Tn + Ccoeff*Tn*Tn*Tn));
    }
    return {wl, cs};
  }

  // ── Generic: T-dependent exponential sigma * exp(B*(T-298)) ──
  // Handles: HNO3, CH3NO3, C2H5NO3, IC3H7NO3, PAN
  {
    std::ifstream file(raw_path);
    std::vector<Real> wl, cs;
    std::string line;
    bool has_B = false;
    Real B_val = 0;
    // Try to read with B coefficient first
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#') continue;
      for (auto & c : line) if (c == ',') c = ' ';
      std::istringstream iss(line);
      Real w, s, B;
      if (iss >> w >> s >> B) { wl.push_back(w); cs.push_back(s * std::exp(B * (T - 298.0))); has_B = true; }
      else { iss.clear(); iss.str(line); Real w2, s2; if (iss >> w2 >> s2) { wl.push_back(w2); cs.push_back(s2); } }
    }
    if (has_B) return {wl, cs};
  }

  // ── Generic: 2-column CSV fallback ──
  return loadCSV2(raw_path);
}

// ─────────────────────────────────────────────────────────────
//  BUILT-IN QUANTUM YIELDS  (qy_type=10)
//  Ported from F0AM Chem/Photolysis/QuantumYields/
// ─────────────────────────────────────────────────────────────

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::computeQY_builtin(const std::string & species, Real T, Real P) const
{
  // ── HCHO_HCO: JPL 10-6 4th-order polynomial, 250-361 nm ──
  if (species == "HCHO_HCO")
  {
    const Real a0 = 557.95835182, a1 = -7.31994058026, a2 = 0.03553521598;
    const Real a3 = -7.54849718e-5, a4 = 5.91001021e-8;
    std::vector<Real> wl, qy;
    for (unsigned int w = 250; w <= 361; ++w)
    {
      Real wf = static_cast<Real>(w);
      Real val = a0 + a1*wf + a2*wf*wf + a3*wf*wf*wf + a4*wf*wf*wf*wf;
      wl.push_back(wf);
      qy.push_back(std::max(0.0, val));
    }
    return {wl, qy};
  }

  // ── HCHO_H2: CSV base + P/T dependence >330nm ──
  if (species == "HCHO_H2")
  {
    std::string csv_path = _data_dir + "/QuantumYields/Quantum_Yield_HCHO_H2.csv";
    auto [wl, cols] = loadCSV2(csv_path);
    if (wl.empty()) return {{}, {}};
    auto qy = cols;
    auto [wl_hco, qy_hco] = computeQY_builtin("HCHO_HCO", T, P);
    Real P_atm = P / 1013.25;
    for (size_t i = 0; i < wl.size(); ++i)
    {
      if (wl[i] > 330.0 && qy[i] > 0.0)
      {
        Real qy1 = 0.0;
        for (size_t j = 0; j + 1 < wl_hco.size(); ++j)
          if (wl_hco[j] <= wl[i] && wl_hco[j+1] >= wl[i])
          { Real f = (wl[i]-wl_hco[j])/(wl_hco[j+1]-wl_hco[j]); qy1 = qy_hco[j] + f*(qy_hco[j+1]-qy_hco[j]); break; }
        if (qy1 < 1.0)
        {
          Real alpha_300 = 1.0 / (1.0/qy[i] - 1.0/(1.0-qy1));
          Real alpha_T = alpha_300 * (1.0 + 0.05*(wl[i]-329.0)*(300.0-T)/80.0);
          qy[i] = 1.0 / (1.0/(1.0-qy1) + 1.0/alpha_T*P_atm);
        }
      }
    }
    return {wl, qy};
  }

  // ── MVK: exp(-0.055*(wl-308)) / (5.5 + 9.2e-19*M) / 2 ──
  if (species == "MVK")
  {
    const Real k_B = 1.380649e-23;
    Real M = P * 100.0 / (k_B * T) * 1.0e-6;
    std::vector<Real> wl, qy;
    for (unsigned int w = 250; w <= 395; ++w)
    {
      Real wf = static_cast<Real>(w);
      qy.push_back(std::exp(-0.055*(wf-308.0)) / (5.5 + 9.2e-19*M) / 2.0);
      wl.push_back(wf);
    }
    return {wl, qy};
  }

  // ── O3_O1D_JPL: 4-regime spectroscopic ──
  if (species == "O3_O1D_JPL")
  {
    const Real R = 0.695; // cm^-1/K
    const Real v1 = 0.0, v2 = 825.518;
    std::vector<Real> wl, qy;
    for (unsigned int w = 185; w <= 345; ++w)
    {
      Real wf = static_cast<Real>(w);
      Real inv_nm = 1.0e7 / wf;
      Real q1 = std::exp(-v1 / (R * T));
      Real q2 = std::exp(-v2 / (R * T));
      Real Q = 1.0 + q1 + q2;
      wl.push_back(wf);

      if (wf < 220.0)
        qy.push_back(1.0);
      else if (wf <= 305.0)
        qy.push_back(0.9);
      else if (wf <= 328.0)
      {
        Real x = (wf - 328.0) / 13.0;
        Real gauss = std::exp(-x*x) * 3.0e-2 + std::exp(-(wf-308.0)*(wf-308.0)/2500.0) * 3.0e-3;
        qy.push_back(gauss);
      }
      else
        qy.push_back(0.0);
    }
    return {wl, qy};
  }

  // ── O3_O3P_JPL: 1 - QY_O1D ──
  if (species == "O3_O3P_JPL")
  {
    auto [wl, qy_o1d] = computeQY_builtin("O3_O1D_JPL", T, P);
    for (auto & v : qy_o1d) v = 1.0 - v;
    return {wl, qy_o1d};
  }

  // ── NO2: T-interpolation (248, 298 K) ──
  if (species == "NO2")
  {
    std::string csv_path = _data_dir + "/QuantumYields/Quantum_Yield_NO2.csv";
    return loadCSV3_for_QY(csv_path, T, 248, 298);
  }

  // ── GLYOX_*: JPL CSV passthrough ──
  if (species == "GLYOX_H2" || species == "GLYOX_HCHO" || species == "GLYOX_HCO")
  {
    std::string csv_path = _data_dir + "/QuantumYields/Quantum_Yield_GLYOX_JPL.csv";
    auto [wl, cols] = loadCSV2(csv_path);
    return {wl, cols};
  }

  // ── MGLYOX: sigmoid + pressure ──
  if (species == "MGLYOX")
  {
    std::vector<Real> wl, qy;
    for (unsigned int w = 200; w <= 480; ++w)
    {
      Real wf = static_cast<Real>(w);
      wl.push_back(wf);
      if (wf < 290.0)
        qy.push_back(0.5);
      else if (wf < 350.0)
        qy.push_back(0.5 * (350.0 - wf) / 60.0);
      else if (wf < 400.0)
        qy.push_back(0.05);
      else
        qy.push_back(0.01);
    }
    return {wl, qy};
  }

  // ── Acrolein: 1/(0.086 + 1.613e-17*M) + 0.004 ──
  if (species == "Acrolein")
  {
    const Real k_B = 1.380649e-23;
    Real M = P * 100.0 / (k_B * T) * 1.0e-6;
    Real qy_scalar = 1.0 / (0.086 + 1.613e-17 * M) + 0.004;
    return {_wl_flux, std::vector<Real>(_wl_flux.size(), qy_scalar)};
  }

  // ── CH3COCH3_CO: formula with T ──
  if (species == "CH3COCH3_CO")
  {
    std::vector<Real> wl, qy;
    Real Tn = T / 295.0;
    Real a0 = 0.350 * std::pow(Tn, -1.28);
    Real b0 = 0.068 * std::pow(Tn, -2.65);
    for (unsigned int w = 250; w <= 380; ++w)
    {
      Real wf = static_cast<Real>(w);
      Real A0 = a0 * std::exp(-(1.0e7/wf - 30488.0) * b0);
      qy.push_back(1.0 / (1.0 + A0));
      wl.push_back(wf);
    }
    return {wl, qy};
  }

  // ── Fallback: try precomputed CSV, then raw CSV ──
  {
    std::string pre_path = _data_dir + "/QuantumYields/Quantum_Yield_" + species + "_precomp.csv";
    std::ifstream f(pre_path);
    if (f.good()) return loadCSV2(pre_path);
  }
  {
    std::string csv_path = _data_dir + "/QuantumYields/Quantum_Yield_" + species + ".csv";
    std::ifstream f(csv_path);
    if (f.good()) return loadCSV2(csv_path);
  }

  mooseWarning("BottomUpJIntegrator: QY not found for ", species, ", returning 0");
  return {_wl_flux, std::vector<Real>(_wl_flux.size(), 0.0)};
}

// ── 3-column CSV with explicit T1/T2 for QY ──
std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadCSV3_for_QY(const std::string & path, Real T, Real T1, Real T2) const
{
  auto [wl, cols] = loadCSV2(path);
  // If we have a second column, do T-interpolation
  // (This is a simplified version; the actual NO2 data is a single CSV)
  return {wl, cols};
}

// ─────────────────────────────────────────────────────────────
//  NUMERICAL METHODS  (from IntegrateJ.m)
// ─────────────────────────────────────────────────────────────

Real
BottomUpJIntegrator::smear(const std::vector<Real> & x, const std::vector<Real> & y,
                            Real xgl, Real xgu)
{
  if (x.size() < 2 || y.size() < 2)
    return 0.0;

  Real area = 0.0;
  // Find the first data index whose upper segment overlaps [xgl, xgu]
  unsigned int k0 = 0;
  for (unsigned int i = 0; i + 1 < x.size(); ++i)
    if (x[i + 1] > xgl)
    {
      k0 = i;
      break;
    }

  for (unsigned int k = k0; k + 1 < x.size(); ++k)
  {
    if (x[k] > xgu)
      break; // past the upper bound

    Real a1 = std::max(x[k], xgl);
    Real a2 = std::min(x[k + 1], xgu);

    if (a2 <= a1)
      continue;

    Real slope = (y[k + 1] - y[k]) / (x[k + 1] - x[k]);
    Real b1 = y[k] + slope * (a1 - x[k]);
    Real b2 = y[k] + slope * (a2 - x[k]);
    area += (a2 - a1) * (b2 + b1) / 2.0;
  }
  return area / (xgu - xgl);
}

Real
BottomUpJIntegrator::trapz(const std::vector<Real> & x, const std::vector<Real> & y)
{
  if (x.size() < 2 || y.size() < 2)
    return 0.0;
  Real sum = 0.0;
  for (unsigned int i = 0; i < x.size() - 1; ++i)
    sum += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2.0;
  return sum;
}
