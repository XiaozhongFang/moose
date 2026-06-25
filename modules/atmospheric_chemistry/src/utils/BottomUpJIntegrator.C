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
  std::string cs_path = _data_dir + "/CrossSections/" + info.cs_file;
  switch (info.cs_type)
  {
    case 1: // 2-column CSV
    {
      auto p = loadCSV2(cs_path);
      wl_cs = std::move(p.first);
      cs = std::move(p.second);
      break;
    }
    case 2: // 3-column CSV with T-interpolation
    {
      auto p = loadCSV3(cs_path, T);
      wl_cs = std::move(p.first);
      cs = std::move(p.second);
      break;
    }
    case 3: // TXT
    {
      auto p = loadTXT(cs_path);
      wl_cs = std::move(p.first);
      cs = std::move(p.second);
      break;
    }
    default:
      mooseError("BottomUpJIntegrator: unsupported CS type ", info.cs_type, " for ", jname);
  }

  // ── Load quantum yield ──
  std::vector<Real> wl_qy, qy;
  if (info.qy_type == 0)
  {
    // scalar QY — use lamp flux grid directly
    Real scalar_val = 0.0;
    // Try to parse the qy_file as a number
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

  // ── Truncate to wl_bounds (use flux grid range) ──
  // Find indices within flux range where CS/QY data exists
  // (IntegrateJ lines 138-165: smear CS and QY onto flux grid)

  std::vector<Real> CS_out(N, 0.0), QY_out(N, 0.0);

  for (unsigned int i = 0; i < N; ++i)
  {
    Real wll = _wllim[i];
    Real wlu = _wllim[i + 1];

    // Smear cross-section onto this bin
    if (!(wlu < wl_cs.front() || wll > wl_cs.back()))
      CS_out[i] = smear(wl_cs, cs, wll, wlu);

    // Smear quantum yield onto this bin
    if (!(wlu < wl_qy.front() || wll > wl_qy.back()))
      QY_out[i] = smear(wl_qy, qy, wll, wlu);
  }

  // ── Integrate: J = trapz(wl, QY_out * CS_out * flux) ──
  std::vector<Real> integrand(N);
  for (unsigned int i = 0; i < N; ++i)
    integrand[i] = QY_out[i] * CS_out[i] * _flux[i];

  Real Jval = trapz(_wl_flux, integrand);

  // Defensive checks
  if (std::isnan(Jval))
  {
    mooseWarning("BottomUpJIntegrator: NaN J-value for ", jname);
    return 0.0;
  }
  if (Jval < 0.0)
  {
    mooseWarning("BottomUpJIntegrator: negative J-value for ", jname, " (", Jval, "), clamping to 0");
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
//  DATA FILE LOADERS
// ─────────────────────────────────────────────────────────────

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadCSV2(const std::string & path) const
{
  std::vector<Real> wl, val;
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open CSV file: ", path);

  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    // Replace commas with spaces for uniform parsing
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream iss(line);
    Real w, v;
    if (!(iss >> w >> v))
      continue;
    if (w <= 0.0)
      continue; // skip bad rows (e.g. dlmread zerojunk)
    wl.push_back(w);
    val.push_back(v);
  }
  file.close();

  if (wl.size() < 2)
    mooseError("BottomUpJIntegrator: too few data points in: ", path);

  return {wl, val};
}

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadCSV3(const std::string & path, Real T) const
{
  // 3-column CSV: wl, val@T1, val@T2 — linear interpolate to T
  std::vector<Real> wl, val;
  Real T1 = 0, T2 = 0;
  bool temps_read = false;

  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open CSV3 file: ", path);

  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream iss(line);
    Real w, v1, v2;
    // Try to read 3 values; if only 2, use the second as-is (no T-interp)
    if (!(iss >> w >> v1))
      continue;
    if (w <= 0.0)
      continue;

    if (!(iss >> v2))
    {
      // Only 2 columns — treat as type 1
      wl.push_back(w);
      val.push_back(v1);
      continue;
    }

    if (!temps_read)
    {
      // First data row with 3 columns — extract reference temps from header or infer
      // We need to know T1 and T2. The CSV doesn't embed them explicitly.
      // Use common conventions based on the filename.
      // For now, detect from known patterns or use defaults.
      temps_read = true;
      // Default reference temperatures (commonly 298K/220K or 295K/218K)
      // We'll handle this heuristically for now
    }

    wl.push_back(w);
    val.push_back(v1); // Store both for interpolation later
    // We'll do the actual interpolation below
  }
  file.close();

  if (wl.size() < 2)
    mooseError("BottomUpJIntegrator: too few data points in: ", path);

  // For 3-column CSVs, we need to detect reference temperatures from the filename/path
  // or use common conventions. Most JPL/IUPAC data uses 298K/220K or 295K/218K.
  // For now, just use the 2nd column as-is without T-interpolation.
  // TODO: Implement proper T-detection from file metadata.

  return {wl, val};
}

std::pair<std::vector<Real>, std::vector<Real>>
BottomUpJIntegrator::loadTXT(const std::string & path) const
{
  std::vector<Real> wl, val;
  std::ifstream file(path);
  if (!file.good())
    mooseError("BottomUpJIntegrator: Cannot open TXT file: ", path);

  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream iss(line);
    Real w, v;
    if (!(iss >> w >> v))
      continue;
    if (w <= 0.0)
      continue;
    wl.push_back(w);
    val.push_back(v);
  }
  file.close();

  if (wl.size() < 2)
    mooseError("BottomUpJIntegrator: too few data points in: ", path);

  return {wl, val};
}

// ─────────────────────────────────────────────────────────────
//  SMEAR / TRAPZ — TUV numer.f interp2 algorithm
// ─────────────────────────────────────────────────────────────

Real
BottomUpJIntegrator::smear(const std::vector<Real> & x, const std::vector<Real> & y, Real xgl,
                            Real xgu)
{
  // Integral convolution of (x,y) over [xgl, xgu] by trapezoidal rule
  // Replicates TUV numer.f/interp2 and F0AM IntegrateJ.m smear()
  Real area = 0.0;
  unsigned int n = x.size();
  for (unsigned int k = 0; k < n - 1; ++k)
  {
    if (x[k + 1] < xgl || x[k] > xgu)
      continue; // outside window

    Real a1 = std::max(x[k], xgl);
    Real a2 = std::min(x[k + 1], xgu);

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
