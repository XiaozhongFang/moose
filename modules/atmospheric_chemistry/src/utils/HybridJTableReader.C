//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HybridJTableReader.h"
#include "MooseUtils.h"
#include <fstream>
#include <algorithm>
#include <cmath>

HybridJTableReader::HybridJTableReader(const std::string & dir)
{
  // Read index
  std::string idx_path = dir + "/index.txt";
  std::ifstream idx(idx_path);
  if (!idx.good())
    mooseError("HybridJTableReader: Cannot open index file: ", idx_path);

  idx >> _ns;  // consume "ns="
  char c; idx >> c >> c;  // skip "="
  idx >> _ns >> _na >> _no >> _nh;

  // Load axes
  loadAxis(_sza_axis, dir + "/axis_sza.dat");
  loadAxis(_alb_axis, dir + "/axis_albedo.dat");
  loadAxis(_o3_axis, dir + "/axis_o3col.dat");
  loadAxis(_alt_axis, dir + "/axis_alt.dat");

  // Load J-value tables
  std::string line;
  while (std::getline(idx, line))
  {
    if (line.empty() || line[0] == '#') continue;
    std::string path = dir + "/" + line;
    std::string fname = line.substr(line.find('/')+1);
    // Extract J name: "table_J4.dat" → "J4"
    auto dot = fname.find(".dat");
    std::string jname = fname.substr(6, dot - 6);

    std::ifstream file(path, std::ios::binary);
    if (!file.good())
      mooseError("HybridJTableReader: Cannot open table file: ", path);

    std::vector<Real> & tbl = _tables[jname];
    unsigned int n = _ns * _na * _no * _nh;
    tbl.resize(n);
    file.read(reinterpret_cast<char *>(tbl.data()), n * sizeof(Real));
    file.close();
  }
}

void
HybridJTableReader::loadAxis(std::vector<Real> & axis, const std::string & filepath)
{
  std::ifstream file(filepath, std::ios::binary);
  if (!file.good())
    mooseError("HybridJTableReader: Cannot open axis file: ", filepath);
  file.seekg(0, std::ios::end);
  unsigned int n = file.tellg() / sizeof(Real);
  axis.resize(n);
  file.seekg(0);
  file.read(reinterpret_cast<char *>(axis.data()), n * sizeof(Real));
}

Real
HybridJTableReader::interpolate(const std::string & jname,
                                Real sza, Real albedo, Real o3col, Real altitude) const
{
  if (sza >= 90.0) return 0.0;
  auto it = _tables.find(jname);
  if (it == _tables.end())
  {
    mooseWarning("HybridJTableReader: J-value '", jname, "' not found, returning 0");
    return 0.0;
  }
  return interp4d(it->second.data(), sza, albedo, o3col, altitude);
}

Real
HybridJTableReader::interp4d(const Real * data,
                              Real sza, Real alb, Real o3, Real alt) const
{
  // Find indices for each axis (clamped to limits)
  struct Idx3 { int i0, i1; Real t; };
  auto find_idx = [](const std::vector<Real> & ax, Real v) -> Idx3 {
    if (v <= ax.front()) return {0, 0, 0.0};
    if (v >= ax.back())  return {(int)ax.size()-2, (int)ax.size()-2, 1.0};
    auto it = std::lower_bound(ax.begin(), ax.end(), v);
    int i = it - ax.begin() - 1;
    return {i, i+1, (v - ax[i])/(ax[i+1] - ax[i])};
  };

  Idx3 si = find_idx(_sza_axis, sza);
  Idx3 ai = find_idx(_alb_axis, alb);
  Idx3 oi = find_idx(_o3_axis, o3);
  Idx3 hi = find_idx(_alt_axis, alt);

  // 4D linear interpolation on MATLAB col-major grid
  auto idx = [&](int s, int a, int o, int h) {
    return s + a * (int)_ns + o * (int)_ns * (int)_na + h * (int)_ns * (int)_na * (int)_no;
  };

  // Fetch all 16 corner values of the 4D hypercube
  // A,B,C,D = corner indices: 0=lo, 1=hi
  auto c = [&](int A, int B, int C, int D) {
    return data[idx(A ? si.i1 : si.i0,  B ? ai.i1 : ai.i0,
                    C ? oi.i1 : oi.i0,  D ? hi.i1 : hi.i0)];
  };

  // 4D linear interpolation = successive 1D lininterp across each dim
  // Step 1: interp altitude (3D → 2D): 8 corners → 4
  Real d00 = c(0,0,0,0) + hi.t * (c(0,0,0,1) - c(0,0,0,0));
  Real d01 = c(0,0,1,0) + hi.t * (c(0,0,1,1) - c(0,0,1,0));
  Real d10 = c(0,1,0,0) + hi.t * (c(0,1,0,1) - c(0,1,0,0));
  Real d11 = c(0,1,1,0) + hi.t * (c(0,1,1,1) - c(0,1,1,0));
  Real e00 = c(1,0,0,0) + hi.t * (c(1,0,0,1) - c(1,0,0,0));
  Real e01 = c(1,0,1,0) + hi.t * (c(1,0,1,1) - c(1,0,1,0));
  Real e10 = c(1,1,0,0) + hi.t * (c(1,1,0,1) - c(1,1,0,0));
  Real e11 = c(1,1,1,0) + hi.t * (c(1,1,1,1) - c(1,1,1,0));

  // Step 2: interp O3 (4 → 2): 4 corners → 2
  Real f00 = d00 + oi.t * (d01 - d00);
  Real f10 = d10 + oi.t * (d11 - d10);
  Real g00 = e00 + oi.t * (e01 - e00);
  Real g10 = e10 + oi.t * (e11 - e10);

  // Step 3: interp albedo (2 → 1): 2 corners → 1
  Real h0 = f00 + ai.t * (f10 - f00);
  Real h1 = g00 + ai.t * (g10 - g00);

  // Step 4: interp SZA (final)
  return h0 + si.t * (h1 - h0);
}
