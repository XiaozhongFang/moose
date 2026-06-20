//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
#include "Moose.h"
#include <vector>
#include <string>
#include <map>

/**
 * Reader + 4D interpolator for F0AM Hybrid J-value lookup tables.
 *
 * Grid: (SZA, albedo, O3column, altitude) — 27×6×11×23
 * Interpolation: 4D linear on log10(J) space.
 *
 * Data format: one binary file per J-value (table_J1.dat ... table_Jn36.dat),
 * axis values in axis_sza.dat / axis_albedo.dat / axis_o3col.dat / axis_alt.dat,
 * metadata in index.txt.
 *
 * Usage:
 *   HybridJTableReader reader(dir);
 *   Real Jval = reader.interpolate("J4", sza_deg, albedo, o3col_du, alt_m);
 */
class HybridJTableReader
{
public:
  HybridJTableReader(const std::string & data_dir);

  /** 4D linear interpolation of a single J-value. Returns 0.0 for SZA >= 90. */
  Real interpolate(const std::string & jname,
                   Real sza, Real albedo, Real o3col, Real altitude) const;

  bool hasJValue(const std::string & jname) const { return _tables.count(jname); }

private:
  void loadAxis(std::vector<Real> & axis, const std::string & filepath);

  /// 4D linear interpolation on regular grid
  Real interp4d(const Real * data, Real sza, Real alb, Real o3, Real alt) const;

  std::vector<Real> _sza_axis, _alb_axis, _o3_axis, _alt_axis;
  unsigned int _ns, _na, _no, _nh;
  std::map<std::string, std::vector<Real>> _tables;
};
