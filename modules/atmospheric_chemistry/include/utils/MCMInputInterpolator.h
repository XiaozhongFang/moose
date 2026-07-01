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
#include <algorithm>
#include <cmath>

/**
 * F0AM-style input interpolation utility.
 *
 * Reads time-value pairs and interpolates at arbitrary query times.
 * Supports linear interpolation (default) and nearest-neighbor.
 *
 * Reference: F0AM Tools/InputInterp.m
 */
class MCMInputInterpolator
{
public:
  enum Method { LINEAR, NEAREST };

  MCMInputInterpolator() : _method(LINEAR) {}

  /** Set interpolation method. */
  void setMethod(Method m) { _method = m; }

  /** Load data from two vectors: time points and corresponding values. */
  void setData(const std::vector<Real> & times, const std::vector<Real> & values)
  {
    _times = times;
    _values = values;
  }

  /** Add a single (time, value) pair. */
  void addPoint(Real t, Real v)
  {
    _times.push_back(t);
    _values.push_back(v);
  }

  /** Sort data by time (required before use if points added out of order). */
  void sort()
  {
    // Pair sort
    std::vector<std::pair<Real, Real>> pairs;
    for (size_t i = 0; i < _times.size(); ++i)
      pairs.emplace_back(_times[i], _values[i]);
    std::sort(pairs.begin(), pairs.end());
    _times.clear();
    _values.clear();
    for (auto & p : pairs)
    {
      _times.push_back(p.first);
      _values.push_back(p.second);
    }
  }

  /** Interpolate at query time t. */
  Real interpolate(Real t) const
  {
    if (_times.empty()) return 0.0;
    if (_times.size() == 1) return _values[0];

    // Before first point: use first value (constant extrapolation)
    if (t <= _times.front()) return _values.front();
    // After last point: use last value
    if (t >= _times.back()) return _values.back();

    // Binary search for interval
    auto it = std::upper_bound(_times.begin(), _times.end(), t);
    size_t i = it - _times.begin() - 1;

    if (_method == NEAREST)
    {
      Real dt_left = t - _times[i];
      Real dt_right = _times[i + 1] - t;
      return (dt_left < dt_right) ? _values[i] : _values[i + 1];
    }

    // LINEAR: weighted average
    Real dt = _times[i + 1] - _times[i];
    if (std::abs(dt) < 1.0e-30) return _values[i];
    Real w = (t - _times[i]) / dt;
    return _values[i] * (1.0 - w) + _values[i + 1] * w;
  }

  /** Clear all data. */
  void clear()
  {
    _times.clear();
    _values.clear();
  }

  /** Whether data has been loaded. */
  bool empty() const { return _times.empty(); }

  /** Number of data points. */
  size_t size() const { return _times.size(); }

private:
  std::vector<Real> _times;
  std::vector<Real> _values;
  Method _method;
};
