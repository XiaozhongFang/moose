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
#include <functional>
#include <algorithm>

/**
 * F0AM-style EKMA (Empirical Kinetics Modeling Approach) scanner.
 *
 * Runs a box model over a 2D grid of VOC and NOx scaling factors,
 * recording O3 peak concentration for each (VOC_factor, NOx_factor) pair.
 *
 * Reference: F0AM Setups/Examples/ExampleSetup_EKMA.m
 *
 * Usage:
 *   EKMAScanner scanner;
 *   scanner.setVOCFactors({0.5, 1.0, 2.0});
 *   scanner.setNOxFactors({0.5, 1.0, 2.0});
 *   scanner.setRunCallback([](Real voc, Real nox) -> Real {
 *     // Run box model with scaled VOC/NOx, return O3 peak
 *     return peak_O3;
 *   });
 *   auto grid = scanner.scan();  // 3x3 matrix
 */
class EKMAScanner
{
public:
  EKMAScanner() = default;

  /** Set VOC scaling factors (e.g. 0.1:0.1:2.0 → 20 values). */
  void setVOCFactors(const std::vector<Real> & factors) { _voc_factors = factors; }

  /** Set NOx scaling factors. */
  void setNOxFactors(const std::vector<Real> & factors) { _nox_factors = factors; }

  /** Set the 'run box model' callback.
   *  Callback receives (voc_factor, nox_factor) and returns O3 peak. */
  void setRunCallback(std::function<Real(Real, Real)> callback) { _callback = std::move(callback); }

  /** Run the full EKMA scan. Returns (NOx_idx x VOC_idx) matrix. */
  std::vector<std::vector<Real>> scan()
  {
    _o3_grid.clear();
    _o3_grid.resize(_nox_factors.size(),
                     std::vector<Real>(_voc_factors.size(), 0.0));

    for (size_t j = 0; j < _nox_factors.size(); ++j)
      for (size_t i = 0; i < _voc_factors.size(); ++i)
        if (_callback)
          _o3_grid[j][i] = _callback(_voc_factors[i], _nox_factors[j]);

    return _o3_grid;
  }

  /** Get the grid (NOx x VOC) of O3 peak values. */
  const std::vector<std::vector<Real>> & ozoneGrid() const { return _o3_grid; }

  /** Number of grid points. */
  size_t nPoints() const { return _voc_factors.size() * _nox_factors.size(); }

  /** Generate default linearly-spaced factors: start:step:end. */
  static std::vector<Real> linspace(Real start, Real end, unsigned int n)
  {
    std::vector<Real> result;
    if (n < 2) { result.push_back(start); return result; }
    Real step = (end - start) / (n - 1);
    for (unsigned int i = 0; i < n; ++i)
      result.push_back(start + i * step);
    return result;
  }

private:
  std::vector<Real> _voc_factors;
  std::vector<Real> _nox_factors;
  std::function<Real(Real, Real)> _callback;
  std::vector<std::vector<Real>> _o3_grid;
};
