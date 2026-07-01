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
#include <map>
#include <cmath>

/**
 * F0AM-style J-value calibration (jcorr).
 *
 * When observed J-values are available for one or more photolysis reactions,
 * computes correction factors and applies them to all photolysis rates.
 *
 * Algorithm (F0AM jcorr_all):
 *   1. For each constrained J: jcorr_i = J_obs_i / J_param_i
 *   2. Average across all constrained J → global jcorr
 *   3. Apply: J_used_i = J_param_i * jcorr_mean
 *
 * Reference: F0AM_ModelCore.m ~ lines 113-117
 */
class JCalibrator
{
public:
  JCalibrator() : _global_factor(1.0) {}

  /** Add a constrained (observed) J value for a specific photolysis reaction. */
  void addObservedJ(unsigned int j_number, Real observed_value)
  {
    _observed_J[j_number] = observed_value;
  }

  /** Clear all observed J data. */
  void clearObserved() { _observed_J.clear(); }

  /**
   * Compute calibration factors given parameterized J values.
   * Must be called after photolysis update and before applying correction.
   */
  void calibrate(const std::map<unsigned int, Real> & param_J_values)
  {
    _per_J_factor.clear();
    Real sum_factor = 0.0;
    size_t n = 0;

    for (auto & [jnum, obs] : _observed_J)
    {
      auto it = param_J_values.find(jnum);
      if (it != param_J_values.end() && std::abs(it->second) > 1e-30)
      {
        Real f = obs / it->second;
        _per_J_factor[jnum] = f;
        sum_factor += f;
        n++;
      }
    }

    _global_factor = (n > 0) ? sum_factor / n : 1.0;
  }

  /** Get global calibration factor (mean of all per-J corrections). */
  Real globalFactor() const { return _global_factor; }

  /** Get per-J calibration factor. Returns global_factor if no per-J data. */
  Real factorForJ(unsigned int j_number) const
  {
    auto it = _per_J_factor.find(j_number);
    return (it != _per_J_factor.end()) ? it->second : _global_factor;
  }

  /**
   * Apply calibration to a set of J values (in-place).
   * Each J_value[i] is multiplied by the corresponding calibration factor.
   */
  void applyTo(std::map<unsigned int, Real> & J_values) const
  {
    for (auto & [jnum, val] : J_values)
      val *= factorForJ(jnum);
  }

  /** Whether any observed J data has been provided. */
  bool hasObservedData() const { return !_observed_J.empty(); }

  /** Number of observed J constraints. */
  unsigned int nObserved() const { return _observed_J.size(); }

private:
  std::map<unsigned int, Real> _observed_J;
  std::map<unsigned int, Real> _per_J_factor;
  Real _global_factor;
};
