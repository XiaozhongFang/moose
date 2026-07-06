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

/**
 * Dilution term for box-model chemistry.
 *
 * Two modes:
 *   1. First-order: dC/dt -= kdil * (C - C_bkgd)
 *   2. Gaussian dispersion (F0AM tgauss): dC/dt -= 1/(tgauss + 2*t) * (C - C_bkgd)
 */
class DilutionHandler
{
public:
  DilutionHandler() = default;

  void setFirstOrder(Real kdil, const std::vector<Real> & conc_bkgd);
  void setGaussian(Real tgauss, const std::vector<Real> & conc_bkgd, Real t_start = 0.0);

  /// Apply dilution term to dC/dt in-place.
  void apply(Real t, const std::vector<Real> & C, std::vector<Real> & dC) const;

  bool active() const { return _active; }

private:
  bool _active = false;
  bool _use_gaussian = false;
  Real _kdil = 0.0;
  Real _tgauss = 0.0;
  Real _t_start = 0.0;
  std::vector<Real> _conc_bkgd;
};
