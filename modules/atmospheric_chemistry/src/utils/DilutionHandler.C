//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DilutionHandler.h"

void
DilutionHandler::setFirstOrder(Real kdil, const std::vector<Real> & conc_bkgd)
{
  _kdil = kdil;
  _conc_bkgd = conc_bkgd;
  _use_gaussian = false;
  _active = true;
}

void
DilutionHandler::setGaussian(Real tgauss, const std::vector<Real> & conc_bkgd, Real t_start)
{
  _tgauss = tgauss;
  _conc_bkgd = conc_bkgd;
  _t_start = t_start;
  _use_gaussian = true;
  _active = true;
}

void
DilutionHandler::apply(Real t, const std::vector<Real> & C, std::vector<Real> & dC) const
{
  if (!_active)
    return;

  if (_use_gaussian)
  {
    // F0AM Gaussian dispersion: dilrate = -1/(tgauss + 2*(t + t_start))
    Real denom = _tgauss + 2.0 * (t + _t_start);
    if (denom > 1.0e-30)
    {
      Real rate = 1.0 / denom;
      for (size_t i = 0; i < C.size() && i < dC.size(); ++i)
        dC[i] -= rate * (C[i] - _conc_bkgd[i]);
    }
  }
  else
  {
    // First-order dilution (AtChem2 DILUTE)
    for (size_t i = 0; i < C.size() && i < dC.size(); ++i)
      dC[i] -= _kdil * (C[i] - _conc_bkgd[i]);
  }
}
