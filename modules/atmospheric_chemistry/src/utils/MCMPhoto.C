//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMPhoto.h"

#include <cmath>

MCMPhoto::MCMPhoto(const std::vector<Real> & cl,
                   const std::vector<Real> & cmm,
                   const std::vector<Real> & cnn,
                   const std::vector<unsigned int> & j_numbers)
  : _cl(cl), _cmm(cmm), _cnn(cnn), _j_numbers(j_numbers)
{
}

void
MCMPhoto::computeJ(Real cos_sza,
                   Real sec_sza,
                   const PhysParams & params,
                   std::vector<Real> & j_out) const
{
  const size_t n = _cl.size();
  j_out.assign(n, 0.0);
  const Real jfac = params.jfac;

  // MCM SZA empirical formula:
  //   J = CL * cos(theta)^CMM * exp(-CNN * sec(theta)) * JFAC
  if (cos_sza > 1.0e-10)
  {
    for (size_t i = 0; i < n; ++i)
      j_out[i] = _cl[i] * std::pow(cos_sza, _cmm[i]) * std::exp(-_cnn[i] * sec_sza) * jfac;
  }
  // else cos_sza <= threshold → J = 0 (no sunlight)
}
