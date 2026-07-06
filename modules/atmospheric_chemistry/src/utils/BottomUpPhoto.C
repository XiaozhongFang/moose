//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BottomUpPhoto.h"

#include <cmath>
#include <cstdlib>

BottomUpPhoto::BottomUpPhoto(const std::string & data_dir,
                             const std::string & flux_file,
                             const std::vector<unsigned int> & j_numbers)
  : _integrator(std::make_unique<BottomUpJIntegrator>(data_dir)),
    _j_numbers(j_numbers),
    _cached_T(0.0),
    _cached_P(0.0),
    _cache_valid(false)
{
  // Load lamp/actinic flux file
  _integrator->loadLampFlux(flux_file);

  // Load reaction mapping (bottomup_jmap.dat)
  _integrator->loadReactionMap("bottomup_jmap.dat");
}

void
BottomUpPhoto::computeJ(Real /*cos_sza*/,
                        Real /*sec_sza*/,
                        const PhysParams & params,
                        std::vector<Real> & j_out) const
{
  const size_t n = _j_numbers.size();
  j_out.assign(n, 0.0);
  const Real jfac = params.jfac;
  const Real T_cur = params.temperature;
  const Real P_cur = params.pressure > 0.0 ? params.pressure : 1013.25;

  // Re-compute J-values only when T or P changes significantly.
  if (!_cache_valid ||
      std::abs(T_cur - _cached_T) > 1.0e-6 ||
      std::abs(P_cur - _cached_P) > 1.0e-6)
  {
    _cached_j = _integrator->computeAllJ(T_cur, P_cur);
    _cached_T = T_cur;
    _cached_P = P_cur;
    _cache_valid = true;
  }

  // Map the cached results to output in J-number order.
  for (size_t i = 0; i < n; ++i)
  {
    const std::string jname = "J" + std::to_string(_j_numbers[i]);
    auto it = _cached_j.find(jname);
    if (it != _cached_j.end())
      j_out[i] = it->second * jfac;
  }
}
