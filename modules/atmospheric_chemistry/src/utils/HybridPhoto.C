//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HybridPhoto.h"

#include <cmath>

HybridPhoto::HybridPhoto(const std::string & table_dir,
                         const std::vector<unsigned int> & j_numbers)
  : _reader(std::make_unique<HybridJTableReader>(table_dir)),
    _j_numbers(j_numbers)
{
}

void
HybridPhoto::computeJ(Real cos_sza,
                      Real /*sec_sza*/,
                      const PhysParams & params,
                      std::vector<Real> & j_out) const
{
  const size_t n = _j_numbers.size();
  j_out.assign(n, 0.0);
  const Real jfac = params.jfac;
  const Real albedo = params.albedo;
  const Real o3column = params.o3column;
  const Real altitude = params.altitude;

  // F0AM Hybrid scheme: log10(J) lookup table; SZA >= 90° → returns -∞ (→ 0.0).
  if (cos_sza > 1.0e-10)
  {
    const Real sza_deg = std::acos(cos_sza) * (180.0 / 3.14159265358979323846);
    for (size_t i = 0; i < n; ++i)
    {
      // Build J-name from j number: J<N>
      const std::string jname = "J" + std::to_string(_j_numbers[i]);
      const Real log10J = _reader->interpolate(jname, sza_deg, albedo, o3column, altitude);
      // Convert back to linear J = 10^log10J * JFAC
      j_out[i] = std::pow(10.0, log10J) * jfac;
    }
  }
  // else cos_sza <= threshold → J = 0 (SZA >= 90°, table returns 0)
}
