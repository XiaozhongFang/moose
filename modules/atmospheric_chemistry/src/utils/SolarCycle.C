//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolarCycle.h"
#include <cmath>

SolarCycle::SolarCycle(Real latitude, Real longitude, int day, int month, int year)
  : _lat(latitude), _lon(longitude), _day(day), _month(month), _year(year)
{
  // Compute day of year
  static const unsigned int days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  unsigned int doy = 0;
  for (int m = 0; m < month - 1; ++m)
    doy += days_in_month[m];
  doy += day;

  // Solar declination (radians) — Madronich 1993
  _declination = 23.45 * M_PI / 180.0 *
      sin(2.0 * M_PI * (284.0 + static_cast<Real>(doy)) / 365.0);

  // Equation of time (minutes)
  Real B = 2.0 * M_PI * (static_cast<Real>(doy) - 81.0) / 364.0;
  _eot = 9.87 * sin(2.0 * B) - 7.53 * cos(B) - 1.5 * sin(B);

  // Precompute sin(lat)*sin(dec) and cos(lat)*cos(dec) — time-invariant
  Real latrad = _lat * M_PI / 180.0;
  _sinld = sin(latrad) * sin(_declination);
  _cosld = cos(latrad) * cos(_declination);
}

Real
SolarCycle::computeCosSZA(Real seconds) const
{
  // Solar time (hours)
  Real solar_time = seconds / 3600.0 + _eot / 60.0;

  // Hour angle (radians)
  Real hour_angle = (solar_time / 12.0 - 1.0) * M_PI;

  // Cosine of solar zenith angle
  Real cos_sza = _sinld + _cosld * cos(hour_angle);

  return std::max(cos_sza, 0.0);
}
