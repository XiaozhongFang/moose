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

/**
 * Solar cycle parameters (Madronich 1993).
 *
 * Computes solar zenith angle, declination, equation of time
 * for given latitude/longitude and day of year.
 */
class SolarCycle
{
public:
  SolarCycle(Real latitude, Real longitude, int day, int month, int year);

  /// Compute cos(SZA) at given seconds since midnight.
  Real computeCosSZA(Real seconds) const;

  Real getDeclination() const { return _declination; }
  Real getSinLD() const { return _sinld; }
  Real getCosLD() const { return _cosld; }
  Real getEQT() const { return _eot; }
  Real getLatitude() const { return _lat; }
  Real getLongitude() const { return _lon; }

private:
  Real _lat, _lon;
  Real _declination, _eot;
  Real _sinld, _cosld;
  int _day, _month, _year;
};
