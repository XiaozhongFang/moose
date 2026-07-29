//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998BenchmarkUtils.h"
#include "MooseError.h"

#include "libmesh/utility.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace
{
constexpr Real mas1998_pi = 3.1415926535897932384626433832795029;
constexpr Real earth_radius_m = 6.378e6;
constexpr Real rotation_period_s = 14.0 * 24.0 * 3600.0;
constexpr Real beta = mas1998_pi / 4.0;

Real
degreesToRadians(const Real degrees)
{
  return degrees * mas1998_pi / 180.0;
}

Real
standardAtmosphereRelativeNumberDensity(const Real z_km)
{
  const Real z_m = std::max(z_km, 0.0) * 1000.0;

  constexpr Real g = 9.80665;
  constexpr Real r_air = 287.05287;

  Real base_z = 0.0;
  Real base_t = 288.15;
  Real base_p = 101325.0;

  const auto advance_layer =
      [&](const Real top_z, const Real lapse, Real & z0, Real & t0, Real & p0)
  {
    const Real dz = top_z - z0;
    if (std::abs(lapse) < 1.0e-14)
      p0 *= std::exp(-g * dz / (r_air * t0));
    else
    {
      const Real t1 = t0 + lapse * dz;
      p0 *= std::pow(t0 / t1, g / (r_air * lapse));
      t0 = t1;
    }
    z0 = top_z;
  };

  const auto eval_layer = [&](const Real lapse, const Real dz, const Real t0, const Real p0)
  {
    if (std::abs(lapse) < 1.0e-14)
      return std::make_pair(t0, p0 * std::exp(-g * dz / (r_air * t0)));

    const Real t = t0 + lapse * dz;
    return std::make_pair(t, p0 * std::pow(t0 / t, g / (r_air * lapse)));
  };

  if (z_m <= 11000.0)
  {
    const auto [t, p] = eval_layer(-0.0065, z_m, base_t, base_p);
    return (p / t) / (base_p / base_t);
  }

  advance_layer(11000.0, -0.0065, base_z, base_t, base_p);
  if (z_m <= 20000.0)
  {
    const auto [t, p] = eval_layer(0.0, z_m - base_z, base_t, base_p);
    return (p / t) / (101325.0 / 288.15);
  }

  advance_layer(20000.0, 0.0, base_z, base_t, base_p);
  if (z_m <= 32000.0)
  {
    const auto [t, p] = eval_layer(0.001, z_m - base_z, base_t, base_p);
    return (p / t) / (101325.0 / 288.15);
  }

  advance_layer(32000.0, 0.001, base_z, base_t, base_p);
  const auto [t, p] = eval_layer(0.0028, z_m - base_z, base_t, base_p);
  return (p / t) / (101325.0 / 288.15);
}

const std::map<std::string, Real> &
groundConcentrations()
{
  static const std::map<std::string, Real> values = {
      {"M", 2.55e19},    {"H2O", 2.55e17}, {"CO", 2.55e12},   {"O2", 5.3295e18}, {"HNO3", 2.55e9},
      {"HO2NO2", 1.0e2}, {"HNO2", 1.0e2},  {"H2O2", 1.0e2},   {"O3", 7.65e11},   {"HCHO", 1.0e2},
      {"CH3OOH", 1.0e2}, {"CH3O2", 1.0e2}, {"CH4", 4.335e13}, {"NO", 1.0e2},     {"NO2", 5.1e9},
      {"NO3", 1.0e2},    {"OH", 1.0e2},    {"HO2", 1.0e2},    {"N2O5", 1.0e2},   {"O1D", 0.0},
      {"O3P", 0.0},      {"O", 0.0}};
  return values;
}
} // namespace

namespace MAS1998
{
Real
verticalDiffusivity(const Real z_km)
{
  if (z_km <= 15.0)
    return 30.0;
  if (z_km <= 17.5)
    return 0.2;
  if (z_km <= 20.0)
    return 0.2 + 0.32 * (z_km - 17.5);
  return std::pow(10.0, 0.05 * z_km - 1.0);
}

Real
airNumberDensity(const Real z_km, const Real ground_number_density)
{
  return ground_number_density * standardAtmosphereRelativeNumberDensity(z_km);
}

Real
groundConcentration(const std::string & species)
{
  const auto & values = groundConcentrations();
  const auto it = values.find(species);
  if (it == values.end())
    mooseError("MAS1998: unsupported species '", species, "'");

  return it->second;
}

bool
hasCylinderInitialCondition(const std::string & species)
{
  return species == "HNO3" || species == "NO";
}

Real
cylinderMaximumConcentration(const std::string & species)
{
  if (species == "HNO3")
    return 4.0e9;
  if (species == "NO")
    return 1.0e9;

  return groundConcentration(species);
}

Real
solidBodyWind(const Real lambda_prime_degrees,
              const Real phi_prime_degrees,
              const unsigned int component)
{
  const Real lambda = degreesToRadians(lambda_prime_degrees + 180.0);
  const Real phi = degreesToRadians(phi_prime_degrees);
  const Real kappa = earth_radius_m / rotation_period_s;

  if (component == 0)
    return 2.0 * mas1998_pi * kappa *
           (std::cos(beta) * std::cos(phi) + std::sin(beta) * std::sin(phi) * std::cos(lambda));
  if (component == 1)
    return -2.0 * mas1998_pi * kappa * std::sin(beta) * std::sin(lambda);

  mooseError("MAS1998: wind component must be 0 (longitude/eastward) or 1 (latitude/northward)");
}

Real
angularDistanceDegrees(const Real lambda_prime_1,
                       const Real phi_prime_1,
                       const Real lambda_prime_2,
                       const Real phi_prime_2)
{
  const Real lambda_1 = degreesToRadians(lambda_prime_1);
  const Real phi_1 = degreesToRadians(phi_prime_1);
  const Real lambda_2 = degreesToRadians(lambda_prime_2);
  const Real phi_2 = degreesToRadians(phi_prime_2);
  const Real arg = std::sin(phi_1) * std::sin(phi_2) +
                   std::cos(phi_1) * std::cos(phi_2) * std::cos(lambda_1 - lambda_2);

  return std::acos(std::max(-1.0, std::min(1.0, arg))) * 180.0 / mas1998_pi;
}
} // namespace MAS1998
