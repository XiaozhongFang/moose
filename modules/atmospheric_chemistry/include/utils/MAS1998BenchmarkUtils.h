//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "MooseTypes.h"

#include <string>

namespace MAS1998
{
Real verticalDiffusivity(Real z_km);
Real airNumberDensity(Real z_km, Real ground_number_density);
Real groundConcentration(const std::string & species);
bool hasCylinderInitialCondition(const std::string & species);
Real cylinderMaximumConcentration(const std::string & species);
Real solidBodyWind(Real lambda_prime_degrees, Real phi_prime_degrees, unsigned int component);
Real angularDistanceDegrees(Real lambda_prime_1,
                            Real phi_prime_1,
                            Real lambda_prime_2,
                            Real phi_prime_2);
} // namespace MAS1998
