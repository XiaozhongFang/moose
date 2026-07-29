//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998AirNumberDensityFunction.h"
#include "MAS1998BenchmarkUtils.h"

registerMooseObject("AtmosphericChemistryApp", MAS1998AirNumberDensityFunction);

namespace
{
unsigned int
airNumberDensityHeightComponent(const MooseEnum & coordinate)
{
  if (coordinate == "x")
    return 0;
  if (coordinate == "y")
    return 1;
  return 2;
}
} // namespace

InputParameters
MAS1998AirNumberDensityFunction::validParams()
{
  InputParameters params = Function::validParams();
  MooseEnum coord("x y z", "z");
  params.addParam<MooseEnum>(
      "height_coordinate", coord, "Coordinate containing the height z in km.");
  params.addParam<Real>(
      "ground_number_density", 2.55e19, "Ground-level air number density in molecules/cm^3.");
  params.addClassDescription(
      "US Standard Atmosphere number-density profile normalized to the MAS1998 ground density.");
  return params;
}

MAS1998AirNumberDensityFunction::MAS1998AirNumberDensityFunction(const InputParameters & parameters)
  : Function(parameters),
    _height_component(airNumberDensityHeightComponent(getParam<MooseEnum>("height_coordinate"))),
    _ground_number_density(getParam<Real>("ground_number_density"))
{
}

Real
MAS1998AirNumberDensityFunction::value(Real /*t*/, const Point & p) const
{
  return MAS1998::airNumberDensity(p(_height_component), _ground_number_density);
}
