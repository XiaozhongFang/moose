//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998VerticalDiffusivityFunction.h"
#include "MAS1998BenchmarkUtils.h"

registerMooseObject("AtmosphericChemistryApp", MAS1998VerticalDiffusivityFunction);

namespace
{
unsigned int
verticalDiffusivityHeightComponent(const MooseEnum & coordinate)
{
  if (coordinate == "x")
    return 0;
  if (coordinate == "y")
    return 1;
  return 2;
}
} // namespace

InputParameters
MAS1998VerticalDiffusivityFunction::validParams()
{
  InputParameters params = Function::validParams();
  MooseEnum coord("x y z", "z");
  params.addParam<MooseEnum>(
      "height_coordinate", coord, "Coordinate containing the height z in km.");
  params.addClassDescription(
      "Vertical turbulent diffusivity K(z) from Spee et al. (1998), Eq. (4.2), in m^2/s.");
  return params;
}

MAS1998VerticalDiffusivityFunction::MAS1998VerticalDiffusivityFunction(
    const InputParameters & parameters)
  : Function(parameters),
    _height_component(verticalDiffusivityHeightComponent(getParam<MooseEnum>("height_coordinate")))
{
}

Real
MAS1998VerticalDiffusivityFunction::value(Real /*t*/, const Point & p) const
{
  return MAS1998::verticalDiffusivity(p(_height_component));
}
