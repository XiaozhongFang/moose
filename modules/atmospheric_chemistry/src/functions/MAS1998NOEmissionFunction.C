//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998NOEmissionFunction.h"

registerMooseObject("AtmosphericChemistryApp", MAS1998NOEmissionFunction);

namespace
{
unsigned int
noEmissionHeightComponent(const MooseEnum & coordinate)
{
  if (coordinate == "x")
    return 0;
  if (coordinate == "y")
    return 1;
  return 2;
}
} // namespace

InputParameters
MAS1998NOEmissionFunction::validParams()
{
  InputParameters params = Function::validParams();
  MooseEnum coord("x y z", "z");
  params.addParam<MooseEnum>(
      "height_coordinate", coord, "Coordinate containing the height z in km.");
  params.addParam<Real>("ground_layer_top", 0.65, "Top of the lowest MAS1998 vertical cell in km.");
  params.addParam<Real>("rate", 1.0e4, "NO emission rate in molecules/(cm^3 s).");
  params.addClassDescription("MAS1998 lowest-layer NO emission source function.");
  return params;
}

MAS1998NOEmissionFunction::MAS1998NOEmissionFunction(const InputParameters & parameters)
  : Function(parameters),
    _height_component(noEmissionHeightComponent(getParam<MooseEnum>("height_coordinate"))),
    _ground_layer_top(getParam<Real>("ground_layer_top")),
    _rate(getParam<Real>("rate"))
{
}

Real
MAS1998NOEmissionFunction::value(Real /*t*/, const Point & p) const
{
  return p(_height_component) <= _ground_layer_top ? _rate : 0.0;
}
