//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AtmosphericSphericalFVTimeDerivative.h"

#include <cmath>

registerADMooseObject("AtmosphericChemistryApp", AtmosphericSphericalFVTimeDerivative);

namespace
{
constexpr Real spherical_time_pi = 3.1415926535897932384626433832795029;

Real
sphericalTimeAngleScale(const MooseEnum & units)
{
  return units == "degrees" ? spherical_time_pi / 180.0 : 1.0;
}
} // namespace

InputParameters
AtmosphericSphericalFVTimeDerivative::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addParam<Real>("radius", 6.378e6, "Sphere radius in meters.");
  MooseEnum coordinate_units("degrees radians", "degrees");
  params.addParam<MooseEnum>(
      "coordinate_units", coordinate_units, "Longitude/latitude coordinate units.");
  params.addParam<Real>(
      "out_of_plane_scale",
      1.0,
      "Scale from the mesh out-of-plane coordinate unit to meters. Use 1000 for km in a 3D "
      "longitude-latitude-height mesh.");
  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";
  params.addClassDescription(
      "Finite-volume time derivative weighted by spherical longitude-latitude cell area.");
  return params;
}

AtmosphericSphericalFVTimeDerivative::AtmosphericSphericalFVTimeDerivative(
    const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _u_dot(_var.adUDot()),
    _radius(getParam<Real>("radius")),
    _angle_scale(sphericalTimeAngleScale(getParam<MooseEnum>("coordinate_units"))),
    _out_of_plane_scale(getParam<Real>("out_of_plane_scale"))
{
  _var.requireQpComputations();

  if (_radius <= 0.0)
    paramError("radius", "The sphere radius must be positive.");
  if (_out_of_plane_scale <= 0.0)
    paramError("out_of_plane_scale", "The out-of-plane scale must be positive.");
}

ADReal
AtmosphericSphericalFVTimeDerivative::computeQpResidual()
{
  const Real phi = _q_point[_qp](1) * _angle_scale;
  return _radius * _radius * std::cos(phi) * _angle_scale * _angle_scale *
         _out_of_plane_scale * _u_dot[_qp];
}
