//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998SpeciesIC.h"
#include "MAS1998BenchmarkUtils.h"

registerMooseObject("AtmosphericChemistryApp", MAS1998SpeciesIC);

namespace
{
unsigned int
speciesICHeightComponent(const MooseEnum & coordinate)
{
  if (coordinate == "x")
    return 0;
  if (coordinate == "y")
    return 1;
  return 2;
}
} // namespace

InputParameters
MAS1998SpeciesIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addRequiredParam<std::string>("species", "MAS1998 species name from Table 1.");
  MooseEnum coord("x y z", "z");
  params.addParam<MooseEnum>(
      "height_coordinate", coord, "Coordinate containing the height z in km.");
  params.addParam<Real>(
      "ground_number_density", 2.55e19, "Ground-level air number density in molecules/cm^3.");
  params.addParam<Real>(
      "plume_radius_degrees",
      20.0,
      "Angular radius of the HNO3/NO cylinder in degrees. The paper figure suggests about 20 "
      "degrees; calibrate this against the original benchmark data if recovered.");
  params.addParam<Real>(
      "plume_center_longitude_degrees", 0.0, "Cylinder center longitude lambda' in degrees.");
  params.addParam<Real>(
      "plume_center_latitude_degrees", 0.0, "Cylinder center latitude phi' in degrees.");
  params.addClassDescription(
      "MAS1998 benchmark initial condition in molecules/cm^3 using paper lambda'/phi' degrees and "
      "height in km.");
  return params;
}

MAS1998SpeciesIC::MAS1998SpeciesIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _species(getParam<std::string>("species")),
    _height_component(speciesICHeightComponent(getParam<MooseEnum>("height_coordinate"))),
    _ground_number_density(getParam<Real>("ground_number_density")),
    _plume_radius_degrees(getParam<Real>("plume_radius_degrees")),
    _plume_center_longitude_degrees(getParam<Real>("plume_center_longitude_degrees")),
    _plume_center_latitude_degrees(getParam<Real>("plume_center_latitude_degrees"))
{
  if (_plume_radius_degrees < 0.0)
    paramError("plume_radius_degrees", "The plume radius must be non-negative.");

  MAS1998::groundConcentration(_species);
}

Real
MAS1998SpeciesIC::value(const Point & p)
{
  Real concentration = MAS1998::groundConcentration(_species);
  if (MAS1998::hasCylinderInitialCondition(_species))
  {
    const Real distance = MAS1998::angularDistanceDegrees(
        p(0), p(1), _plume_center_longitude_degrees, _plume_center_latitude_degrees);
    if (distance <= _plume_radius_degrees)
      concentration = MAS1998::cylinderMaximumConcentration(_species);
  }

  const Real density_ratio =
      MAS1998::airNumberDensity(p(_height_component), _ground_number_density) /
      MAS1998::airNumberDensity(0.0, _ground_number_density);

  return concentration * density_ratio;
}
