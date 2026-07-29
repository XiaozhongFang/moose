//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "InitialCondition.h"

class MAS1998SpeciesIC : public InitialCondition
{
public:
  static InputParameters validParams();

  MAS1998SpeciesIC(const InputParameters & parameters);

  virtual Real value(const Point & p) override;

private:
  const std::string _species;
  const unsigned int _height_component;
  const Real _ground_number_density;
  const Real _plume_radius_degrees;
  const Real _plume_center_longitude_degrees;
  const Real _plume_center_latitude_degrees;
};
