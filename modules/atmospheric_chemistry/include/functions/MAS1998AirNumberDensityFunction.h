//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "Function.h"

class MAS1998AirNumberDensityFunction : public Function
{
public:
  static InputParameters validParams();

  MAS1998AirNumberDensityFunction(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;
  virtual RealGradient gradient(Real t, const Point & p) const override;

private:
  const unsigned int _height_component;
  const Real _ground_number_density;
};
