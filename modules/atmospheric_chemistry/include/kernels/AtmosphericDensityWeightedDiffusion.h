//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "Kernel.h"

class Function;

/**
 * One-dimensional turbulent diffusion in density-weighted atmospheric form:
 *
 *   d/dx_i (rho K d(c / rho)/dx_i)
 */
class AtmosphericDensityWeightedDiffusion : public Kernel
{
public:
  static InputParameters validParams();

  AtmosphericDensityWeightedDiffusion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:
  const Function & _diffusivity;
  const Function & _density;
  const unsigned int _component;
  const Real _coordinate_scale;
};
