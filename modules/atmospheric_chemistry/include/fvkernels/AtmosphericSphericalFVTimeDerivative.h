//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "FVElementalKernel.h"

/**
 * FV time derivative weighted by spherical shell surface area:
 *
 *   a^2 cos(phi) dc/dt
 *
 * for longitude-latitude coordinates.
 */
class AtmosphericSphericalFVTimeDerivative : public FVElementalKernel
{
public:
  static InputParameters validParams();

  AtmosphericSphericalFVTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  const ADVariableValue & _u_dot;
  const Real _radius;
  const Real _angle_scale;
  const Real _out_of_plane_scale;
};
