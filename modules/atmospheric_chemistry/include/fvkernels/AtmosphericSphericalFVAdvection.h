//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "FVFluxKernel.h"

class Function;
namespace libMesh
{
class Elem;
}

/**
 * Horizontal finite-volume advection on a longitude-latitude spherical shell.
 */
class AtmosphericSphericalFVAdvection : public FVFluxKernel
{
public:
  static InputParameters validParams();

  AtmosphericSphericalFVAdvection(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

private:
  enum class FluxScheme
  {
    Framework,
    MAS1998Limited,
    MAS1998Unlimited
  };

  ADReal advectedFaceValue(Real normal_flux) const;
  ADReal mas1998FaceValue(Real normal_flux, unsigned int axis) const;
  Real coordinateDelta(const Point & target, const Point & origin, unsigned int axis) const;
  const libMesh::Elem * farUpwindElem(const libMesh::Elem & upwind_elem,
                                      const libMesh::Elem & downwind_elem,
                                      unsigned int axis) const;
  Real limiterValue(Real theta) const;

  const Function & _u_wind;
  const Function & _v_wind;
  const Real _radius;
  const Real _angle_scale;
  const Real _out_of_plane_scale;
  const Real _longitude_period;
  Moose::FV::InterpMethod _advected_interp_method;
  const FluxScheme _flux_scheme;
};
