//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AtmosphericDensityWeightedDiffusion.h"
#include "Function.h"

registerMooseObject("AtmosphericChemistryApp", AtmosphericDensityWeightedDiffusion);

namespace
{
unsigned int
coordinateComponent(const MooseEnum & coordinate)
{
  if (coordinate == "x")
    return 0;
  if (coordinate == "y")
    return 1;
  return 2;
}
} // namespace

InputParameters
AtmosphericDensityWeightedDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<FunctionName>("diffusivity", "Turbulent diffusivity K.");
  params.addRequiredParam<FunctionName>("density", "Air number-density function rho.");
  MooseEnum component("x y z", "z");
  params.addParam<MooseEnum>(
      "component", component, "Coordinate direction for the density-weighted diffusion operator.");
  params.addClassDescription(
      "Density-weighted atmospheric turbulent diffusion d_i(rho K d_i(c / rho)).");
  return params;
}

AtmosphericDensityWeightedDiffusion::AtmosphericDensityWeightedDiffusion(
    const InputParameters & parameters)
  : Kernel(parameters),
    _diffusivity(getFunction("diffusivity")),
    _density(getFunction("density")),
    _component(coordinateComponent(getParam<MooseEnum>("component")))
{
}

Real
AtmosphericDensityWeightedDiffusion::computeQpResidual()
{
  const Point & p = _q_point[_qp];
  const Real rho = _density.value(_t, p);
  if (rho <= 0.0)
    mooseError("AtmosphericDensityWeightedDiffusion requires a positive density.");

  const Real grad_rho_over_rho = _density.gradient(_t, p)(_component) / rho;
  const Real flux =
      _diffusivity.value(_t, p) * (_grad_u[_qp](_component) - _u[_qp] * grad_rho_over_rho);

  return _grad_test[_i][_qp](_component) * flux;
}

Real
AtmosphericDensityWeightedDiffusion::computeQpJacobian()
{
  const Point & p = _q_point[_qp];
  const Real rho = _density.value(_t, p);
  if (rho <= 0.0)
    mooseError("AtmosphericDensityWeightedDiffusion requires a positive density.");

  const Real grad_rho_over_rho = _density.gradient(_t, p)(_component) / rho;

  return _diffusivity.value(_t, p) * _grad_test[_i][_qp](_component) *
         (_grad_phi[_j][_qp](_component) - _phi[_j][_qp] * grad_rho_over_rho);
}
