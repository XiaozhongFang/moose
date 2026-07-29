//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MAS1998SolidBodyWindComponent.h"
#include "MAS1998BenchmarkUtils.h"

registerMooseObject("AtmosphericChemistryApp", MAS1998SolidBodyWindComponent);

namespace
{
unsigned int
windComponent(const MooseEnum & component)
{
  return component == "u" ? 0 : 1;
}
} // namespace

InputParameters
MAS1998SolidBodyWindComponent::validParams()
{
  InputParameters params = Function::validParams();
  MooseEnum component("u v", "u");
  params.addParam<MooseEnum>("component", component, "Wind component: u or v from Eq. (4.1).");
  params.addClassDescription(
      "MAS1998 solid-body rotation wind component in m/s using lambda' and phi' coordinates in "
      "degrees.");
  return params;
}

MAS1998SolidBodyWindComponent::MAS1998SolidBodyWindComponent(const InputParameters & parameters)
  : Function(parameters), _component(windComponent(getParam<MooseEnum>("component")))
{
}

Real
MAS1998SolidBodyWindComponent::value(Real /*t*/, const Point & p) const
{
  return MAS1998::solidBodyWind(p(0), p(1), _component);
}
