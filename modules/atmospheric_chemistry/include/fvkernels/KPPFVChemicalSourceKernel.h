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
 * FV chemical source term backed by KPP-generated RHS and dense Jacobian
 * material properties.
 */
class KPPFVChemicalSourceKernel : public FVElementalKernel
{
public:
  static InputParameters validParams();

  KPPFVChemicalSourceKernel(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

private:
  Real volumeWeight() const;

  const unsigned int _species_index;
  const Real _unit_conversion;
  const bool _spherical_volume_weight;
  const Real _radius;
  const Real _angle_scale;
  const Real _out_of_plane_scale;

  const MaterialProperty<std::vector<Real>> & _kpp_rhs;
  const MaterialProperty<std::vector<Real>> & _kpp_jacobian_dense;

  std::vector<const Moose::Functor<ADReal> *> _species_functors;
};
