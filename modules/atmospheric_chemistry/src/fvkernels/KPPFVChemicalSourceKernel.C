//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "KPPFVChemicalSourceKernel.h"
#include "metaphysicl/raw_type.h"

#include <cmath>

registerADMooseObject("AtmosphericChemistryApp", KPPFVChemicalSourceKernel);

namespace
{
constexpr Real kpp_fv_source_pi = 3.1415926535897932384626433832795029;

Real
kppFVSourceAngleScale(const MooseEnum & units)
{
  return units == "degrees" ? kpp_fv_source_pi / 180.0 : 1.0;
}
} // namespace

InputParameters
KPPFVChemicalSourceKernel::validParams()
{
  InputParameters params = FVElementalKernel::validParams();

  params.addRequiredParam<unsigned int>("species_index",
                                        "Row index of this species in KPP VAR order.");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "all_species", "All chemical species variables/functors in KPP VAR order.");
  params.addParam<Real>("unit_conversion", 1.0,
                        "Variable-unit to molec/cm^3 conversion factor. Use M/1e9 for ppb.");
  params.addParam<bool>(
      "spherical_volume_weight",
      false,
      "Apply longitude-latitude spherical volume weights to the source term.");
  params.addParam<Real>("radius", 6.378e6, "Sphere radius in meters.");
  MooseEnum coordinate_units("degrees radians", "degrees");
  params.addParam<MooseEnum>(
      "coordinate_units", coordinate_units, "Longitude/latitude coordinate units.");
  params.addParam<Real>(
      "out_of_plane_scale",
      1.0,
      "Scale from the mesh out-of-plane coordinate unit to meters. Use 1000 for km in a 3D "
      "longitude-latitude-height mesh.");

  params.addClassDescription(
      "Finite-volume KPP chemical source term with KPP analytical Jacobian linearization.");
  return params;
}

KPPFVChemicalSourceKernel::KPPFVChemicalSourceKernel(const InputParameters & params)
  : FVElementalKernel(params),
    _species_index(getParam<unsigned int>("species_index")),
    _unit_conversion(getParam<Real>("unit_conversion")),
    _spherical_volume_weight(getParam<bool>("spherical_volume_weight")),
    _radius(getParam<Real>("radius")),
    _angle_scale(kppFVSourceAngleScale(getParam<MooseEnum>("coordinate_units"))),
    _out_of_plane_scale(getParam<Real>("out_of_plane_scale")),
    _kpp_rhs(getMaterialProperty<std::vector<Real>>("kpp_rhs")),
    _kpp_jacobian_dense(getMaterialProperty<std::vector<Real>>("kpp_jacobian_dense"))
{
  const auto & species_names = getParam<std::vector<MooseFunctorName>>("all_species");
  if (_species_index >= species_names.size())
    paramError("species_index",
               "The species index must be smaller than the all_species vector size.");

  _species_functors.reserve(species_names.size());
  for (const auto & species_name : species_names)
    _species_functors.push_back(&getFunctor<ADReal>(species_name));

  if (_unit_conversion <= 0.0)
    paramError("unit_conversion", "The unit conversion must be positive.");
  if (_radius <= 0.0)
    paramError("radius", "The sphere radius must be positive.");
  if (_out_of_plane_scale <= 0.0)
    paramError("out_of_plane_scale", "The out-of-plane scale must be positive.");
}

Real
KPPFVChemicalSourceKernel::volumeWeight() const
{
  if (!_spherical_volume_weight)
    return 1.0;

  const Real phi = _q_point[_qp](1) * _angle_scale;
  return _radius * _radius * std::cos(phi) * _angle_scale * _angle_scale *
         _out_of_plane_scale;
}

ADReal
KPPFVChemicalSourceKernel::computeQpResidual()
{
  const auto & rhs = _kpp_rhs[_qp];
  const auto & jac = _kpp_jacobian_dense[_qp];
  const auto n_species = _species_functors.size();
  const auto required = n_species * n_species;

  if (rhs.size() <= _species_index)
    mooseError("KPPFVChemicalSourceKernel: kpp_rhs has size ",
               rhs.size(),
               " but species_index is ",
               _species_index);
  if (jac.size() < required)
    mooseError("KPPFVChemicalSourceKernel: kpp_jacobian_dense has size ",
               jac.size(),
               " but expected at least ",
               required);

  ADReal source = rhs[_species_index] / _unit_conversion;
  const auto elem = makeElemArg(_current_elem);
  const auto state = determineState();
  for (const auto col : make_range(n_species))
  {
    const ADReal c = (*_species_functors[col])(elem, state);
    source += jac[_species_index * n_species + col] * (c - MetaPhysicL::raw_value(c));
  }

  return -volumeWeight() * source;
}
