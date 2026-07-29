//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "AtmosphericSphericalFVAdvection.h"
#include "FEProblemBase.h"
#include "Function.h"
#include "Steady.h"
#include "metaphysicl/raw_type.h"

#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"

#include <algorithm>
#include <cmath>

registerADMooseObject("AtmosphericChemistryApp", AtmosphericSphericalFVAdvection);

namespace
{
constexpr Real spherical_advection_pi = 3.1415926535897932384626433832795029;

Real
sphericalAdvectionAngleScale(const MooseEnum & units)
{
  return units == "degrees" ? spherical_advection_pi / 180.0 : 1.0;
}

Real
sphericalAdvectionLongitudePeriod(const MooseEnum & units, const Real longitude_period)
{
  if (longitude_period > 0.0)
    return longitude_period;

  return units == "degrees" ? 360.0 : 2.0 * spherical_advection_pi;
}
} // namespace

InputParameters
AtmosphericSphericalFVAdvection::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addRequiredParam<FunctionName>("u_wind", "Eastward/zonal wind function in m/s.");
  params.addRequiredParam<FunctionName>("v_wind", "Northward/meridional wind function in m/s.");
  params.addParam<Real>("radius", 6.378e6, "Sphere radius in meters.");
  MooseEnum coordinate_units("degrees radians", "degrees");
  params.addParam<MooseEnum>(
      "coordinate_units", coordinate_units, "Longitude/latitude coordinate units.");
  params.addParam<Real>(
      "out_of_plane_scale",
      1.0,
      "Scale from the mesh out-of-plane coordinate unit to meters. Use 1000 for km in a 3D "
      "longitude-latitude-height mesh.");
  params.addParam<Real>("longitude_period",
                        0.0,
                        "Longitude period in mesh coordinate units. A non-positive value uses 360 "
                        "degrees or 2*pi radians based on coordinate_units.");
  params += Moose::FV::advectedInterpolationParameter();
  MooseEnum flux_scheme("framework mas1998_limited mas1998_unlimited", "framework");
  params.addParam<MooseEnum>(
      "flux_scheme",
      flux_scheme,
      "Face value scheme. 'framework' uses advected_interp_method. The MAS1998 schemes use the "
      "third-order upwind reconstruction from the benchmark, with or without the published "
      "positivity limiter.");

  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC |
          Moose::RelationshipManagerType::COUPLING,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      { FVRelationshipManagerInterface::setRMParamsAdvection(obj_params, rm_params, 2); });

  params.addClassDescription(
      "Finite-volume conservative horizontal advection on a longitude-latitude spherical shell.");
  return params;
}

AtmosphericSphericalFVAdvection::AtmosphericSphericalFVAdvection(const InputParameters & params)
  : FVFluxKernel(params),
    _u_wind(getFunction("u_wind")),
    _v_wind(getFunction("v_wind")),
    _radius(getParam<Real>("radius")),
    _angle_scale(sphericalAdvectionAngleScale(getParam<MooseEnum>("coordinate_units"))),
    _out_of_plane_scale(getParam<Real>("out_of_plane_scale")),
    _longitude_period(sphericalAdvectionLongitudePeriod(getParam<MooseEnum>("coordinate_units"),
                                                        getParam<Real>("longitude_period"))),
    _flux_scheme(getParam<MooseEnum>("flux_scheme") == "mas1998_limited"
                     ? FluxScheme::MAS1998Limited
                     : (getParam<MooseEnum>("flux_scheme") == "mas1998_unlimited"
                            ? FluxScheme::MAS1998Unlimited
                            : FluxScheme::Framework))
{
  const bool need_more_ghosting =
      Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
  if (need_more_ghosting && _tid == 0)
    getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")
        ->setErrorOnJacobianNonzeroReallocation(false);

  if (dynamic_cast<Steady *>(_app.getExecutioner()))
  {
    const MooseEnum not_available_with_steady("sou min_mod vanLeer quick venkatakrishnan");
    const std::string chosen_scheme =
        static_cast<std::string>(getParam<MooseEnum>("advected_interp_method"));
    if (not_available_with_steady.find(chosen_scheme) != not_available_with_steady.items().end())
      paramError("advected_interp_method",
                 "The given advected interpolation cannot be used with steady-state runs!");
  }

  if (_radius <= 0.0)
    paramError("radius", "The sphere radius must be positive.");
  if (_out_of_plane_scale <= 0.0)
    paramError("out_of_plane_scale", "The out-of-plane scale must be positive.");
  if (_longitude_period <= 0.0)
    paramError("longitude_period", "The longitude period must be positive.");
}

ADReal
AtmosphericSphericalFVAdvection::advectedFaceValue(const Real normal_flux) const
{
  const auto state = determineState();
  const auto & limiter_time = _subproblem.isTransient()
                                  ? Moose::StateArg(1, Moose::SolutionIterationType::Time)
                                  : Moose::StateArg(1, Moose::SolutionIterationType::Nonlinear);
  const bool elem_is_upwind = normal_flux >= 0.0;
  const auto face = makeFace(*_face_info,
                             Moose::FV::limiterType(_advected_interp_method),
                             elem_is_upwind,
                             false,
                             &limiter_time);

  return _var(face, state);
}

Real
AtmosphericSphericalFVAdvection::coordinateDelta(const Point & target,
                                                 const Point & origin,
                                                 const unsigned int axis) const
{
  Real delta = target(axis) - origin(axis);
  if (axis == 0)
  {
    const Real half_period = 0.5 * _longitude_period;
    if (delta > half_period)
      delta -= _longitude_period;
    else if (delta < -half_period)
      delta += _longitude_period;
  }

  return delta;
}

const libMesh::Elem *
AtmosphericSphericalFVAdvection::farUpwindElem(const libMesh::Elem & upwind_elem,
                                               const libMesh::Elem & downwind_elem,
                                               const unsigned int axis) const
{
  const auto upwind_centroid = upwind_elem.vertex_average();
  const auto downwind_centroid = downwind_elem.vertex_average();
  const Real upstream_direction = coordinateDelta(upwind_centroid, downwind_centroid, axis);
  const libMesh::Elem * best_elem = nullptr;
  Real best_projection = 0.0;

  for (const auto side : make_range(upwind_elem.n_sides()))
  {
    const auto * const candidate = upwind_elem.neighbor_ptr(side);
    if (!candidate || candidate == libMesh::remote_elem || candidate == &downwind_elem ||
        !candidate->active())
      continue;

    const auto candidate_centroid = candidate->vertex_average();
    const Real projection =
        coordinateDelta(candidate_centroid, upwind_centroid, axis) * upstream_direction;
    if (projection > best_projection)
    {
      best_projection = projection;
      best_elem = candidate;
    }
  }

  return best_elem;
}

Real
AtmosphericSphericalFVAdvection::limiterValue(const Real theta) const
{
  const Real third_order = 1.0 / 3.0 + theta / 6.0;
  if (_flux_scheme == FluxScheme::MAS1998Unlimited)
    return third_order;

  return std::max(0.0, std::min(std::min(1.0, theta), third_order));
}

ADReal
AtmosphericSphericalFVAdvection::mas1998FaceValue(const Real normal_flux,
                                                  const unsigned int axis) const
{
  const auto state = determineState();
  const bool elem_is_upwind = normal_flux >= 0.0;
  const auto * const upwind_elem = elem_is_upwind ? _face_info->elemPtr() : _face_info->neighborPtr();
  const auto * const downwind_elem =
      elem_is_upwind ? _face_info->neighborPtr() : _face_info->elemPtr();

  if (!upwind_elem || !downwind_elem)
    return advectedFaceValue(normal_flux);

  const auto upwind = makeElemArg(upwind_elem);
  const auto downwind = makeElemArg(downwind_elem);
  const ADReal c_upwind = _var(upwind, state);
  const ADReal c_downwind = _var(downwind, state);
  const auto * const far_upwind_elem = farUpwindElem(*upwind_elem, *downwind_elem, axis);
  if (!far_upwind_elem)
    return c_upwind;

  const auto far_upwind = makeElemArg(far_upwind_elem);
  const ADReal c_far_upwind = _var(far_upwind, state);
  const Real denominator = MetaPhysicL::raw_value(c_downwind - c_upwind);
  if (std::abs(denominator) <= 1e-30)
    return c_upwind;

  const Real theta = MetaPhysicL::raw_value(c_upwind - c_far_upwind) / denominator;
  const Real phi = limiterValue(theta);
  return c_upwind + phi * (c_downwind - c_upwind);
}

ADReal
AtmosphericSphericalFVAdvection::computeQpResidual()
{
  const Point & p = _face_info->faceCentroid();
  const Real phi = p(1) * _angle_scale;
  const Real cos_phi = std::cos(phi);

  if (std::abs(_normal(0)) >= std::abs(_normal(1)) && std::abs(_normal(0)) > 0.5)
  {
    const Real normal_flux =
        _normal(0) * _radius * _angle_scale * _out_of_plane_scale * _u_wind.value(_t, p);
    const ADReal c_face = _flux_scheme == FluxScheme::Framework ? advectedFaceValue(normal_flux)
                                                                : mas1998FaceValue(normal_flux, 0);
    return normal_flux * c_face;
  }

  if (std::abs(_normal(1)) > 0.5)
  {
    const Real normal_flux = _normal(1) * _radius * _angle_scale * cos_phi *
                             _out_of_plane_scale * _v_wind.value(_t, p);
    const ADReal c_face = _flux_scheme == FluxScheme::Framework ? advectedFaceValue(normal_flux)
                                                                : mas1998FaceValue(normal_flux, 1);
    return normal_flux * c_face;
  }

  return 0.0;
}
