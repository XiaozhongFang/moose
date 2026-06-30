#include "SurfaceStabilizationKernel.h"

registerMooseObject("CutFEMApp", SurfaceStabilizationKernel);

InputParameters
SurfaceStabilizationKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();

  params.addParam<Real>("gamma", 1.0,
    "Stabilization parameter gamma in [0, 1].");
  params.addParam<unsigned int>("k", 1,
    "Derivative order (k=1 for 2nd order PDE).");
  params.addParam<Real>("c_Gamma", 1.0,
    "Surface penalty constant (O(1)).");

  params.addClassDescription(
    "Surface ghost penalty stabilization for CutFEM.\n\n"
    "Penalizes normal derivatives on the interface Γ_h:\n"
    "  s_{h,Γ}(w,v) = Σ c'_j · h^{2(j-1+γ)} ∫_Γ D_n^j w · D_n^j v dS\n\n"
    "References:\n"
    "  - Larson & Zahedi (2020)");

  return params;
}

SurfaceStabilizationKernel::SurfaceStabilizationKernel(const InputParameters & params)
  : InterfaceKernel(params),
    _gamma(getParam<Real>("gamma")),
    _k(getParam<unsigned int>("k")),
    _c_Gamma(getParam<Real>("c_Gamma"))
{
  if (_gamma < 0.0 || _gamma > 1.0)
    paramError("gamma", "Must be in [0, 1]. Got " + std::to_string(_gamma));
  if (_k < 1)
    paramError("k", "Must be >= 1. Got " + std::to_string(_k));
  if (_c_Gamma <= 0.0)
    paramError("c_Gamma", "Must be positive. Got " + std::to_string(_c_Gamma));

  _exponent = 2.0 * (static_cast<Real>(_k) - 1.0 + _gamma);
}

Real
SurfaceStabilizationKernel::computePenaltyCoeff()
{
  unsigned int dim = _current_elem->dim();
  Real h_F = std::pow(_current_elem->volume(), 1.0 / static_cast<Real>(dim));
  return _c_Gamma * std::pow(h_F, _exponent);
}

Real
SurfaceStabilizationKernel::computeQpResidual(Moose::DGResidualType type)
{
  // Penalize normal derivative: D_n u * D_n v
  // D_n u = ∇u·n, D_n v = ∇v·n (same on both sides for surface Γ)
  Real penalty = computePenaltyCoeff();
  Real u_n = _grad_u[_qp] * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
    {
      Real v_n = _grad_test[_i][_qp] * _normals[_qp];
      return penalty * u_n * v_n * _JxW[_qp];
    }
    case Moose::Neighbor:
    {
      Real v_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      return penalty * u_n * v_n * _JxW[_qp];
    }
    default:
      mooseError("Unknown DGResidualType");
  }
}

Real
SurfaceStabilizationKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  Real penalty = computePenaltyCoeff();

  // Jacobian: ∂(penalty * D_n u * D_n v) / ∂u = penalty * D_n φ * D_n v
  // Since D_n is linear, this is same as residual with φ replacing u

  switch (type)
  {
    case Moose::ElementElement:
    {
      Real phi_n = _grad_phi[_j][_qp] * _normals[_qp];
      Real v_n = _grad_test[_i][_qp] * _normals[_qp];
      return penalty * phi_n * v_n * _JxW[_qp];
    }
    case Moose::ElementNeighbor:
    {
      Real phi_n = _grad_phi[_j][_qp] * _normals[_qp];
      Real v_n = _grad_test[_i][_qp] * _normals[_qp];
      return penalty * phi_n * v_n * _JxW[_qp];
    }
    case Moose::NeighborElement:
    {
      Real phi_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      Real v_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      return penalty * phi_n * v_n * _JxW[_qp];
    }
    case Moose::NeighborNeighbor:
    {
      Real phi_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      Real v_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      return penalty * phi_n * v_n * _JxW[_qp];
    }
    default:
      mooseError("Unknown DGJacobianType");
  }
}
