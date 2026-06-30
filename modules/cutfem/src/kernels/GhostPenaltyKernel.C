#include "GhostPenaltyKernel.h"

registerMooseObject("CutFEMApp", GhostPenaltyKernel);

InputParameters
GhostPenaltyKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();

  params.addParam<Real>(
      "gamma",
      1.0,
      "Stabilization parameter gamma in [0, 1]. 1.0 (default) gives maximum stabilization.");
  params.addParam<unsigned int>(
      "k",
      1,
      "Derivative order (k=1 for 2nd order PDE, k=2 for 4th order). Must be >= 1.");
  params.addParam<Real>("c_F",
                        1.0,
                        "Face penalty constant. Controls stabilization strength. Default is 1.0.");

  params.addClassDescription(
      "Ghost Penalty stabilization kernel for CutFEM.\n\n"
      "Implements face-based gradient jump penalty:\n"
      "  s_h(w,v) = sum_j c_F * h^{2(j-1+gamma)} \\int_F [D_n^j w] . [D_n^j v] dS\n"
      "where [D_n^j v] is the jump of the j-th normal derivative across face F.\n"
      "\n"
      "References:\n"
      "  - Burman, Hansbo, Larson (2015) CutFEM method\n"
      "  - Larson & Zahedi (2020) High-order CutFEM stabilization");

  return params;
}

GhostPenaltyKernel::GhostPenaltyKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _gamma(getParam<Real>("gamma")),
    _k(getParam<unsigned int>("k")),
    _c_F(getParam<Real>("c_F"))
{
  // Validate parameters
  if (_gamma < 0.0 || _gamma > 1.0)
    paramError("gamma", "Must be in range [0.0, 1.0]. Received: " + std::to_string(_gamma));

  if (_k < 1)
    paramError("k", "Must be >= 1. Received: " + std::to_string(_k));

  if (_c_F <= 0.0)
    paramError("c_F", "Must be positive. Received: " + std::to_string(_c_F));

  // Compute the exponent for the penalty coefficient:
  // exponent = 2(k - 1 + gamma)
  _exponent = 2.0 * (static_cast<Real>(_k) - 1.0 + _gamma);
}

Real
GhostPenaltyKernel::computeQpResidual(Moose::DGResidualType type)
{
  // Compute normal component of element gradient
  Real u_grad_n_element = _grad_u[_qp] * _normals[_qp];

  // Compute normal component of neighbor gradient
  Real u_grad_n_neighbor = _grad_neighbor_value[_qp] * _normals[_qp];

  // Jump of normal gradient: [u_n] = u_n^+ - u_n^-
  Real u_jump = u_grad_n_element - u_grad_n_neighbor;

  // Compute penalty coefficient
  Real penalty = computePenaltyCoeff();

  // Residual: sign * penalty * [u_n] * (test_gradient·n) * JxW
  // Element side: positive sign, using element test function
  // Neighbor side: negative sign (normal flips), using neighbor test function
  switch (type)
  {
    case Moose::Element:
    {
      Real test_grad_n = _grad_test[_i][_qp] * _normals[_qp];
      return penalty * u_jump * test_grad_n * _JxW[_qp];
    }
    case Moose::Neighbor:
    {
      Real test_grad_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      return -penalty * u_jump * test_grad_n * _JxW[_qp];
    }
    default:
      mooseError("Unknown DGResidualType in GhostPenaltyKernel");
  }
}

Real
GhostPenaltyKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  // Compute penalty coefficient
  Real penalty = computePenaltyCoeff();

  // Compute normal component of shape function gradient on element side
  Real phi_grad_n = _grad_phi[_j][_qp] * _normals[_qp];

  Real jac = 0.0;

  switch (type)
  {
    case Moose::ElementElement:
    {
      Real test_grad_n = _grad_test[_i][_qp] * _normals[_qp];
      jac = penalty * phi_grad_n * test_grad_n * _JxW[_qp];
      break;
    }

    case Moose::ElementNeighbor:
    {
      // ∂r⁺/∂u⁻ = -penalty * (∇φ⁻·n) * (∇ψ⁺·n)
      Real phi_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      Real test_n = _grad_test[_i][_qp] * _normals[_qp];
      jac = -penalty * phi_n * test_n * _JxW[_qp];
      break;
    }

    case Moose::NeighborElement:
    {
      // ∂r⁻/∂u⁺ = -penalty * (∇φ⁺·n) * (∇ψ⁻·n)
      Real test_neighbor_grad_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      jac = -penalty * phi_grad_n * test_neighbor_grad_n * _JxW[_qp];
      break;
    }

    case Moose::NeighborNeighbor:
    {
      // ∂r⁻/∂u⁻ = +penalty * (∇φ⁻·n) * (∇ψ⁻·n)
      Real phi_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      Real test_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      jac = penalty * phi_n * test_n * _JxW[_qp];
      break;
    }

    default:
      mooseError("Unknown Jacobian type in GhostPenaltyKernel");
  }

  return jac;
}

Real
GhostPenaltyKernel::computePenaltyCoeff()
{
  // Compute characteristic face length
  // h = sqrt(|F|) where |F| is the area of the face
  // Approximate from element volume: h_F ≈ (elem_volume)^{1/d}
  unsigned int dim = _current_elem->dim();
  Real elem_volume = _current_elem->volume();
  Real h_F = std::pow(elem_volume, 1.0 / static_cast<Real>(dim));

  // Compute h_F^{2(k-1+gamma)}
  Real h_F_power = std::pow(h_F, _exponent);

  // Penalty coefficient: c_F * h_F^{exponent}
  Real penalty = _c_F * h_F_power;

  return penalty;
}
