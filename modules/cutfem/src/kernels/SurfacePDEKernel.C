#include "SurfacePDEKernel.h"

registerMooseObject("CutFEMApp", SurfacePDEKernel);

InputParameters
SurfacePDEKernel::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription(
      "Laplace-Beltrami surface PDE kernel.\n\n"
      "Implements the tangential gradient weak form:\n"
      "  ∫_Γ (P ∇u) · (P ∇v) dS\n"
      "where P = I - n⊗n projects onto the tangent space.\n\n"
      "References:\n"
      "  - Olshanskii, Reusken, Grande (2009)\n"
      "  - Larson & Zahedi (2020)");
  return params;
}

SurfacePDEKernel::SurfacePDEKernel(const InputParameters & parameters)
  : InterfaceKernel(parameters)
{
}

Real
SurfacePDEKernel::computeQpResidual(Moose::DGResidualType type)
{
  // Tangential gradient of solution:
  // P ∇u = ∇u - (n·∇u) n
  Real u_grad_dot_n = _grad_u[_qp] * _normals[_qp];
  RealVectorValue grad_u_tang = _grad_u[_qp] - u_grad_dot_n * _normals[_qp];

  switch (type)
  {
    case Moose::Element:
    {
      // Tangential gradient of test function on element side
      Real test_grad_dot_n = _grad_test[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test[_i][_qp] - test_grad_dot_n * _normals[_qp];

      return (grad_u_tang * grad_test_tang) * _JxW[_qp];
    }
    case Moose::Neighbor:
    {
      // Tangential gradient of test function on neighbor side
      Real test_grad_dot_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test_neighbor[_i][_qp] -
                                       test_grad_dot_n * _normals[_qp];

      return (grad_u_tang * grad_test_tang) * _JxW[_qp];
    }
    default:
      mooseError("Unknown DGResidualType in SurfacePDEKernel");
  }
}

Real
SurfacePDEKernel::computeQpJacobian(Moose::DGJacobianType type)
{
  switch (type)
  {
    case Moose::ElementElement:
    {
      Real phi_grad_dot_n = _grad_phi[_j][_qp] * _normals[_qp];
      RealVectorValue grad_phi_tang = _grad_phi[_j][_qp] - phi_grad_dot_n * _normals[_qp];

      Real test_grad_dot_n = _grad_test[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test[_i][_qp] - test_grad_dot_n * _normals[_qp];

      return (grad_phi_tang * grad_test_tang) * _JxW[_qp];
    }
    case Moose::ElementNeighbor:
    {
      Real phi_grad_dot_n = _grad_phi[_j][_qp] * _normals[_qp];
      RealVectorValue grad_phi_tang = _grad_phi[_j][_qp] - phi_grad_dot_n * _normals[_qp];

      // Neighbor test function uses the same n (tangential projection is sign-independent)
      Real test_grad_dot_n = _grad_test[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test[_i][_qp] - test_grad_dot_n * _normals[_qp];

      return (grad_phi_tang * grad_test_tang) * _JxW[_qp];
    }
    case Moose::NeighborElement:
    {
      Real phi_grad_dot_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      RealVectorValue grad_phi_tang = _grad_phi_neighbor[_j][_qp] -
                                      phi_grad_dot_n * _normals[_qp];

      Real test_grad_dot_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test_neighbor[_i][_qp] -
                                       test_grad_dot_n * _normals[_qp];

      return (grad_phi_tang * grad_test_tang) * _JxW[_qp];
    }
    case Moose::NeighborNeighbor:
    {
      Real phi_grad_dot_n = _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      RealVectorValue grad_phi_tang = _grad_phi_neighbor[_j][_qp] -
                                      phi_grad_dot_n * _normals[_qp];

      Real test_grad_dot_n = _grad_test_neighbor[_i][_qp] * _normals[_qp];
      RealVectorValue grad_test_tang = _grad_test_neighbor[_i][_qp] -
                                       test_grad_dot_n * _normals[_qp];

      return (grad_phi_tang * grad_test_tang) * _JxW[_qp];
    }
    default:
      mooseError("Unknown DGJacobianType in SurfacePDEKernel");
  }
}
