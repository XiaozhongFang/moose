#include "CutFEMCombinedKernel.h"
#include "CutCellQuadratureUserObject.h"

registerMooseObject("CutFEMApp", CutFEMCombinedKernel);

InputParameters
CutFEMCombinedKernel::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<UserObjectName>(
      "cut_cell_quadrature",
      "CutCellQuadratureUserObject providing cut cell quadrature data");

  params.addParam<Real>("source", 0.0,
                        "Body force source term value (constant)");

  params.addClassDescription(
      "Combined diffusion + source kernel with cut cell integration.\n\n"
      "Solves -∇·(∇u) = source on a LevelSet-defined domain.\n"
      "Cut elements integrate over φ<0 using sub-triangle Gauss quadrature.\n"
      "Uncut elements use standard quadrature.\n\n"
      "Single-kernel design avoids multi-kernel quadrature conflicts.");

  return params;
}

CutFEMCombinedKernel::CutFEMCombinedKernel(const InputParameters & params)
  : Kernel(params),
    _cut_quad(getUserObject<CutCellQuadratureUserObject>("cut_cell_quadrature")),
    _source_value(getParam<Real>("source"))
{
}

void
CutFEMCombinedKernel::computeResidual()
{
  dof_id_type eid = _current_elem->id();

  if (!_cut_quad.isCutElement(eid))
  {
    // Standard integration for uncut elements
    // _source_value is 0 → only diffusion term contributes
    Kernel::computeResidual();
    return;
  }

  const auto & qpts = _cut_quad.getQuadraturePoints(eid);
  const auto & qwts = _cut_quad.getQuadratureWeights(eid);
  if (qpts.empty())
    return;

  _assembly.reinitAtPhysical(_current_elem, qpts);
  _var.prepare();
  _assembly.prepare();

  for (_qp = 0; _qp < static_cast<int>(qpts.size()); ++_qp)
  {
    for (_i = 0; _i < _test.size(); ++_i)
    {
      // Diffusion term: ∇φ_i · ∇u
      _local_re(_i) += (_grad_test[_i][_qp] * _grad_u[_qp]) * qwts[_qp];
      // Source term: φ_i * source (if non-zero)
      if (_source_value != 0.0)
        _local_re(_i) -= _test[_i][_qp] * _source_value * qwts[_qp];
    }
  }

  _fe_problem.reinitElem(_current_elem, _tid);
}

void
CutFEMCombinedKernel::computeJacobian()
{
  dof_id_type eid = _current_elem->id();

  if (!_cut_quad.isCutElement(eid))
  {
    Kernel::computeJacobian();
    return;
  }

  const auto & qpts = _cut_quad.getQuadraturePoints(eid);
  const auto & qwts = _cut_quad.getQuadratureWeights(eid);
  if (qpts.empty())
    return;

  _assembly.reinitAtPhysical(_current_elem, qpts);
  _var.prepare();
  _assembly.prepare();

  for (_qp = 0; _qp < static_cast<int>(qpts.size()); ++_qp)
  {
    for (_i = 0; _i < _test.size(); ++_i)
    {
      for (_j = 0; _j < _phi.size(); ++_j)
      {
        // Diffusion Jacobian: ∇φ_i · ∇φ_j
        _local_ke(_i, _j) += (_grad_test[_i][_qp] * _grad_phi[_j][_qp]) * qwts[_qp];
        // Source Jacobian: φ_i * φ_j (if source is non-zero and depends on u)
        // For constant source, no Jacobian contribution
      }
    }
  }

  _fe_problem.reinitElem(_current_elem, _tid);
}

Real
CutFEMCombinedKernel::computeQpResidual()
{
  Real r = _grad_test[_i][_qp] * _grad_u[_qp];
  if (_source_value != 0.0)
    r -= _test[_i][_qp] * _source_value;
  return r;
}

Real
CutFEMCombinedKernel::computeQpJacobian()
{
  return _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
