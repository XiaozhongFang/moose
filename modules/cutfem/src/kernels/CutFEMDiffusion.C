#include "CutFEMDiffusion.h"
#include "CutCellQuadratureUserObject.h"

registerMooseObject("CutFEMApp", CutFEMDiffusion);

InputParameters
CutFEMDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<UserObjectName>(
      "cut_cell_quadrature",
      "CutCellQuadratureUserObject providing cut cell quadrature data");

  params.addClassDescription(
      "Diffusion kernel with cut cell integration for CutFEM.\n\n"
      "On elements cut by a LevelSet interface, reinitializes FE data\n"
      "at sub-triangle Gauss quadrature points via Assembly::reinitAtPhysical.\n"
      "On uncut elements, uses standard quadrature.\n\n"
      "Note: Only one CutFEM kernel per variable is supported. For problems\n"
      "with multiple terms, use CutFEMCombinedKernel or add terms directly.");

  return params;
}

CutFEMDiffusion::CutFEMDiffusion(const InputParameters & params)
  : Kernel(params),
    _cut_quad(getUserObject<CutCellQuadratureUserObject>("cut_cell_quadrature"))
{
}

void
CutFEMDiffusion::computeResidual()
{
  dof_id_type eid = _current_elem->id();
  if (!_cut_quad.isCutElement(eid))
    return Kernel::computeResidual();

  const auto & qpts = _cut_quad.getQuadraturePoints(eid);
  const auto & qwts = _cut_quad.getQuadratureWeights(eid);
  if (qpts.empty())
    return;

  // Reinit FE at custom physical points; weight qpts via qwts, not _JxW
  _assembly.reinitAtPhysical(_current_elem, qpts);
  _var.prepare();
  _assembly.prepare();

  for (_qp = 0; _qp < static_cast<int>(qpts.size()); ++_qp)
    for (_i = 0; _i < _test.size(); ++_i)
      _local_re(_i) += (_grad_test[_i][_qp] * _grad_u[_qp]) * qwts[_qp];

  // Restore standard quadrature for subsequent kernels
  _fe_problem.reinitElem(_current_elem, _tid);
}

void
CutFEMDiffusion::computeJacobian()
{
  dof_id_type eid = _current_elem->id();
  if (!_cut_quad.isCutElement(eid))
    return Kernel::computeJacobian();

  const auto & qpts = _cut_quad.getQuadraturePoints(eid);
  const auto & qwts = _cut_quad.getQuadratureWeights(eid);
  if (qpts.empty())
    return;

  _assembly.reinitAtPhysical(_current_elem, qpts);
  _var.prepare();
  _assembly.prepare();

  for (_qp = 0; _qp < static_cast<int>(qpts.size()); ++_qp)
    for (_i = 0; _i < _test.size(); ++_i)
      for (_j = 0; _j < _phi.size(); ++_j)
        _local_ke(_i, _j) += (_grad_test[_i][_qp] * _grad_phi[_j][_qp]) * qwts[_qp];

  _fe_problem.reinitElem(_current_elem, _tid);
}

Real
CutFEMDiffusion::computeQpResidual()
{
  return _grad_test[_i][_qp] * _grad_u[_qp];
}

Real
CutFEMDiffusion::computeQpJacobian()
{
  return _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}
