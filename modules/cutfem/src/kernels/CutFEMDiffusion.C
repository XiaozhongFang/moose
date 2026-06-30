#include "CutFEMDiffusion.h"
#include "CutCellQuadratureUserObject.h"

registerMooseObject("CutFEMApp", CutFEMDiffusion);

InputParameters
CutFEMDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addRequiredParam<UserObjectName>(
      "cut_cell_quadrature",
      "Name of the CutCellQuadratureUserObject providing cut cell quadrature data");

  params.addClassDescription(
      "Diffusion kernel with cut cell integration for CutFEM.\n\n"
      "On cut elements, uses _fe_problem.reinitElemPhys() to evaluate\n"
      "FE data at sub-triangle Gauss quadrature points. On uncut elements,\n"
      "uses standard element quadrature.\n\n"
      "Requires a CutCellQuadratureUserObject with a LevelSet function\n"
      "defining the physical domain (φ < 0).");

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
  dof_id_type elem_id = _current_elem->id();

  if (!_cut_quad.isCutElement(elem_id))
  {
    // Uncut element: standard integration
    Kernel::computeResidual();
    return;
  }

  // Cut element: use sub-triangle Gauss quadrature
  const auto & qpts = _cut_quad.getQuadraturePoints(elem_id);
  const auto & qwts = _cut_quad.getQuadratureWeights(elem_id);

  if (qpts.empty())
    return;

  // Reinit FE at custom physical quadrature points
  // This updates _grad_u, _grad_test, _grad_phi, _q_point, _JxW, etc.
  _fe_problem.reinitElemPhys(_current_elem, qpts, _tid);

  // Now use standard qp loop — FE data is at custom points
  for (_qp = 0; _qp < (int)qpts.size(); ++_qp)
  {
    // _grad_u[_qp] is now at the custom quadrature point qpts[_qp]
    for (_i = 0; _i < _test.size(); ++_i)
      _local_re(_i) += computeQpResidual();
  }

  // Restore standard quadrature for subsequent operations
  _fe_problem.reinitElem(_current_elem, _tid);
}

void
CutFEMDiffusion::computeJacobian()
{
  dof_id_type elem_id = _current_elem->id();

  if (!_cut_quad.isCutElement(elem_id))
  {
    Kernel::computeJacobian();
    return;
  }

  const auto & qpts = _cut_quad.getQuadraturePoints(elem_id);
  const auto & qwts = _cut_quad.getQuadratureWeights(elem_id);

  if (qpts.empty())
    return;

  _fe_problem.reinitElemPhys(_current_elem, qpts, _tid);

  for (_qp = 0; _qp < (int)qpts.size(); ++_qp)
    for (_i = 0; _i < _test.size(); ++_i)
      for (_j = 0; _j < _phi.size(); ++_j)
        _local_ke(_i, _j) += computeQpJacobian();

  _fe_problem.reinitElem(_current_elem, _tid);
}

Real
CutFEMDiffusion::computeQpResidual()
{
  // Use custom quadrature weight (not _JxW which is 1.0 after reinitElemPhys)
  dof_id_type eid = _current_elem->id();
  Real w = _cut_quad.isCutElement(eid)
               ? _cut_quad.getQuadratureWeights(eid)[_qp]
               : _JxW[_qp];
  return _grad_test[_i][_qp] * _grad_u[_qp] * w;
}

Real
CutFEMDiffusion::computeQpJacobian()
{
  dof_id_type eid = _current_elem->id();
  Real w = _cut_quad.isCutElement(eid)
               ? _cut_quad.getQuadratureWeights(eid)[_qp]
               : _JxW[_qp];
  return _grad_test[_i][_qp] * _grad_phi[_j][_qp] * w;
}
