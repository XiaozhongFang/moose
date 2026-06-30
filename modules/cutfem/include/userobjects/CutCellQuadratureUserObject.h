#pragma once

#include "GeneralUserObject.h"
#include "MooseMesh.h"
#include "Function.h"

/**
 * Cut Cell Quadrature UserObject for CutFEM (Phase 2)
 *
 * Detects elements cut by a Level Set interface by evaluating the level set
 * function at element nodes and checking for sign changes.
 *
 * Phase 2 status:
 *   - Cut detection: implemented
 *   - Marching Cubes: stub (WP5)
 *   - Sub-element quadrature: stub (WP6)
 *
 * Usage:
 * \code
 * [UserObjects]
 *   [cut_cell_quad]
 *     type = CutCellQuadratureUserObject
 *     level_set_function = level_set_func
 *     quadrature_order = 4
 *     execute_on = 'INITIAL'
 *   []
 * []
 * \endcode
 */
class CutCellQuadratureUserObject : public GeneralUserObject
{
public:
  static InputParameters validParams();

  CutCellQuadratureUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & y) override;
  virtual void finalize() override;

  /// Check if element is cut by interface (by element ID)
  bool isCutElement(dof_id_type elem_id) const;

  /// Get quadrature points for cut element
  const std::vector<Point> & getQuadraturePoints(dof_id_type elem_id) const;

  /// Get quadrature weights for cut element
  const std::vector<Real> & getQuadratureWeights(dof_id_type elem_id) const;

  /// Element classification
  enum class ElementType : int
  {
    NO_INTERSECTION = 0,
    CUT_BY_INTERFACE = 2
  };

  /// Get element classification
  ElementType getElementType(dof_id_type elem_id) const;

  /// Total cut elements found
  unsigned int getTotalCutElements() const;

protected:
  MooseMesh & _mesh;

  unsigned int _quad_order;
  Real _level_set_tolerance;

  std::map<dof_id_type, std::vector<Point>> _quad_points_map;
  std::map<dof_id_type, std::vector<Real>> _quad_weights_map;
  std::set<dof_id_type> _cut_elements;
  std::map<dof_id_type, ElementType> _element_type_map;
  unsigned int _total_cut_elements = 0;

  const Function & _level_set_func;
};
