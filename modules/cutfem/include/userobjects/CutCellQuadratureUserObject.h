#pragma once

#include "GeneralUserObject.h"
#include "MooseMesh.h"

/**
 * Cut Cell Quadrature UserObject for CutFEM (Phase 2)
 *
 * Automatically detects elements cut by a Level Set interface and computes
 * proper quadrature points and weights for integration on cut cell portions.
 *
 * This is the core component enabling CutFEM's ability to handle integration
 * on non-conforming cuts. Phase 1 uses Ghost Penalty on faces. Phase 2 will
 * use this UO for volume integration on cut elements.
 *
 * Algorithm overview:
 * 1. For each element, evaluate the Level Set function at nodes
 * 2. Detect cut: if sign changes across nodes → element intersected
 * 3. Find edge-intersection points using bisection/linear interpolation
 * 4. Subdivide element into sub-elements using Marching Cubes
 * 5. Apply standard Gauss quadrature to each sub-element
 * 6. Cache results for efficient access during assembly
 *
 * References:
 * - Saye (2015, 2022): High-order quadrature on level sets
 * - ngsxfem library: Reference implementation of cut cell algorithms
 * - Burman et al. (2018): Isoparametric mapping strategy
 *
 * Usage:
 * \code
 * [UserObjects]
 *   [cut_cell_quadrature]
 *     type = CutCellQuadratureUserObject
 *     level_set_variable = phi
 *     quadrature_order = 4           # Gauss quadrature order (1-8)
 *     level_set_tolerance = 1e-12    # Tolerance for detecting cuts
 *     execute_on = 'INITIAL LINEAR'  # Update during solve
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

  // Public query interface
  
  /**
   * Check if element is intersected by level set interface (by element ID)
   *
   * @param elem_id: element ID
   * @return true if element is cut by interface
   */
  bool isCutElement(dof_id_type elem_id) const;

  /**
   * Get quadrature points for cut element
   *
   * @param elem_id: element ID
   * @return reference to vector of quadrature points
   */
  const std::vector<Point> & getQuadraturePoints(dof_id_type elem_id) const;

  /**
   * Get quadrature weights for cut element
   *
   * @param elem_id: element ID
   * @return reference to vector of quadrature weights
   */
  const std::vector<Real> & getQuadratureWeights(dof_id_type elem_id) const;

  /// Element classification enum
  enum class ElementType : int
  {
    NO_INTERSECTION = 0,    ///< Element not cut by interface
    FULL_ELEMENT = 1,       ///< Element fully inside or outside
    CUT_BY_INTERFACE = 2    ///< Element intersected by interface
  };

  /**
   * Get element classification
   *
   * @param elem_id: element ID
   * @return ElementType classification
   */
  ElementType getElementType(dof_id_type elem_id) const;

  /// Get total number of cut elements detected
  unsigned int getTotalCutElements() const;

  /// Get total number of sub-elements generated
  unsigned int getTotalSubElements() const;

protected:
  /// Reference to mesh
  MooseMesh & _mesh;

  /// Level set variable name for interface representation: φ(x)
  AuxVariableName _var_name;

  /// Quadrature rule order (1-8)
  unsigned int _quad_order;

  /// Tolerance for level set evaluation
  Real _level_set_tolerance;

  /// Virtual interface width for stabilization
  Real _interface_width;

  /// Storage: map from element ID to quadrature points
  std::map<dof_id_type, std::vector<Point>> _quad_points_map;

  /// Storage: map from element ID to quadrature weights
  std::map<dof_id_type, std::vector<Real>> _quad_weights_map;

  /// Set of cut element IDs
  std::set<dof_id_type> _cut_elements;

  /// Element type classification map
  std::map<dof_id_type, ElementType> _element_type_map;

  /// Counters for statistics
  unsigned int _total_cut_elements = 0;
  unsigned int _total_sub_elements = 0;

  // Private implementation methods

  /**
   * Marching Cubes algorithm: subdivide element cut by interface (Phase 2)
   *
   * Takes an element and the level set values at its nodes, performs
   * sub-element subdivision to separate positive and negative regions.
   *
   * @param elem: element to subdivide
   * @param phi_at_nodes: level set values at element nodes
   * @return vector of sub-element quadrature points
   *
   * Reference: Saye (2015) quadrature algorithm
   */
  std::vector<Point> marchingCubes(const Elem * elem,
                                    const std::vector<Real> & phi_at_nodes);

  /**
   * Get quadrature weights for given points
   *
   * @param elem_id: element ID
   * @param quad_points: vector of quadrature points
   * @return vector of corresponding weights
   */
  std::vector<Real> getQuadratureWeights(dof_id_type elem_id,
                                         const std::vector<Point> & quad_points) const;
};
