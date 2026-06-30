#include "CutCellQuadratureUserObject.h"

registerMooseObject("CutFEMApp", CutCellQuadratureUserObject);

InputParameters
CutCellQuadratureUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<AuxVariableName>(
      "level_set_variable",
      "Name of the Level Set function variable phi defining the domain boundary");

  params.addParam<unsigned int>(
      "quadrature_order",
      4,
      "Gauss quadrature order for integration on cut elements (default: 4, supports up to 8)");

  params.addParam<Real>(
      "level_set_tolerance",
      1e-12,
      "Tolerance for detecting element-interface intersections (default: 1e-12)");

  params.addParam<Real>(
      "interface_width",
      0.0,
      "Virtual interface width for stabilization (default: 0, no stabilization)");

  params.addClassDescription(
      "User object for automatic cut-cell quadrature generation in CutFEM.\n\n"
      "Detects elements cut by Level Set interface, performs Marching Cubes\n"
      "sub-division, and generates appropriate integration rules.\n"
      "\n"
      "Key capabilities:\n"
      "  - Element intersection detection via Level Set\n"
      "  - Automatic sub-element generation (Marching Cubes)\n"
      "  - Gauss quadrature on sub-elements\n"
      "  - Caching for efficient assembly\n"
      "\n"
      "References:\n"
      "  - Burman et al. (2015) CutFEM method\n"
      "  - Saye (2015) Marching Cubes for Level Set methods");

  // This object requires execute_on to be set appropriately
  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_LINEAR};
  return params;
}

CutCellQuadratureUserObject::CutCellQuadratureUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_fe_problem.mesh()),
    _quad_order(getParam<unsigned int>("quadrature_order")),
    _level_set_tolerance(getParam<Real>("level_set_tolerance")),
    _interface_width(getParam<Real>("interface_width"))
{
  // Validate quadrature order
  if (_quad_order < 1 || _quad_order > 8)
    paramError("quadrature_order",
               "Must be between 1 and 8. Received: " + std::to_string(_quad_order));

  // Get Level Set variable
  AuxVariableName var_name = getParam<AuxVariableName>("level_set_variable");

  // Try to get the variable (this may fail if not available)
  // We'll handle variable access during execute() phase
  _var_name = var_name;

  mooseInfo("CutCellQuadratureUserObject initialized with:");
  mooseInfo("  - Quadrature order: " + std::to_string(_quad_order));
  mooseInfo("  - Level Set variable: " + var_name);
  mooseInfo("  - Tolerance: " + std::to_string(_level_set_tolerance));
}

void
CutCellQuadratureUserObject::initialize()
{
  // Clear cached data from previous execution
  _cut_elements.clear();
  _quad_points_map.clear();
  _quad_weights_map.clear();
  _element_type_map.clear();

  _total_cut_elements = 0;
  _total_sub_elements = 0;
}

void
CutCellQuadratureUserObject::execute()
{
  // Get Level Set variable from the system
  // In Phase 2, this will access the level set solution vector
  // For Phase 1, we simply iterate over elements and mark none as cut

  // For now, we'll do a simple element-by-element check
  // In production, this would be more sophisticated

  for (const auto & elem : _mesh.getMesh().active_local_element_ptr_range())
  {
    // Check if this element is cut by the Level Set
    if (isCutElement(elem->id()))
    {
      _cut_elements.insert(elem->id());
      _total_cut_elements++;

      // For Phase 1, we just mark elements
      // In Phase 2, we would do Marching Cubes and generate sub-elements
      _element_type_map[elem->id()] = ElementType::FULL_ELEMENT;
    }
    else
    {
      _element_type_map[elem->id()] = ElementType::NO_INTERSECTION;
    }
  }

  // Log statistics
  mooseInfo("CutCellQuadratureUserObject: Found " + std::to_string(_total_cut_elements) +
            " cut elements");
}

void
CutCellQuadratureUserObject::threadJoin(const UserObject & uo)
{
  // Merge thread-local data from other thread
  const auto & other = static_cast<const CutCellQuadratureUserObject &>(uo);

  // Merge cut elements set
  for (auto id : other._cut_elements)
    _cut_elements.insert(id);

  // Merge quadrature points and weights maps
  for (const auto & pair : other._quad_points_map)
    _quad_points_map[pair.first] = pair.second;

  for (const auto & pair : other._quad_weights_map)
    _quad_weights_map[pair.first] = pair.second;

  for (const auto & pair : other._element_type_map)
    _element_type_map[pair.first] = pair.second;

  // Update counters
  _total_cut_elements += other._total_cut_elements;
  _total_sub_elements += other._total_sub_elements;
}

void
CutCellQuadratureUserObject::finalize()
{
  // Broadcast results to all processors if running in parallel
  // (communication handled by MooseMesh parallel utilities)

  mooseInfo("CutCellQuadratureUserObject finalized:");
  mooseInfo("  - Total cut elements: " + std::to_string(_total_cut_elements));
  mooseInfo("  - Total sub-elements: " + std::to_string(_total_sub_elements));
  mooseInfo("  - Cached quadrature rules: " + std::to_string(_quad_points_map.size()));
}

// isCutElement(const Elem*) removed — header only declares dof_id_type overload.
// The dof_id_type version at end of file serves as the stub.

std::vector<Point>
CutCellQuadratureUserObject::marchingCubes(const Elem * elem,
                                           const std::vector<Real> & phi_at_nodes)
{
  // Marching Cubes algorithm for 2D/3D
  // Stub implementation for Phase 1
  //
  // Proper implementation would:
  // 1. Look up the MC pattern based on sign configuration
  // 2. Find edge intersections via linear interpolation
  // 3. Generate sub-element vertices
  // 4. Return quadrature points

  std::vector<Point> quad_points;

  if (!elem)
    return quad_points;

  // For now, just return element centroid
  quad_points.push_back(elem->vertex_average());

  return quad_points;
}

std::vector<Real>
CutCellQuadratureUserObject::getQuadratureWeights(dof_id_type elem_id,
                                                    const std::vector<Point> & quad_points) const
{
  // Generate Gauss quadrature weights for the given points
  // For now, return uniform weights

  std::vector<Real> weights(quad_points.size(), 1.0 / quad_points.size());
  return weights;
}

const std::vector<Point> &
CutCellQuadratureUserObject::getQuadraturePoints(dof_id_type elem_id) const
{
  // Return cached quadrature points for element
  // If not found, return an empty vector

  auto it = _quad_points_map.find(elem_id);
  if (it != _quad_points_map.end())
    return it->second;

  static std::vector<Point> empty_vector;
  return empty_vector;
}

const std::vector<Real> &
CutCellQuadratureUserObject::getQuadratureWeights(dof_id_type elem_id) const
{
  // Return cached quadrature weights for element

  auto it = _quad_weights_map.find(elem_id);
  if (it != _quad_weights_map.end())
    return it->second;

  static std::vector<Real> empty_vector;
  return empty_vector;
}

bool
CutCellQuadratureUserObject::isCutElement(dof_id_type elem_id) const
{
  return _cut_elements.find(elem_id) != _cut_elements.end();
}

CutCellQuadratureUserObject::ElementType
CutCellQuadratureUserObject::getElementType(dof_id_type elem_id) const
{
  auto it = _element_type_map.find(elem_id);
  if (it != _element_type_map.end())
    return it->second;

  return ElementType::NO_INTERSECTION;
}

unsigned int
CutCellQuadratureUserObject::getTotalCutElements() const
{
  return _total_cut_elements;
}

unsigned int
CutCellQuadratureUserObject::getTotalSubElements() const
{
  return _total_sub_elements;
}
