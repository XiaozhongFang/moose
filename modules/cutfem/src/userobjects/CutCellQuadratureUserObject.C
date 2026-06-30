#include "CutCellQuadratureUserObject.h"
#include "MarchingCubes.h"

registerMooseObject("CutFEMApp", CutCellQuadratureUserObject);

InputParameters
CutCellQuadratureUserObject::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addRequiredParam<FunctionName>(
      "level_set_function",
      "Level Set function φ(x) defining the domain boundary.\n"
      "Evaluated at element nodes to detect cuts by sign change.");

  params.addParam<unsigned int>(
      "quadrature_order", 4,
      "Gauss quadrature order for integration on cut elements (1-8)");

  params.addParam<Real>(
      "level_set_tolerance", 1e-12,
      "Tolerance for detecting element-interface intersections");

  params.addClassDescription(
      "Detects elements cut by a Level Set interface.\n"
      "Phase 2: cut detection implemented; Marching Cubes and sub-element\n"
      "quadrature are stubs (WP5/WP6).\n\n"
      "References:\n"
      "  - Burman et al. (2015) CutFEM method\n"
      "  - Saye (2015) High-order quadrature on level sets");

  params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_BEGIN};
  return params;
}

CutCellQuadratureUserObject::CutCellQuadratureUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _mesh(_fe_problem.mesh()),
    _quad_order(getParam<unsigned int>("quadrature_order")),
    _level_set_tolerance(getParam<Real>("level_set_tolerance")),
    _level_set_func(getFunction("level_set_function"))
{
  if (_quad_order < 1 || _quad_order > 8)
    paramError("quadrature_order",
               "Must be between 1 and 8. Got: " + std::to_string(_quad_order));
}

void
CutCellQuadratureUserObject::initialize()
{
  _cut_elements.clear();
  _quad_points_map.clear();
  _quad_weights_map.clear();
  _element_type_map.clear();
  _sub_triangles.clear();
  _total_cut_elements = 0;
  _total_sub_elements = 0;
}

void
CutCellQuadratureUserObject::execute()
{
  for (const auto & elem : _mesh.getMesh().active_local_element_ptr_range())
  {
    // Evaluate level set function at element vertices
    unsigned int n_vert = elem->n_vertices();
    std::vector<Real> phi_at_nodes(n_vert);

    for (unsigned int n = 0; n < n_vert; ++n)
    {
      const Node & node = elem->node_ref(n);
      phi_at_nodes[n] = _level_set_func.value(_t, node);
    }

    // Detect cut by sign change
    auto [min_it, max_it] = std::minmax_element(phi_at_nodes.begin(), phi_at_nodes.end());
    Real phi_min = *min_it;
    Real phi_max = *max_it;

    if (phi_min * phi_max < -_level_set_tolerance ||
        std::abs(phi_min) < _level_set_tolerance ||
        std::abs(phi_max) < _level_set_tolerance)
    {
      _cut_elements.insert(elem->id());
      _total_cut_elements++;
      _element_type_map[elem->id()] = ElementType::CUT_BY_INTERFACE;

      // Generate sub-triangles via Marching Cubes
      if (elem->n_vertices() == 4)
      {
        auto triangles = MarchingCubes2D::triangulateCutElement(elem, phi_at_nodes);
        _sub_triangles[elem->id()] = triangles;
        _total_sub_elements += triangles.size();

        // Generate Gauss quadrature points on sub-triangles
        std::vector<Point> qpts;
        std::vector<Real> qwts;
        MarchingCubes2D::triangleGaussQuadrature(triangles, _quad_order, qpts, qwts);
        _quad_points_map[elem->id()] = qpts;
        _quad_weights_map[elem->id()] = qwts;
      }
    }
    else
      _element_type_map[elem->id()] = ElementType::NO_INTERSECTION;
  }
}

void
CutCellQuadratureUserObject::threadJoin(const UserObject & uo)
{
  const auto & other = static_cast<const CutCellQuadratureUserObject &>(uo);
  for (auto id : other._cut_elements)
    _cut_elements.insert(id);
  for (const auto & pair : other._quad_points_map)
    _quad_points_map[pair.first] = pair.second;
  for (const auto & pair : other._quad_weights_map)
    _quad_weights_map[pair.first] = pair.second;
  for (const auto & pair : other._element_type_map)
    _element_type_map[pair.first] = pair.second;
  for (const auto & pair : other._sub_triangles)
    _sub_triangles[pair.first] = pair.second;
  _total_cut_elements += other._total_cut_elements;
  _total_sub_elements += other._total_sub_elements;
}

void
CutCellQuadratureUserObject::finalize()
{
  _mesh.comm().set_union(_cut_elements);

  // Build cut element ID string for test verification
  std::string ids;
  for (auto id : _cut_elements)
    ids += std::to_string(id) + " ";

  // Verify quadrature: weight sum should equal total sub-triangle area
  Real total_weight_sum = 0.0;
  Real total_sub_area = 0.0;
  for (const auto & pair : _sub_triangles)
  {
    total_sub_area += MarchingCubes2D::totalArea(pair.second);
    auto w_it = _quad_weights_map.find(pair.first);
    if (w_it != _quad_weights_map.end())
      for (auto w : w_it->second)
        total_weight_sum += w;
  }
  Real area_diff = std::abs(total_weight_sum - total_sub_area);

  mooseInfo("CutCellQuadratureUserObject: " +
            std::to_string(_total_cut_elements) + " cut, " +
            std::to_string(_total_sub_elements) + " tris, " +
            "area=" + std::to_string(total_sub_area) + ", " +
            "weights=" + std::to_string(total_weight_sum) + ", " +
            "diff=" + std::to_string(area_diff) +
            (ids.empty() ? "" : " IDs: " + ids));
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

const std::vector<Point> &
CutCellQuadratureUserObject::getQuadraturePoints(dof_id_type elem_id) const
{
  auto it = _quad_points_map.find(elem_id);
  if (it != _quad_points_map.end())
    return it->second;
  static std::vector<Point> empty;
  return empty;
}

const std::vector<Real> &
CutCellQuadratureUserObject::getQuadratureWeights(dof_id_type elem_id) const
{
  auto it = _quad_weights_map.find(elem_id);
  if (it != _quad_weights_map.end())
    return it->second;
  static std::vector<Real> empty;
  return empty;
}

unsigned int
CutCellQuadratureUserObject::getTotalCutElements() const
{
  return _total_cut_elements;
}
