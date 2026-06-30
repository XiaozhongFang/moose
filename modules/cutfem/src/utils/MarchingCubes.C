#include "MarchingCubes.h"
#include <cmath>
#include <algorithm>

// Case table for 2D quadrilateral Marching Squares.
// Bitmask: bit 0=v0, bit 1=v1, bit 2=v2, bit 3=v3
// Vertex indices 0-3, edge intersection indices: 4(v0-v1), 5(v1-v2), 6(v2-v3), 7(v3-v0)

const int MarchingCubes2D::_num_tris[16] = {
    0, // case  0 (0000): all outside
    1, // case  1 (0001): v0 inside
    1, // case  2 (0010): v1 inside
    2, // case  3 (0011): v0,v1 inside
    1, // case  4 (0100): v2 inside
    2, // case  5 (0101): v0,v2 inside (opposite corners)
    2, // case  6 (0110): v1,v2 inside
    3, // case  7 (0111): v0,v1,v2 inside (v3 outside)
    1, // case  8 (1000): v3 inside
    2, // case  9 (1001): v0,v3 inside
    2, // case 10 (1010): v1,v3 inside (opposite corners)
    3, // case 11 (1011): v0,v1,v3 inside (v2 outside)
    2, // case 12 (1100): v2,v3 inside
    3, // case 13 (1101): v0,v2,v3 inside (v1 outside)
    3, // case 14 (1110): v1,v2,v3 inside (v0 outside)
    2, // case 15 (1111): all inside
};

const int MarchingCubes2D::_case_table[16][3][3] = {
    // Case 0 (0000): all outside
    {{0,0,0},{0,0,0},{0,0,0}},
    // Case 1 (0001): v0 inside — tri [v0, p01, p30]
    {{0,4,7},{0,0,0},{0,0,0}},
    // Case 2 (0010): v1 inside — tri [v1, p12, p01]
    {{1,5,4},{0,0,0},{0,0,0}},
    // Case 3 (0011): v0,v1 inside — tri [v0,v1,p12], [v0,p12,p30]
    {{0,1,5},{0,5,7},{0,0,0}},
    // Case 4 (0100): v2 inside — tri [v2, p23, p12]
    {{2,6,5},{0,0,0},{0,0,0}},
    // Case 5 (0101): v0,v2 opposite — tris [v0,p01,p30], [v2,p12,p23]
    {{0,4,7},{2,5,6},{0,0,0}},
    // Case 6 (0110): v1,v2 inside — tris [v1,v2,p23], [v1,p23,p01]
    {{1,2,6},{1,6,4},{0,0,0}},
    // Case 7 (0111): v0,v1,v2 inside — fan from v0: [v0,v1,v2],[v0,v2,p23],[v0,p23,p30]
    {{0,1,2},{0,2,6},{0,6,7}},
    // Case 8 (1000): v3 inside — tri [v3, p30, p23]
    {{3,7,6},{0,0,0},{0,0,0}},
    // Case 9 (1001): v0,v3 inside — tris [v0,p01,p23], [v0,p23,v3]
    {{0,4,6},{0,6,3},{0,0,0}},
    // Case 10 (1010): v1,v3 opposite — tris [v1,p01,p12], [v3,p23,p30]
    {{1,4,5},{3,6,7},{0,0,0}},
    // Case 11 (1011): v0,v1,v3 inside — fan: [v0,v1,p12],[v0,p12,p23],[v0,p23,v3]
    {{0,1,5},{0,5,6},{0,6,3}},
    // Case 12 (1100): v2,v3 inside — tris [v2,v3,p30], [v2,p30,p12]
    {{2,3,7},{2,7,5},{0,0,0}},
    // Case 13 (1101): v0,v2,v3 inside — fan from v0: [v0,p01,p12],[v0,p12,v2],[v0,v2,v3]
    {{0,4,5},{0,5,2},{0,2,3}},
    // Case 14 (1110): v1,v2,v3 inside — fan from v1: [v1,v2,v3],[v1,v3,p30],[v1,p30,p01]
    {{1,2,3},{1,3,7},{1,7,4}},
    // Case 15 (1111): all inside — quad as 2 tris [v0,v1,v3], [v1,v2,v3]
    {{0,1,3},{1,2,3},{0,0,0}},
};

libMesh::Point
MarchingCubes2D::edgeIntersection(const libMesh::Point & va,
                                  const libMesh::Point & vb,
                                  libMesh::Real phi_a,
                                  libMesh::Real phi_b)
{
  libMesh::Real lambda = phi_a / (phi_a - phi_b);
  return va + lambda * (vb - va);
}

libMesh::Real
MarchingCubes2D::totalArea(
    const std::vector<std::vector<libMesh::Point>> & triangles)
{
  libMesh::Real area = 0.0;
  for (const auto & tri : triangles)
  {
    if (tri.size() < 3)
      continue;
    libMesh::Real v1x = tri[1](0) - tri[0](0);
    libMesh::Real v1y = tri[1](1) - tri[0](1);
    libMesh::Real v2x = tri[2](0) - tri[0](0);
    libMesh::Real v2y = tri[2](1) - tri[0](1);
    area += 0.5 * std::abs(v1x * v2y - v1y * v2x);
  }
  return area;
}

std::vector<std::vector<libMesh::Point>>
MarchingCubes2D::triangulateCutElement(
    const libMesh::Elem * elem,
    const std::vector<libMesh::Real> & phi_nodes)
{
  std::vector<std::vector<libMesh::Point>> triangles;

  if (!elem || elem->n_vertices() != 4 || phi_nodes.size() < 4)
    return triangles;

  // Build bitmask: bit i = 1 if vertex i is inside (phi < 0)
  int case_idx = 0;
  for (unsigned int i = 0; i < 4; ++i)
    if (phi_nodes[i] < 0.0)
      case_idx |= (1 << i);

  int num_tris = _num_tris[case_idx];

  // Collect all 8 potential points: 4 vertices + 4 edge intersections
  libMesh::Point pts[8];
  for (unsigned int i = 0; i < 4; ++i)
    pts[i] = elem->point(i);

  if ((phi_nodes[0] < 0) != (phi_nodes[1] < 0))
    pts[4] = edgeIntersection(pts[0], pts[1], phi_nodes[0], phi_nodes[1]);
  if ((phi_nodes[1] < 0) != (phi_nodes[2] < 0))
    pts[5] = edgeIntersection(pts[1], pts[2], phi_nodes[1], phi_nodes[2]);
  if ((phi_nodes[2] < 0) != (phi_nodes[3] < 0))
    pts[6] = edgeIntersection(pts[2], pts[3], phi_nodes[2], phi_nodes[3]);
  if ((phi_nodes[3] < 0) != (phi_nodes[0] < 0))
    pts[7] = edgeIntersection(pts[3], pts[0], phi_nodes[3], phi_nodes[0]);

  // Build triangles from lookup table
  for (int t = 0; t < num_tris; ++t)
  {
    int i0 = _case_table[case_idx][t][0];
    int i1 = _case_table[case_idx][t][1];
    int i2 = _case_table[case_idx][t][2];
    triangles.push_back({pts[i0], pts[i1], pts[i2]});
  }

  return triangles;
}

void
MarchingCubes2D::triangleGaussQuadrature(
    const std::vector<std::vector<libMesh::Point>> & triangles,
    unsigned int order,
    std::vector<libMesh::Point> & points,
    std::vector<libMesh::Real> & weights)
{
  points.clear();
  weights.clear();

  // Barycentric coordinates and weights for triangle Gauss quadrature.
  // Format: each entry is {ξ1, ξ2, ξ3, w_barycentric}
  // Physical point = ξ1*a + ξ2*b + ξ3*c, weight = area * w_barycentric
  struct BaryPoint { libMesh::Real xi1, xi2, xi3, w; };

  std::vector<BaryPoint> ref_pts;

  switch (order)
  {
    case 1:
      ref_pts = {{1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0}};
      break;
    case 2:
      ref_pts = {{0.5, 0.5, 0.0, 1.0/3.0},
                 {0.5, 0.0, 0.5, 1.0/3.0},
                 {0.0, 0.5, 0.5, 1.0/3.0}};
      break;
    case 3:
      ref_pts = {{1.0/3.0, 1.0/3.0, 1.0/3.0, -27.0/48.0},
                 {0.6, 0.2, 0.2, 25.0/48.0},
                 {0.2, 0.6, 0.2, 25.0/48.0},
                 {0.2, 0.2, 0.6, 25.0/48.0}};
      break;
    case 4:
      ref_pts = {{0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322},
                 {0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322},
                 {0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322},
                 {0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011},
                 {0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011},
                 {0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011}};
      break;
    default:
      // Fall back to order 1
      ref_pts = {{1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0}};
      break;
  }

  for (const auto & tri : triangles)
  {
    if (tri.size() < 3)
      continue;

    const auto & a = tri[0];
    const auto & b = tri[1];
    const auto & c = tri[2];

    // Compute triangle area: 0.5 * |(b-a) × (c-a)|
    libMesh::Real v1x = b(0) - a(0), v1y = b(1) - a(1);
    libMesh::Real v2x = c(0) - a(0), v2y = c(1) - a(1);
    libMesh::Real area = 0.5 * std::abs(v1x * v2y - v1y * v2x);

    for (const auto & rp : ref_pts)
    {
      // Physical point = ξ1*a + ξ2*b + ξ3*c
      libMesh::Point p(rp.xi1 * a(0) + rp.xi2 * b(0) + rp.xi3 * c(0),
                       rp.xi1 * a(1) + rp.xi2 * b(1) + rp.xi3 * c(1),
                       0.0);
      points.push_back(p);
      weights.push_back(area * rp.w);
    }
  }
}
