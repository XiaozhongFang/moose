#pragma once

#include "libmesh/point.h"
#include "libmesh/elem.h"
#include <vector>

/**
 * Marching Cubes algorithm for 2D quadrilateral elements.
 *
 * Given a level set function evaluated at element vertices, generates
 * sub-triangles covering the region where φ < 0 (inside the domain).
 *
 * Edge and vertex numbering (CCW):
 *
 *   v3 ---e2--- v2
 *   |           |
 *   e3          e1
 *   |           |
 *   v0 ---e0--- v1
 *
 * Sub-triangle vertex indices:
 *   0-3 = original element vertices (v0..v3)
 *   4   = intersection on edge 0 (v0-v1)
 *   5   = intersection on edge 1 (v1-v2)
 *   6   = intersection on edge 2 (v2-v3)
 *   7   = intersection on edge 3 (v3-v0)
 */
class MarchingCubes2D
{
public:
  /**
   * Compute sub-triangles for a cut quadrilateral element.
   */
  static std::vector<std::vector<libMesh::Point>>
  triangulateCutElement(const libMesh::Elem * elem,
                        const std::vector<libMesh::Real> & phi_nodes);

  /**
   * Compute the total area of sub-triangles.
   */
  static libMesh::Real
  totalArea(const std::vector<std::vector<libMesh::Point>> & triangles);

  /**
   * Generate Gauss quadrature points and weights for a set of sub-triangles.
   *
   * @param triangles  Sub-triangles from triangulateCutElement
   * @param order      Gauss quadrature order (1-4)
   * @param points     Output: quadrature points in physical coordinates
   * @param weights    Output: quadrature weights
   */
  static void
  triangleGaussQuadrature(const std::vector<std::vector<libMesh::Point>> & triangles,
                          unsigned int order,
                          std::vector<libMesh::Point> & points,
                          std::vector<libMesh::Real> & weights);

  /**
   * Get the interface segments (line segments on φ=0) within a cut element.
   * For a quadrilateral, the interface is the segment(s) connecting edge intersections.
   *
   * @param elem       The element
   * @param phi_nodes  Level set values at vertices
   * @return           List of interface segments, each = {p_start, p_end}
   */
  static std::vector<std::pair<libMesh::Point, libMesh::Point>>
  interfaceSegments(const libMesh::Elem * elem,
                    const std::vector<libMesh::Real> & phi_nodes);

  /**
   * Generate 1D Gauss quadrature along interface segments.
   *
   * @param segments   Interface segments from interfaceSegments()
   * @param order      1D Gauss quadrature order (1-4)
   * @param points     Output: quadrature points in physical coordinates
   * @param weights    Output: quadrature weights (includes segment length Jacobian)
   */
  static void
  interfaceGaussQuadrature(
      const std::vector<std::pair<libMesh::Point, libMesh::Point>> & segments,
      unsigned int order,
      std::vector<libMesh::Point> & points,
      std::vector<libMesh::Real> & weights);

private:
  static libMesh::Point
  edgeIntersection(const libMesh::Point & va, const libMesh::Point & vb,
                   libMesh::Real phi_a, libMesh::Real phi_b);

  /// Lookup table: for each of the 16 cases, up to 3 triangles, each with 3 vertex indices
  static const int _case_table[16][3][3];
  static const int _num_tris[16];
};
