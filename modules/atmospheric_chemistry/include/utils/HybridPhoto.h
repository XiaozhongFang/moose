//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PhotolysisEngine.h"
#include "HybridJTableReader.h"

#include <memory>
#include <string>
#include <vector>

/**
 * F0AM Hybrid photolysis engine (4D TUV lookup table).
 *
 * Uses HybridJTableReader to interpolate J-values on a 4D grid
 * of (SZA, albedo, O3 column, altitude).  The table stores
 * log10(J); the result is converted back to linear space and
 * scaled by JFAC from PhysParams.
 *
 * Required PhysParams fields:
 *   - jfac       JFAC scaling factor
 *   - albedo     Surface albedo [0–1]
 *   - o3column   Ozone column [Dobson Units]
 *   - altitude   Altitude [m]
 */
class HybridPhoto : public PhotolysisEngine
{
public:
  /**
   * Construct from a directory containing F0AM TUV tables.
   *
   * @param table_dir  Path to the directory with table_J*.dat files
   *                   (relative to the mechanism root, as in the existing
   *                   hybrid_table_dir convention).
   * @param j_numbers  1-based J<N> numbers that this engine will compute.
   */
  HybridPhoto(const std::string & table_dir,
              const std::vector<unsigned int> & j_numbers);

  /**
   * Compute photolysis rates via 4D table lookup.
   *
   * @param cos_sza  Cosine of SZA (converted to degrees for the table)
   * @param sec_sza  Secant of SZA (not used by Hybrid scheme)
   * @param params   Physical parameters (jfac, albedo, o3column, altitude)
   * @param j_out    Output vector, resized to nJ()
   */
  void computeJ(Real cos_sza,
                Real sec_sza,
                const PhysParams & params,
                std::vector<Real> & j_out) const override;

  /// Number of photolysis rates.
  unsigned int nJ() const override { return _j_numbers.size(); }

private:
  std::unique_ptr<HybridJTableReader> _reader;
  std::vector<unsigned int> _j_numbers;
};
