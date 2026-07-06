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
#include "BottomUpJIntegrator.h"

#include <memory>
#include <map>
#include <string>
#include <vector>

/**
 * Bottom-Up photolysis engine — cross-section × quantum-yield × lamp-flux
 * integration following the F0AM IntegrateJ.m / J_BottomUp.m algorithm.
 *
 * Uses BottomUpJIntegrator to compute J-values from spectral data files.
 * This scheme depends on temperature and pressure (from PhysParams) for
 * cross-section interpolation and temperature-dependent quantum yields.
 *
 * Required PhysParams fields:
 *   - jfac       JFAC scaling factor
 *   - temperature  [K]
 *   - pressure     [mbar] (0 = use default air density)
 */
class BottomUpPhoto : public PhotolysisEngine
{
public:
  /**
   * Construct from a data directory and flux file.
   *
   * @param data_dir  Path to the data directory containing CrossSections/
   *                  and bottomup_jmap.dat.
   * @param flux_file Lamp/actinic flux filename (relative to data_dir).
   * @param j_numbers 1-based J<N> numbers that this engine will compute.
   */
  BottomUpPhoto(const std::string & data_dir,
                const std::string & flux_file,
                const std::vector<unsigned int> & j_numbers);

  /**
   * Compute photolysis rates via spectral integration.
   *
   * Uses cached results from BottomUpJIntegrator — the integrator
   * re-computes only when temperature or pressure changes (detected
   * via internal caching in BottomUpJIntegrator).
   *
   * @param cos_sza  Not used by BottomUp scheme (pass any value)
   * @param sec_sza  Not used by BottomUp scheme (pass any value)
   * @param params   Physical parameters (jfac, temperature, pressure)
   * @param j_out    Output vector, resized to nJ()
   */
  void computeJ(Real cos_sza,
                Real sec_sza,
                const PhysParams & params,
                std::vector<Real> & j_out) const override;

  /// Number of photolysis rates.
  unsigned int nJ() const override { return _j_numbers.size(); }

private:
  std::unique_ptr<BottomUpJIntegrator> _integrator;
  std::vector<unsigned int> _j_numbers;

  /// Cache for computed J-values keyed by J name.
  mutable std::map<std::string, Real> _cached_j;
  mutable Real _cached_T;
  mutable Real _cached_P;
  mutable bool _cache_valid;
};
