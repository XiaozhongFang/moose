//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"
#include "IMechanism.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Parameters for constructing a PhotolysisEngine via the factory.
 *
 * Bundles all data needed by any of the three photolysis schemes so
 * the factory can dispatch to the correct subclass constructor.
 */
struct PhotolysisParams
{
  // ── MCM_SZA ──
  std::vector<Real> j_CL;                    ///< CL coefficients
  std::vector<Real> j_CMM;                   ///< CMM exponents
  std::vector<Real> j_CNN;                   ///< CNN coefficients
  std::vector<unsigned int> j_numbers;       ///< 1-based J<N> numbers

  // ── HYBRID ──
  std::string hybrid_table_dir;              ///< Directory with F0AM TUV tables
  std::vector<unsigned int> hybrid_j_numbers;///< J numbers for hybrid scheme

  // ── BOTTOMUP ──
  std::string bottomup_data_dir;             ///< Data directory (CrossSections/, etc.)
  std::string bottomup_flux_file;            ///< Lamp/actinic flux filename
};

/**
 * Abstract photolysis rate engine.
 *
 * Three implementations provide J-values for different use cases:
 *   MCMPhoto      — empirical CL·cos^CMM·exp(-CNN·sec) parameterization
 *   HybridPhoto   — F0AM TUV 4D lookup table interpolation
 *   BottomUpPhoto — cross-section × quantum-yield × lamp-flux integration
 */
class PhotolysisEngine
{
public:
  enum Method { MCM_SZA, HYBRID, BOTTOMUP };

  virtual ~PhotolysisEngine() = default;

  /// Compute photolysis rates J[0..nJ-1].
  virtual void computeJ(Real cos_sza,
                         Real sec_sza,
                         const PhysParams & params,
                         std::vector<Real> & j_out) const = 0;

  /// Number of photolysis rates.
  virtual unsigned int nJ() const = 0;

  /// Factory: create engine from method enum + construction parameters.
  static std::unique_ptr<PhotolysisEngine> create(Method method,
                                                    const PhotolysisParams & pparams);
};
