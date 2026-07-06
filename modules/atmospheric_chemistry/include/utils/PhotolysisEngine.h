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
#include <vector>

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

  /// Factory: create engine from MooseEnum.
  static std::unique_ptr<PhotolysisEngine> create(Method method);
};
