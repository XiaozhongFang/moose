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

#include <vector>

/**
 * MCM empirical SZA photolysis engine.
 *
 * Computes J-values using the standard MCM parameterization:
 *   J[N] = CL[N] * cos(θ)^CMM[N] * exp(-CNN[N] * sec(θ)) * JFAC
 *
 * where CL, CMM, and CNN are empirical coefficients read from the
 * photolysis-rates file accompanying the MCM mechanism, and JFAC
 * is a user-specified scaling factor from PhysParams.
 */
class MCMPhoto : public PhotolysisEngine
{
public:
  MCMPhoto(const std::vector<Real> & cl,
           const std::vector<Real> & cmm,
           const std::vector<Real> & cnn,
           const std::vector<unsigned int> & j_numbers);

  /**
   * Compute photolysis rates into j_out.
   *
   * @param cos_sza  Cosine of the solar zenith angle
   * @param sec_sza  Secant of SZA (1/cos_sza, capped to avoid Inf)
   * @param params   Physical parameters (uses jfac)
   * @param j_out    Output vector, resized to nJ()
   */
  void computeJ(Real cos_sza,
                Real sec_sza,
                const PhysParams & params,
                std::vector<Real> & j_out) const override;

  /// Number of photolysis rates (= size of CL/CMM/CNN arrays).
  unsigned int nJ() const override { return _cl.size(); }

private:
  std::vector<Real> _cl;
  std::vector<Real> _cmm;
  std::vector<Real> _cnn;
  std::vector<unsigned int> _j_numbers;
};
