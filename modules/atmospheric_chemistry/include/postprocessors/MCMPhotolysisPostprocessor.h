//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"
#include "MCMBoxModel.h"

/**
 * Postprocessor that outputs photolysis rate values (J-values) from the
 * MCMBoxModel.  Equivalent to AtChem2's photolysisRates.output.
 *
 * Supports two modes:
 *   - output_all = true:  outputs ALL J values as a vector (for CSV column output)
 *   - output_all = false: outputs a SINGLE named J value (e.g., J1, J4)
 */
class MCMPhotolysisPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  MCMPhotolysisPostprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual Real getValue() const override;

protected:
  /// Reference to the MCMBoxModel UserObject
  const MCMBoxModel & _box_model;

  /// If true, output all J values (this PP represents the first one)
  bool _output_all;

  /// Specific J index to output (1-based J number, e.g., 1 for J1)
  unsigned int _j_number;

  /// Cached value
  Real _j_value;

  /// Compute cosx for SZA-based J calculation
  Real _cosx;
};
