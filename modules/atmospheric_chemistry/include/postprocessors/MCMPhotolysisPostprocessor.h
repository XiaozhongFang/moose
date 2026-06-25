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

// Forward declaration for coupled-mode material support
class MCMRatesMaterial;

/**
 * Postprocessor that outputs photolysis rate values (J-values).
 * Supports both box mode (MCMBoxModel) and coupled mode (MCMRatesMaterial).
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
  /// Box-mode: reference to the MCMBoxModel UserObject (null if coupled mode)
  const MCMBoxModel * _box_model;

  /// Coupled-mode: J-value material property (from MCMRatesMaterial)
  const MaterialProperty<std::vector<Real>> * _j_material_prop;

  /// Coupled-mode: J-number-to-index mapping (from MCMRatesMaterial)
  const MaterialProperty<std::vector<unsigned int>> * _j_number_list_prop;

  /// If true, output all J values (this PP represents the first one)
  bool _output_all;

  /// Specific J index to output (1-based J number, e.g., 1 for J1)
  unsigned int _j_number;

  /// Cached value
  Real _j_value;

  /// Compute cosx for SZA-based J calculation
  Real _cosx;
};
