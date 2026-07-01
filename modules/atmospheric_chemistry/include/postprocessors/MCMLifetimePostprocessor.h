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
 * Postprocessor that computes the chemical lifetime of a species.
 * tau = C / loss_rate  (F0AM lifetime.m equivalent)
 *
 * Formula: tau = C[s] / speciesLossRate(s)  [seconds]
 * If loss rate is near zero, returns a large sentinel value (1e12 s).
 */
class MCMLifetimePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MCMLifetimePostprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const MCMBoxModel & _box_model;
  const unsigned int _species_index;
  const std::vector<VariableName> _species_vars;
  mutable std::vector<Real> _C_buffer;
  const std::vector<Real> & _buildC() const;
};
