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
 * Postprocessor that outputs the rate of a single chemical reaction.
 * Rate = k[r] * prod(C[reactants])  (molec/cm³/s or ppb/s)
 *
 * Analogous to AtChem2 productionRates.output / lossRates.output
 * and F0AM ExtractRates.m.
 */
class MCMReactionRatePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MCMReactionRatePostprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const MCMBoxModel & _box_model;
  const unsigned int _reaction_index;
  const std::vector<VariableName> _species_vars;
  mutable std::vector<Real> _C_buffer;
  const std::vector<Real> & _buildC() const;
};
