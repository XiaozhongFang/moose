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

class MooseVariableScalar;

/**
 * Postprocessor that computes the total RO2 (peroxy radical) concentration
 * from the MCMBoxModel, equivalent to AtChem2's RO2 output in
 * environmentVariables.output.
 */
class MCMRO2Postprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();

  MCMRO2Postprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual Real getValue() const override;

protected:
  /// Reference to the MCMBoxModel UserObject
  const MCMBoxModel & _box_model;

  /// Pointers to species ScalarVariables for reading concentrations
  std::vector<const MooseVariableScalar *> _species_vars;

  /// Cached RO2 sum value
  Real _ro2_value;
};
