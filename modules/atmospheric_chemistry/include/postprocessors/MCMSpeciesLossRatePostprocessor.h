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
 * Postprocessor that outputs the total loss rate for a chemical species.
 * Loss = sum over reactions consuming the species (negative stoichiometry)
 * in molec/cm³/s.
 *
 * Analogous to AtChem2 lossRates.output.
 */
class MCMSpeciesLossRatePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MCMSpeciesLossRatePostprocessor(const InputParameters & params);

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

/**
 * Postprocessor that outputs the total production rate for a chemical species.
 * Production = sum over reactions producing the species (positive stoichiometry)
 * in molec/cm³/s.
 *
 * Analogous to AtChem2 productionRates.output.
 */
class MCMSpeciesProductionRatePostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MCMSpeciesProductionRatePostprocessor(const InputParameters & params);

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
