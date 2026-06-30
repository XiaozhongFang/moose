//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * AuxKernel that computes RO2 = sum of peroxy radical species concentrations.
 * Used in coupled (FEM) mode to expose the RO2 diagnostic variable.
 *
 * The list of RO2 species is provided by AtmosphericChemistryAction from the
 * mechanism parser (same list used by MCMRatesMaterial internally).
 */
class MCMRO2Aux : public AuxKernel
{
public:
  static InputParameters validParams();
  MCMRO2Aux(const InputParameters & params);

protected:
  virtual Real computeValue() override;

private:
  /// Coupled variable values for each RO2 species
  std::vector<const VariableValue *> _ro2_vals;
};
