//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralVectorPostprocessor.h"
#include "MCMBoxModel.h"

/**
 * VectorPostprocessor that outputs the list of RO2 (peroxy radical) species
 * detected by MCMBoxModel from the parsed mechanism file.
 *
 * Outputs one count vector and one vector per RO2 species name:
 *   - ro2_count (1 element): total number of detected RO2 species
 *   - <species name> (1 element): 1 if that species was detected as RO2
 *
 * Used in CSVDiff tests to validate the detected RO2 species by name. Missing
 * or extra RO2 species appear as CSV header differences.
 */
class MCMRO2ListPostprocessor : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  MCMRO2ListPostprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;

protected:
  /// Reference to the MCMBoxModel UserObject
  const MCMBoxModel & _box_model;

  /// Vector: RO2 species count (single element)
  VectorPostprocessorValue & _ro2_count;
  /// Vectors keyed by detected RO2 species names.
  std::vector<VectorPostprocessorValue *> _ro2_species_flags;
};
