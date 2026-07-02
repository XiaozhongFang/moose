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
 * Outputs two vectors:
 *   - ro2_count (1 element): total number of detected RO2 species
 *   - ro2_species (N elements): 0-based species indices for each RO2
 *
 * Used in CSVDiff tests to validate that the RO2 detection is correct.
 * The gold CSV matches the VPP output format:
 *   ro2_count,ro2_species
 *   117,<index_0>
 *   0,<index_1>
 *   ...
 *
 * CSVDiff compares columns common to both the gold and test output.
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
  /// Vector: RO2 species indices
  VectorPostprocessorValue & _ro2_species;
};
