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
 * Postprocessor that outputs solar parameters from MCMBoxModel (box mode)
 * or MCMRatesMaterial (coupled mode, via material properties solar_*).
 * Mirrors AtChem2's zenith_data_mod variables:
 *   cosx, secx, lha, sinld, cosld, eqtime, latitude, longitude, declination
 */
class MCMSolarPostprocessor : public GeneralPostprocessor
{
public:
  static InputParameters validParams();
  MCMSolarPostprocessor(const InputParameters & params);

  virtual void initialize() override {}
  virtual void execute() override;
  virtual Real getValue() const override;

protected:
  /// Box-mode: reference to the MCMBoxModel UserObject (null if coupled mode)
  const MCMBoxModel * _box_model;

  /// Coupled-mode: material property for the selected solar parameter (null if box mode)
  const MaterialProperty<Real> * _solar_prop;

  /// Which solar parameter to output
  enum ParamType { COSX, SECX, LHA, SINLD, COSLD, EQTIME, LAT, LON, DEC } _param;
  Real _value;
};
