//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "KPPGeneratedMechanism.h"
#include "BottomUpJIntegrator.h"

#include <map>
#include <memory>

/**
 * Material that evaluates a KPP-generated mechanism RHS and sparse analytical
 * Jacobian at each quadrature point.
 */
class KPPMechanismMaterial : public Material
{
public:
  static InputParameters validParams();

  KPPMechanismMaterial(const InputParameters & params);

protected:
  virtual void computeQpProperties() override;

private:
  Real airDensity() const;
  void applyPhotolysisGlobals();

  KPPGeneratedMechanism _mechanism;

  const Real _temperature;
  const Real _air_density;
  const Real _water_vapor;
  const Real _pressure;
  const Real _jfac;
  const bool _roof_open;
  const MooseEnum _units;
  const MooseEnum _photolysis_scheme;
  const std::string _lamp_flux_file;
  const std::string _bottomup_data_dir;

  std::unique_ptr<BottomUpJIntegrator> _bottomup_integrator;
  std::map<std::string, Real> _cached_bottomup_j;
  Real _cached_bottomup_temperature;
  Real _cached_bottomup_pressure;
  bool _bottomup_j_valid;

  unsigned int _n_species;
  std::vector<const VariableValue *> _species_vals;

  MaterialProperty<std::vector<Real>> & _kpp_rhs;
  MaterialProperty<std::vector<Real>> & _kpp_jacobian_dense;
};
