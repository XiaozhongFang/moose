//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"
#include "MechanismLoader.h"

#include <string>
#include <vector>
#include <map>

/**
 * Action for FEM transport + chemistry (coupled mode) atmospheric chemistry.
 *
 * Creates MooseVariableFE (family=LAGRANGE) for each species,
 * MCMRatesMaterial or KPPMechanismMaterial for rate evaluation, and
 * TimeDerivative + source kernel for each species.
 *
 * Suitable for spatially-resolved simulations (5--50 species).
 *
 * Replaces [AtmosphericChemistry] mode=coupled.
 */
class AtmosphericChemistryCoupledAction : public Action
{
public:
  static InputParameters validParams();

  AtmosphericChemistryCoupledAction(const InputParameters & params);

  virtual void act() override;

protected:
  /// Build the reactant index matrix for the Material
  std::vector<std::vector<Real>> buildReactantMatrix() const;
  /// Derive the KPP shared-library path from the mechanism file.
  std::string kppLibraryPath(const std::string & mechanism_file) const;

  /// Create FE variables for each chemical species
  void actAddVariable();
  /// Create MCMRatesMaterial
  void actAddMaterial();
  /// Create TimeDerivative + ChemicalSourceKernel for each species
  void actAddKernel();
  /// Whether species advection kernels should be created
  bool hasAdvection() const;
  /// Whether species diffusion kernels should be created
  bool hasDiffusion() const;
  /// Whether density-weighted species diffusion kernels should be created
  bool hasDensityWeightedDiffusion() const;
  /// Whether species variables and kernels should use finite volume objects
  bool useFV() const;
  /// Whether FV time/source kernels should use longitude-latitude spherical weights
  bool useSphericalFVTime() const;
  /// Whether spherical FV horizontal advection kernels should be created
  bool hasSphericalFVAdvection() const;
  /// Create optional transport kernels for a species
  void addTransportKernels(const std::string & species_name);
  /// Create optional FV transport kernels for a species
  void addFVTransportKernels(const std::string & species_name);

  /// Parsed species names in mechanism order
  std::vector<std::string> _species;
  /// RO2 peroxy-radical species
  std::vector<std::string> _ro2_species;
  /// Full mechanism data (from MechanismLoader)
  MechanismData _mech_data;
  /// Chemical solver/backend selector.
  std::string _chem_solver;
  /// Whether this Action is using a KPP-generated mechanism.
  bool _use_kpp;
  /// Whether RO2 diagnostic variable should be created
  bool _ro2_diagnostic_enabled;
};
