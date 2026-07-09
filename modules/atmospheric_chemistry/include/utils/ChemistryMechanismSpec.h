//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MechanismLoader.h"

#include <string>
#include <vector>
#include <set>

/**
 * Shared configuration and parsed mechanism data for atmospheric chemistry Actions.
 *
 * Encapsulates the result of parsing a mechanism file (species list, RO2 species,
 * reaction rates, J-values) so both AtmosphericChemistryBoxAction and
 * AtmosphericChemistryCoupledAction can share the same loading logic without
 * duplicating the MechanismLoader or parseKPPSpecies() code.
 *
 * Usage:
 *   ChemistryMechanismSpec spec(_mechanism_file, _chem_solver, _app,
 *                               getParam<...>(...), ...);
 *   spec.species()        → std::vector<std::string>
 *   spec.ro2Species()     → std::vector<std::string>
 *   spec.mechanismData()  → MechanismData &
 */
class ChemistryMechanismSpec
{
public:
  /**
   * Parse a mechanism file (MCM FACSIMILE .fac or KPP .kpp) and extract
   * all species names, reaction data, and photolysis parameters.
   *
   * @param mechanism_file   Path to the mechanism file (.fac or .kpp)
   * @param chem_solver      The chem_solver parameter string (from Action params)
   * @param mcm_version      MCM version string (for .fac files)
   * @param photo_path       Path to MCM photolysis-rates file
   * @param peroxy_path      Path to MCM peroxy-radicals file
   * @param input_files      MOOSE input file paths (for relative resolution)
   */
  ChemistryMechanismSpec(const std::string & mechanism_file,
                         const std::string & chem_solver,
                         const std::string & mcm_version,
                         const std::string & photo_path,
                         const std::string & peroxy_path,
                         const std::vector<std::string> & input_files);

  ///@{
  /// Accessors
  const std::vector<std::string> & species() const { return _species; }
  const std::vector<std::string> & ro2Species() const { return _ro2_species; }
  const MechanismData & mechanismData() const { return _mech_data; }
  const std::set<std::string> & baseVariables() const { return _base_vars; }
  bool isKPP() const { return _is_kpp; }
  ///@}

private:
  /// Parse KPP .spc files to extract species names (recursive #INCLUDE/#MODEL)
  std::vector<std::string> parseKPPSpecies(const std::string & kpp_file) const;
  /// Parse species names from KPP-generated Parameters.h when available.
  std::vector<std::string> parseKPPGeneratedSpecies(const std::string & kpp_file) const;

  /// Parsed mechanism data (populated for .fac files)
  MechanismData _mech_data;

  /// Species names in mechanism order
  std::vector<std::string> _species;

  /// RO2 peroxy-radical species (from .fac files)
  std::vector<std::string> _ro2_species;

  /// Base/environment variable names referenced by rate expressions
  std::set<std::string> _base_vars;

  /// Whether the mechanism is KPP format
  bool _is_kpp;
};
