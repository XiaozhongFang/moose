//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

/**
 * Self-contained result of loading a chemical mechanism from a .fac file.
 *
 * Holds all data extracted from the mechanism, including the resolved
 * photolysis file path and the full set of photolysis J<N> parameters
 * from the photolysis-rates file.  Designed to replace the ad-hoc data
 * loading and path-resolution logic embedded in AtmosphericChemistryAction.
 */
struct MechanismData
{
  /// A single chemical reaction from the mechanism
  struct Reaction
  {
    std::string rate_expression;                               ///< Rate coefficient expression string
    std::vector<std::pair<Real, std::string>> reactants;       ///< (stoichiometric_coeff, species_name)
    std::vector<std::pair<Real, std::string>> products;        ///< (stoichiometric_coeff, species_name)
  };

  // ── Core mechanism data ────────────────────────────────────────────
  std::vector<std::string> species;
  std::vector<Reaction> reactions;
  std::vector<std::vector<Real>> stoichiometric_matrix;         ///< stoichiometry[species_idx][reaction_idx]

  // ── Rate coefficients (topologically sorted) ───────────────────────
  std::map<std::string, std::string> rate_coefficients;         ///< name → expression (raw)
  std::map<std::string, std::string> converted_coefficients;    ///< name → expression (converted for fparser)
  std::vector<std::string> eval_order;                          ///< coefficient names in topological order
  std::vector<std::string> reaction_rate_expressions;           ///< per-reaction rate expression string
  std::set<std::string> coefficient_names;                      ///< all coefficient names
  std::set<std::string> base_variables;                         ///< TEMP/M/O2/N2/H2O + J<N> names

  // ── RO2 diagnostic ─────────────────────────────────────────────────
  std::vector<std::string> ro2_species;

  // ── Photolysis parameters (mechanism-referenced J<N> only) ────────
  std::vector<unsigned int> j_numbers;
  std::vector<Real> j_CL;                    ///< CL coefficient per J<N>
  std::vector<Real> j_CMM;                   ///< CMM coefficient per J<N>
  std::vector<Real> j_CNN;                   ///< CNN coefficient per J<N>

  // ── Resolved photolysis file path ──────────────────────────────────
  std::string resolved_photo_path;           ///< Absolute or input-relative photo file path

  // ── Full photolysis parameter set (all J<N> entries from file) ────
  // Used by coupled mode for MCMPhotolysisPostprocessor output.
  std::vector<unsigned int> j_numbers_all;
  std::vector<Real> j_cl_values;             ///< CL per entry
  std::vector<Real> j_cmm_values;            ///< CMM per entry
  std::vector<Real> j_cnn_values;            ///< CNN per entry
};

/**
 * Standalone mechanism loader.
 *
 * Parses an MCM Facsimile (.fac) mechanism file, resolves the photolysis
 * parameter file path, and returns a self-contained MechanismData struct.
 * Encapsulates the MCMFacsimileParser call, path resolution, and full
 * photolysis file re-read — replacing the duplicated logic previously
 * embedded in AtmosphericChemistryAction.
 *
 * This class has no MOOSE object dependencies: it only needs file paths
 * and string data.  It can be used independently of any Action or
 * UserObject.
 */
class MechanismLoader
{
public:
  /// Reaction type alias for convenience
  using Reaction = MechanismData::Reaction;

  /**
   * Load and parse a complete chemical mechanism.
   *
   * @param mechanism_file  Path to the .fac mechanism file (relative to
   *                        the working directory or an input-file search path).
   * @param photo_path      User-specified photolysis parameter file path.
   * @param mcm_version     MCM version string (e.g. "v3.3.1").
   * @param peroxy_path     Path to the peroxy-radicals reference file.
   * @param input_file_dirs Directories of the input files being processed,
   *                        used for relative path resolution of the photolysis file.
   * @return Fully populated MechanismData.
   */
  static MechanismData load(const std::string & mechanism_file,
                            const std::string & photo_path,
                            const std::string & mcm_version,
                            const std::string & peroxy_path,
                            const std::vector<std::string> & input_file_dirs);

private:
  /**
   * Resolve the photolysis parameter file path.
   *
   * Search order:
   *   1. Try photo_path as-is (relative to working directory).
   *   2. Try relative to each input file directory.
   *   3. Try relative to the mechanism file's directory.
   *
   * @param photo_path       User-provided photolysis file path.
   * @param input_file_dirs  Input file directories for resolution.
   * @param mech_file_dir    Directory of the .fac mechanism file.
   * @return Resolved path, or the original photo_path if none found.
   */
  static std::string resolvePhotoPath(const std::string & photo_path,
                                       const std::vector<std::string> & input_file_dirs,
                                       const std::string & mech_file_dir);

  /**
   * Read all J<N> entries from a photolysis-rates parameter file.
   *
   * Parses the MCM-format photolysis file and extracts every J<N>
   * entry's (CL, CMM, CNN) parameters, regardless of whether they are
   * referenced by the mechanism.
   *
   * @param resolved_path  Resolved path to the photolysis file.
   * @param[out] j_numbers     Vector of J<N> numbers.
   * @param[out] j_cl          Vector of CL values.
   * @param[out] j_cmm         Vector of CMM values.
   * @param[out] j_cnn         Vector of CNN values.
   */
  static void loadFullPhotolysisSet(const std::string & resolved_path,
                                     std::vector<unsigned int> & j_numbers,
                                     std::vector<Real> & j_cl,
                                     std::vector<Real> & j_cmm,
                                     std::vector<Real> & j_cnn);
};
