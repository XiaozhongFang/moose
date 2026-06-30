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
#include <vector>
#include <string>
#include <map>

/**
 * Bottom-Up photolysis integrator — computes J-values by integrating
 *   J = ∫ QY(λ) · CS(λ) · F(λ) dλ
 * over the lamp/actinic flux wavelength grid, following the F0AM
 * IntegrateJ.m / J_BottomUp.m algorithm.
 *
 * ### Supported data formats
 * | Type | Columns | Description                            |
 * |------|---------|----------------------------------------|
 * | 0    | —       | Scalar (quantum yield only)            |
 * | 1    | 2       | CSV: wl [nm], value                    |
 * | 2    | 3       | CSV: wl, val@T1, val@T2 (linear T-interp) |
 * | 3    | 2       | TXT: space/tab-separated wl, value     |
 *
 * ### Mapping file format (bottomup_jmap.dat)
 * One line per reaction:
 *   JNAME  CS_FILE  CS_TYPE  QY_FILE  QY_TYPE
 * Entries with CS_FILE = "SCALAR" are ignored (valid only for QY).
 *
 * ### Usage
 * @code
 *   BottomUpJIntegrator integrator(data_dir);
 *   integrator.loadLampFlux("ExampleLightFlux.txt");
 *   integrator.loadReactionMap("bottomup_jmap.dat");
 *   Real J4 = integrator.computeJ("J4", 298.0, 1013.0);
 *   // or compute all at once:
 *   auto allJ = integrator.computeAllJ(298.0, 1013.0);
 * @endcode
 */
class BottomUpJIntegrator
{
public:
  BottomUpJIntegrator(const std::string & data_dir);

  /// Load lamp/actinic flux file (2-column: wl[nm], flux[photons/cm²/s/nm])
  void loadLampFlux(const std::string & flux_file);

  /// Load reaction mapping from a data file
  void loadReactionMap(const std::string & map_file);

  /// Register a single reaction's CS and QY sources manually
  void addReaction(const std::string & jname, const std::string & cs_file, int cs_type,
                   const std::string & qy_file, int qy_type);

  /// Compute J-value for one reaction at given T, P
  Real computeJ(const std::string & jname, Real T, Real P) const;

  /// Register a built-in cross-section formula (cs_type=10)
  void addBuiltinCS(const std::string & jname, const std::string & species,
                    const std::string & qy_file, int qy_type);
  /// Register a built-in quantum yield formula (qy_type=10)
  void addBuiltinQY(const std::string & jname, const std::string & cs_file, int cs_type,
                    const std::string & species);

  /// Compute all registered J-values
  std::map<std::string, Real> computeAllJ(Real T, Real P) const;

  bool hasJValue(const std::string & jname) const { return _reactions.count(jname) > 0; }

  unsigned int nJValues() const { return _reactions.size(); }

private:
  /// Load a 2-column data file → (wl, value)
  std::pair<std::vector<Real>, std::vector<Real>> loadCSV2(const std::string & path) const;

  /// Load a 3-column temperature-dependent CSV → (wl, value_interpolated_at_T)
  std::pair<std::vector<Real>, std::vector<Real>> loadCSV3(const std::string & path, Real T) const;

  /// Load a TXT file (space/tab separated, 2 columns) → (wl, value)
  std::pair<std::vector<Real>, std::vector<Real>> loadTXT(const std::string & path) const;

  /// Smear / convolve data onto a wavelength bin [xgl, xgu] (TUV numer.f interp2 algorithm)
  static Real smear(const std::vector<Real> & x, const std::vector<Real> & y, Real xgl, Real xgu);

  /// Trapezoidal integration
  static Real trapz(const std::vector<Real> & x, const std::vector<Real> & y);

  /// Built-in cross-section computation (called when cs_type==10)
  std::pair<std::vector<Real>, std::vector<Real>> computeCS_builtin(
      const std::string & species, Real T, Real P) const;
  /// Built-in quantum yield computation (called when qy_type==10)
  std::pair<std::vector<Real>, std::vector<Real>> computeQY_builtin(
      const std::string & species, Real T, Real P) const;
  /// 3-column CSV with explicit T1/T2 (for QY data without precomputed headers)
  std::pair<std::vector<Real>, std::vector<Real>> loadCSV3_for_QY(
      const std::string & path, Real T, Real T1, Real T2) const;

  std::string _data_dir;

  /// Lamp flux wavelength grid (nm) — this is the COMMON grid for all integrations
  std::vector<Real> _wl_flux;
  /// Lamp flux values [photons/cm²/s/nm]
  std::vector<Real> _flux;

  /// Flux wavelength bin boundaries: wllim[i] = wl[i] - 0.5*dwl, wllim[i+1] = wl[i] + 0.5*dwl
  std::vector<Real> _wllim;
  /// Bin widths
  std::vector<Real> _dwl;

  struct ReactionInfo
  {
    std::string cs_file; ///< Cross-section data file (relative to data_dir/CrossSections/)
    int cs_type;         ///< 1=CSV2, 2=CSV3(T-interp), 3=TXT
    std::string qy_file; ///< Quantum yield data file
    int qy_type;         ///< 0=scalar, 1=CSV2, 2=CSV3, 3=TXT
  };
  std::map<std::string, ReactionInfo> _reactions;
};
