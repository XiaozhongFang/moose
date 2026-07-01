//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMBoxModel.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "pcrecpp.h"
#include <algorithm>
#include <regex>
#include <sstream>
#include <utility>

registerMooseObject("AtmosphericChemistryApp", MCMBoxModel);

// ---- StoichMatrix::build ----------------------------------------------------

void
StoichMatrix::build(const ParsedMechanism & mech, Format fmt)
{
  format = fmt;
  nSpecies = mech.species.size();
  nReactions = mech.reactions.size();

  // Clear all format vectors unconditionally — prevents stale data if build()
  // is ever called more than once (e.g. mechanism reload with different format).
  csr_cols.clear(); csr_vals.clear(); csr_row_ptr.clear();
  lhs_species.clear(); lhs_coeff.clear(); lhs_row_ptr.clear();
  rhs_species.clear(); rhs_coeff.clear(); rhs_row_ptr.clear();
  dense.clear();
  csc_rows.clear(); csc_c_vals.clear(); csc_col_ptr.clear();

  // Build name -> index map (needed by COO format for reactant/product lookup)
  std::unordered_map<std::string, unsigned int> name_to_idx;
  for (unsigned int i = 0; i < nSpecies; ++i)
    name_to_idx[mech.species[i]] = i;

  switch (format)
  {
    case CSR:
    {
      csr_row_ptr.resize(nReactions + 1);
      csr_row_ptr[0] = 0;
      csr_cols.clear();
      csr_vals.clear();
      for (unsigned int r = 0; r < nReactions; ++r)
      {
        for (unsigned int s = 0; s < nSpecies; ++s)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csr_cols.push_back((int)s);
            csr_vals.push_back(val);
          }
        }
        csr_row_ptr[r + 1] = csr_cols.size();
      }
      break;
    }
    case COO:
    {
      lhs_row_ptr.resize(nReactions + 1);
      rhs_row_ptr.resize(nReactions + 1);
      lhs_row_ptr[0] = rhs_row_ptr[0] = 0;
      lhs_species.clear(); lhs_coeff.clear();
      rhs_species.clear(); rhs_coeff.clear();

      for (unsigned int r = 0; r < nReactions; ++r)
      {
        const auto & rx = mech.reactions[r];
        for (const auto & [coeff, name] : rx.reactants)
        {
          auto it = name_to_idx.find(name);
          if (it != name_to_idx.end())
          {
            lhs_species.push_back((int)it->second);
            lhs_coeff.push_back(coeff);
          }
        }
        lhs_row_ptr[r + 1] = lhs_species.size();

        for (const auto & [coeff, name] : rx.products)
        {
          auto it = name_to_idx.find(name);
          if (it != name_to_idx.end())
          {
            rhs_species.push_back((int)it->second);
            rhs_coeff.push_back(coeff);
          }
        }
        rhs_row_ptr[r + 1] = rhs_species.size();
      }
      break;
    }
    case DENSE:
    {
      dense.assign(nReactions, std::vector<Real>(nSpecies, 0.0));
      for (unsigned int s = 0; s < nSpecies; ++s)
        for (unsigned int r = 0; r < nReactions; ++r)
          dense[r][s] = mech.stoichiometry[s][r];
      break;
    }
    case CSC:
    {
      // Build species-major (column) storage from dense stoichiometry.
      // Simultaneously populate CSR fields as a row-iteration forward index.
      csc_col_ptr.resize(nSpecies + 1);
      csc_col_ptr[0] = 0;
      csc_rows.clear();
      csc_c_vals.clear();

      csr_row_ptr.resize(nReactions + 1);
      csr_row_ptr[0] = 0;
      csr_cols.clear();
      csr_vals.clear();

      // Two-pass: first count nonzeros per species, then fill
      for (unsigned int s = 0; s < nSpecies; ++s)
      {
        for (unsigned int r = 0; r < nReactions; ++r)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csc_rows.push_back((int)r);
            csc_c_vals.push_back(val);
          }
        }
        csc_col_ptr[s + 1] = csc_rows.size();
      }

      // Build CSR forward index for row iteration
      for (unsigned int r = 0; r < nReactions; ++r)
      {
        for (unsigned int s = 0; s < nSpecies; ++s)
        {
          Real val = mech.stoichiometry[s][r];
          if (std::abs(val) > 1e-30)
          {
            csr_cols.push_back((int)s);
            csr_vals.push_back(val);
          }
        }
        csr_row_ptr[r + 1] = csr_cols.size();
      }
      break;
    }
  }
}

// ---- MCMBoxModel ------------------------------------------------------------

InputParameters
MCMBoxModel::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += FunctionParserUtils<false>::validParams();
  params.addParam<std::string>(
      "mechanism_file", "", "Path to MCM Facsimile (.fac) mechanism file for auto-parsing");
  params.addParam<std::string>(
      "photolysis_file", "", "Path to MCM photolysis-rates data file (e.g. mcm_photolysis_rates_v3.3.1.dat)");
  params.addParam<Real>("temperature", 298.15, "Ambient temperature (K)");
  params.addParam<Real>("air_density", 2.46e19,
      "Air number density (molecules/cm^3).  If press is set (>0), computed "
      "dynamically from press/temp using ideal gas law (AtChem2 calcAirDensity).");
  params.addParam<Real>("water_vapor", 2.46e17,
      "Background water vapor (molecules/cm^3).  If rh is set (>=0), computed "
      "dynamically from rh/temp/press using Vaisala 2013 formula (AtChem2 convertRHtoH2O).");
  params.addParam<Real>("press", 0.0,
      "Pressure (mbar).  If >0, air_density is computed via ideal gas law: "
      "M = 1e-6 * NA/R * (press*100 / temp).  Default 0 = use air_density directly.");
  params.addParam<Real>("rh", -1.0,
      "Relative humidity (%).  If >=0, water_vapor is computed via Vaisala 2013. "
      "Default -1 = use water_vapor directly (AtChem2 sentinel convention).");
  params.addParam<Real>("dilute", 0.0,
      "Dilution rate coefficient kdil (/s).  If >0, dC/dt -= kdil*(C-C_bkgd).");
  params.addParam<Real>("blheight", 0.0,
      "Boundary layer height (m) — informational, used with dilute.");
  params.addParam<Real>("latitude", 51.51, "Latitude (deg N)");
  params.addParam<Real>("longitude", 0.13, "Longitude (deg E)");
  params.addParam<unsigned int>("day", 21, "Day of month");
  params.addParam<unsigned int>("month", 6, "Month");
  params.addParam<unsigned int>("year", 2010, "Year");
  MooseEnum units_enum("molec_cm3 ppb", "molec_cm3");
  params.addParam<MooseEnum>("units", units_enum,
      "Concentration units for input/output: 'molec_cm3' (molecules/cm³, default) "
      "or 'ppb' (parts per billion by volume). When 'ppb', the ChemistryODEKernel "
      "automatically converts between ppb and molec/cm³ using the air density.");

  params.addParam<Real>("jfac", 1.0, "JFAC scaling factor");
  params.addParam<Real>("default_ic", 0.0,
      "Default initial concentration (molec/cm³) for species without explicit ICs.\n"
      "Set to a small positive value (e.g. 1e6) to prevent Jacobian singularities\n"
      "when many species start at zero. 0.0 = no default (original F0AM behavior).");
  params.addParam<bool>("use_limiting_reagent", false,
      "Enable F0AM-style limiting-reagent formulation for RO2+RO2 termination "
      "reactions: rate = k * min([A],[B])² instead of k * [A] * [B]. "
      "Default false (standard MCM chemistry). Set true for F0AM-compatible "
      "RO2 termination or when comparing against F0AM reference outputs.");
  MooseEnum stoich_fmt("CSR COO DENSE CSC", "CSR");
  params.addParam<MooseEnum>(
      "stoich_format", stoich_fmt,
      "Stoichiometric matrix storage format.  CSR = compressed sparse row (PETSc AIJ-compatible, "
      "HPC default).  COO = AtChem2-style split reactant/product.  DENSE = dense 2D array "
      "(best for < ~50 species).  CSC = compressed sparse column (species-major; enables "
      "column queries 'which reactions involve species X?').  Analogous to PETSc -mat_type.");

  MooseEnum photo_scheme("MCM_SZA HYBRID BOTTOMUP", "MCM_SZA");
  params.addParam<MooseEnum>("photolysis_scheme", photo_scheme,
      "Photolysis scheme: MCM_SZA (empirical SZA formula), "
      "HYBRID (F0AM TUV 4D lookup table interpolation), or "
      "BOTTOMUP (F0AM lab-chamber cross-section × quantum-yield × lamp-flux integration)");

  params.addParam<std::string>("hybrid_table_dir", "",
      "Directory containing F0AM Hybrid J table files. Required if photolysis_scheme=HYBRID.");

  params.addParam<std::string>("lamp_flux_file", "",
      "Path to lamp/actinic flux file (relative to bottomup_data_dir). "
      "Required if photolysis_scheme=BOTTOMUP.");

  params.addParam<std::string>("bottomup_data_dir",
      "../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup",
      "Directory containing BottomUp photolysis data files "
      "(CrossSections/, QuantumYields/, bottomup_jmap.dat, lamp flux file).");

  params.addParam<Real>("albedo", 0.1, "Surface albedo (0-1), used by HYBRID scheme");
  params.addParam<Real>("o3column", 350.0, "O3 column in Dobson Units, used by HYBRID scheme");
  params.addParam<Real>("altitude", 0.0, "Altitude in meters, used by HYBRID scheme");

  MooseEnum integrator_enum("moose petsc_ts", "moose");
  params.addParam<MooseEnum>("integrator", integrator_enum,
      "ODE integrator for box mode: 'moose' (default, through MOOSE's Newton solver) "
      "or 'petsc_ts' (PETSc TS, bypasses MOOSE solver for chemistry)");

  MooseEnum ts_type_enum("bdf arkimex sundials", "bdf");
  params.addParam<MooseEnum>("petsc_ts_type", ts_type_enum,
      "PETSc TS integrator type: 'bdf' (default), 'arkimex', or 'sundials' (CVODE).");
  params.addParam<Real>("petsc_ts_rtol", 1e-6,
      "Relative tolerance for PETSc TS adaptive integrator.");
  params.addParam<Real>("petsc_ts_atol", 1e-10,
      "Absolute tolerance for PETSc TS adaptive integrator.");

  params.addClassDescription(
      "Centralized box model UserObject for atmospheric chemistry ODE systems.");
  return params;
}

MCMBoxModel::MCMBoxModel(const InputParameters & params)
  : GeneralUserObject(params),
    FunctionParserUtils<false>(params),
    _n_species(0), _n_reactions(0),
    _photolysis_method(MCM_SZA),
    _units_ppb(getParam<MooseEnum>("units") == "ppb"),
    _ppb_to_molec(1.0),  // default 1.0 (no conversion); set dynamically in evaluateCoefficients
    _lat(getParam<Real>("latitude")),
    _lon(getParam<Real>("longitude")),
    _day((int)getParam<unsigned int>("day")),
    _month((int)getParam<unsigned int>("month")),
    _year((int)getParam<unsigned int>("year")),
    _kdil(0.0),
    _tgauss(0.0),
    _t_start_dil(0.0),
    _use_gaussian(false),
    _roof_open(true),
    _j_index_start(0),
    _temperature(getParam<Real>("temperature")),
    _air_density(getParam<Real>("air_density")),
    _water_vapor(getParam<Real>("water_vapor")),
    _press(getParam<Real>("press")),
    _rh(getParam<Real>("rh")),
    _blheight(getParam<Real>("blheight")),
    _jfac(getParam<Real>("jfac")),
    _t(0.0),
    _dirty(true),
    _cached_bottomup_T(0.0),
    _cached_bottomup_P(0.0),
    _bottomup_j_valid(false),
    _use_petsc_ts(getParam<MooseEnum>("integrator") == "petsc_ts")
{
  // Load Hybrid table reader if photolysis scheme is HYBRID
  auto scheme = getParam<MooseEnum>("photolysis_scheme");
  if (scheme == "HYBRID")
  {
    std::string dir = getParam<std::string>("hybrid_table_dir");
    if (dir.empty())
      mooseError("MCMBoxModel: hybrid_table_dir is required when photolysis_scheme=HYBRID");
    if (!dir.empty() && dir[0] == '/')
      mooseError("MCMBoxModel: hybrid_table_dir must be relative, got absolute: ", dir);
    enableHybridPhotolysis(dir);
  }
  else if (scheme == "BOTTOMUP")
  {
    _photolysis_method = BOTTOMUP;
    std::string data_dir = getParam<std::string>("bottomup_data_dir");
    std::string flux_file = getParam<std::string>("lamp_flux_file");
    if (flux_file.empty())
      mooseError("MCMBoxModel: lamp_flux_file is required when photolysis_scheme=BOTTOMUP");
    loadBottomUpData(data_dir, flux_file);
  }
}

MCMBoxModel::~MCMBoxModel()
{
  // Clean up PETSc TS objects (safe if never initialized — pointers are null)
  if (_ts_J) { MatDestroy(&_ts_J); _ts_J = nullptr; }
  if (_ts_X) { VecDestroy(&_ts_X); _ts_X = nullptr; }
  if (_ts)   { TSDestroy(&_ts);   _ts = nullptr; }
}

void
MCMBoxModel::initialize()
{
  // Parse .fac file only once
  if (_n_species > 0) return;

  std::string mech_file = getParam<std::string>("mechanism_file");
  if (!mech_file.empty())
  {
    std::string photo_file = getParam<std::string>("photolysis_file");
    MCMFacsimileParser parser;
    ParsedMechanism mech = parser.parse(mech_file, photo_file);
    bool use_lr = getParam<bool>("use_limiting_reagent");
    loadMechanism(mech, use_lr);
    _console << "MCMBoxModel: Loaded " << _n_species << " species, "
             << _n_reactions << " reactions from " << mech_file << std::endl;

    // Apply default IC (if > 0) to uninitialized species.
    // This prevents Jacobian singularities when many species start at zero.
    // 1.0e6 molec/cm³ is chemically negligible (~4e-5 ppb).
    Real default_ic = getParam<Real>("default_ic");
    if (default_ic > 0.0)
    {
      for (const auto & sp_name : _species_names)
      {
        MooseVariableScalar & sv = _subproblem.getScalarVariable(0, sp_name);
        if (sv.sln()[0] == 0.0)
          sv.setValue(0, default_ic);
      }
      _console << "MCMBoxModel: Set default_ic = " << default_ic
               << " molec/cm³ for uninitialized species" << std::endl;
    }

    // Pre-compute rate coefficients (_k) and time-invariant solar parameters
    // (DEC, sinld, cosld, eqtime) so the first time step starts with correct
    // thermal chemistry even if the solver skips the initial residual evaluation.
    // Without this, _k stays at the default 1.0 until the first getDCdt() call,
    // which may not happen in the expected order during the first nonlinear solve.
    if (!_coeff_parsers.empty())
    {
      evaluateCoefficients();
      // Reset time-varying solar params to safe defaults — per-timestep
      // evaluateCoefficients() will overwrite them during the first solve.
      _solar_cosx = 0.0;
      _solar_secx = 1.0e2;
      _solar_lha = 0.0;
    }

    // Set up dilution if dilute > 0
    Real dilute = getParam<Real>("dilute");
    if (dilute > 0.0)
    {
      _kdil = dilute;
      _conc_bkgd.assign(_n_species, 0.0);
      _console << "MCMBoxModel: Dilution enabled, kdil = " << dilute << " /s" << std::endl;
    }

    // Initialize PETSc TS integrator if box mode uses standalone TS integration
    if (_use_petsc_ts)
      setupPETScTS();
  }
}

void
MCMBoxModel::computeDCdt(const std::vector<Real> & C, std::vector<Real> & dC) const
{
  // Lazy initialization: the framework may call computeDCdt() (via getDCdt())
  // before GeneralUserObject::initialize().  If the mechanism hasn't been
  // loaded yet, force initialization now so chemistry runs on the first step.
  if (_n_species == 0)
    const_cast<MCMBoxModel*>(this)->initialize();

  dC.assign(_n_species, 0.0);

  if (_n_reactions == 0 || _n_species == 0)
    return;

  // Step 1: compute reactant products G[i] = C[iG[i][0]] * C[iG[i][1]] * C[iG[i][2]]
  // Sentinel -1 means "no reactant" → concentration = 1.0 (F0AM padding convention).
  // Per.14: _scratch_G / _scratch_rates are mutable members pre-allocated once,
  // avoiding heap allocation on every dC/dt call (~272 KiB per call for full MCM).
  _scratch_G.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real c0 = (_iG[r][0] >= 0) ? C[_iG[r][0]] : 1.0;
    Real c1 = (_iG[r][1] >= 0) ? C[_iG[r][1]] : 1.0;
    Real c2 = (_iG[r][2] >= 0) ? C[_iG[r][2]] : 1.0;
    _scratch_G[r] = c0 * c1 * c2;
  }

  // Limiting-reagent override (F0AM RO2 termination reactions):
  // rate = k * min([RO2_i], [RO2_j]) instead of k * [RO2_i] * [RO2_j]
  // For RO2+RO2, the reaction is limited by whichever RO2 concentration is smaller.
  // OFF by default — opt-in via use_limiting_reagent=true to preserve gold CSV compatibility.
  if (_use_limiting_reagent && !_limiting_reagent.empty())
  {
    for (unsigned int r = 0; r < _n_reactions; ++r)
    {
      if (!_limiting_reagent[r])
        continue;
      const int i0 = _iG[r][0], i1 = _iG[r][1];
      if (i0 >= 0 && i1 >= 0)
      {
        Real c0 = C[i0], c1 = C[i1];
        Real minc = std::min(c0, c1);
        _scratch_G[r] = minc * minc;  // [min]² replaces [A]*[B]
      }
    }
  }

  // Step 2: rates = k .* G
  _scratch_rates.assign(_n_reactions, 0.0);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _scratch_rates[r] = _k[r] * _scratch_G[r];

  // Step 3: dC = f^T * rates  — format-agnostic iteration via StoichMatrix.
  // forEachInRow dispatches based on _stoich.format (CSR or COO); the lambda
  // is fully inlined, so there is zero virtual-call overhead.
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real rate = _scratch_rates[r];
    _stoich.forEachInRow(r, [&](int s, Real coeff) {
      dC[s] += coeff * rate;
    });
  }

  // Apply dilution if kdil > 0 (AtChem2 DILUTE parameter)
  if (_kdil > 0.0 && !_conc_bkgd.empty())
    for (unsigned int i = 0; i < _n_species; ++i)
      dC[i] -= _kdil * (C[i] - _conc_bkgd[i]);
}

void
MCMBoxModel::computeJacobianTriplets(
    const std::vector<Real> & C,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  J.clear();

  if (_n_reactions == 0 || _n_species == 0)
    return;

  // F0AM-style 3-term product Jacobian — format-agnostic via StoichMatrix.
  //   rate_r = k_r * C[i0] * C[i1] * C[i2]
  //   d(rate_r)/dC[j] = k_r * sum_{k where iG[k]==j} prod_{m≠k} C[iG[m]]
  //   J[s][j] += f[r][s] * d(rate_r)/dC[j]

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    // Sentinel -1 → concentration = 1.0, no derivative contribution
    const Real c0 = (i0 >= 0) ? C[i0] : 1.0;
    const Real c1 = (i1 >= 0) ? C[i1] : 1.0;
    const Real c2 = (i2 >= 0) ? C[i2] : 1.0;
    const Real k = _k[r];

    // Helper: for reactant j with derivative drate, emit J(s,j) for every
    // species s that participates in reaction r.
    auto emit_contrib = [&](unsigned int j, Real drate) {
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        J.emplace_back((unsigned int)s, j, drate * coeff);
      });
    };

    // Limiting-reagent (LR) Jacobian:
    // For LR reactions: rate = k * min([A],[B])^2  (F0AM RO2 termination)
    //   ∂rate/∂C_min = 2*k*min([A],[B])
    //   ∂rate/∂C_other = 0
    // OFF by default — opt-in via use_limiting_reagent=true.
    if (_use_limiting_reagent && !_limiting_reagent.empty() && r < _limiting_reagent.size() && _limiting_reagent[r])
    {
      if (i0 >= 0 && i1 >= 0)
      {
        Real minc = std::min(c0, c1);
        unsigned int min_idx = (c0 <= c1) ? (unsigned int)i0 : (unsigned int)i1;
        Real drate = 2.0 * k * minc;
        emit_contrib(min_idx, drate);
        // The non-limiting reactant has zero derivative contribution
      }
      continue; // Skip standard Jacobian for LR reactions
    }

    // Contribution from C[i0]: dr/dC[i0] = k * c1 * c2
    if (i0 >= 0) emit_contrib((unsigned int)i0, k * c1 * c2);
    // Contribution from C[i1]: dr/dC[i1] = k * c0 * c2
    if (i1 >= 0 && i1 != i0) emit_contrib((unsigned int)i1, k * c0 * c2);
    // Contribution from C[i2]: dr/dC[i2] = k * c0 * c1
    if (i2 >= 0 && i2 != i0 && i2 != i1) emit_contrib((unsigned int)i2, k * c0 * c1);
  }
}

// -- Cached single-species interface --

Real
MCMBoxModel::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  // _dirty flag is set by markDirty() in EVERY reinit() call (610x per residual),
  // but the first getDCdt() call rebuilds the cache and clears _dirty.
  // Subsequent calls (remaining 609 kernels) find _dirty=false and return cached.
  if ((_dirty) || _cached_dC.size() != _n_species)
  {
    _cached_C = C;  // store for evaluateCoefficients species lookup (deep copy)
    if (!_coeff_parsers.empty())
      const_cast<MCMBoxModel*>(this)->evaluateCoefficients();
    _cached_dC.resize(_n_species);
    computeDCdt(C, _cached_dC);
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_dC[idx] : 0.0;
}

Real
MCMBoxModel::getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const
{
  if (_dirty || _cached_diag_J.size() != _n_species)
  {
    _cached_C = C;
    _buildJacobianCache();
    _dirty = false;
  }
  return (idx < _n_species) ? _cached_diag_J[idx] : 0.0;
}

Real
MCMBoxModel::getJacobianOffDiagonal(unsigned int i, unsigned int j, const std::vector<Real> & C) const
{
  if (_dirty || _cached_od_row_ptr.size() != _n_species + 1)
  {
    _cached_C = C;
    _buildJacobianCache();
    _dirty = false;
  }
  // CSR binary search: row i columns are sorted at [_od_row_ptr[i] .. _od_row_ptr[i+1])
  size_t lo = _cached_od_row_ptr[i];
  size_t hi = _cached_od_row_ptr[i + 1];
  while (lo < hi)
  {
    size_t mid = lo + (hi - lo) / 2;
    unsigned int col = _cached_od_cols[mid];
    if (col == j)
      return _cached_od_vals[mid];
    else if (col < j)
      lo = mid + 1;
    else
      hi = mid;
  }
  return 0.0;
}

void
MCMBoxModel::_buildJacobianCache() const
{
  _cached_diag_J.assign(_n_species, 0.0);

  if (_n_reactions == 0 || _n_species == 0)
  {
    _cached_od_row_ptr.assign(_n_species + 1, 0);
    _cached_od_cols.clear();
    _cached_od_vals.clear();
    return;
  }

  // Per-row temporary accumulator: collects (col, val) pairs
  std::vector<std::vector<std::pair<unsigned int, Real>>> temp(_n_species);

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const int i0 = _iG[r][0], i1 = _iG[r][1], i2 = _iG[r][2];
    // Sentinel -1 → concentration = 1.0
    const Real c0 = (i0 >= 0) ? _cached_C[i0] : 1.0;
    const Real c1 = (i1 >= 0) ? _cached_C[i1] : 1.0;
    const Real c2 = (i2 >= 0) ? _cached_C[i2] : 1.0;
    const Real k = _k[r];

    // Helper: accumulate Jacobian contribution (s, j) += val
    auto accum = [&](unsigned int j, Real drate) {
      if (std::abs(drate) < 1e-30) return;
      _stoich.forEachInRow(r, [&](int s, Real coeff) {
        Real val = drate * coeff;
        unsigned int us = (unsigned int)s;
        if (j == us)
          _cached_diag_J[us] += val;
        else
          temp[us].emplace_back(j, val);
      });
    };

    // Limiting-reagent (LR) Jacobian: rate = k * min([A],[B])²
    // Only the limiting species contributes: ∂rate/∂C_min = 2*k*min([A],[B])
    // OFF by default — opt-in via use_limiting_reagent=true.
    bool is_lr = (_use_limiting_reagent && !_limiting_reagent.empty() && r < _limiting_reagent.size() && _limiting_reagent[r]);
    if (is_lr && i0 >= 0 && i1 >= 0)
    {
      Real minc = std::min(c0, c1);
      unsigned int min_idx = (c0 <= c1) ? (unsigned int)i0 : (unsigned int)i1;
      Real drate = 2.0 * k * minc;
      accum(min_idx, drate);
      continue; // Skip standard Jacobian for LR reactions
    }

    // dr/dC[i0] = k * c1 * c2
    accum((unsigned int)i0, k * c1 * c2);
    // dr/dC[i1] = k * c0 * c2
    if (i1 != i0) accum((unsigned int)i1, k * c0 * c2);
    // dr/dC[i2] = k * c0 * c1
    if (i2 != i0 && i2 != i1) accum((unsigned int)i2, k * c0 * c1);
  }

  // Flatten temp to CSR: sort each row, merge duplicates
  _cached_od_row_ptr.resize(_n_species + 1);
  _cached_od_row_ptr[0] = 0;
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    auto & row = temp[i];
    if (row.size() <= 1)
    {
      if (row.size() == 1)
        _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i] + 1;
      else
        _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i];
      continue;
    }
    std::sort(row.begin(), row.end());
    size_t w = 0;
    for (size_t k = 1; k < row.size(); ++k)
    {
      if (row[w].first == row[k].first)
        row[w].second += row[k].second;
      else
        row[++w] = row[k];
    }
    row.resize(w + 1);
    _cached_od_row_ptr[i + 1] = _cached_od_row_ptr[i] + row.size();
  }

  // Copy to flat CSR arrays
  size_t nnz = _cached_od_row_ptr[_n_species];
  _cached_od_cols.resize(nnz);
  _cached_od_vals.resize(nnz);
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    size_t base = _cached_od_row_ptr[i];
    for (size_t k = 0; k < temp[i].size(); ++k)
    {
      _cached_od_cols[base + k] = temp[i][k].first;
      _cached_od_vals[base + k] = temp[i][k].second;
    }
  }
}

void
MCMBoxModel::loadMechanism(const ParsedMechanism & mech, bool use_limiting_reagent)
{
  _n_species = mech.species.size();
  _n_reactions = mech.reactions.size();
  _species_names = mech.species;
  _reaction_names = mech.reaction_names;

  // Build RO2 species index list from parser's explicit ro2_species
  _ro2_indices.clear();
  for (const auto & ro2_name : mech.ro2_species)
  {
    auto it = std::find(_species_names.begin(), _species_names.end(), ro2_name);
    if (it != _species_names.end())
      _ro2_indices.push_back((unsigned int)(it - _species_names.begin()));
  }

  // Build stoichiometric matrix in the selected format (parameter "stoich_format").
  // CSR: compact, PETSc AIJ-compatible, optimal for HPC.
  // COO: AtChem2-style split reactant/product — enables loss/production diagnostics.
  StoichMatrix::Format fmt = StoichMatrix::CSR;
  if (isParamValid("stoich_format"))
  {
    MooseEnum fmt_enum = getParam<MooseEnum>("stoich_format");
    if (fmt_enum == "COO")        fmt = StoichMatrix::COO;
    else if (fmt_enum == "DENSE") fmt = StoichMatrix::DENSE;
    else if (fmt_enum == "CSC")   fmt = StoichMatrix::CSC;
  }
  _stoich.build(mech, fmt);

  // Copy reactant indices — convert vector<vector<int>> to flat array<int,3>
  // Per.16: single contiguous allocation vs 17k independent heap allocations.
  _iG.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    const auto & src = mech.reactant_indices[r];
    _iG[r] = {src[0], src[1], src[2]};
  }

  // Load limiting-reagent flags (F0AM RO2 termination reactions)
  // Opt-in via use_limiting_reagent=true (default false — preserves gold CSV compat).
  _use_limiting_reagent = use_limiting_reagent;
  _limiting_reagent = mech.is_limiting_reagent;
  _limiting_reactant = mech.limiting_reactant;

  // --- Evaluate rate coefficients (fparser for complex expressions) ---
  _k.assign(_n_reactions, 1.0);

  // Try simple numeric parse first; fall back to fparser
  std::map<std::string, Real> coeff_map;
  bool need_fparser = false;

  for (unsigned int i = 0; i < mech.coefficient_names.size(); ++i)
  {
    const std::string & expr = mech.coefficient_expressions[i];
    try
    {
      Real val = std::stod(expr);
      coeff_map[mech.coefficient_names[i]] = val;
    }
    catch (...)
    {
      need_fparser = true;
      break;
    }
  }

  if (!need_fparser)
  {
    // Simple case: all coefficients are numeric (e.g. tutorial_5sp)
    for (unsigned int r = 0; r < _n_reactions; ++r)
    {
      const std::string & expr = mech.reactions[r].rate_expression;
      auto it = coeff_map.find(expr);
      if (it != coeff_map.end())
        _k[r] = it->second;
      else
        try { _k[r] = std::stod(expr); } catch (...) {}
    }
  }
  else
  {
    // Complex case: use fparser (see evaluateCoefficients)
    _console << "MCMBoxModel: " << mech.coefficient_names.size()
             << " coefficients — using fparser for complex expressions" << std::endl;
    setupFparser(mech);
  }

  _dirty = true;
}

void
MCMBoxModel::setupFparser(const ParsedMechanism & mech)
{
  auto coeff_exprs = mech.coefficient_expressions;
  std::vector<std::string> coeff_names = mech.coefficient_names;
  std::vector<std::string> rxn_exprs(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rxn_exprs[r] = mech.reactions[r].rate_expression;

  // Convert J<N> to PHOTOJN (fparser doesn't allow < >)
  auto replace_j = [](std::string & s) {
    std::string result;
    for (size_t i = 0; i < s.size(); )
    {
      if (i + 1 < s.size() && s[i] == 'J' && s[i+1] == '<')
      {
        std::string num; i += 2;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { num += s[i]; i++; }
        result += "PHOTOJ" + num;
        if (i < s.size() && s[i] == '>') i++;
      }
      else { result += s[i]; i++; }
    }
    s = result;
  };
  for (auto & e : coeff_exprs) replace_j(e);
  for (auto & e : rxn_exprs) replace_j(e);

  // Detect photolysis J numbers
  std::set<int> j_numbers;
  pcrecpp::RE re_j("PHOTOJ([0-9]+)");
  for (auto & e : coeff_exprs)
  { int jnum; pcrecpp::StringPiece sp(e); while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum); }
  for (auto & e : rxn_exprs)
  { int jnum; pcrecpp::StringPiece sp(e); while (re_j.FindAndConsume(&sp, &jnum)) j_numbers.insert(jnum); }

  // Build name-to-index mapping for ALL possible variables.
  // This is the authoritative index registry used by both evaluateCoefficients
  // and the per-parser variable indirection arrays.
  _name_to_index["TEMP"] = 0; _name_to_index["M"] = 1;
  _name_to_index["O2"] = 2; _name_to_index["N2"] = 3;
  _name_to_index["H2O"] = 4;

  for (unsigned int i = 0; i < coeff_names.size(); ++i)
    _name_to_index[coeff_names[i]] = 5 + i;
  for (unsigned int i = 0; i < _n_species; ++i)
    _name_to_index[_species_names[i]] = 5 + coeff_names.size() + i;

  bool has_ro2_in_species = false;
  for (auto & s : _species_names)
    if (s == "RO2") { has_ro2_in_species = true; break; }
  if (!has_ro2_in_species)
    _name_to_index["RO2"] = 5 + coeff_names.size() + _n_species;
  unsigned int n_extra_vars = has_ro2_in_species ? 0 : 1;
  _j_index_start = 5 + coeff_names.size() + _n_species + n_extra_vars;

  for (auto & n : j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(n);
    _name_to_index[jname] = _j_index_start + _name_to_index.size() - _j_index_start;
  }

  unsigned int n_j = j_numbers.size();
  _func_params.resize(5 + coeff_names.size() + _n_species + n_extra_vars + n_j, 0.0);

  // ── Parse coefficients using ParseAndDeduceVariables ──
  // Each parser auto-discovers its own (small) variable list, avoiding the
  // giant Vars string that breaks fparser for mechanisms with >~400 species.
  _coeff_parsers.resize(coeff_names.size());
  _coeff_var_indices.resize(coeff_names.size());
  _coeff_local_params.resize(coeff_names.size());

  for (unsigned int i = 0; i < coeff_names.size(); ++i)
  {
    _coeff_parsers[i] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_coeff_parsers[i]);

    std::string discoveredVars;
    int nFound = 0;
    if (_coeff_parsers[i]->ParseAndDeduceVariables(coeff_exprs[i], discoveredVars, &nFound) >= 0)
      mooseError("MCMBoxModel: Bad coefficient '", coeff_names[i], "': ", coeff_exprs[i], "\n",
                 _coeff_parsers[i]->ErrorMsg());

    // Build per-parser index into _func_params
    auto & indices = _coeff_var_indices[i];
    indices.reserve(nFound);
    std::istringstream vstream(discoveredVars);
    std::string vname;
    while (std::getline(vstream, vname, ','))
    {
      auto it = _name_to_index.find(vname);
      if (it != _name_to_index.end())
        indices.push_back(it->second);
      else
      {
        mooseWarning("MCMBoxModel: coefficient '", coeff_names[i],
                     "' references undefined variable '", vname,
                     "'. The .fac file may be missing MCM standard rate "
                     "constant definitions. Use scripts/extract_mcm_k.py to "
                     "inject them from MCMv331_K.m. Auto-registering with value 0.");
        unsigned int idx = _func_params.size();
        _func_params.push_back(0.0);
        _name_to_index[vname] = idx;
        indices.push_back(idx);
      }
    }
    _coeff_local_params[i].resize(indices.size());

    if (!_disable_fpoptimizer)
      _coeff_parsers[i]->Optimize();
  }

  // ── Parse reaction expressions using ParseAndDeduceVariables ──
  _reaction_parsers.resize(_n_reactions);
  _reaction_var_indices.resize(_n_reactions);
  _reaction_local_params.resize(_n_reactions);

  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    _reaction_parsers[r] = std::make_shared<SymFunction>();
    setParserFeatureFlags(_reaction_parsers[r]);

    std::string discoveredVars;
    int nFound = 0;
    if (_reaction_parsers[r]->ParseAndDeduceVariables(rxn_exprs[r], discoveredVars, &nFound) >= 0)
      mooseError("MCMBoxModel: Bad reaction ", r, ": ", rxn_exprs[r], "\n",
                 _reaction_parsers[r]->ErrorMsg());

    auto & indices = _reaction_var_indices[r];
    indices.reserve(nFound);
    std::istringstream vstream(discoveredVars);
    std::string vname;
    while (std::getline(vstream, vname, ','))
    {
      auto it = _name_to_index.find(vname);
      if (it != _name_to_index.end())
        indices.push_back(it->second);
      else
      {
        mooseWarning("MCMBoxModel: reaction ", r,
                     " references undefined variable '", vname,
                     "'. The .fac file may be missing MCM standard rate "
                     "constant definitions. Use scripts/extract_mcm_k.py to "
                     "inject them from MCMv331_K.m. Auto-registering with value 0.");
        unsigned int idx = _func_params.size();
        _func_params.push_back(0.0);
        _name_to_index[vname] = idx;
        indices.push_back(idx);
      }
    }
    _reaction_local_params[r].resize(indices.size());

    if (!_disable_fpoptimizer)
      _reaction_parsers[r]->Optimize();
  }

  // Store photolysis parameters for SZA-based J calculation.
  // When the photolysis file is empty (BOTTOMUP mode), use the locally-detected
  // J numbers from expression scanning instead.
  if (!mech.j_numbers.empty())
  {
    _j_numbers = mech.j_numbers;
    _j_CL_vals = mech.j_CL;
    _j_CMM_vals = mech.j_CMM;
    _j_CNN_vals = mech.j_CNN;
  }
  else
  {
    _j_numbers.assign(j_numbers.begin(), j_numbers.end());
    _j_CL_vals.assign(_j_numbers.size(), 0.0);
    _j_CMM_vals.assign(_j_numbers.size(), 1.0);
    _j_CNN_vals.assign(_j_numbers.size(), 0.0);
  }
  _n_j_vars = _j_numbers.size();

  // Pre-compute J photo indices into _func_params (Per.14 — avoids string+map in evaluateCoefficients).
  // MUST be after _j_numbers assignment above (built from mech.j_numbers, not the local j_numbers from
  // expression scanning above, which may detect a different set).
  _j_photo_indices.clear();
  _j_photo_indices.reserve(_j_numbers.size());
  for (auto & jn : _j_numbers)
  {
    std::string jname = "PHOTOJ" + std::to_string(jn);
    auto it = _name_to_index.find(jname);
    if (it != _name_to_index.end())
      _j_photo_indices.push_back(it->second);
    else
      _j_photo_indices.push_back((unsigned int)-1); // sentinel: skip in evaluateCoefficients
  }
  if (_n_j_vars > 0)
    _console << "MCMBoxModel: " << _n_j_vars << " photolysis J values loaded" << std::endl;

  // Build fast pre-compiled handlers for coefficient and reaction expressions.
  // These bypass fparser tree traversal (AST walk) for common patterns.
  // Non-matching expressions get nullptr → fall back to fparser.
  _coeff_fast.resize(_coeff_parsers.size());
  for (unsigned int i = 0; i < _coeff_parsers.size(); ++i)
    _coeff_fast[i] = compileFastHandler(coeff_exprs[i], _coeff_var_indices[i]);

  _reaction_fast.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    _reaction_fast[r] = compileFastHandler(rxn_exprs[r], _reaction_var_indices[r]);
}

MCMBoxModel::FastHandler
MCMBoxModel::compileFastHandler(const std::string & expr,
                                 const std::vector<unsigned int> & var_indices) const
{
  // Pattern 1: Single variable reference — most reaction expressions (e.g. "KMT01").
  // Just read _func_params at the pre-computed index.
  if (var_indices.size() == 1 && std::regex_match(expr, std::regex("^[A-Za-z_][A-Za-z0-9_]*$")))
  {
    unsigned int idx = var_indices[0];
    return [idx](const std::vector<Real> & p) { return p[idx]; };
  }

  // Pattern 2: Simple numeric constant.
  // Patterns like "2.0e-11", "1.234"
  {
    char * end = nullptr;
    double val = std::strtod(expr.c_str(), &end);
    if (end && *end == '\0')
      return [val](const std::vector<Real> &) { return val; };
  }

  // Pattern 3: Simple Arrhenius: A*exp(B/TEMP)
  // e.g. "1.0e-11*exp(-550/TEMP)", "2.54e-12*exp(360/TEMP)"
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) { return A * std::exp(B / p[0]); };
    }
  }

  // Pattern 4: Modified Arrhenius: A*(TEMP/300)^B*exp(C/TEMP)
  // e.g. "1.44e-13*(TEMP/300)^4*exp(825/TEMP)"
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*\(TEMP/300\)\^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      double C = std::stod(m[3].str());
      return [A, B, C](const std::vector<Real> & p) {
        return A * std::pow(p[0] / 300.0, B) * std::exp(C / p[0]);
      };
    }
  }

  // Pattern 5: Power temp: A*(TEMP/300)^B
  // e.g. "1.0e-11*(TEMP/300)^2"
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*\(TEMP/300\)\^([0-9.eE+\-]+)$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) {
        return A * std::pow(p[0] / 300.0, B);
      };
    }
  }

  // Pattern 6: Arrhenius with M: A*exp(B/TEMP)*M
  // e.g. "3.28e-28*M*(TEMP/300)^-6.87" – handled by pattern 7
  // This is just A*exp(B/TEMP)*M → p[1] is M
  {
    std::smatch m;
    if (std::regex_match(expr, m, std::regex(
        R"(^([0-9.eE+\-]+)\*exp\(([0-9.eE+\-]+)/TEMP\)\*M$)")))
    {
      double A = std::stod(m[1].str());
      double B = std::stod(m[2].str());
      return [A, B](const std::vector<Real> & p) {
        return A * std::exp(B / p[0]) * p[1];
      };
    }
  }

  // No pattern matched — fall back to fparser
  return nullptr;
}

void
MCMBoxModel::evaluateCoefficients()
{
  if (_coeff_parsers.empty()) return;

  // Compute M (air density) dynamically from press/temp if press > 0
  Real M_val = _air_density;
  if (_press > 0.0)
  {
    // AtChem2 calcAirDensity: M = 1e-6 * NA/R * (press_pa / temp)
    // press_mbar → press_pa: * 100
    // NA = 6.02214129e23, R = 8.3144621
    constexpr Real NA_over_R = 6.02214129e23 / 8.3144621;
    M_val = 1.0e-6 * NA_over_R * (_press * 100.0 / _temperature);
  }

  // Update ppb conversion factor from the actual air density used
  _ppb_to_molec = (_units_ppb) ? M_val / 1.0e9 : 1.0;

  // Compute H2O dynamically from rh/temp/press if rh >= 0
  Real H2O_val = _water_vapor;
  if (_rh >= 0.0)
  {
    // Vaisala 2013 convertRHtoH2O (AtChem2 atmosphereFunctions.f90)
    Real temp_c = _temperature - 273.15;
    // Eq.6: water vapour saturation pressure (mbar)
    Real wvp = (_rh / 100.0) * 6.116441 * std::pow(10.0, (7.591386 * temp_c) / (temp_c + 240.7263));
    // Eq.18: volume of water vapour per volume of dry air (parts per unit)
    Real press_mbar = (_press > 0.0) ? _press : 1013.25;
    Real h2o_ppu = wvp / (press_mbar - wvp);
    // Convert ppu to molecule/cm³
    H2O_val = h2o_ppu * M_val;
  }

  _func_params[0] = _temperature;
  _func_params[1] = M_val;
  _func_params[2] = 0.21 * M_val;  // O2 volume fraction
  _func_params[3] = 0.78 * M_val;  // N2 volume fraction
  _func_params[4] = H2O_val;

  // Compute photolysis J values into _func_params for fparser evaluation.
  // MCM_SZA: J = l * cosx^m * exp(-n * secx) * JFAC * ROOF
  // BOTTOMUP: J = ∫ QY*CS*Flux dλ  (pre-computed by BottomUpJIntegrator)
  if (_photolysis_method == BOTTOMUP && _bottomup_integrator)
  {
    // BottomUp J-values depend only on T and P (constant for chamber).
    // Cache to avoid recomputing full numerical integration at every Newton iteration
    // — matches F0AM behavior where J-values are computed once at startup.
    const Real T_cur = _temperature;
    const Real P_cur = _press > 0 ? _press : 1013.25;
    if (!_bottomup_j_valid ||
        std::abs(T_cur - _cached_bottomup_T) > 1.0e-6 ||
        std::abs(P_cur - _cached_bottomup_P) > 1.0e-6)
    {
      _cached_bottomup_j = _bottomup_integrator->computeAllJ(T_cur, P_cur);
      _cached_bottomup_T = T_cur;
      _cached_bottomup_P = P_cur;
      _bottomup_j_valid = true;
    }
    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    for (size_t i = 0; i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      std::string jname = "J" + std::to_string(jn);
      auto it = _cached_bottomup_j.find(jname);
      _func_params[idx] = (it != _cached_bottomup_j.end()) ? it->second * _jfac * roof_factor : 0.0;
    }
  }
  else if (!_j_CL_vals.empty())
  {
    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    Real cosx = calculateCosSZA(_t);
    Real secx = _solar_secx;  // cached by calculateCosSZA

    // Per.14: use pre-computed _j_photo_indices instead of string+map lookup.
    // Sentinel value (unsigned int)-1 means J number was not registered; skip it.
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      if (cosx > 1.0e-10)
        _func_params[idx] = _j_CL_vals[i] * std::pow(cosx, _j_CMM_vals[i])
                          * std::exp(-_j_CNN_vals[i] * secx) * _jfac * roof_factor;
      else
        _func_params[idx] = 0.0;
    }
  }

  // Auto-calibrate J values from observed data (F0AM jcorr)
  if (_jcalibrator && _jcalibrator->hasObservedData())
  {
    // Collect parameterized J values (without _jfac/roof_factor — those are
    // separate corrections applied before calibration)
    std::map<unsigned int, Real> param_J;
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      param_J[jn] = _func_params[idx];
    }
    _jcalibrator->calibrate(param_J);
    _jcalibrator->applyTo(param_J);
    // Write calibrated J values back to _func_params
    for (size_t i = 0; i < _j_CL_vals.size() && i < _j_photo_indices.size(); ++i)
    {
      unsigned int idx = _j_photo_indices[i];
      if (idx == (unsigned int)-1) continue;
      unsigned int jn = (_j_numbers.size() > i) ? _j_numbers[i] : (unsigned int)(i + 1);
      auto it = param_J.find(jn);
      if (it != param_J.end())
        _func_params[idx] = it->second;
    }
  }

  // Per.19: before evaluating any parser, populate ALL _func_params slots.
  // Species and RO2 must be set before coefficient evaluation because some
  // coefficient expressions reference species concentrations directly
  // (e.g. "KMT06" references "C5H8").
  for (unsigned int s = 0; s < _n_species; ++s)
  {
    auto it = _name_to_index.find(_species_names[s]);
    if (it != _name_to_index.end())
      _func_params[it->second] = (_cached_C.empty() ? 0.0 :
          (s < _cached_C.size() ? _cached_C[s] : 0.0));
  }
  auto it_ro2 = _name_to_index.find("RO2");
  if (it_ro2 != _name_to_index.end())
  {
    Real ro2_sum = 0.0;
    for (auto idx : _ro2_indices)
      ro2_sum += (_cached_C.empty() ? 0.0 :
          (idx < _cached_C.size() ? _cached_C[idx] : 0.0));
    _func_params[it_ro2->second] = ro2_sum;
  }

  // Evaluate coefficients in topological order.
  // Use fast pre-compiled handlers when available; fall back to fparser for
  // complex expressions (fall-off formulas, RO2-dependent, etc.).
  unsigned int n_coeff = _coeff_parsers.size();
  for (unsigned int i = 0; i < n_coeff; ++i)
  {
    Real val;
    if (_coeff_fast[i])
      val = _coeff_fast[i](_func_params);
    else
    {
      const auto & indices = _coeff_var_indices[i];
      auto & local = _coeff_local_params[i];
      for (size_t j = 0; j < indices.size(); ++j)
        local[j] = _func_params[indices[j]];
      val = evaluate(_coeff_parsers[i], local);
    }
    if (std::isnan(val) || std::isinf(val)) val = 0.0;
    _func_params[5 + i] = val;
  }

  // Evaluate reaction rate expressions → _k
  // ~93% of reaction expressions are single-variable references (e.g. "KMT01")
  // handled by the fast path; complex expressions fall back to fparser.
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real val;
    if (_reaction_fast[r])
      val = _reaction_fast[r](_func_params);
    else
    {
      const auto & indices = _reaction_var_indices[r];
      auto & local = _reaction_local_params[r];
      for (size_t j = 0; j < indices.size(); ++j)
        local[j] = _func_params[indices[j]];
      val = evaluate(_reaction_parsers[r], local);
    }
    _k[r] = (std::isnan(val) || std::isinf(val)) ? 0.0 : val;
  }
}

unsigned int
MCMBoxModel::computeDayOfYear() const
{
  unsigned int days_in_months[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0)
    days_in_months[1] = 29;
  // AtChem2 convention: Jan 1 = 0, Jan 2 = 1, ... (0-based day of year).
  // This matches AtChem2's date_mod::calcDayOfYear:
  //   result = sum(monthList(1:month-1)) + day - 1
  unsigned int doy = 0;
  for (unsigned int m = 0; m < (unsigned int)(_month - 1); ++m)
    doy += days_in_months[m];
  doy += _day - 1;
  return doy;
}

Real
MCMBoxModel::calculateCosSZA(Real t) const
{
  constexpr Real pi = 3.14159265358979323846;
  unsigned int doy = computeDayOfYear();
  unsigned int days_in_year =
      ((_year % 4 == 0 && _year % 100 != 0) || _year % 400 == 0) ? 366 : 365;
  Real theta = 2.0 * pi * (Real)doy / (Real)days_in_year;
  Real dec = 0.006918 - 0.399912 * cos(theta) + 0.070257 * sin(theta) -
             0.006758 * cos(2.0 * theta) + 0.000907 * sin(2.0 * theta) -
             0.002697 * cos(3.0 * theta) + 0.001480 * sin(3.0 * theta);
  Real eqt = 0.000075 + 0.001868 * cos(theta) - 0.032077 * sin(theta) -
             0.014615 * cos(2.0 * theta) - 0.040849 * sin(2.0 * theta);
  Real current_frac_hour = std::fmod(t / 3600.0, 24.0);
  Real lat_rad = _lat * pi / 180.0;
  Real lha = pi * ((current_frac_hour / 12.0) - (1.0 + _lon / 180.0)) + eqt;
  Real cosx = cos(lha) * cos(lat_rad) * cos(dec) + sin(lat_rad) * sin(dec);

  // Cache all solar parameters for MCMSolarPostprocessor
  _declination = dec;  // store for DEC postprocessor
  _eot = eqt;
  _solar_lha = lha;
  _solar_sinld = sin(lat_rad) * sin(dec);
  _solar_cosld = cos(lat_rad) * cos(dec);
  _solar_eqt = eqt;

  if (cosx <= 0.0) cosx = 0.0;
  Real secx = (cosx > 1.0e-10) ? (1.0 / cosx) : 1.0e2;
  _solar_cosx = cosx;
  _solar_secx = secx;
  return cosx;
}

// -- Constrained species --
void
MCMBoxModel::setConstrainedSpecies(const std::vector<std::string> & names)
{
  _constrained_set.clear();
  for (auto & name : names)
    for (unsigned int i = 0; i < _n_species; ++i)
      if (_species_names[i] == name)
        _constrained_set.insert(i);
  _constrained_values.assign(_constrained_set.size(), 0.0);
}

void
MCMBoxModel::updateConstrainedValues(const std::vector<Real> & values)
{
  if (values.size() == _constrained_set.size())
    _constrained_values = values;
}

void
MCMBoxModel::computeDCdtFull(const std::vector<Real> & C_full, std::vector<Real> & dC) const
{
  computeDCdt(C_full, dC);
  // Zero out constrained species derivatives
  for (auto idx : _constrained_set)
    dC[idx] = 0.0;
}

// -- Reaction rate queries --
Real
MCMBoxModel::reactionRate(unsigned int r, const std::vector<Real> & C) const
{
  if (r >= _n_reactions || C.size() != _n_species) return 0.0;
  Real c0 = (_iG[r][0] >= 0) ? C[_iG[r][0]] : 1.0;
  Real c1 = (_iG[r][1] >= 0) ? C[_iG[r][1]] : 1.0;
  Real c2 = (_iG[r][2] >= 0) ? C[_iG[r][2]] : 1.0;
  return _k[r] * c0 * c1 * c2;
}

Real
MCMBoxModel::speciesReactionRate(unsigned int s, unsigned int r, const std::vector<Real> & C) const
{
  return _stoich.get(r, s) * reactionRate(r, C);
}

void
MCMBoxModel::allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const
{
  rates.resize(_n_reactions);
  for (unsigned int r = 0; r < _n_reactions; ++r)
    rates[r] = reactionRate(r, C);
}

Real
MCMBoxModel::speciesLossRate(unsigned int s, const std::vector<Real> & C) const
{
  Real total = 0.0;
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real coeff = _stoich.get(r, s);
    if (coeff < 0.0)
      total += (-coeff) * reactionRate(r, C);
  }
  return total;
}

Real
MCMBoxModel::speciesProductionRate(unsigned int s, const std::vector<Real> & C) const
{
  Real total = 0.0;
  for (unsigned int r = 0; r < _n_reactions; ++r)
  {
    Real coeff = _stoich.get(r, s);
    if (coeff > 0.0)
      total += coeff * reactionRate(r, C);
  }
  return total;
}

// -- Photolysis --
void
MCMBoxModel::enableHybridPhotolysis(const std::string & table_dir)
{
  _photolysis_method = HYBRID;
  _hybrid_reader = std::make_unique<HybridJTableReader>(table_dir);
}

// -- BottomUp photolysis (F0AM chamber mode) --
void
MCMBoxModel::loadBottomUpData(const std::string & data_dir, const std::string & flux_file)
{
  _bottomup_integrator = std::make_unique<BottomUpJIntegrator>(data_dir);
  _bottomup_integrator->loadLampFlux(flux_file);
  _bottomup_integrator->loadReactionMap("bottomup_jmap.dat");
}

// -- Solar cycle (Madronich 1993) --
void
MCMBoxModel::setSolarCycle(Real lat, Real lon, int day, int month, int year)
{
  _lat = lat * M_PI / 180.0;
  _lon = lon * M_PI / 180.0;
  _day = day;
  _month = month;
  _year = year;

  unsigned int doy = computeDayOfYear();
  unsigned int days_in_year =
      ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) ? 366 : 365;
  Real theta = 2.0 * M_PI * (Real)doy / (Real)days_in_year;
  Real B = 2.0 * M_PI * ((Real)doy - 1.0) / (Real)days_in_year;

  _declination = 0.006918 - 0.399912 * std::cos(theta) + 0.070257 * std::sin(theta) -
                 0.006758 * std::cos(2.0 * theta) + 0.000907 * std::sin(2.0 * theta) -
                 0.002697 * std::cos(3.0 * theta) + 0.001480 * std::sin(3.0 * theta);
  _eot = 0.165 * std::sin(2.0 * B) - 0.126 * std::cos(B) - 0.025 * std::sin(B);
}

Real
MCMBoxModel::cosSZA(Real seconds) const
{
  Real hour = seconds / 3600.0;
  Real ha = (hour - 12.0 + _eot + _lon*12.0/M_PI) * M_PI / 12.0;
  return std::sin(_lat)*std::sin(_declination)
       + std::cos(_lat)*std::cos(_declination)*std::cos(ha);
}

// -- Dilution --
void
MCMBoxModel::setDilution(Real kdil, const std::vector<Real> & bg)
{ _kdil = kdil; _conc_bkgd = bg; }

Real
MCMBoxModel::getRO2Sum(const std::vector<Real> & C) const
{
  Real sum = 0.0;
  for (auto idx : _ro2_indices)
    if (idx < C.size())
      sum += C[idx];
  return sum;
}

Real
MCMBoxModel::getJValue(unsigned int j_number) const
{
  // First try the fparser parameter buffer (fast path for J numbers referenced
  // in reaction/coefficient expressions — these are pre-computed).
  std::string jname = "PHOTOJ" + std::to_string(j_number);
  auto it = _name_to_index.find(jname);
  if (it != _name_to_index.end() && it->second < _func_params.size())
    return _func_params[it->second];

  // Fallback: compute J directly from MCM SZA formula for J numbers that have
  // photolysis parameters (CL/CMM/CNN) but aren't referenced in any reaction.
  // These are in the photolysis-rates file and AtChem2 outputs them, but they
  // don't participate in chemistry.
  for (size_t i = 0; i < _j_numbers.size(); ++i)
  {
    if (_j_numbers[i] == j_number)
    {
      const Real roof_factor = _roof_open ? 1.0 : 0.0;
      Real cosx = _solar_cosx;  // cached during evaluateCoefficients
      if (cosx > 1.0e-10)
        return _j_CL_vals[i] * std::pow(cosx, _j_CMM_vals[i])
               * std::exp(-_j_CNN_vals[i] / cosx) * _jfac * roof_factor;
      return 0.0;
    }
  }
  return 0.0;
}

void
MCMBoxModel::setGaussianDispersion(Real tgauss, const std::vector<Real> & conc_bkgd,
                                    Real t_start)
{
  _tgauss = tgauss;
  _conc_bkgd = conc_bkgd;
  _t_start_dil = t_start;
  _use_gaussian = true;
  _kdil = 0.0; // disable first-order dilution
  _console << "MCMBoxModel: Gaussian dispersion enabled, tgauss=" << tgauss << "s"
           << std::endl;
}

void
MCMBoxModel::computeDCdtWithDilution(const std::vector<Real> & C, std::vector<Real> & dC) const
{
  computeDCdt(C, dC);
  if (_use_gaussian)
  {
    // F0AM Gaussian dispersion: dilrate = -1/(tgauss + 2*(t+t_start))
    Real denom = _tgauss + 2.0 * (_t + _t_start_dil);
    Real dilrate = (denom > 1.0e-30) ? 1.0 / denom : 1.0e30;
    for (unsigned int i = 0; i < _n_species; ++i)
      dC[i] -= dilrate * (C[i] - (_conc_bkgd.empty() ? 0.0 : _conc_bkgd[i]));
  }
  else
  {
    // First-order dilution (AtChem2 DILUTE)
    for (unsigned int i = 0; i < _n_species; ++i)
      dC[i] -= _kdil * (C[i] - (_conc_bkgd.empty() ? 0.0 : _conc_bkgd[i]));
  }
}

// -- Convergence --
Real
MCMBoxModel::checkConvergence(const std::vector<Real> & Cp, const std::vector<Real> & Cc,
                               const std::vector<std::string> & cs) const
{
  Real mx = 0.0;
  for (unsigned int i = 0; i < _n_species; ++i)
  {
    if (!cs.empty()) {
      bool ok = false;
      for (auto & s : cs) if (s == _species_names[i]) { ok = true; break; }
      if (!ok) continue;
    }
    Real d = std::abs(Cc[i] - Cp[i]) / std::max(std::abs(Cp[i]), 1e-30);
    if (d > mx) mx = d;
  }
  return mx;
}

// ===== PETSc TS standalone integrator =====

void
MCMBoxModel::execute()
{
  // PETSc TS mode: integrate chemistry at TIMESTEP_END, bypassing MOOSE's solver.
  // No kernels (ODETimeDerivative / ChemistryODEKernel) exist in this mode —
  // MOOSE's solve is a no-op (residual=0), and this execute() sets the solution.
  if (!_use_petsc_ts || _n_species == 0)
    return;

  // Get current and previous times from FEProblemBase
  FEProblemBase & fe_problem = static_cast<FEProblemBase &>(_subproblem);
  PetscReal t_end = fe_problem.time();
  PetscReal dt = fe_problem.dt();
  PetscReal t_start = t_end - dt;

  // Skip initial call at t=0 (no step yet)
  if (t_start < 0.0 || dt <= 0.0)
    return;

  // Build concentration vector from ScalarVariable values (source of truth
  // for ICs and persistent state — the NL solution vector may not be
  // populated until after the first solve).
  NonlinearSystemBase & nl = fe_problem.getNonlinearSystemBase(0);
  {
    PetscScalar *x_arr;
    VecGetArray(_ts_X, &x_arr);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      x_arr[i] = sv.sln()[0];
    }
    VecRestoreArray(_ts_X, &x_arr);
  }

  // Run PETSc TS integration from t_start to t_end
  runPETScStep(t_start, t_end);

  // Write TS solution to BOTH the libMesh solution and the variable cache.
  // The CSV output reads from sln() → _u.  sv.setValue() updates _u.
  // Also update sys.solution for consistency (used by next step's predictor).
  {
    PetscScalar *x_arr;
    VecGetArray(_ts_X, &x_arr);
    NumericVector<Number> & sys_sol = *nl.system().solution;
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      dof_id_type dof = sv.dofIndices()[0];
      // Write to ALL locations — sln() reads from _u, CSV writes from sln()
      sv.setValue(0, x_arr[i]);    // _u[i] = val (for sln() → output)
      sys_sol.set(dof, x_arr[i]);  // sys.solution (for next predictor)
      nl.solution().set(dof, x_arr[i]);  // nl.solution (backup)
    }
    sys_sol.close();
    nl.solution().close();
    VecRestoreArray(_ts_X, &x_arr);
    *const_cast<NumericVector<Number>*>(nl.currentSolution()) = nl.solution();
  }

  // Diagnostics removed — verify using CSV output comparison

  _console << "MCMBoxModel: TS step t=[" << t_start << "," << t_end
           << "] dt=" << dt << " completed" << std::endl;
}

void
MCMBoxModel::setupPETScTS()
{
  PetscErrorCode ierr;

  // Load TS parameters from input
  _ts_type = std::string(getParam<MooseEnum>("petsc_ts_type"));
  _ts_rtol = getParam<Real>("petsc_ts_rtol");
  _ts_atol = getParam<Real>("petsc_ts_atol");

  // Helper macro for error checking without PetscCheck (which returns a value)
  // and without requiring MPI_Comm (use PETSC_COMM_SELF for sequential TS).
#define PETSC_TRY(expr) do { ierr = (expr); if (ierr) mooseError("PETSc error ", ierr, " at ", __FILE__, ":", __LINE__); } while(0)

  PETSC_TRY(TSCreate(PETSC_COMM_SELF, &_ts));
  PETSC_TRY(TSSetProblemType(_ts, TS_NONLINEAR));
  PETSC_TRY(TSSetRHSFunction(_ts, nullptr, tsRHSFunction, this));
  PETSC_TRY(VecCreateSeq(PETSC_COMM_SELF, _n_species, &_ts_X));

  // Use dense matrix for Jacobian — 610×610 ~ 3MB, avoids sparse preallocation issues.
  // PETSc's BDF/ARKIMEX with SuperLU_DIST or LU handles dense efficiently at this size.
  PetscInt n = static_cast<PetscInt>(_n_species);
  PETSC_TRY(MatCreateSeqDense(PETSC_COMM_SELF, n, n, nullptr, &_ts_J));
  PETSC_TRY(MatSetFromOptions(_ts_J));

  PETSC_TRY(TSSetRHSJacobian(_ts, _ts_J, _ts_J, tsRHSJacobian, this));
  PETSC_TRY(TSSetType(_ts, _ts_type.c_str()));

  // For ARKIMEX, configure as fully implicit with 4th order scheme
  {
    TSType actual_type;
    PETSC_TRY(TSGetType(_ts, &actual_type));
    if (std::string(actual_type) == "arkimex")
    {
      PETSC_TRY(TSARKIMEXSetFullyImplicit(_ts, PETSC_TRUE));
      PETSC_TRY(TSARKIMEXSetType(_ts, TSARKIMEX4));
    }
  }

  PETSC_TRY(TSSetTolerances(_ts, _ts_atol, nullptr, _ts_rtol, nullptr));
  // Allow unlimited SNES failures (retry with smaller step), matching extchem.c
  PETSC_TRY(TSSetMaxSNESFailures(_ts, -1));
  PETSC_TRY(TSSetSolution(_ts, _ts_X));
  PETSC_TRY(TSSetFromOptions(_ts));

#undef PETSC_TRY

  _console << "MCMBoxModel: PETSc TS initialized (" << _n_species << " species, type="
           << _ts_type << ", rtol=" << _ts_rtol << ", atol=" << _ts_atol << ")" << std::endl;
}

void
MCMBoxModel::runPETScStep(PetscReal t0, PetscReal t1)
{
  PetscErrorCode ierr;

#define PETSC_TRY(expr) do { ierr = (expr); if (ierr) mooseError("PETSc error ", ierr, " at ", __FILE__, ":", __LINE__); } while(0)

  PETSC_TRY(TSSetTime(_ts, t0));
  PETSC_TRY(TSSetMaxTime(_ts, t1));

  // Set initial step size (very small for stiff chemistry)
  PetscReal step0 = std::max((t1 - t0) * 1.0e-6, 1.0e-10);
  PETSC_TRY(TSSetTimeStep(_ts, step0));
  // Also set adaptive step limits
  TSAdapt adapt;
  PETSC_TRY(TSGetAdapt(_ts, &adapt));
  PETSC_TRY(TSAdaptSetStepLimits(adapt, 1.0e-12, (t1 - t0)));
  PETSC_TRY(TSSetSolution(_ts, _ts_X));

  // Evaluate rate coefficients at the midpoint of this interval.
  _t = 0.5 * (t0 + t1);
  if (!_coeff_parsers.empty())
    evaluateCoefficients();

  // NOTE: _ts_X was filled from sv.sln() in execute() above.
  // _k is evaluated at interval midpoint for constant photolysis.

  // Run the integrator
  PETSC_TRY(TSSolve(_ts, _ts_X));

  // Get integration statistics
  PetscInt steps;
  TSConvergedReason reason;
  PETSC_TRY(TSGetStepNumber(_ts, &steps));
  PETSC_TRY(TSGetConvergedReason(_ts, &reason));

#undef PETSC_TRY

  // Diagnostic: print key species after TS integration
  {
    const PetscScalar *dbg_arr;
    VecGetArrayRead(_ts_X, &dbg_arr);
    auto c5h8 = _name_to_index.find("C5H8");
    auto no2 = _name_to_index.find("NO2");
    auto h2o2 = _name_to_index.find("H2O2");
    auto o3 = _name_to_index.find("O3");
    auto oh = _name_to_index.find("OH");
    _console << "MCMBoxModel: TS [" << t0 << "," << t1 << "] C5H8="
             << (c5h8 != _name_to_index.end() ? dbg_arr[c5h8->second] : -1.0)
             << " NO2=" << (no2 != _name_to_index.end() ? dbg_arr[no2->second] : -1.0)
             << " O3=" << (o3 != _name_to_index.end() ? dbg_arr[o3->second] : -1.0)
             << " OH=" << (oh != _name_to_index.end() ? dbg_arr[oh->second] : -1.0)
             << std::endl;
    VecRestoreArrayRead(_ts_X, &dbg_arr);
  }

  _console << "MCMBoxModel: TS step [" << t0 << "," << t1 << "] "
           << steps << " internal steps, reason=" << reason << std::endl;
}

// static
PetscErrorCode
MCMBoxModel::tsRHSFunction(TS ts, PetscReal t, Vec C, Vec F, void *ctx)
{
  MCMBoxModel *model = static_cast<MCMBoxModel *>(ctx);
  const PetscScalar *c_arr;
  PetscScalar *f_arr;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(C, &c_arr));
  PetscCall(VecGetArray(F, &f_arr));

  // NOTE: evaluateCoefficients() is NOT called here — it's called once per
  // MOOSE timestep (in runPETScStep) with the midpoint time.  Calling it
  // on every TS internal step would modify shared state (_func_params, _k)
  // in ways that break the TS integrator's assumptions.  For the short
  // MOOSE timesteps (dt=100s), SZA-driven photolysis changes are negligible.

  // Build std::vector around PETSc array (no copy for 610 elements - acceptable)
  std::vector<Real> C_vec(c_arr, c_arr + model->_n_species);
  std::vector<Real> dC_vec;
  model->computeDCdt(C_vec, dC_vec);

  // NOTE: computeDCdt uses _k from the most recent evaluateCoefficients()
  // call in runPETScStep().  _k is constant throughout the TS interval.

  for (unsigned int i = 0; i < model->_n_species; ++i)
    f_arr[i] = dC_vec[i];

  PetscCall(VecRestoreArrayRead(C, &c_arr));
  PetscCall(VecRestoreArray(F, &f_arr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// static
PetscErrorCode
MCMBoxModel::tsRHSJacobian(TS ts, PetscReal t, Vec C, Mat Amat, Mat Pmat, void *ctx)
{
  MCMBoxModel *model = static_cast<MCMBoxModel *>(ctx);
  const PetscScalar *c_arr;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(C, &c_arr));

  // Build concentration vector
  std::vector<Real> C_vec(c_arr, c_arr + model->_n_species);

  // Compute Jacobian triplets
  std::vector<std::tuple<unsigned int, unsigned int, Real>> J;
  model->computeJacobianTriplets(C_vec, J);

  // Insert into PETSc matrix — use ADD_VALUES because multiple reactions
  // contribute to the same (row,col) pair, and the triplets from
  // computeJacobianTriplets are per-reaction contributions.
  PetscCall(MatZeroEntries(Pmat));
  for (const auto & [row, col, val] : J)
  {
    PetscInt irow = static_cast<PetscInt>(row);
    PetscInt icol = static_cast<PetscInt>(col);
    PetscScalar pval = val;
    PetscCall(MatSetValues(Pmat, 1, &irow, 1, &icol, &pval, ADD_VALUES));
  }

  PetscCall(MatAssemblyBegin(Pmat, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Pmat, MAT_FINAL_ASSEMBLY));
  if (Amat != Pmat)
  {
    PetscCall(MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY));
  }

  PetscCall(VecRestoreArrayRead(C, &c_arr));
  PetscFunctionReturn(PETSC_SUCCESS);
}
