//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMBoxModel.h"
#include "MCMRuntimeMechanism.h"
#include "MCMFacsimileParser.h"
#include "BottomUpJIntegrator.h"
#include "MooseVariableScalar.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "pcrecpp.h"
#ifdef KPP_ENABLED
#include "KppBoxIntegrator.h"
#endif
#include <algorithm>
#include <cmath>
#include <regex>
#include <sstream>
#include <utility>

registerMooseObject("AtmosphericChemistryApp", MCMBoxModel);

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
  params.addParam<bool>("roof_open", true,
      "Roof (chamber cover) open. false forces photolysis rates to zero.");
  params.addParam<Real>("default_ic", 0.0,
      "Default initial concentration (molec/cm³) for species without explicit ICs.\n"
      "Set to a small positive value (e.g. 1e6) to prevent Jacobian singularities\n"
      "when many species start at zero. 0.0 = no default (original F0AM behavior).");
  params.addParam<bool>("use_limiting_reagent", false,
      "Enable F0AM-style limiting-reagent formulation for RO2+RO2 termination "
      "reactions: rate = k * min([RO2_i], [RO2]) instead of k * [RO2_i] * [RO2]. "
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

  MooseEnum integrator_enum("moose_implicit petsc_ts", "moose_implicit");
  params.addParam<MooseEnum>("box_solver_mode", integrator_enum,
      "Box ODE solver mode: 'moose_implicit' (default, through MOOSE's Newton solver) "
      "or 'petsc_ts' (PETSc TS, bypasses MOOSE solver for chemistry)");

  MooseEnum ts_type_enum("bdf arkimex eimex rosw mimex beuler cn rk theta ssp sundials", "bdf");
  params.addParam<MooseEnum>("solver_type", ts_type_enum,
      "ODE solver type: 'bdf' (default), 'arkimex', or 'sundials' (CVODE).");
  params.addParam<Real>("solver_rtol", 1e-6,
      "Relative tolerance for the ODE solver's adaptive integrator.");
  params.addParam<Real>("solver_atol", 1e-10,
      "Absolute tolerance for the ODE solver's adaptive integrator.");

  // --- Parameters forwarded from AtmosphericChemistryAction ---
  MooseEnum solver_enum(
      "moose_implicit petsc_ts sundials kpp_rosenbrock kpp_sdirk kpp_runge_kutta",
      "petsc_ts");
  params.addParam<MooseEnum>("chem_solver", solver_enum,
      "Chemical solver backend.");
  MooseEnum mech_fmt("MCM_FACSIMILE KPP", "MCM_FACSIMILE");
  params.addParam<MooseEnum>("mechanism_format", mech_fmt,
      "Mechanism format (forwarded from AtmosphericChemistryAction).");

  params.addClassDescription(
      "Centralized box model UserObject for atmospheric chemistry ODE systems.");
  return params;
}

MCMBoxModel::MCMBoxModel(const InputParameters & params)
  : GeneralUserObject(params),
    FunctionParserUtils<false>(params),
    _mechanism(nullptr),
    _n_species(0), _n_reactions(0),
    _units_ppb(getParam<MooseEnum>("units") == "ppb"),
    _ppb_to_molec(1.0),
    _lat(getParam<Real>("latitude")),
    _lon(getParam<Real>("longitude")),
    _day((int)getParam<unsigned int>("day")),
    _month((int)getParam<unsigned int>("month")),
    _year((int)getParam<unsigned int>("year")),
    _kdil(0.0),
    _tgauss(0.0),
    _t_start_dil(0.0),
    _use_gaussian(false),
    _roof_open(getParam<bool>("roof_open")),
    _jfac(getParam<Real>("jfac")),
    _t(0.0),
    _bottomup_j_valid(false),
    _use_box_solver(false)  // derived in constructor body
{
  // ---- Derive solver selection from chem_solver parameter ----
  _chem_solver = std::string(getParam<MooseEnum>("chem_solver"));
  _use_sundials = (_chem_solver == "sundials");
  _use_kpp = (_chem_solver == "kpp_rosenbrock" ||
              _chem_solver == "kpp_sdirk" ||
              _chem_solver == "kpp_runge_kutta");
  _use_petsc_ts = (_chem_solver == "petsc_ts");

  // All non-moose_implicit solvers are self-driven modes
  if (_use_sundials || _use_kpp || _use_petsc_ts)
    _use_box_solver = true;

  // KPP solvers need mechanism_format=KPP (enforced by Action, double-check here)
  if (_use_kpp)
  {
#ifndef KPP_ENABLED
    mooseError("MCMBoxModel: chem_solver=", _chem_solver,
               " requires KPP support (KPP_ENABLED). "
               "Rebuild the atmospheric_chemistry module with KPP support.");
#endif
  }

  // Read solver type (only meaningful for petsc_ts)
  _solver_type = std::string(getParam<MooseEnum>("solver_type"));

  // Read solver tolerances ONCE here so all solver paths
  // share the same parameter source.
  _solver_rtol = getParam<Real>("solver_rtol");
  _solver_atol = getParam<Real>("solver_atol");

  // Store photolysis settings for use when mechanism is created in initialize()
  // (The mechanism is not yet created — we cache the config and apply it later)
  _photolysis_scheme = std::string(getParam<MooseEnum>("photolysis_scheme"));
  if (_photolysis_scheme == "HYBRID")
  {
    _hybrid_table_dir = getParam<std::string>("hybrid_table_dir");
    if (_hybrid_table_dir.empty())
      mooseError("MCMBoxModel: hybrid_table_dir is required when photolysis_scheme=HYBRID");
    if (!_hybrid_table_dir.empty() && _hybrid_table_dir[0] == '/')
      mooseError("MCMBoxModel: hybrid_table_dir must be relative, got absolute: ", _hybrid_table_dir);
  }
  else if (_photolysis_scheme == "BOTTOMUP")
  {
    _bottomup_data_dir = getParam<std::string>("bottomup_data_dir");
    _lamp_flux_file = getParam<std::string>("lamp_flux_file");
    if (_lamp_flux_file.empty())
      mooseError("MCMBoxModel: lamp_flux_file is required when photolysis_scheme=BOTTOMUP");
  }

  // Create the box integrator strategy — must be done in the constructor
  // because ChemistryODEKernel accesses getIntegrator() during its own
  // construction (before initialize() is called).
  if (_use_sundials)
    _integrator = std::make_unique<SundialsBoxIntegrator>(
        *this, _app, _solver_rtol, _solver_atol);
  else if (_use_kpp)
  {
#ifdef KPP_ENABLED
    std::string mech_path = getParam<std::string>("mechanism_file");
    _integrator = std::make_unique<KppBoxIntegrator>(
        _app, mech_path, _solver_rtol, _solver_atol, _chem_solver);
#else
    mooseError("MCMBoxModel: chem_solver=", _chem_solver,
               " requires KPP support. Rebuild the atmospheric_chemistry module with "
               "KPP support.");
#endif
  }
  else if (_use_box_solver)
    _integrator = std::make_unique<PetscTSIntegrator>(*this);
  else
    _integrator = std::make_unique<MooseImplicitIntegrator>(*this);
}

MCMBoxModel::~MCMBoxModel()
{
  // Clean up PETSc TS objects (safe if never initialized — pointers are null).
  // The sundials path never creates PETSc TS/Vec/Mat, so these will all be null
  // in that case and the cleanup is a no-op.
  if (_ts_J) { CHKERRABORT(PETSC_COMM_SELF, MatDestroy(&_ts_J)); _ts_J = nullptr; }
  if (_ts_X) { CHKERRABORT(PETSC_COMM_SELF, VecDestroy(&_ts_X)); _ts_X = nullptr; }
  if (_ts)   { CHKERRABORT(PETSC_COMM_SELF, TSDestroy(&_ts));   _ts = nullptr; }
}

void
MCMBoxModel::initialize()
{
  // Parse .fac file only once
  if (_n_species > 0) return;

  std::string mech_file = getParam<std::string>("mechanism_file");

  // KPP mode: mechanism is handled by KppBoxIntegrator (.so).
  // Get species info from the ScalarVariables set up by the Action.
  if (_use_kpp)
  {
    _console << "MCMBoxModel: KPP mode — using " << mech_file << std::endl;
#ifdef KPP_ENABLED
    const auto * kpp_ptr = static_cast<const KppBoxIntegrator *>(_integrator.get());
    _species_names = kpp_ptr->speciesNames();
#endif
    _n_species = _species_names.size();
    _n_reactions = 0;
    _console << "MCMBoxModel: Loaded " << _n_species << " species from KPP mechanism"
             << std::endl;

    if (_photolysis_scheme == "BOTTOMUP")
    {
      _kpp_bottomup_integrator = std::make_unique<BottomUpJIntegrator>(_bottomup_data_dir);
      _kpp_bottomup_integrator->loadLampFlux(_lamp_flux_file);
      _kpp_bottomup_integrator->loadReactionMap("bottomup_jmap.dat");
      _kpp_bottomup_j_valid = false;
    }

    // Apply default IC (if > 0) to uninitialized species.
    {
      Real default_ic = getParam<Real>("default_ic");
      if (default_ic > 0.0)
        for (const auto & sp_name : _species_names)
        {
          MooseVariableScalar & sv = _subproblem.getScalarVariable(0, sp_name);
          if (sv.sln()[0] == 0.0)
            sv.setValue(0, default_ic);
        }
    }

    // Save initial concentrations for self-driven integrators.
    _initial_conc.resize(_n_species, 0.0);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      _initial_conc[i] = sv.sln()[0];
    }
    return;
  }

  if (!mech_file.empty())
  {
    std::string photo_file = getParam<std::string>("photolysis_file");
    _console << "MCMBoxModel: parsing " << mech_file << "..." << std::endl;
    MCMFacsimileParser parser;
    ParsedMechanism mech = parser.parse(mech_file, photo_file);
    _console << "MCMBoxModel: parsing complete." << std::endl;
    bool use_lr = getParam<bool>("use_limiting_reagent");

    // Create MCMRuntimeMechanism — all chemical computation is delegated to it
    {
      MooseEnum stoich_fmt_enum = getParam<MooseEnum>("stoich_format");
      StoichMatrix::Format stoich_fmt = StoichMatrix::CSR;
      if (stoich_fmt_enum == "COO")        stoich_fmt = StoichMatrix::COO;
      else if (stoich_fmt_enum == "DENSE") stoich_fmt = StoichMatrix::DENSE;
      else if (stoich_fmt_enum == "CSC")   stoich_fmt = StoichMatrix::CSC;

      _mechanism = std::make_unique<MCMRuntimeMechanism>(
          mech, use_lr, stoich_fmt,
          !getParam<bool>("disable_fpoptimizer") ? false : true,
          false, true, 3);
    }

    // Configure photolysis on the mechanism
    {
      auto * mech_ptr = static_cast<MCMRuntimeMechanism*>(_mechanism.get());
      mech_ptr->setSolarParams(_lat, _lon, _day, _month, _year);
      mech_ptr->setJFac(_jfac);
      mech_ptr->setRoofOpen(_roof_open);
      if (_photolysis_scheme == "HYBRID")
        mech_ptr->enableHybridPhotolysis(_hybrid_table_dir);
      else if (_photolysis_scheme == "BOTTOMUP")
        mech_ptr->loadBottomUpData(_bottomup_data_dir, _lamp_flux_file);
    }

    // Cache metadata for fast query
    _n_species = _mechanism->nSpecies();
    _n_reactions = _mechanism->nReactions();
    _species_names = _mechanism->speciesNames();
    _console << "MCMBoxModel: Loaded " << _n_species << " species, "
             << _n_reactions << " reactions from " << mech_file << std::endl;

    // Apply default IC (if > 0) to uninitialized species.
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

    // Save initial concentrations for self-driven integrators.
    // The MOOSE Newton solve may clear the solution (u_old = 0 for the
    // first step), so we must reload the ICs before running the TS.
    _initial_conc.resize(_n_species, 0.0);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      _initial_conc[i] = sv.sln()[0];
    }

    // Sync physical params to mechanism and pre-evaluate coefficients
    {
      _console << "MCMBoxModel: evaluating coefficients..." << std::endl;
      PhysParams p;
      p.temperature = getParam<Real>("temperature");
      p.air_density = getParam<Real>("air_density");
      p.water_vapor = getParam<Real>("water_vapor");
      p.pressure = getParam<Real>("press");
      p.rh = getParam<Real>("rh");
      p.blheight = getParam<Real>("blheight");
      p.jfac = _jfac;
      p.latitude = _lat;
      p.longitude = _lon;
      _mechanism->setCurrentTime(_t);
      _mechanism->updateParams(p);
      _mechanism->markDirty();
      // Trigger initial coefficient evaluation
      auto * mech_ptr = static_cast<MCMRuntimeMechanism*>(_mechanism.get());
      mech_ptr->evaluateCoefficients();
      // Reset time-varying solar params to safe defaults
      _solar_cosx = 0.0;
      _solar_secx = 1.0e2;
      _solar_lha = 0.0;
      _console << "MCMBoxModel: coefficients evaluated." << std::endl;
    }

    // Set up dilution if dilute > 0
    Real dilute = getParam<Real>("dilute");
    if (dilute > 0.0)
    {
      _kdil = dilute;
      _conc_bkgd.assign(_n_species, 0.0);
      _console << "MCMBoxModel: Dilution enabled, kdil = " << dilute << " /s" << std::endl;
    }

    // Initialize the per-step solver.  For sundials, the SUNDIALS context is
    // built on the fly inside SundialsBoxIntegrator::solveSundialsCVODE() per
    // time step — no persistent setup needed here.
    if (_use_sundials)
    {
      _console << "MCMBoxModel: SUNDIALS direct solver will be used ("
               << _n_species << " species, rtol=" << _solver_rtol
               << ", atol=" << _solver_atol << ")" << std::endl;
    }
    else if (_use_box_solver)
    {
      _console << "MCMBoxModel: setting up PETSc TS..." << std::endl;
      setupPETScTS();
      _console << "MCMBoxModel: PETSc TS setup complete." << std::endl;
    }
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

  if (!_mechanism)
  {
    dC.assign(_n_species, 0.0);
    return;
  }

  // The mechanism evaluates concentration-dependent rate coefficients for C.
  auto * mech = static_cast<MCMRuntimeMechanism*>(_mechanism.get());
  mech->computeDCdt(C, dC);

  // Apply dilution if kdil > 0 (AtChem2 DILUTE parameter) — box-model specific
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
  if (_n_reactions == 0 || _n_species == 0 || !_mechanism)
    return;

  auto * mech = static_cast<MCMRuntimeMechanism*>(_mechanism.get());
  mech->computeJacobianTriplets(C, J);
}

void
MCMBoxModel::computeJacobianCSRValues(const std::vector<Real> & C, std::vector<Real> & values) const
{
  values.clear();
  if (_n_reactions == 0 || _n_species == 0 || !_mechanism)
    return;

  auto * mech = static_cast<MCMRuntimeMechanism*>(_mechanism.get());
  mech->computeJacobianCSRValues(C, values);
}

// -- Cached single-species interface --

Real
MCMBoxModel::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  // Lazy initialization: the framework may call getDCdt() (via computeResidual)
  // before GeneralUserObject::initialize().  If the mechanism hasn't been
  // loaded yet, force initialization now so chemistry runs on the first step.
  if (_n_species == 0)
    const_cast<MCMBoxModel*>(this)->initialize();
  if (!_mechanism) return 0.0;
  return _mechanism->getDCdt(idx, C);
}

Real
MCMBoxModel::getJacobianDiagonal(unsigned int idx, const std::vector<Real> & C) const
{
  if (!_mechanism) return 0.0;
  return _mechanism->getJacobianDiagonal(idx, C);
}

Real
MCMBoxModel::getJacobianOffDiagonal(unsigned int i, unsigned int j, const std::vector<Real> & C) const
{
  if (!_mechanism) return 0.0;
  return _mechanism->getJacobianOffDiagonal(i, j, C);
}

// _buildJacobianCache is now handled by MCMRuntimeMechanism

// loadMechanism is now handled by MCMRuntimeMechanism constructor

// setupFparser is now handled by MCMRuntimeMechanism

// compileFastHandler is now handled by MCMRuntimeMechanism

// evaluateCoefficients is now handled by MCMRuntimeMechanism

// computeDayOfYear and calculateCosSZA are now handled by MCMRuntimeMechanism

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

// -- Reaction rate queries (delegated to mechanism) --
Real
MCMBoxModel::reactionRate(unsigned int r, const std::vector<Real> & C) const
{
  return _mechanism ? _mechanism->reactionRate(r, C) : 0.0;
}

Real
MCMBoxModel::speciesReactionRate(unsigned int s, unsigned int r, const std::vector<Real> & C) const
{
  return _mechanism ? _mechanism->speciesReactionRate(s, r, C) : 0.0;
}

void
MCMBoxModel::allReactionRates(const std::vector<Real> & C, std::vector<Real> & rates) const
{
  if (_mechanism)
    _mechanism->allReactionRates(C, rates);
  else
    rates.clear();
}

Real
MCMBoxModel::speciesLossRate(unsigned int s, const std::vector<Real> & C) const
{
  return _mechanism ? _mechanism->speciesLossRate(s, C) : 0.0;
}

Real
MCMBoxModel::speciesProductionRate(unsigned int s, const std::vector<Real> & C) const
{
  return _mechanism ? _mechanism->speciesProductionRate(s, C) : 0.0;
}

// -- Photolysis (delegated to mechanism) --
void
MCMBoxModel::enableHybridPhotolysis(const std::string & table_dir)
{
  if (_mechanism)
    static_cast<MCMRuntimeMechanism*>(_mechanism.get())->enableHybridPhotolysis(table_dir);
  _photolysis_scheme = "HYBRID";
  _hybrid_table_dir = table_dir;
}

// -- BottomUp photolysis (F0AM chamber mode) --
void
MCMBoxModel::loadBottomUpData(const std::string & data_dir, const std::string & flux_file)
{
  if (_mechanism)
    static_cast<MCMRuntimeMechanism*>(_mechanism.get())->loadBottomUpData(data_dir, flux_file);
  _photolysis_scheme = "BOTTOMUP";
  _bottomup_data_dir = data_dir;
  _lamp_flux_file = flux_file;
}

// -- Solar cycle (Madronich 1993) --
void
MCMBoxModel::setSolarCycle(Real lat, Real lon, int day, int month, int year)
{
  _lat = lat;
  _lon = lon;
  _day = day;
  _month = month;
  _year = year;
  if (_mechanism)
    static_cast<MCMRuntimeMechanism*>(_mechanism.get())->setSolarParams(lat, lon, day, month, year);
}

Real
MCMBoxModel::cosSZA(Real seconds) const
{
  return _mechanism ? static_cast<MCMRuntimeMechanism*>(_mechanism.get())->cosSZA(seconds) : 0.0;
}

// -- Dilution --
void
MCMBoxModel::setDilution(Real kdil, const std::vector<Real> & bg)
{ _kdil = kdil; _conc_bkgd = bg; }

Real
MCMBoxModel::getRO2Sum(const std::vector<Real> & C) const
{
  return _mechanism ? _mechanism->getRO2Sum(C) : 0.0;
}

Real
MCMBoxModel::getJValue(unsigned int j_number) const
{
  return _mechanism ? _mechanism->getJValue(j_number) : 0.0;
}

const std::vector<std::string> &
MCMBoxModel::getRO2Species() const
{
  if (_mechanism)
    return static_cast<MCMRuntimeMechanism*>(_mechanism.get())->ro2Species();
  return _ro2_species_names;
}

const std::vector<unsigned int> &
MCMBoxModel::getRO2Indices() const
{
  if (_mechanism)
    return static_cast<MCMRuntimeMechanism*>(_mechanism.get())->ro2Indices();
  return _ro2_indices;
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

PhysParams
MCMBoxModel::currentPhysParams() const
{
  PhysParams p;
  p.temperature = getParam<Real>("temperature");
  p.air_density = getParam<Real>("air_density");
  p.water_vapor = getParam<Real>("water_vapor");
  p.pressure = getParam<Real>("press");
  p.rh = getParam<Real>("rh");
  p.blheight = getParam<Real>("blheight");
  p.jfac = _jfac;
  p.latitude = _lat;
  p.longitude = _lon;
  return p;
}

std::map<std::string, Real>
MCMBoxModel::kppGlobalValues() const
{
  const PhysParams p = currentPhysParams();

  Real air_density = p.air_density;
  if (p.pressure > 0.0)
  {
    constexpr Real NA_over_R = 6.02214129e23 / 8.3144621;
    air_density = 1.0e-6 * NA_over_R * (p.pressure * 100.0 / p.temperature);
  }

  Real h2o = p.water_vapor;
  if (p.rh >= 0.0)
  {
    const Real temp_c = p.temperature - 273.15;
    const Real wvp =
        (p.rh / 100.0) * 6.116441 * std::pow(10.0, (7.591386 * temp_c) / (temp_c + 240.7263));
    const Real press_mbar = (p.pressure > 0.0) ? p.pressure : 1013.25;
    const Real h2o_ppu = wvp / (press_mbar - wvp);
    h2o = h2o_ppu * air_density;
  }

  std::map<std::string, Real> globals;
  globals["TEMP"] = p.temperature;
  globals["AIR"] = air_density;
  globals["M"] = air_density;
  globals["O2"] = 0.21 * air_density;
  globals["N2"] = 0.78 * air_density;
  globals["H2O"] = h2o;

  if (_photolysis_scheme == "BOTTOMUP" && _kpp_bottomup_integrator)
  {
    const Real T_cur = p.temperature;
    const Real P_cur = p.pressure > 0.0 ? p.pressure : 1013.25;
    if (!_kpp_bottomup_j_valid ||
        std::abs(T_cur - _kpp_cached_bottomup_T) > 1.0e-6 ||
        std::abs(P_cur - _kpp_cached_bottomup_P) > 1.0e-6)
    {
      _kpp_cached_bottomup_j = _kpp_bottomup_integrator->computeAllJ(T_cur, P_cur);
      _kpp_cached_bottomup_T = T_cur;
      _kpp_cached_bottomup_P = P_cur;
      _kpp_bottomup_j_valid = true;
    }

    const Real roof_factor = _roof_open ? 1.0 : 0.0;
    for (const auto & [jname, value] : _kpp_cached_bottomup_j)
      globals[jname] = value * _jfac * roof_factor;
  }

  return globals;
}

// ===== PETSc TS standalone integrator =====

void
MCMBoxModel::execute()
{
  // Self-driven mode (PETSc TS or SUNDIALS): integrate chemistry at
  // TIMESTEP_END, bypassing MOOSE's solver.  ChemistryODEKernel residuals are 0
  // and this execute() sets the solution directly.
  // MOOSE-driven mode: no-op — ChemistryODEKernel provides residuals/Jacobians
  // and the Newton solver handles the integration.
  //
  // Lazy initialization: GeneralUserObject default execute_on = TIMESTEP_END,
  // so initialize() may not have been called before the first execute().
  // Force initialization now so self-driven integrators have the mechanism
  // available on the first step.
  if (_n_species == 0)
    const_cast<MCMBoxModel*>(this)->initialize();

  if (!_integrator->selfDriven())
    return;

  // Still no mechanism?  No chemistry to run (e.g. test without mechanism_file).
  if (_n_species == 0)
    return;

  // Get current and previous times from FEProblemBase
  FEProblemBase & fe_problem = static_cast<FEProblemBase &>(_subproblem);
  Real t_end = fe_problem.time();
  Real dt = fe_problem.dt();
  Real t_start = t_end - dt;

  // Skip initial call at t=0 (no step yet)
  if (t_start < 0.0 || dt <= 0.0)
    return;

  // ---- SUNDIALS direct path (solver_type = sundials) ----
  if (_use_sundials)
  {
    // Build std::vector<Real> concentration from ScalarVariable values.
    std::vector<Real> C(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      C[i] = sv.sln()[0];
    }

    // Evaluate rate coefficients at the step midpoint so time-invariant
    // photolysis parameters (BottomUp) match the behavior of the PETSc path.
    _t = 0.5 * (t_start + t_end);
    _mechanism->setCurrentTime(_t);
    if (_mechanism)
      static_cast<MCMRuntimeMechanism*>(_mechanism.get())->evaluateCoefficients();
    _mechanism->markDirty();

    // Call SUNDIALS CVODE wrapper.
    SundialsBoxIntegrator * sundials_ptr =
        static_cast<SundialsBoxIntegrator *>(_integrator.get());
    sundials_ptr->solveSundialsCVODE(t_start, t_end, C);

    // Write integrated state back to all ScalarVariable storage locations.
    NonlinearSystemBase & nl = fe_problem.getNonlinearSystemBase(0);
    NumericVector<Number> & sys_sol = *nl.system().solution;
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      dof_id_type dof = sv.dofIndices()[0];
      sv.setValue(0, C[i]);
      sys_sol.set(dof, C[i]);
      nl.solution().set(dof, C[i]);
    }
    sys_sol.close();
    nl.solution().close();
    *const_cast<NumericVector<Number> *>(nl.currentSolution()) = nl.solution();

    _console << "MCMBoxModel: SUNDIALS step t=[" << t_start << "," << t_end
             << "] dt=" << dt << " completed" << std::endl;
    return;
  }

  // ---- KPP direct path (chem_solver = kpp_*) ----
  if (_use_kpp)
  {
#ifdef KPP_ENABLED
    std::vector<Real> C(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      C[i] = sv.sln()[0];
      // Fallback: reload initial concentration if Newton solve cleared it
      if (C[i] == 0.0 && i < _initial_conc.size() && _initial_conc[i] > 0.0)
        C[i] = _initial_conc[i];
    }
    _console << "MCMBoxModel: starting KPP step t=[" << t_start << "," << t_end
             << "] nspec=" << _n_species
             << " C[0]=" << (C.size() > 0 ? C[0] : 0.0)
             << std::endl;

    // Call KPP INTEGRATE via KppBoxIntegrator.
    KppBoxIntegrator * kpp_ptr =
        static_cast<KppBoxIntegrator *>(_integrator.get());
    const auto globals = kppGlobalValues();
    kpp_ptr->solve(t_start, t_end, C, globals);

    // Write integrated state back to all ScalarVariable storage locations.
    NonlinearSystemBase & nl = fe_problem.getNonlinearSystemBase(0);
    NumericVector<Number> & sys_sol = *nl.system().solution;
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      dof_id_type dof = sv.dofIndices()[0];
      sv.setValue(0, C[i]);
      sys_sol.set(dof, C[i]);
      nl.solution().set(dof, C[i]);
    }
    sys_sol.close();
    nl.solution().close();
    *const_cast<NumericVector<Number> *>(nl.currentSolution()) = nl.solution();

    _console << "MCMBoxModel: KPP step t=[" << t_start << "," << t_end
             << "] dt=" << dt << " completed" << std::endl;
    return;
#else
    // KPP_ENABLED not defined — this path should never be reached because
    // the constructor validates _use_kpp requires KPP_ENABLED.
    mooseError("MCMBoxModel: KPP path reached but KPP_ENABLED not defined.");
#endif
  }

  // ---- PETSc TS path (solver_type = bdf/arkimex/eimex/...) ----
  // Build concentration vector from ScalarVariable values (source of truth
  // for ICs and persistent state — the NL solution vector may not be
  // populated until after the first solve).
  NonlinearSystemBase & nl = fe_problem.getNonlinearSystemBase(0);
  {
    PetscScalar *x_arr;
    CHKERRABORT(PETSC_COMM_SELF, VecGetArray(_ts_X, &x_arr));
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, _species_names[i]);
      x_arr[i] = sv.sln()[0];
      // Fallback: if the MOOSE Newton solve cleared the solution (u_old=0
      // for the first step), reload the initial concentration.
      if (x_arr[i] == 0.0 && i < _initial_conc.size() && _initial_conc[i] > 0.0)
        x_arr[i] = _initial_conc[i];
    }
    CHKERRABORT(PETSC_COMM_SELF, VecRestoreArray(_ts_X, &x_arr));
  }

  // Run PETSc TS integration from t_start to t_end
  runPETScStep(t_start, t_end);

  // Write TS solution to BOTH the libMesh solution and the variable cache.
  // The CSV output reads from sln() → _u.  sv.setValue() updates _u.
  // Also update sys.solution for consistency (used by next step's predictor).
  {
    PetscScalar *x_arr;
    CHKERRABORT(PETSC_COMM_SELF, VecGetArray(_ts_X, &x_arr));
    std::vector<Real> C_final(x_arr, x_arr + _n_species);
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
    if (_subproblem.hasScalarVariable("RO2"))
    {
      MooseVariableScalar & sv = _subproblem.getScalarVariable(0, "RO2");
      dof_id_type dof = sv.dofIndices()[0];
      const Real ro2_sum = _mechanism->getRO2Sum(C_final);
      sv.setValue(0, ro2_sum);
      sys_sol.set(dof, ro2_sum);
      nl.solution().set(dof, ro2_sum);
    }
    sys_sol.close();
    nl.solution().close();
    CHKERRABORT(PETSC_COMM_SELF, VecRestoreArray(_ts_X, &x_arr));
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

  // ODE solver parameters (_solver_type/_rtol/_atol) are read once in the
  // constructor so both PETSc TS and SUNDIALS paths share the same source.
  // (Kept here: the values are valid across all TS solver types.)

  // Helper macro for error checking without PetscCheck (which returns a value)
  // and without requiring MPI_Comm (use PETSC_COMM_SELF for sequential TS).
#define PETSC_TRY(expr) do { ierr = (expr); if (ierr) mooseError("PETSc error ", ierr, " at ", __FILE__, ":", __LINE__); } while(0)

  PETSC_TRY(TSCreate(PETSC_COMM_SELF, &_ts));
  PETSC_TRY(TSSetProblemType(_ts, TS_NONLINEAR));
  // NOTE: TSSetOptionsPrefix is deliberately NOT used here.
  //
  // Setting an options prefix would prevent the TS's internal linear solver
  // (KSP/PC) from reading the Executioner's petsc_options_iname settings
  // (e.g. -pc_type lu -pc_factor_shift_type NONZERO), which are essential
  // for convergence on stiff chemical systems.
  //
  // The TS type is set explicitly via TSSetType() and tolerances via
  // TSSetTolerances(), so the primary risk of executioner-side -ts_*
  // options leaking in is mitigated by code-level overrides.  The user
  // controls the chemical solver via the chem_solver/chem_solver_type
  // parameters in the [AtmosphericChemistry] block.
  PETSC_TRY(TSSetRHSFunction(_ts, nullptr, tsRHSFunction, this));
  PETSC_TRY(VecCreateSeq(PETSC_COMM_SELF, _n_species, &_ts_X));

  // Preallocate the PETSc matrix from the fixed chemical reaction graph.
  // The chamber mechanism is sparse: each reaction only couples products and
  // reactants from that reaction.  A dense 610x610 matrix spent most of the
  // runtime zeroing and factorizing structural zeros.
  PetscInt n = static_cast<PetscInt>(_n_species);
  const auto * runtime_mech = dynamic_cast<MCMRuntimeMechanism *>(_mechanism.get());
  if (!runtime_mech)
    mooseError("MCMBoxModel: PETSc TS sparse Jacobian requires MCMRuntimeMechanism.");

  const auto & jac_row_ptr = runtime_mech->jacobianRowPtr();
  const auto & jac_cols = runtime_mech->jacobianCols();
  if (jac_row_ptr.size() != _n_species + 1)
    mooseError("MCMBoxModel: invalid sparse Jacobian row pointer size.");

  std::vector<PetscInt> row_nnz(_n_species, 0);
  for (const auto i : make_range(_n_species))
    row_nnz[i] = static_cast<PetscInt>(jac_row_ptr[i + 1] - jac_row_ptr[i]);

  _ts_jac_cols.resize(jac_cols.size());
  for (const auto i : index_range(jac_cols))
    _ts_jac_cols[i] = static_cast<PetscInt>(jac_cols[i]);
  _ts_jac_values_real.resize(jac_cols.size(), 0.0);
  _ts_jac_values.resize(jac_cols.size(), 0.0);

  _console << "MCMBoxModel: Jacobian matrix: " << n << "x" << n
           << " sparse AIJ, nnz=" << jac_cols.size()
           << " (fill=" << (100.0 * (Real)jac_cols.size() / ((Real)n * (Real)n))
           << "%)" << std::endl;
  PETSC_TRY(MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 0, row_nnz.data(), &_ts_J));
  PETSC_TRY(MatSetFromOptions(_ts_J));

  PETSC_TRY(TSSetRHSJacobian(_ts, _ts_J, _ts_J, tsRHSJacobian, this));
  PETSC_TRY(TSSetType(_ts, _solver_type.c_str()));

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

  PETSC_TRY(TSSetTolerances(_ts, _solver_atol, nullptr, _solver_rtol, nullptr));
  // Allow unlimited SNES failures (retry with smaller step), matching extchem.c
  PETSC_TRY(TSSetMaxSNESFailures(_ts, -1));
  PETSC_TRY(TSSetSolution(_ts, _ts_X));
  PETSC_TRY(TSSetFromOptions(_ts));
  // Re-assert after TSSetFromOptions, which may reset matrix options via
  // the TS's internal option processing (especially under TSType sundials).
  PETSC_TRY(MatSetOption(_ts_J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

#undef PETSC_TRY

  _console << "MCMBoxModel: PETSc TS initialized (" << _n_species << " species, type="
           << _solver_type << ", rtol=" << _solver_rtol << ", atol=" << _solver_atol << ")" << std::endl;
}

void
MCMBoxModel::runPETScStep(PetscReal t0, PetscReal t1)
{
  PetscErrorCode ierr;

#define PETSC_TRY(expr) do { ierr = (expr); if (ierr) mooseError("PETSc error ", ierr, " at ", __FILE__, ":", __LINE__); } while(0)

  PETSC_TRY(TSSetTime(_ts, t0));
  PETSC_TRY(TSSetMaxTime(_ts, t1));

  // First chemistry interval starts conservatively because the chamber IC has
  // many zero species.  Later intervals reuse PETSc's accepted step size so
  // each output interval does not repeat the same tiny-step startup ramp.
  const PetscReal max_step = t1 - t0;
  PetscReal step0 = std::max(max_step * 1.0e-6, 1.0e-10);
  if (_ts_last_dt > 0.0)
    step0 = std::min(std::max(_ts_last_dt, 1.0e-10), max_step);
  PETSC_TRY(TSSetTimeStep(_ts, step0));
  // Also set adaptive step limits
  TSAdapt adapt;
  PETSC_TRY(TSGetAdapt(_ts, &adapt));
  PETSC_TRY(TSAdaptSetStepLimits(adapt, 1.0e-12, (t1 - t0)));
  PETSC_TRY(TSSetSolution(_ts, _ts_X));

  // Use midpoint time for time-dependent photolysis during this interval.
  // Concentration-dependent coefficients are evaluated from the current TS
  // state in the RHS/Jacobian callbacks.
  _t = 0.5 * (t0 + t1);
  _mechanism->setCurrentTime(_t);

  // Run the integrator
  PETSC_TRY(TSSolve(_ts, _ts_X));

  // Get integration statistics
  PetscInt steps;
  TSConvergedReason reason;
  PETSC_TRY(TSGetStepNumber(_ts, &steps));
  PETSC_TRY(TSGetConvergedReason(_ts, &reason));
  PETSC_TRY(TSGetTimeStep(_ts, &_ts_last_dt));

#undef PETSC_TRY

  // Diagnostic: print key species after TS integration
  {
    const PetscScalar *dbg_arr;
    CHKERRABORT(PETSC_COMM_SELF, VecGetArrayRead(_ts_X, &dbg_arr));
    const auto & sp_names = _mechanism->speciesNames();
    auto lookup = [&](const std::string & name) -> PetscInt {
      for (size_t i = 0; i < sp_names.size(); ++i)
        if (sp_names[i] == name) return static_cast<PetscInt>(i);
      return -1;
    };
    PetscInt c5h8 = lookup("C5H8"), no2 = lookup("NO2");
    PetscInt o3 = lookup("O3"), oh = lookup("OH");
    _console << "MCMBoxModel: TS [" << t0 << "," << t1 << "] C5H8="
             << (c5h8 >= 0 ? dbg_arr[c5h8] : -1.0)
             << " NO2=" << (no2 >= 0 ? dbg_arr[no2] : -1.0)
             << " O3=" << (o3 >= 0 ? dbg_arr[o3] : -1.0)
             << " OH=" << (oh >= 0 ? dbg_arr[oh] : -1.0)
             << std::endl;
    CHKERRABORT(PETSC_COMM_SELF, VecRestoreArrayRead(_ts_X, &dbg_arr));
  }

  _console << "MCMBoxModel: TS step [" << t0 << "," << t1 << "] "
           << steps << " internal steps, reason=" << reason << std::endl;
}

// static
PetscErrorCode
MCMBoxModel::tsRHSFunction(TS /*ts*/, PetscReal /*t*/, Vec C, Vec F, void *ctx)
{
  MCMBoxModel *model = static_cast<MCMBoxModel *>(ctx);
  const PetscScalar *c_arr;
  PetscScalar *f_arr;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(C, &c_arr));
  PetscCall(VecGetArray(F, &f_arr));

  // Build std::vector around PETSc array (no copy for 610 elements - acceptable)
  std::vector<Real> C_vec(c_arr, c_arr + model->_n_species);
  std::vector<Real> dC_vec;
  model->computeDCdt(C_vec, dC_vec);

  for (unsigned int i = 0; i < model->_n_species; ++i)
    f_arr[i] = dC_vec[i];

  PetscCall(VecRestoreArrayRead(C, &c_arr));
  PetscCall(VecRestoreArray(F, &f_arr));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// static
PetscErrorCode
MCMBoxModel::tsRHSJacobian(TS /*ts*/, PetscReal /*t*/, Vec C, Mat Amat, Mat Pmat, void *ctx)
{
  MCMBoxModel *model = static_cast<MCMBoxModel *>(ctx);
  const PetscScalar *c_arr;

  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(C, &c_arr));

  // Build concentration vector
  std::vector<Real> C_vec(c_arr, c_arr + model->_n_species);

  model->computeJacobianCSRValues(C_vec, model->_ts_jac_values_real);
  if (model->_ts_jac_values.size() != model->_ts_jac_values_real.size())
    model->_ts_jac_values.resize(model->_ts_jac_values_real.size());
  for (const auto i : index_range(model->_ts_jac_values_real))
    model->_ts_jac_values[i] = static_cast<PetscScalar>(model->_ts_jac_values_real[i]);

  const auto * runtime_mech = static_cast<MCMRuntimeMechanism *>(model->_mechanism.get());
  const auto & row_ptr = runtime_mech->jacobianRowPtr();

  PetscCall(MatZeroEntries(Pmat));
  for (const auto row : make_range(model->_n_species))
  {
    const PetscInt irow = static_cast<PetscInt>(row);
    const auto begin = row_ptr[row];
    const auto ncols = row_ptr[row + 1] - begin;
    if (ncols == 0)
      continue;
    PetscCall(MatSetValues(Pmat,
                           1,
                           &irow,
                           static_cast<PetscInt>(ncols),
                           model->_ts_jac_cols.data() + begin,
                           model->_ts_jac_values.data() + begin,
                           INSERT_VALUES));
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

// ===== SUNDIALS direct solver dispatcher =====

void
MCMBoxModel::solveSundialsCVODEWrapper(Real t0, Real t1, std::vector<Real> & C)
{
  // This public entry forwards to the SundialsBoxIntegrator member via the
  // BoxIntegrator unique_ptr.  Called only when _use_sundials is true.
  SundialsBoxIntegrator * sundials_ptr =
      static_cast<SundialsBoxIntegrator *>(_integrator.get());
  if (!sundials_ptr)
    mooseError("MCMBoxModel: _use_sundials is true but _integrator is not a "
               "SundialsBoxIntegrator — this should not happen.");
  sundials_ptr->solveSundialsCVODE(t0, t1, C);
}
