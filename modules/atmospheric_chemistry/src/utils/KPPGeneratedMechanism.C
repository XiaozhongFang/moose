//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KPPGeneratedMechanism.h"

#include "Moose.h"

#include <dlfcn.h>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cmath>

KPPGeneratedMechanism::KPPGeneratedMechanism(const std::string & lib_path)
  : _lib_handle(nullptr),
    _fun(nullptr),
    _jac(nullptr),
    _init(nullptr),
    _update_rconst(nullptr),
    _kpp_C(nullptr),
    _kpp_VAR(nullptr),
    _kpp_FIX(nullptr),
    _kpp_RCONST(nullptr),
    _lu_irow(nullptr),
    _lu_icol(nullptr),
    _lu_crow(nullptr),
    _n_species(0),
    _n_variable(0),
    _n_reactions(0),
    _jac_nnz(0),
    _roof_open(true),
    _jfac(1.0),
    _t(0.0),
    _n_j_vals(0)
{
  // Open the KPP shared library
  _lib_handle = dlopen(lib_path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!_lib_handle)
    mooseError("KPPGeneratedMechanism: cannot load KPP shared library \"",
               lib_path, "\": ", dlerror());

  // Resolve KPP C API function symbols
  _fun  = reinterpret_cast<KppFunFn>(dlsym(_lib_handle, "Fun"));
  _jac  = reinterpret_cast<KppJacFn>(dlsym(_lib_handle, "Jac_SP"));
  _init = reinterpret_cast<KppInitFn>(dlsym(_lib_handle, "Initialize"));

  if (!_fun || !_jac || !_init)
  {
    std::string msg = "KPPGeneratedMechanism: failed to resolve KPP API symbols "
                      "(Fun, Jac_SP, Initialize) from \"";
    msg += lib_path + "\"";
    dlclose(_lib_handle);
    _lib_handle = nullptr;
    mooseError(msg);
  }

  // Resolve KPP dimension accessors (injected by kpp/build/Makefile)
  using GetIntFn = int (*)(void);
  auto kpp_get_nspec_fn  = reinterpret_cast<GetIntFn>(dlsym(_lib_handle, "kpp_get_nspec"));
  auto kpp_get_nvar_fn   = reinterpret_cast<GetIntFn>(dlsym(_lib_handle, "kpp_get_nvar"));
  auto kpp_get_nreact_fn = reinterpret_cast<GetIntFn>(dlsym(_lib_handle, "kpp_get_nreact"));

  if (!kpp_get_nspec_fn)
  {
    std::string msg = "KPPGeneratedMechanism: failed to resolve KPP dimension";
    msg += " function 'kpp_get_nspec' from \"" + lib_path + "\"";
    dlclose(_lib_handle);
    _lib_handle = nullptr;
    mooseError(msg);
  }

  _n_species   = static_cast<unsigned int>(kpp_get_nspec_fn());
  _n_variable  = kpp_get_nvar_fn ? static_cast<unsigned int>(kpp_get_nvar_fn()) : _n_species;
  _n_reactions = kpp_get_nreact_fn ? static_cast<unsigned int>(kpp_get_nreact_fn()) : 0;

  // Resolve KPP global variable pointers
  _kpp_C      = reinterpret_cast<double *>(dlsym(_lib_handle, "C"));
  _kpp_VAR    = reinterpret_cast<double *>(dlsym(_lib_handle, "VAR"));
  _kpp_FIX    = reinterpret_cast<double *>(dlsym(_lib_handle, "FIX"));
  _kpp_RCONST = reinterpret_cast<double *>(dlsym(_lib_handle, "RCONST"));

  if (!_kpp_C)
  {
    std::string msg = "KPPGeneratedMechanism: failed to resolve KPP global ";
    msg += "array 'C' from \"" + lib_path + "\"";
    dlclose(_lib_handle);
    _lib_handle = nullptr;
    mooseError(msg);
  }

  // Resolve sparse Jacobian structure (LU_*) — these are integral to
  // KPP and should always be present in a compiled KPP library.
  _lu_irow = reinterpret_cast<int *>(dlsym(_lib_handle, "LU_IROW"));
  _lu_icol = reinterpret_cast<int *>(dlsym(_lib_handle, "LU_ICOL"));
  _lu_crow = reinterpret_cast<int *>(dlsym(_lib_handle, "LU_CROW"));

  if (!_lu_irow || !_lu_icol || !_lu_crow)
    mooseError("KPPGeneratedMechanism: failed to resolve sparse Jacobian "
               "structure (LU_IROW, LU_ICOL, LU_CROW) from \"",
               lib_path, "\". Rebuild the KPP library.");

  // NNZ = LU_CROW[NVAR] (total non-zero entries in the Jacobian)
  _jac_nnz = _lu_crow[_n_variable];

  // Attempt to read species names from KPP's SPC_NAMES global.
  // KPP may export SPC_NAMES as char* SPC_NAMES[NSPEC][32] or similar.
  // If not available, fall back to generic names.
  char ** kpp_spc_names = reinterpret_cast<char **>(dlsym(_lib_handle, "SPC_NAMES"));
  if (kpp_spc_names)
  {
    _species_names.reserve(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      if (kpp_spc_names[i])
        _species_names.emplace_back(kpp_spc_names[i]);
      else
        _species_names.emplace_back("spc_" + std::to_string(i));
    }
  }
  else
  {
    _species_names.reserve(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
      _species_names.emplace_back("species_" + std::to_string(i));
  }

  // Initialize KPP global state (sets up FIX, RCONST, etc.)
  _init();

  // Count J-values by scanning RCONST names (heuristic)
  // KPP typically names rate constants KJ_<N> or J_<N> for photolysis.
  // For now, default to 0 — updateParams will re-count.
  _n_j_vals = 0;
}

KPPGeneratedMechanism::~KPPGeneratedMechanism()
{
  if (_lib_handle)
  {
    dlclose(_lib_handle);
    _lib_handle = nullptr;
  }
}

void
KPPGeneratedMechanism::updateParams(const PhysParams & params)
{
  // Update the KPP global RCONST array with new rate coefficients.
  //
  // KPP-generated mechanisms typically provide an Update_RCONST() function
  // that re-evaluates rate coefficients from current TEMP, AIR, H2O, etc.
  // If the mechanism exports Update_RCONST, we call it after setting globals.
  //
  // If Update_RCONST is not available (e.g. older generation), the user
  // must ensure the .so's global arrays are updated externally.

  using UpdateRCONSTFn = void (*)(void);
  auto update_rconst = reinterpret_cast<UpdateRCONSTFn>(
      dlsym(_lib_handle, "Update_RCONST"));

  if (update_rconst)
  {
    // Set KPP global physical parameters if the symbols are available
    double * kpp_temp = reinterpret_cast<double *>(dlsym(_lib_handle, "TEMP"));
    double * kpp_air  = reinterpret_cast<double *>(dlsym(_lib_handle, "AIR"));

    if (kpp_temp)
      *kpp_temp = static_cast<double>(params.temperature);

    // AIR density: either from params.air_density or computed from pressure
    double air_dens = static_cast<double>(params.air_density);
    if (params.pressure > 0.0)
    {
      // n/V = P/(kB*T) in molec/cm³
      // pressure in mbar → Pa = mbar * 100
      // kB = 1.380649e-23 J/K
      // 1 m³ = 1e6 cm³
      air_dens = (params.pressure * 100.0) /
                 (1.380649e-23 * params.temperature) * 1e-6;
    }
    if (kpp_air)
      *kpp_air = air_dens;

    // H2O: from water_vapor or computed from RH
    double * kpp_h2o = reinterpret_cast<double *>(dlsym(_lib_handle, "H2O"));
    if (kpp_h2o)
    {
      double h2o = static_cast<double>(params.water_vapor);
      if (params.rh >= 0.0)
      {
        // Saturation vapor pressure (Clausius-Clapeyron approximation)
        double es = 6.112 * std::exp(17.67 * (params.temperature - 273.15) /
                                      (params.temperature - 29.65));
        double p_h2o = (params.rh / 100.0) * es * 100.0; // Pa
        h2o = p_h2o / (1.380649e-23 * params.temperature) * 1e-6; // molec/cm³
      }
      *kpp_h2o = h2o;
    }

    // Call Update_RCONST to re-evaluate all rate coefficients
    update_rconst();
  }
  else
  {
    // Fallback: set RCONST directly if no Update_RCONST is available.
    // In this case the user must populate RCONST externally.
    // For now just warn once (not every call).
    static bool warned = false;
    if (!warned)
    {
      mooseDoOnce(mooseWarning(
          "KPPGeneratedMechanism: KPP library does not export Update_RCONST. "
          "Rate coefficients will NOT be automatically updated."));
      warned = true;
    }
  }
}

void
KPPGeneratedMechanism::computeRHS(Real /*t*/,
                                   const std::vector<Real> & C,
                                   const PhysParams & /*params*/,
                                   std::vector<Real> & dC_dt) const
{
  unsigned int n = _n_species;
  if (C.size() < n || dC_dt.size() < n)
    mooseError("KPPGeneratedMechanism::computeRHS: concentration vector size mismatch");

  // Copy concentrations into KPP's global C[] array
  for (unsigned int i = 0; i < n; ++i)
    _kpp_C[i] = static_cast<double>(C[i]);

  // Call KPP's Fun() which writes into the global arrays and Ydot.
  // Fun signature: Fun(Y, FIX, RCONST, Ydot)
  // KPP uses the global C[] as Y internally; we pass it explicitly.
  std::vector<double> Ydot(n, 0.0);
  _fun(_kpp_C, _kpp_FIX, _kpp_RCONST, Ydot.data());

  // Copy result back to MOOSE Real vector
  for (unsigned int i = 0; i < n; ++i)
    dC_dt[i] = static_cast<Real>(Ydot[i]);
}

void
KPPGeneratedMechanism::computeJacobian(
    Real /*t*/,
    const std::vector<Real> & C,
    const PhysParams & /*params*/,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  unsigned int n = _n_species;
  if (C.size() < n)
    mooseError("KPPGeneratedMechanism::computeJacobian: concentration vector size mismatch");

  // Copy concentrations into KPP's global C[] array
  for (unsigned int i = 0; i < n; ++i)
    _kpp_C[i] = static_cast<double>(C[i]);

  // KPP's Jac_SP writes into global JVS[] sparse Jacobian storage.
  // The KPP sparse Jacobian format uses:
  //   JVS[LU_NONZERO] — non-zero Jacobian values
  //   LU_IROW[LU_NONZERO] — row indices (1-based)
  //   LU_ICOL[LU_NONZERO] — column indices (1-based)
  //
  // We resolve these globals at runtime (they may not exist in all KPP versions).

  // Resolve KPP sparse Jacobian function (Jac_SP) and LU structure.
  // Jac_SP(Y, FIX, RCONST, JVS) evaluates the sparse Jacobian and writes
  // to the caller-provided JVS array.  LU_IROW/LU_ICOL give the structure.
  //
  // NNZ = LU_CROW[NVAR] — pre-computed in constructor as _jac_nnz.

  if (!_jac || !_lu_irow || !_lu_icol || _jac_nnz == 0)
  {
    mooseDoOnce(mooseWarning(
        "KPPGeneratedMechanism::computeJacobian: KPP Jacobian not available. "
        "Returning empty Jacobian."));
    J.clear();
    return;
  }

  // Call KPP's Jac_SP() with our own JVS buffer
  std::vector<double> jvs(_jac_nnz);
  _jac(_kpp_C, _kpp_FIX, _kpp_RCONST, jvs.data());

  // Convert sparse Jacobian to triplet format
  int nnz = static_cast<int>(_jac_nnz);
  J.clear();
  J.reserve(static_cast<std::size_t>(nnz));
  for (int k = 0; k < nnz; ++k)
  {
    // KPP uses 1-based indexing; convert to 0-based
    unsigned int row = static_cast<unsigned int>(_lu_irow[k] - 1);
    unsigned int col = static_cast<unsigned int>(_lu_icol[k] - 1);
    Real val = static_cast<Real>(jvs[k]);
    J.emplace_back(row, col, val);
  }
}

SpeciesRates
KPPGeneratedMechanism::computeSpeciesRates(
    Real /*t*/,
    const std::vector<Real> & /*C*/,
    const PhysParams & /*params*/) const
{
  // KPP does not natively provide per-species production/loss rates
  // separate from the net RHS.  This would require iterating over
  // reactions and their stoichiometry — information not easily
  // extracted from KPP globals without the mechanism's equation-by-equation
  // data.
  //
  // For the scaffold, return zero vectors.

  SpeciesRates rates;
  rates.production.assign(_n_species, 0.0);
  rates.loss.assign(_n_species, 0.0);
  return rates;
}

// ===== Single-species accessors (delegate to computeRHS / computeJacobian) =====

Real
KPPGeneratedMechanism::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  // Compute full RHS and return the requested component.
  // Inefficient for repeated calls — caches could be added later.
  PhysParams dummy;
  std::vector<Real> dC_dt(_n_species);
  computeRHS(_t, C, dummy, dC_dt);
  return dC_dt[idx];
}

Real
KPPGeneratedMechanism::getJacobianDiagonal(unsigned int idx,
                                            const std::vector<Real> & C) const
{
  // Compute full Jacobian and extract diagonal element
  PhysParams dummy;
  std::vector<std::tuple<unsigned int, unsigned int, Real>> J;
  computeJacobian(_t, C, dummy, J);
  for (const auto & entry : J)
  {
    if (std::get<0>(entry) == idx && std::get<1>(entry) == idx)
      return std::get<2>(entry);
  }
  return 0.0;
}

Real
KPPGeneratedMechanism::getJacobianOffDiagonal(unsigned int i,
                                               unsigned int j,
                                               const std::vector<Real> & C) const
{
  // Compute full Jacobian and extract the requested off-diagonal element
  PhysParams dummy;
  std::vector<std::tuple<unsigned int, unsigned int, Real>> J;
  computeJacobian(_t, C, dummy, J);
  for (const auto & entry : J)
  {
    if (std::get<0>(entry) == i && std::get<1>(entry) == j)
      return std::get<2>(entry);
  }
  return 0.0;
}

// ===== Photolysis =====

Real
KPPGeneratedMechanism::getJValue(unsigned int /*j_number*/) const
{
  // KPP does not expose photolysis J-values through a standard API.
  // If the mechanism uses photolysis, they are embedded in RCONST entries.
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::getJValue: not implemented for KPP mechanisms."));
  return 0.0;
}

// ===== Reaction diagnostics =====

Real
KPPGeneratedMechanism::getRO2Sum(const std::vector<Real> & /*C*/) const
{
  // RO2 sum requires knowing which species are peroxy radicals.
  // This metadata is not available from KPP globals alone.
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::getRO2Sum: RO2 species indices not available "
      "from KPP library.  Returning 0."));
  return 0.0;
}

Real
KPPGeneratedMechanism::reactionRate(unsigned int /*r*/,
                                     const std::vector<Real> & /*C*/) const
{
  // KPP does not provide a per-reaction rate query API.
  // Would require access to the mechanism's stoichiometry data.
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::reactionRate: not implemented."));
  return 0.0;
}

Real
KPPGeneratedMechanism::speciesReactionRate(unsigned int /*s*/,
                                            unsigned int /*r*/,
                                            const std::vector<Real> & /*C*/) const
{
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::speciesReactionRate: not implemented."));
  return 0.0;
}

Real
KPPGeneratedMechanism::speciesLossRate(unsigned int /*s*/,
                                        const std::vector<Real> & /*C*/) const
{
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::speciesLossRate: not implemented."));
  return 0.0;
}

Real
KPPGeneratedMechanism::speciesProductionRate(unsigned int /*s*/,
                                              const std::vector<Real> & /*C*/) const
{
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::speciesProductionRate: not implemented."));
  return 0.0;
}

void
KPPGeneratedMechanism::allReactionRates(const std::vector<Real> & /*C*/,
                                         std::vector<Real> & rates) const
{
  mooseDoOnce(mooseWarning(
      "KPPGeneratedMechanism::allReactionRates: not implemented."));
  rates.assign(_n_reactions, 0.0);
}
