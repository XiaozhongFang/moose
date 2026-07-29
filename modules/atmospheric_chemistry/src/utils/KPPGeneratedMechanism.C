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
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>

namespace
{
std::recursive_mutex kpp_generated_mechanism_mutex;

bool
samePhysParams(const PhysParams & a, const PhysParams & b)
{
  return a.temperature == b.temperature && a.air_density == b.air_density &&
         a.water_vapor == b.water_vapor && a.pressure == b.pressure && a.rh == b.rh &&
         a.jfac == b.jfac && a.latitude == b.latitude && a.longitude == b.longitude &&
         a.albedo == b.albedo && a.o3column == b.o3column && a.altitude == b.altitude &&
         a.cos_sza == b.cos_sza && a.blheight == b.blheight && a.j_vals == b.j_vals;
}
}

KPPGeneratedMechanism::KPPGeneratedMechanism(const std::string & lib_path)
  : _lib_handle(nullptr),
    _fun(nullptr),
    _jac(nullptr),
    _init(nullptr),
    _update_rconst(nullptr),
    _kpp_C(nullptr),
    _kpp_VAR_ptr(nullptr),
    _kpp_FIX_ptr(nullptr),
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
    _n_j_vals(0),
    _cached_time(std::numeric_limits<Real>::quiet_NaN()),
    _rhs_valid(false),
    _jacobian_valid(false),
    _cached_params_valid(false)
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
  _update_rconst = reinterpret_cast<KppVoidFn>(dlsym(_lib_handle, "Update_RCONST"));

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
  _kpp_C       = reinterpret_cast<double *>(dlsym(_lib_handle, "C"));
  _kpp_VAR_ptr = reinterpret_cast<double **>(dlsym(_lib_handle, "VAR"));
  _kpp_FIX_ptr = reinterpret_cast<double **>(dlsym(_lib_handle, "FIX"));
  _kpp_RCONST  = reinterpret_cast<double *>(dlsym(_lib_handle, "RCONST"));

  if (!_kpp_C || !_kpp_VAR_ptr || !_kpp_FIX_ptr || !_kpp_RCONST)
  {
    std::string msg = "KPPGeneratedMechanism: failed to resolve KPP global ";
    msg += "arrays (C, VAR, FIX, RCONST) from \"" + lib_path + "\"";
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

  // NNZ = LU_CROW[NVAR] (total non-zero entries in the Jacobian).
  // KPP-generated C stores LU_IROW/LU_ICOL as 0-based indices.
  _jac_nnz = _lu_crow[_n_variable];

  // Attempt to read species names from KPP's SPC_NAMES global.
  // KPP may export SPC_NAMES as char* SPC_NAMES[NSPEC][32] or similar.
  // If not available, fall back to generic names.
  char ** kpp_spc_names = reinterpret_cast<char **>(dlsym(_lib_handle, "SPC_NAMES"));
  if (kpp_spc_names)
  {
    _all_species_names.reserve(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
    {
      if (kpp_spc_names[i])
        _all_species_names.emplace_back(kpp_spc_names[i]);
      else
        _all_species_names.emplace_back("spc_" + std::to_string(i));
    }
  }
  else
  {
    _all_species_names.reserve(_n_species);
    for (unsigned int i = 0; i < _n_species; ++i)
      _all_species_names.emplace_back("species_" + std::to_string(i));
  }
  _species_names.assign(_all_species_names.begin(),
                        _all_species_names.begin() + static_cast<std::ptrdiff_t>(_n_variable));

  // Initialize KPP global state (sets up FIX, RCONST, etc.)
  _init();

  // Count J-values by scanning RCONST names (heuristic)
  // KPP typically names rate constants KJ_<N> or J_<N> for photolysis.
  // For now, default to 0 — updateParams will re-count.
  _n_j_vals = 0;
  _cached_C.assign(_n_variable, 0.0);
  _cached_rhs.assign(_n_variable, 0.0);
}

KPPGeneratedMechanism::~KPPGeneratedMechanism()
{
  if (_lib_handle)
  {
    dlclose(_lib_handle);
    _lib_handle = nullptr;
  }
}

bool
KPPGeneratedMechanism::setGlobal(const std::string & name, Real value)
{
  std::lock_guard<std::recursive_mutex> lock(kpp_generated_mechanism_mutex);

  double * value_ptr = reinterpret_cast<double *>(dlsym(_lib_handle, name.c_str()));
  if (!value_ptr)
    return false;

  const double new_value = static_cast<double>(value);
  if (*value_ptr != new_value)
  {
    *value_ptr = new_value;
    markDirty();
  }
  else
    *value_ptr = new_value;

  return true;
}

void
KPPGeneratedMechanism::updateParams(const PhysParams & params)
{
  std::lock_guard<std::recursive_mutex> lock(kpp_generated_mechanism_mutex);

  _rhs_valid = false;
  _jacobian_valid = false;
  _cached_params = params;
  _cached_params_valid = true;
  // Update the KPP global RCONST array with new rate coefficients.
  //
  // KPP-generated mechanisms typically provide an Update_RCONST() function
  // that re-evaluates rate coefficients from current TEMP, AIR, H2O, etc.
  // If the mechanism exports Update_RCONST, we call it after setting globals.
  //
  // If Update_RCONST is not available (e.g. older generation), the user
  // must ensure the .so's global arrays are updated externally.

  if (_update_rconst)
  {
    // Set KPP global physical parameters if the symbols are available.
    auto set_global = [this](const char * name, double value)
    {
      double * value_ptr = reinterpret_cast<double *>(dlsym(_lib_handle, name));
      if (value_ptr)
        *value_ptr = value;
    };

    set_global("TEMP", static_cast<double>(params.temperature));
    set_global("TIME", static_cast<double>(_t));
    set_global("CFACTOR", 1.0);

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
    const double pressure_pa =
        params.pressure > 0.0
            ? static_cast<double>(params.pressure) * 100.0
            : air_dens * 1.0e6 * 1.380649e-23 * static_cast<double>(params.temperature);
    set_global("PRESS", pressure_pa);
    set_global("AIR", air_dens);
    set_global("M", air_dens);
    set_global("O2", 0.21 * air_dens);
    set_global("N2", 0.78 * air_dens);

    // H2O: from water_vapor or computed from RH
    double h2o = static_cast<double>(params.water_vapor);
    if (params.rh >= 0.0)
    {
      // Saturation vapor pressure (Clausius-Clapeyron approximation)
      double es = 6.112 * std::exp(17.67 * (params.temperature - 273.15) /
                                    (params.temperature - 29.65));
      double p_h2o = (params.rh / 100.0) * es * 100.0; // Pa
      h2o = p_h2o / (1.380649e-23 * params.temperature) * 1e-6; // molec/cm³
    }
    set_global("H2O", h2o);

    double * var = _kpp_VAR_ptr ? *_kpp_VAR_ptr : nullptr;
    double * fix = _kpp_FIX_ptr ? *_kpp_FIX_ptr : nullptr;
    if (fix && var && _all_species_names.size() >= _n_species)
    {
      const double o2 = 0.21 * air_dens;
      const double n2 = 0.78 * air_dens;
      for (unsigned int i = _n_variable; i < _n_species; ++i)
      {
        const auto fix_index = i - _n_variable;

        const std::string & name = _all_species_names[i];
        if (name == "M" || name == "AIR")
          fix[fix_index] = air_dens;
        else if (name == "O2")
          fix[fix_index] = o2;
        else if (name == "N2")
          fix[fix_index] = n2;
        else if (name == "H2O")
          fix[fix_index] = h2o;
      }
    }

    auto update_sun = reinterpret_cast<KppVoidFn>(dlsym(_lib_handle, "Update_SUN"));
    if (update_sun)
      update_sun();

    double * kpp_sun = reinterpret_cast<double *>(dlsym(_lib_handle, "SUN"));
    if (kpp_sun)
    {
      if (!_roof_open)
        *kpp_sun = 0.0;
      else
        *kpp_sun *= static_cast<double>(_jfac * params.jfac);
    }

    // Call Update_RCONST to re-evaluate all rate coefficients
    _update_rconst();
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
KPPGeneratedMechanism::computeRHS(Real t,
                                   const std::vector<Real> & C,
                                   const PhysParams & params,
                                   std::vector<Real> & dC_dt) const
{
  std::lock_guard<std::recursive_mutex> lock(kpp_generated_mechanism_mutex);

  unsigned int n = _n_variable;
  if (C.size() < n || dC_dt.size() < n)
    mooseError("KPPGeneratedMechanism::computeRHS: concentration vector size mismatch");

  const bool cache_hit = _rhs_valid && _cached_C.size() == n &&
                         _cached_time == t &&
                         _cached_params_valid && samePhysParams(_cached_params, params) &&
                         std::equal(_cached_C.begin(), _cached_C.end(), C.begin());
  if (cache_hit)
  {
    for (unsigned int i = 0; i < n; ++i)
      dC_dt[i] = _cached_rhs[i];
    return;
  }

  double * var = *_kpp_VAR_ptr;
  double * fix = *_kpp_FIX_ptr;
  for (unsigned int i = 0; i < n; ++i)
    var[i] = static_cast<double>(C[i]);

  _t = t;
  const_cast<KPPGeneratedMechanism *>(this)->updateParams(params);

  std::vector<double> Ydot(n, 0.0);
  _fun(var, fix, _kpp_RCONST, Ydot.data());

  for (unsigned int i = 0; i < n; ++i)
  {
    dC_dt[i] = static_cast<Real>(Ydot[i]);
    _cached_rhs[i] = dC_dt[i];
    _cached_C[i] = C[i];
  }
  _cached_time = t;
  _cached_params = params;
  _cached_params_valid = true;
  _rhs_valid = true;
  _jacobian_valid = false;
}

void
KPPGeneratedMechanism::computeJacobian(
    Real t,
    const std::vector<Real> & C,
    const PhysParams & params,
    std::vector<std::tuple<unsigned int, unsigned int, Real>> & J) const
{
  std::lock_guard<std::recursive_mutex> lock(kpp_generated_mechanism_mutex);

  unsigned int n = _n_variable;
  if (C.size() < n)
    mooseError("KPPGeneratedMechanism::computeJacobian: concentration vector size mismatch");

  const bool cache_hit = _jacobian_valid && _cached_C.size() == n &&
                         _cached_time == t &&
                         _cached_params_valid && samePhysParams(_cached_params, params) &&
                         std::equal(_cached_C.begin(), _cached_C.end(), C.begin());
  if (cache_hit)
  {
    J = _cached_jacobian;
    return;
  }

  double * var = *_kpp_VAR_ptr;
  double * fix = *_kpp_FIX_ptr;
  for (unsigned int i = 0; i < n; ++i)
    var[i] = static_cast<double>(C[i]);

  _t = t;
  const_cast<KPPGeneratedMechanism *>(this)->updateParams(params);

  // KPP's Jac_SP writes into global JVS[] sparse Jacobian storage.
  // The KPP sparse Jacobian format uses:
  //   JVS[LU_NONZERO] — non-zero Jacobian values
  //   LU_IROW[LU_NONZERO] — row indices (0-based in generated C)
  //   LU_ICOL[LU_NONZERO] — column indices (0-based in generated C)
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
  _jac(var, fix, _kpp_RCONST, jvs.data());

  // Convert sparse Jacobian to triplet format
  int nnz = static_cast<int>(_jac_nnz);
  J.clear();
  J.reserve(static_cast<std::size_t>(nnz));
  for (int k = 0; k < nnz; ++k)
  {
    if (_lu_irow[k] < 0 || _lu_icol[k] < 0)
      continue;

    unsigned int row = static_cast<unsigned int>(_lu_irow[k]);
    unsigned int col = static_cast<unsigned int>(_lu_icol[k]);
    if (row >= n || col >= n)
      continue;

    Real val = static_cast<Real>(jvs[k]);
    J.emplace_back(row, col, val);
  }
  _cached_C.assign(C.begin(), C.begin() + static_cast<std::ptrdiff_t>(n));
  _cached_jacobian = J;
  _cached_time = t;
  _cached_params = params;
  _cached_params_valid = true;
  _jacobian_valid = true;
}

SpeciesRates
KPPGeneratedMechanism::computeSpeciesRates(
    Real /*t*/,
    const std::vector<Real> & /*C*/,
    const PhysParams & /*params*/) const
{
  mooseError("KPPGeneratedMechanism::computeSpeciesRates: separated production/loss "
             "rates require stoichiometry metadata that KPP shared libraries do not export.");
}

// ===== Single-species accessors (delegate to computeRHS / computeJacobian) =====

Real
KPPGeneratedMechanism::getDCdt(unsigned int idx, const std::vector<Real> & C) const
{
  // Compute full RHS and return the requested component.
  // Inefficient for repeated calls — caches could be added later.
  PhysParams dummy;
  std::vector<Real> dC_dt(_n_variable);
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
  mooseError("KPPGeneratedMechanism::getJValue: individual photolysis J values "
             "are not exported by KPP shared libraries.");
}

// ===== Reaction diagnostics =====

Real
KPPGeneratedMechanism::getRO2Sum(const std::vector<Real> & /*C*/) const
{
  // RO2 sum requires knowing which species are peroxy radicals.
  // This metadata is not available from KPP globals alone.
  mooseError("KPPGeneratedMechanism::getRO2Sum: RO2 species metadata is not "
             "available from KPP shared libraries.");
}

Real
KPPGeneratedMechanism::reactionRate(unsigned int /*r*/,
                                     const std::vector<Real> & /*C*/) const
{
  mooseError("KPPGeneratedMechanism::reactionRate: per-reaction rates require "
             "stoichiometry metadata that KPP shared libraries do not export.");
}

Real
KPPGeneratedMechanism::speciesReactionRate(unsigned int /*s*/,
                                            unsigned int /*r*/,
                                            const std::vector<Real> & /*C*/) const
{
  mooseError("KPPGeneratedMechanism::speciesReactionRate: per-reaction species "
             "rates require stoichiometry metadata that KPP shared libraries do not export.");
}

Real
KPPGeneratedMechanism::speciesLossRate(unsigned int /*s*/,
                                        const std::vector<Real> & /*C*/) const
{
  mooseError("KPPGeneratedMechanism::speciesLossRate: separated production/loss "
             "rates require stoichiometry metadata that KPP shared libraries do not export.");
}

Real
KPPGeneratedMechanism::speciesProductionRate(unsigned int /*s*/,
                                              const std::vector<Real> & /*C*/) const
{
  mooseError("KPPGeneratedMechanism::speciesProductionRate: separated production/loss "
             "rates require stoichiometry metadata that KPP shared libraries do not export.");
}

void
KPPGeneratedMechanism::allReactionRates(const std::vector<Real> & /*C*/,
                                         std::vector<Real> & rates) const
{
  rates.clear();
  mooseError("KPPGeneratedMechanism::allReactionRates: per-reaction rates require "
             "stoichiometry metadata that KPP shared libraries do not export.");
}
