// kpp/adapter/KppBoxIntegrator.h
#pragma once

#include "BoxIntegrator.h"
#include "ConsoleStreamInterface.h"
#include "KppFortranBridge.h"

#include <map>
#include <string>
#include <vector>

class MooseApp;

/**
 * KPP-generated code integrator for box-mode chemistry.
 *
 * Wraps KPP-generated code via dlopen/dlsym (C kpp_adapter.c bridge).
 * Self-driven mode (selfDriven() == true) — all BoxIntegrator residual/Jacobian
 * methods return 0, and MCMBoxModel::execute() calls solve() to run KPP's
 * Rosenbrock integrator for the full timestep.
 *
 * Available only when the module is compiled with KPP_ENABLED.
 */
class KppBoxIntegrator : public BoxIntegrator, public ConsoleStreamInterface
{
public:
  KppBoxIntegrator(MooseApp & app,
                   const std::string & mech_name,
                   Real rtol = 1e-4,
                   Real atol = 1e-6,
                   const std::string & solver_type = "rosenbrock");

  ~KppBoxIntegrator() override;

  // -- BoxIntegrator interface (self-driven -> all zeros) --
  Real computeResidual(unsigned int /*species_idx*/,
                        const std::vector<Real> & /*C*/) const override { return 0.0; }
  Real computeJacobianDiagonal(unsigned int /*species_idx*/,
                                const std::vector<Real> & /*C*/) const override { return 1.0; }
  Real computeJacobianOffDiagonal(unsigned int /*species_idx*/,
                                   unsigned int /*jvar*/,
                                   const std::vector<Real> & /*C*/) const override { return 0.0; }
  void reinit(Real /*time*/) const override {}
  bool selfDriven() const override { return true; }
  Real ppbToMolec() const override { return 1.0; }

  /// Advance the chemistry from t0 to t1 using KPP INTEGRATE().
  void solve(Real t0,
             Real t1,
             std::vector<Real> & C,
             const std::map<std::string, Real> & globals = {}) const;

  const std::vector<std::string> & speciesNames() const { return _species_names; }

private:
  const Real _rtol;
  const Real _atol;
  const std::string _solver_type;

  void * _lib_handle;

  // Function pointer types for the KPP C bridge API
  using InitFn      = void (*)(void);
  using IntegrateFn = int (*)(double *, double, double, double, double);
  using SetConcFn   = void (*)(double *, int);
  using GetConcFn   = void (*)(double *, int);
  using GetIntFn    = int (*)(void);
  using NameFn      = const char * (*)(int);

  void setGlobal(const std::string & name, Real value) const;

  InitFn      _kpp_init;
  IntegrateFn _kpp_integrate;
  SetConcFn   _kpp_set_conc;
  GetConcFn   _kpp_get_conc;
  GetIntFn    _kpp_get_nvar;
  NameFn      _kpp_get_species_name;
  std::vector<std::string> _species_names;
};
