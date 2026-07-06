// kpp/adapter/KppBoxIntegrator.h
#pragma once

#include "BoxIntegrator.h"
#include "ConsoleStreamInterface.h"
#include "KppFortranBridge.h"

#include <string>
#include <vector>

class MooseApp;

/**
 * KPP-generated code integrator for box-mode chemistry.
 *
 * Wraps KPP-generated Fortran code via dlopen/dlsym.  Self-driven mode
 * (selfDriven() == true) — all BoxIntegrator residual/Jacobian methods
 * return 0, and MCMBoxModel::execute() calls solve() to run KPP's
 * INTEGRATE() subroutine for the full timestep.
 *
 * Available only when the module is compiled with KPP_ENABLED.
 */
class KppBoxIntegrator : public BoxIntegrator, public ConsoleStreamInterface
{
public:
  KppBoxIntegrator(MooseApp & app,
                   Real rtol = 1e-4,
                   Real atol = 1e-6,
                   const std::string & solver_type = "rosenbrock");

  ~KppBoxIntegrator() override;

  // -- BoxIntegrator interface (self-driven -> all zeros) --
  Real computeResidual(unsigned int species_idx,
                        const std::vector<Real> & C) const override { return 0.0; }
  Real computeJacobianDiagonal(unsigned int species_idx,
                                const std::vector<Real> & C) const override { return 0.0; }
  Real computeJacobianOffDiagonal(unsigned int species_idx,
                                   unsigned int jvar,
                                   const std::vector<Real> & C) const override { return 0.0; }
  void reinit(Real time) const override {}
  bool selfDriven() const override { return true; }
  Real ppbToMolec() const override { return 1.0; }

  /// Advance the chemistry from t0 to t1 using KPP INTEGRATE().
  /// @param t0  Start time (s)
  /// @param t1  End time (s)
  /// @param C   On entry: current concentrations; on exit: integrated concentrations
  void solve(Real t0, Real t1, std::vector<Real> & C) const;

private:
  const Real _rtol;
  const Real _atol;
  const std::string _solver_type;

  // dlopen handle and function pointers
  void * _lib_handle;

  // Resolved function pointers (set in constructor)
  using InitFn      = void (*)(void);
  using IntegrateFn = void (*)(double *, double *, int *, double *, double *, int *);
  using GetConcFn   = void (*)(double *, int *);
  using SetConcFn   = void (*)(double *, int *);
  using UpdateRFn   = void (*)(void);

  InitFn      _kpp_init;
  IntegrateFn _kpp_integrate;
  GetConcFn   _kpp_get_conc;
  SetConcFn   _kpp_set_conc;
  UpdateRFn   _kpp_update_rconst;
};
