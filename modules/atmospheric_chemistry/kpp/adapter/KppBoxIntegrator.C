// kpp/adapter/KppBoxIntegrator.C
#include "KppBoxIntegrator.h"
#include "MooseApp.h"

#include <dlfcn.h>
#include <string>

KppBoxIntegrator::KppBoxIntegrator(MooseApp & app,
                                     Real rtol,
                                     Real atol,
                                     const std::string & solver_type)
  : BoxIntegrator(),
    ConsoleStreamInterface(app),
    _rtol(rtol),
    _atol(atol),
    _solver_type(solver_type),
    _lib_handle(nullptr),
    _kpp_init(nullptr),
    _kpp_integrate(nullptr),
    _kpp_get_conc(nullptr),
    _kpp_set_conc(nullptr),
    _kpp_update_rconst(nullptr)
{
  // Open the KPP shared library.
  // The library name follows the convention libkpp_<mech>.so, where <mech>
  // is the mechanism name.  For now the library path must be set via the
  // KPP_LIB environment variable or a default search path.
  const char * kpp_lib_path = std::getenv("KPP_LIB");
  std::string lib_name = kpp_lib_path ? std::string(kpp_lib_path) : "libkpp_generated.so";

  _lib_handle = dlopen(lib_name.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!_lib_handle)
  {
    _console << "KppBoxIntegrator: failed to dlopen " << lib_name
             << ": " << dlerror() << std::endl;
    mooseError("KppBoxIntegrator: cannot load KPP shared library: ", lib_name);
  }

  // Resolve symbols
  _kpp_init = reinterpret_cast<InitFn>(dlsym(_lib_handle, "kpp_init_"));
  _kpp_integrate = reinterpret_cast<IntegrateFn>(dlsym(_lib_handle, "kpp_integrate_"));
  _kpp_get_conc  = reinterpret_cast<GetConcFn>(dlsym(_lib_handle, "kpp_get_conc_"));
  _kpp_set_conc  = reinterpret_cast<SetConcFn>(dlsym(_lib_handle, "kpp_set_conc_"));
  _kpp_update_rconst = reinterpret_cast<UpdateRFn>(dlsym(_lib_handle, "kpp_update_rconst_"));

  if (!_kpp_init || !_kpp_integrate || !_kpp_get_conc || !_kpp_set_conc || !_kpp_update_rconst)
  {
    std::string msg = "KppBoxIntegrator: failed to resolve KPP symbols from ";
    msg += lib_name;
    dlclose(_lib_handle);
    _lib_handle = nullptr;
    mooseError(msg);
  }

  // Initialize KPP global state
  _kpp_init();

  _console << "KppBoxIntegrator: loaded " << lib_name
           << " (rtol=" << _rtol << ", atol=" << _atol
           << ", solver=" << _solver_type << ")" << std::endl;
}

KppBoxIntegrator::~KppBoxIntegrator()
{
  if (_lib_handle)
  {
    dlclose(_lib_handle);
    _lib_handle = nullptr;
  }
}

void
KppBoxIntegrator::solve(Real t0, Real t1, std::vector<Real> & C) const
{
  int n = static_cast<int>(C.size());
  if (n == 0)
    return;

  // Write concentration vector to KPP
  _kpp_set_conc(C.data(), &n);

  // Configure KPP integrator control arrays
  int icntrl[20] = {0};
  double rcntrl[20] = {0.0};
  double rstatus[20] = {0.0};

  // ICNTRL defaults (zero = KPP default behavior):
  //   ICNTRL(1) = 0: autonomous (default)
  //   ICNTRL(2) = 0: stop at tout (default)
  // RCNTRL defaults:
  //   RCNTRL(1) = _rtol (relative tolerance)
  //   RCNTRL(2) = _atol (absolute tolerance)
  //   RCNTRL(3) = STEPMIN (default 0.0 = KPP internal default)
  //   RCNTRL(4) = STEPMAX (default 0.0 = KPP internal default)
  rcntrl[0] = static_cast<double>(_rtol);  // relative tolerance
  rcntrl[1] = static_cast<double>(_atol);  // absolute tolerance

  // Call KPP INTEGRATE
  double tin  = static_cast<double>(t0);
  double tout = static_cast<double>(t1);
  int ierr = 0;

  _kpp_integrate(&tin, &tout, icntrl, rcntrl, rstatus, &ierr);

  if (ierr < 0)
  {
    _console << "KppBoxIntegrator: INTEGRATE failed, ierr=" << ierr
             << " t=[" << t0 << "," << t1 << "]" << std::endl;
  }

  // Read integrated concentration back
  _kpp_get_conc(C.data(), &n);

  // Clamp negative concentrations to zero (numerical artifacts)
  for (int i = 0; i < n; ++i)
    if (C[i] < 0.0 && C[i] > -1.0e-20)
      C[i] = 0.0;

  // Diagnostics
  long int nst = static_cast<long int>(rstatus[0]);  // steps taken
  _console << "KppBoxIntegrator: t=[" << t0 << "," << t1
           << "] nst=" << nst
           << " ierr=" << ierr
           << " | C[0]=" << (n > 0 ? C[0] : 0.0)
           << std::endl;
}
