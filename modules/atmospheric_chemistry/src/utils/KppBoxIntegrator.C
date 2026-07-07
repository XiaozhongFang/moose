// kpp/adapter/KppBoxIntegrator.C
#include "KppBoxIntegrator.h"
#include "MooseApp.h"

#include <dlfcn.h>
#include <string>

KppBoxIntegrator::KppBoxIntegrator(MooseApp & app,
                                     const std::string & mech_name,
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
    _kpp_set_conc(nullptr),
    _kpp_get_conc(nullptr)
{
  // Derive .so path: KPP_LIB env var takes precedence, then auto-discover
  std::string lib_name;
  const char * kpp_lib_env = std::getenv("KPP_LIB");
  if (kpp_lib_env)
  {
    lib_name = kpp_lib_env;
  }
  else if (!mech_name.empty())
  {
    // Convention: <mech_dir>/kpp_build_<mech>/libkpp_<mech>.so
    // mech_name is the full path to the .kpp file — strip extension for the name
    auto slash = mech_name.find_last_of('/');
    std::string dir = (slash != std::string::npos) ? mech_name.substr(0, slash + 1) : "";
    std::string base = mech_name.substr(slash + 1);
    auto dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    lib_name = dir + "kpp_build_" + base + "/libkpp_" + base + ".so";
  }
  else
  {
    lib_name = "libkpp_generated.so";
  }

  _lib_handle = dlopen(lib_name.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!_lib_handle)
  {
    _console << "KppBoxIntegrator: failed to dlopen " << lib_name
             << ": " << dlerror() << std::endl;
    mooseError("KppBoxIntegrator: cannot load KPP shared library: ", lib_name);
  }

  // Resolve symbols (plain C names, no Fortran underscore mangling)
  _kpp_init      = reinterpret_cast<InitFn>(dlsym(_lib_handle, "kpp_init"));
  _kpp_integrate = reinterpret_cast<IntegrateFn>(dlsym(_lib_handle, "kpp_integrate"));
  _kpp_get_conc  = reinterpret_cast<GetConcFn>(dlsym(_lib_handle, "kpp_get_conc"));
  _kpp_set_conc  = reinterpret_cast<SetConcFn>(dlsym(_lib_handle, "kpp_set_conc"));

  if (!_kpp_init || !_kpp_integrate || !_kpp_get_conc || !_kpp_set_conc)
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

  // Allocate KPP-style double array
  std::vector<double> Y(n);
  for (int i = 0; i < n; ++i)
    Y[i] = static_cast<double>(C[i]);

  // Copy into KPP global state (kpp_set_conc uses memcpy into C[NSPEC])
  _kpp_set_conc(Y.data(), n);

  // Call KPP integration wrapper
  int ierr = _kpp_integrate(Y.data(),
                             static_cast<double>(t0),
                             static_cast<double>(t1),
                             static_cast<double>(_rtol),
                             static_cast<double>(_atol));

  if (ierr < 0)
  {
    _console << "KppBoxIntegrator: integration failed, ierr=" << ierr
             << " t=[" << t0 << "," << t1 << "]" << std::endl;
    mooseError("KppBoxIntegrator: Rosenbrock integration failed with ierr=", ierr,
               " at t=[", t0, ", ", t1, "]. "
               "Check mechanism tolerances or initial conditions.");
  }

  // Read integrated concentration back from KPP global state
  _kpp_get_conc(Y.data(), n);

  // Copy back to MOOSE Real vector with negative clamp
  for (int i = 0; i < n; ++i)
  {
    C[i] = static_cast<Real>(Y[i]);
    if (C[i] < 0.0 && C[i] > -1.0e-20)
      C[i] = 0.0;
  }

  _console << "KppBoxIntegrator: t=[" << t0 << "," << t1
           << "] ierr=" << ierr
           << " | C[0]=" << (n > 0 ? C[0] : 0.0)
           << std::endl;
}
