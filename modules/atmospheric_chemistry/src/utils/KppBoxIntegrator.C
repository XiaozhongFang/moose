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
    _kpp_set_conc(nullptr),
    _kpp_get_conc(nullptr)
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
