// kpp/adapter/KppFortranBridge.h
#pragma once

#include "Moose.h"

/// C bridge for KPP-generated code (via kpp_adapter.c).
///
/// KPP generates C code with function names Fun(), Initialize(), Rosenbrock(),
/// etc.  The kpp_adapter.c wrapper wraps these into the simplified API below.
/// KppBoxIntegrator resolves these symbols via dlsym from the compiled .so.
///
/// This header documents the bridge interface.  For the actual adapter
/// implementation see kpp/adapter/kpp_adapter.c.

extern "C" {

/// Initialize KPP global state (called once at construction).
void kpp_init(void);

/// Integrate chemistry from t0 to t1.
/// @param Y       Concentration array (in/out, length = NSPEC)
/// @param t0      Start time (s)
/// @param t1      End time (s)
/// @param rtol    Relative tolerance
/// @param atol    Absolute tolerance
/// @return        0 on success, negative on failure
int kpp_integrate(double Y[], double t0, double t1,
                   double rtol, double atol);

/// Number of KPP variable species.
int kpp_get_nvar(void);

/// KPP species name at a zero-based C[]/VAR[] index.
const char * kpp_get_species_name(int i);

/// Copy external concentrations into KPP global C[] array.
/// @param c   Concentration vector (length n)
/// @param n   Number of species
void kpp_set_conc(double c[], int n);

/// Copy KPP global C[] array to external buffer.
/// @param c   Output buffer (length n)
/// @param n   Number of species
void kpp_get_conc(double c[], int n);

}  // extern "C"
