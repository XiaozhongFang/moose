// kpp/adapter/KppFortranBridge.h
#pragma once

#include "Moose.h"

/// C++ ↔ Fortran bridge for KPP-generated code.
///
/// KPP generates Fortran subroutines with lowercase names and trailing
/// underscore mangling (gfortran convention).  This header provides
/// extern "C" wrappers that the KppBoxIntegrator can call.
///
/// Each generated mechanism exports these symbols.  The KppBoxIntegrator
/// dlopen()s the corresponding .so and resolves them via dlsym().

extern "C" {

/// Initialize KPP global state (called once at construction).
void kpp_init_(void);

/// Main integration entry point.
/// @param[in]     tin     Start time (s)
/// @param[in]     tout    End time (s)
/// @param[in,out] icntrl  Integer control array (ICNTRL, length 20)
/// @param[in,out] rcntrl  Real control array (RCNTRL, length 20)
/// @param[out]    rstatus Real status array (RSTATUS, length 20)
/// @param[out]    ierr    Error flag (0 = success)
void kpp_integrate_(double *tin, double *tout,
                     int *icntrl, double *rcntrl,
                     double *rstatus, int *ierr);

/// Get current concentration array.
/// @param[out] c   Concentration vector (length n)
/// @param[in]  n   Number of species
void kpp_get_conc_(double *c, int *n);

/// Set concentration array from external source.
/// @param[in] c   Concentration vector (length n)
/// @param[in] n   Number of species
void kpp_set_conc_(double *c, int *n);

/// Update rate coefficients (temperature, density, photolysis).
/// Called at the start of each MOOSE timestep before INTEGRATE.
void kpp_update_rconst_(void);

/// Update sunlight intensity (solar zenith angle).
void kpp_update_sun_(void);

}  // extern "C"
