// kpp/adapter/kpp_adapter.c
// Lightweight C wrapper providing a simplified interface for KPP-generated code.
//
// Compiled together with KPP-generated .c files for a specific mechanism.
// The KppBoxIntegrator dlopen()s the resulting .so and resolves these symbols.

#include "small_test_Parameters.h"
#include "small_test_Global.h"
#include "small_test_Sparse.h"

#include <string.h>

// --- KPP C API (generated) ---
extern void Initialize(void);
extern void Fun(double Y[], double FIX[], double RCONST[], double Ydot[]);
extern void Jac_SP(double Y[], double FIX[], double RCONST[], double JVS[]);
extern int Rosenbrock(double Y[], double Tstart, double Tend,
                      void (*ode_Fun)(double, double[], double[]),
                      void (*ode_Jac)(double, double[], double[]),
                      double RPAR[], int IPAR[]);

// --- These functions are called by Rosenbrock internally ---
void FunTemplate(double T, double Y[], double Ydot[]);
void JacTemplate(double T, double Y[], double JVS[]);

// ============================================================
// Bridge API — resolved by KppBoxIntegrator via dlsym
// ============================================================

/// Initialize KPP global state.
void kpp_init(void) {
    Initialize();
}

/// Set concentration array C[NSPEC] from external source.
void kpp_set_conc(double c[], int n) {
    int m = (n < NSPEC) ? n : NSPEC;
    memcpy(C, c, m * sizeof(double));
}

/// Get concentration array from KPP global state.
void kpp_get_conc(double c[], int n) {
    int m = (n < NSPEC) ? n : NSPEC;
    memcpy(c, C, m * sizeof(double));
}

/// Integrate chemistry from t0 to t1.
/// Y[] is concentration array. Updated in place.
/// Returns 0 on success, negative on failure.
int kpp_integrate(double Y[], double t0, double t1,
                   double rtol, double atol) {
    // KPP Rosenbrock uses RPAR/IPAR for configuration
    double RPAR[20];
    int IPAR[20];
    memset(RPAR, 0, sizeof(RPAR));
    memset(IPAR, 0, sizeof(IPAR));

    // Set tolerances via RPAR (Rosenbrock-specific)
    // Note: Rosenbrock's RPAR isn't standard — the tolerances are typically
    // set via global variables ATOL/RTOL in the KPP model.
    // This call uses KPP's default tolerances if ATOL/RTOL aren't set.

    return Rosenbrock(Y, t0, t1, FunTemplate, JacTemplate, RPAR, IPAR);
}
