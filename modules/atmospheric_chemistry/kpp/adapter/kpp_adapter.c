// kpp/adapter/kpp_adapter.c
// Generic C bridge for KPP-generated mechanisms.
// KPP always uses the same function/variable names across ALL mechanisms,
// so this file needs ZERO mechanism-specific includes.
//
// Compiled together with KPP-generated .c files for a specific mechanism.
// KppBoxIntegrator dlopen()s the resulting .so and resolves:
//   kpp_init, kpp_set_conc, kpp_get_conc, kpp_integrate

#include <string.h>
#include <math.h>

// KPP global concentration array (defined in _Global.h by KPP-generated code)
extern double C[];

// KPP functions (names fixed by KPP tool, independent of mechanism)
extern void Initialize(void);
extern int Rosenbrock(double Y[], double Tstart, double Tend,
                      void (*ode_Fun)(double, double[], double[]),
                      void (*ode_Jac)(double, double[], double[]),
                      double RPAR[], int IPAR[]);

// Injected into _Main.c by kpp/build/Makefile: int kpp_get_nspec(void) { return NSPEC; }
extern int kpp_get_nspec(void);

// Templates called by Rosenbrock internally
void FunTemplate(double T, double Y[], double Ydot[]);
void JacTemplate(double T, double Y[], double JVS[]);

// ============================================================
// Bridge API — resolved by KppBoxIntegrator via dlsym
// ============================================================

void kpp_init(void) {
    Initialize();
}

void kpp_set_conc(double c[], int n) {
    int nspec = kpp_get_nspec();
    int m = (n < nspec) ? n : nspec;
    memcpy(C, c, m * sizeof(double));
}

void kpp_get_conc(double c[], int n) {
    int nspec = kpp_get_nspec();
    int m = (n < nspec) ? n : nspec;
    memcpy(c, C, m * sizeof(double));
    // Zero any output elements beyond NSPEC (safety padding)
    for (int i = m; i < n; i++)
        c[i] = 0.0;
}

int kpp_integrate(double Y[], double t0, double t1,
                   double rtol, double atol) {
    double RPAR[20];
    int IPAR[20];
    memset(RPAR, 0, sizeof(RPAR));
    memset(IPAR, 0, sizeof(IPAR));
    return Rosenbrock(Y, t0, t1, FunTemplate, JacTemplate, RPAR, IPAR);
}
