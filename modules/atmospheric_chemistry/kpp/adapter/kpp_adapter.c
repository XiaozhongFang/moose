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

// KPP global state (defined in _Global.h by KPP-generated code)
extern double C[];
extern double ATOL[];
extern double RTOL[];
extern double STEPMIN;
extern double STEPMAX;
extern char * SPC_NAMES[];

// KPP functions (names fixed by KPP tool, independent of mechanism)
extern void Initialize(void);
#ifndef KPP_USE_GENERATED_INTEGRATE
extern int Rosenbrock(double Y[], double Tstart, double Tend,
                      double AbsTol[], double RelTol[],
                      void (*ode_Fun)(double, double[], double[]),
                      void (*ode_Jac)(double, double[], double[]),
                      double RPAR[], int IPAR[]);
#else
extern void INTEGRATE(double TIN, double TOUT);
#endif

// Injected into _Main.c by kpp/build/Makefile: int kpp_get_nspec(void) { return NSPEC; }
extern int kpp_get_nspec(void);
extern int kpp_get_nvar(void);

// Templates called by Rosenbrock internally
void FunTemplate(double T, double Y[], double Ydot[]);
void JacTemplate(double T, double Y[], double JVS[]);

// ============================================================
// Bridge API — resolved by KppBoxIntegrator via dlsym
// ============================================================

void kpp_init(void) {
    Initialize();
}

const char * kpp_get_species_name(int i) {
    int nspec = kpp_get_nspec();
    if (i < 0 || i >= nspec)
        return "";
    return SPC_NAMES[i] ? SPC_NAMES[i] : "";
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
    int nvar = kpp_get_nvar();
    for (int i = 0; i < nvar; ++i) {
        RTOL[i] = rtol;
        ATOL[i] = atol;
    }
    STEPMIN = fmax(fabs(t1 - t0) * 1.0e-6, 1.0e-12);
    STEPMAX = fabs(t1 - t0);
    kpp_set_conc(Y, nvar);

#ifdef KPP_USE_GENERATED_INTEGRATE
    INTEGRATE(t0, t1);
    kpp_get_conc(Y, nvar);
    return 1;
#else
    double RPAR[20];
    int IPAR[20];
    memset(RPAR, 0, sizeof(RPAR));
    memset(IPAR, 0, sizeof(IPAR));
    IPAR[0] = 0;
    IPAR[1] = 1;
    IPAR[3] = 5;
    RPAR[2] = STEPMIN;
    int ierr = Rosenbrock(Y, t0, t1, ATOL, RTOL,
                          FunTemplate, JacTemplate, RPAR, IPAR);
    if (ierr > 0) {
        STEPMIN = RPAR[11];
        kpp_set_conc(Y, nvar);
    }
    return ierr;
#endif
}
