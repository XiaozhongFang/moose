#INCLUDE atoms.kpp

#DEFVAR
O1D    = O;
O3P    = O;
HNO2   = H + N + O + O;
HO2NO2 = H + N + O + O + O + O;
O3     = O + O + O;
HNO3   = H + N + O + O + O;
NO     = N + O;
NO2    = N + O + O;
NO3    = N + O + O + O;
N2O5   = N + N + O + O + O + O + O;
OH     = O + H;
HO2    = H + O + O;
H2O2   = H + H + O + O;
HCHO   = C + H + H + O;
CH3OOH = C + H + H + H + H + O + O;
CH3O2  = C + H + H + H + O + O;
CH4    = C + H + H + H + H;

#DEFFIX
CO2 = C + O + O;
H2O = H + H + O;
CO  = C + O;
H2  = H + H;
O2  = O + O;
M   = IGNORE;

#INITVALUES
ALL_SPEC = 1.0E2;
M        = 2.55E19;
H2O      = 2.55E17;
CO       = 2.55E12;
O2       = 5.3295E18;
HNO3     = 2.55E9;
O3       = 7.65E11;
CH4      = 4.335E13;
NO2      = 5.1E9;
O3P      = 0.0;
O1D      = 0.0;
CFACTOR  = 1.0;
