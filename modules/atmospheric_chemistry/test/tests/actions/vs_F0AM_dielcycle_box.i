# vs F0AM DielCycle — SOAS 2013 campaign diurnal cycle — box mode
# Full reproduction of F0AM ExampleSetup_DielCycle:
#   2908 species, 8797 reactions (MCM v3.3.1 subset)
#   36 constrained species from SOAS observations
#   3×24h spinup, 3600s IntTime, MCM SZA photolysis
#   Centerville, AL: 32.903°N, 87.250°W, 2013-06-30
#   jcorr=0.5, kdil=1/(24*3600) s⁻¹

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/MCMv331_DielExampleChemistry.fac'
  temperature = 296.914
  air_density = 2.44594e19
  water_vapor = 2.08e17
  press = 751.684
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  photolysis_scheme = MCM_SZA
  latitude = 32.903
  longitude = -87.250
  day = 30
  month = 6
  year = 2013
  jfac = 0.5
  dilute = 1.1574074074074073e-5
[]

[ICs]
  # ── Free-evolving species (HoldMe=0 in F0AM) ──
  [O3_ic]    type = ScalarConstantIC  variable = O3    value = 5.273e11  []  # 21.34 ppb
  [OH_ic]    type = ScalarConstantIC  variable = OH    value = 1.726e5   []  # 6.99e-6 ppb
  [NO_ic]    type = ScalarConstantIC  variable = NO    value = 1.144e8   []  # 0.00463 ppb
  [NO2_ic]   type = ScalarConstantIC  variable = NO2   value = 1.651e10  []  # 0.668 ppb
  # ── Key constrained species (HoldMe=1) — initial values from SOAS hour 0 ──
  [H2_ic]    type = ScalarConstantIC  variable = H2    value = 1.358e13  []  # 550 ppb
  [CO_ic]    type = ScalarConstantIC  variable = CO    value = 3.790e12  []  # 153.5 ppb
  [CH4_ic]   type = ScalarConstantIC  variable = CH4   value = 4.371e13  []  # 1770 ppb
  [H2O2_ic]  type = ScalarConstantIC  variable = H2O2  value = 3.506e9   []  # 0.142 ppb
  [C5H8_ic]  type = ScalarConstantIC  variable = C5H8  value = 9.288e10  []  # 3.76 ppb
  [HCHO_ic]  type = ScalarConstantIC  variable = HCHO  value = 5.676e10  []  # 2.30 ppb (approx)
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 259200   # 72h = 3 × 24 × 3600
  dt = 3600           # 1h steps
  l_max_its = 50
  l_tol = 1e-3
  nl_max_its = 8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  [TimeStepper]
    type = ConstantDT
    dt = 3600
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_dielcycle_box'
  time_step_interval = 1
[]
