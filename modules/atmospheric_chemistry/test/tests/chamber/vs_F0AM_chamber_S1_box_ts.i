# Same config but relaxed TS tolerances for speed test
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac'
    temperature = 298.0
    air_density = 2.46e19
    water_vapor = 3.12e17
    press = 1013.0
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = '../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
    jfac = 1.0
    chem_solver = petsc_ts
    chem_solver_type = bdf
    chem_solver_rtol = 1e-3
    chem_solver_atol = 1e-6
  []
[]

[ICs]
  [C5H8_ic]    type = ScalarConstantIC  variable = C5H8  value = 2.47e11   []
  [NO2_ic]     type = ScalarConstantIC  variable = NO2   value = 2.47e9    []
  [H2O2_ic]    type = ScalarConstantIC  variable = H2O2  value = 4.94e12   []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 200
  dt = 100
  l_max_its = 200
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  [TimeStepper]
    type = ConstantDT
    dt = 100
  []
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_chamber_S1_box_ts_fast'
  time_step_interval = 1
[]
