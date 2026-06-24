# vs AtChem2 CVODE — MCM inorganic — box mode
# Full reproduction of AtChem2 model_mcm_inorg:
#   29 species, 71 reactions, 12h diurnal cycle (48 steps × 900s)
#   SZA photolysis, T=288.15K, P=1013.25mbar
#   lat=51.51, lon=0.13, 2010-06-21

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/atchem2_example.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [CH4_ic]  type = ScalarConstantIC  variable = CH4  value = 4.9e13  []
  [CO_ic]   type = ScalarConstantIC  variable = CO   value = 3.6e12  []
  [O3_ic]   type = ScalarConstantIC  variable = O3   value = 5.2e11  []
  [NO2_ic]  type = ScalarConstantIC  variable = NO2  value = 2.4e11  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 43200
  dt = 900
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = ConstantDT
    dt = 900
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_AtChem2_inorg_box'
  time_step_interval = 1
[]
