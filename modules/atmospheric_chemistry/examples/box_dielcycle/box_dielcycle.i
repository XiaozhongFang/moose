# F0AM-style diurnal cycle box model example
# Simulates 24h diurnal cycle with full MCM v3.3.1 subset
# Mechanism: MCMv331_DielExampleChemistry (2908 species, 8797 reactions)
# Photolysis: MCM_SZA with diurnal variation, lat=51.51, lon=0.13

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = 'doc/content/modules/atmospheric_chemistry/database/MCMv331_DielExampleChemistry.fac'
  mcm_photolysis_file = 'doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  temperature = 288.15
  press = 1013.25
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
  units = ppb
  jfac = 1.0
[]

[ICs]
  [O3_ic]   type = ScalarConstantIC  variable = O3   value = 30.0 []
  [NO_ic]   type = ScalarConstantIC  variable = NO   value = 0.1 []
  [NO2_ic]  type = ScalarConstantIC  variable = NO2  value = 0.5 []
  [CO_ic]   type = ScalarConstantIC  variable = CO   value = 150.0 []
  [CH4_ic]  type = ScalarConstantIC  variable = CH4  value = 1900.0 []
  [H2O_ic]  type = ScalarConstantIC  variable = H2O  value = 2.0e17 []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 86400
  dt = 900
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  csv = true
  time_interval = 6
  execute_on = 'TIMESTEP_END'
[]
