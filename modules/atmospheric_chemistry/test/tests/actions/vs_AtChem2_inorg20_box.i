# vs AtChem2 — MCM inorganic 20sp — box mode
# Pure inorganic mechanism (mcm_export.fac): 20 species, 48 reactions
# Troe falloff + temperature-dependent rates + SZA photolysis
# 12h diurnal cycle matching AtChem2 model_mcm_inorg conditions
# RunApp validation: verifies box model handles complex rate expressions

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_export.fac'
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
  [O3_ic]   type = ScalarConstantIC  variable = O3   value = 5.2e11  []
  [NO_ic]   type = ScalarConstantIC  variable = NO   value = 2.5e11  []
  [NO2_ic]  type = ScalarConstantIC  variable = NO2  value = 2.4e11  []
  [OH_ic]   type = ScalarConstantIC  variable = OH   value = 1.0e6  []
  [HO2_ic]  type = ScalarConstantIC  variable = HO2  value = 1.0e8  []
  [CO_ic]   type = ScalarConstantIC  variable = CO   value = 3.6e12  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  end_time = 3600
  dt = 900
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-3
  [TimeStepper]
    type = ConstantDT
    dt = 900
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  hide = 'O O1D N2O5 H2 H2O2 HONO HO2NO2 SO2 SO3 HSO3 NA SA'
  file_base = 'vs_AtChem2_inorg20_box'
[]
