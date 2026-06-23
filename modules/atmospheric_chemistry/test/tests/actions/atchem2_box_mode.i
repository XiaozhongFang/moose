# MOOSE 盒子模式测试（ScalarVariable + ChemistryODEKernel）
# 使用 AtChem2 示例机制（MCM v3.3.1 子集）
# 对标 atchem2_validation.i，但使用 mode=box

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
  # SZA-based photolysis with AtChem2 defaults
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [CH4_ic]
    type = ScalarConstantIC
    variable = CH4
    value = 4.9e13
  []
  [CO_ic]
    type = ScalarConstantIC
    variable = CO
    value = 3.6e12
  []
  [O3_ic]
    type = ScalarConstantIC
    variable = O3
    value = 5.2e11
  []
  [NO2_ic]
    type = ScalarConstantIC
    variable = NO2
    value = 2.4e11
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 900
  end_time = 3600
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

[Postprocessors]
  [O3_val]
    type = ScalarVariable
    variable = O3
  []
  [NO_val]
    type = ScalarVariable
    variable = NO
  []
  [NO2_val]
    type = ScalarVariable
    variable = NO2
  []
  [OH_val]
    type = ScalarVariable
    variable = OH
  []
  [HO2_val]
    type = ScalarVariable
    variable = HO2
  []
  [CH4_val]
    type = ScalarVariable
    variable = CH4
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = atchem2_box_mode
  time_step_interval = 1
[]
