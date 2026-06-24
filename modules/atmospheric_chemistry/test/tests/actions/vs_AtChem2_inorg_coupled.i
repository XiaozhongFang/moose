# vs AtChem2 CVODE — MCM inorganic — coupled mode (FEM 0D transport)
# Same mechanism and conditions as vs_AtChem2_inorg_box.i.
# Both tests shall produce identical results — only the implementation
# differs (ScalarVariable ODE vs FEM with 1-element mesh).
#
# Gold file: vs_AtChem2_inorg_box.csv (shared with box-mode test)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/atchem2_example.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  press = 1013.25
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [CH4_ic]  type = ConstantIC  variable = CH4  value = 4.9e13  []
  [CO_ic]   type = ConstantIC  variable = CO   value = 3.6e12  []
  [O3_ic]   type = ConstantIC  variable = O3   value = 5.2e11  []
  [NO2_ic]  type = ConstantIC  variable = NO2  value = 2.4e11  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 90
  end_time = 43200
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-3
  [TimeStepper]
    type = ConstantDT
    dt = 90
  []
[]

[Postprocessors]
  [O3_val]    type = ElementAverageValue  variable = O3   []
  [NO_val]    type = ElementAverageValue  variable = NO   []
  [NO2_val]   type = ElementAverageValue  variable = NO2  []
  [OH_val]    type = ElementAverageValue  variable = OH   []
  [HO2_val]   type = ElementAverageValue  variable = HO2  []
  [CO_val]    type = ElementAverageValue  variable = CO   []
  [CH4_val]   type = ElementAverageValue  variable = CH4  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_AtChem2_inorg_coupled'
[]
