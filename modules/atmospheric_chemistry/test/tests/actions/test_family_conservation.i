# Test: Family conservation (DAE method)
# Uses atchem2_example.fac mechanism which has NO, NO2, O3, etc.
# Family: NOx = NO2 + NO  (slack = NO2)

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
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  temperature = 298.15
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010

  # Family conservation: NOx = NO2 + NO (slack=NO2, member=NO)
  family_names = 'NOx'
  family_members = 'NO2 NO'
  family_scaling = '1.0 1.0'
[]

[ICs]
  [NO_ic]   type = ScalarConstantIC  variable = NO   value = 1.0e10 []
  [NO2_ic]  type = ScalarConstantIC  variable = NO2  value = 5.0e11 []
  [NO3_ic]  type = ScalarConstantIC  variable = NO3  value = 0 []
  [N2O5_ic] type = ScalarConstantIC  variable = N2O5 value = 0 []
  [O3_ic]   type = ScalarConstantIC  variable = O3   value = 1.0e12 []
  [HONO_ic] type = ScalarConstantIC  variable = HONO value = 0 []
  [HNO3_ic] type = ScalarConstantIC  variable = HNO3 value = 0 []
  [HO2NO2_ic] type = ScalarConstantIC variable = HO2NO2 value = 0 []
  [OH_ic]   type = ScalarConstantIC  variable = OH   value = 1.0e6 []
  [HO2_ic]  type = ScalarConstantIC  variable = HO2  value = 1.0e8 []
  [H2O2_ic] type = ScalarConstantIC  variable = H2O2 value = 1.0e10 []
  [CO_ic]   type = ScalarConstantIC  variable = CO   value = 5.0e13 []
  [CH4_ic]  type = ScalarConstantIC  variable = CH4  value = 5.0e13 []
  [HCHO_ic] type = ScalarConstantIC  variable = HCHO value = 1.0e10 []
  [CH3O2_ic] type = ScalarConstantIC variable = CH3O2 value = 0 []
  [CH3OOH_ic] type = ScalarConstantIC variable = CH3OOH value = 0 []
  [CH3OH_ic] type = ScalarConstantIC variable = CH3OH value = 0 []
  [CH3NO3_ic] type = ScalarConstantIC variable = CH3NO3 value = 0 []
  [O_ic]    type = ScalarConstantIC  variable = O    value = 0 []
  [O1D_ic]  type = ScalarConstantIC  variable = O1D  value = 0 []
  [CL_ic]   type = ScalarConstantIC  variable = CL   value = 0 []
  [SA_ic]  type = ScalarConstantIC  variable = SA  value = 0 []
  [SO3_ic] type = ScalarConstantIC  variable = SO3 value = 0 []
  [SO2_ic] type = ScalarConstantIC  variable = SO2 value = 0 []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 100
  dt = 10
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'TIMESTEP_END'
[]
