# 0D smoke test for the recovered MAS1998 methane mechanism.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = 'chemistry/mas1998_methane.kpp'
    temperature = 288.15
    air_density = 2.55e19
    water_vapor = 2.55e17
    chem_solver = kpp_rosenbrock
  []
[]

[ICs]
  [O1D_ic]    type = ScalarConstantIC variable = O1D    value = 0.0 []
  [CH4_ic]    type = ScalarConstantIC variable = CH4    value = 4.335e13 []
  [HNO2_ic]   type = ScalarConstantIC variable = HNO2   value = 1.0e2 []
  [H2O2_ic]   type = ScalarConstantIC variable = H2O2   value = 1.0e2 []
  [N2O5_ic]   type = ScalarConstantIC variable = N2O5   value = 1.0e2 []
  [HNO3_ic]   type = ScalarConstantIC variable = HNO3   value = 2.55e9 []
  [HO2NO2_ic] type = ScalarConstantIC variable = HO2NO2 value = 1.0e2 []
  [CH3OOH_ic] type = ScalarConstantIC variable = CH3OOH value = 1.0e2 []
  [HCHO_ic]   type = ScalarConstantIC variable = HCHO   value = 1.0e2 []
  [CH3O2_ic]  type = ScalarConstantIC variable = CH3O2  value = 1.0e2 []
  [NO3_ic]    type = ScalarConstantIC variable = NO3    value = 1.0e2 []
  [O3P_ic]    type = ScalarConstantIC variable = O3P    value = 0.0 []
  [NO_ic]     type = ScalarConstantIC variable = NO     value = 1.0e2 []
  [OH_ic]     type = ScalarConstantIC variable = OH     value = 1.0e2 []
  [NO2_ic]    type = ScalarConstantIC variable = NO2    value = 5.1e9 []
  [O3_ic]     type = ScalarConstantIC variable = O3     value = 7.65e11 []
  [HO2_ic]    type = ScalarConstantIC variable = HO2    value = 1.0e2 []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = implicit-euler
  start_time = 0
  end_time = 60
  dt = 60
  nl_max_its = 20
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_max_its = 50
  l_tol = 1e-10
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  checkpoint = false
  csv = false
  exodus = false
[]
