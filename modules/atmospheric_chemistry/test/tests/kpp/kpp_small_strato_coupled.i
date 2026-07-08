# kpp_small_strato_coupled.i -- KPP small_strato coupled-mode smoke test

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'kpp_small_strato/small_strato.kpp'
    temperature = 270.0
    air_density = 8.120e16
    chem_solver = kpp_rosenbrock
  []
[]

[ICs]
  [O1D_ic]
    type = ConstantIC
    variable = O1D
    value = 9.906e1
  []
  [O_ic]
    type = ConstantIC
    variable = O
    value = 6.624e8
  []
  [O3_ic]
    type = ConstantIC
    variable = O3
    value = 5.326e11
  []
  [NO_ic]
    type = ConstantIC
    variable = NO
    value = 8.725e8
  []
  [NO2_ic]
    type = ConstantIC
    variable = NO2
    value = 2.240e8
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = implicit-euler
  start_time = 43200
  end_time = 45000
  dt = 900
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
  exodus = false
[]
