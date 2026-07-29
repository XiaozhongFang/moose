[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -30
    xmax = 30
    ymin = -30
    ymax = 30
    zmin = 0
    zmax = 20
    nx = 2
    ny = 2
    nz = 2
  []
[]

[Variables]
  [HNO3]
  []
  [NO]
  []
  [O3]
  []
[]

[ICs]
  [HNO3_ic]
    type = MAS1998SpeciesIC
    variable = HNO3
    species = HNO3
  []
  [NO_ic]
    type = MAS1998SpeciesIC
    variable = NO
    species = NO
  []
  [O3_ic]
    type = MAS1998SpeciesIC
    variable = O3
    species = O3
  []
[]

[Functions]
  [K_profile]
    type = MAS1998VerticalDiffusivityFunction
  []
  [air_density]
    type = MAS1998AirNumberDensityFunction
  []
  [u_wind]
    type = MAS1998SolidBodyWindComponent
    component = u
  []
  [v_wind]
    type = MAS1998SolidBodyWindComponent
    component = v
  []
  [NO_emission]
    type = MAS1998NOEmissionFunction
  []
[]

[Postprocessors]
  [K_10km]
    type = FunctionValuePostprocessor
    function = K_profile
    point = '0 0 10'
    execute_on = INITIAL
  []
  [K_16km]
    type = FunctionValuePostprocessor
    function = K_profile
    point = '0 0 16'
    execute_on = INITIAL
  []
  [K_18km]
    type = FunctionValuePostprocessor
    function = K_profile
    point = '0 0 18'
    execute_on = INITIAL
  []
  [K_22km]
    type = FunctionValuePostprocessor
    function = K_profile
    point = '0 0 22'
    execute_on = INITIAL
  []
  [rho_0km]
    type = FunctionValuePostprocessor
    function = air_density
    point = '0 0 0'
    execute_on = INITIAL
  []
  [rho_20km]
    type = FunctionValuePostprocessor
    function = air_density
    point = '0 0 20'
    execute_on = INITIAL
  []
  [u_center]
    type = FunctionValuePostprocessor
    function = u_wind
    point = '0 0 0'
    execute_on = INITIAL
  []
  [v_center]
    type = FunctionValuePostprocessor
    function = v_wind
    point = '0 0 0'
    execute_on = INITIAL
  []
  [NO_emit_lowest]
    type = FunctionValuePostprocessor
    function = NO_emission
    point = '0 0 0.3'
    execute_on = INITIAL
  []
  [NO_emit_above]
    type = FunctionValuePostprocessor
    function = NO_emission
    point = '0 0 1.0'
    execute_on = INITIAL
  []
  [HNO3_max]
    type = NodalExtremeValue
    variable = HNO3
    value_type = max
    execute_on = INITIAL
  []
  [NO_max]
    type = NodalExtremeValue
    variable = NO
    value_type = max
    execute_on = INITIAL
  []
  [O3_max]
    type = NodalExtremeValue
    variable = O3
    value_type = max
    execute_on = INITIAL
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Problem]
  solve = false
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = INITIAL
  file_base = mas1998_physics
[]
