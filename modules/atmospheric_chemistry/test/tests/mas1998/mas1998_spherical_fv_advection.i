[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 360
    ymin = -80
    ymax = 80
    nx = 8
    ny = 4
  []
[]

[Variables]
  [tracer]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
  []
[]

[ICs]
  [tracer_ic]
    type = FunctionIC
    variable = tracer
    function = '1 + 0.05*sin(0.017453292519943295*x) + 0.02*sin(0.017453292519943295*y)'
  []
[]

[Functions]
  [u_wind]
    type = MAS1998SolidBodyWindComponent
    component = u
  []
  [v_wind]
    type = MAS1998SolidBodyWindComponent
    component = v
  []
[]

[FVKernels]
  [time]
    type = AtmosphericSphericalFVTimeDerivative
    variable = tracer
    coordinate_units = degrees
  []
  [advection]
    type = AtmosphericSphericalFVAdvection
    variable = tracer
    u_wind = u_wind
    v_wind = v_wind
    coordinate_units = degrees
    longitude_period = 360
    advected_interp_method = upwind
    flux_scheme = mas1998_limited
  []
[]

[Postprocessors]
  [tracer_avg]
    type = ElementAverageValue
    variable = tracer
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = implicit-euler
  num_steps = 1
  dt = 10
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-10
  l_max_its = 50
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  checkpoint = false
  csv = true
  exodus = false
[]
