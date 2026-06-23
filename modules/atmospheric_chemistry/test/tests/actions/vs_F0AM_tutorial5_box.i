# vs F0AM ode15s — tutorial_5sp — box mode (ScalarVariable 0D ODE)
# Reproduces F0AM LearnF0AM_ODE.mlx: 5 species, 6 reactions
# tspan=[0,10]s · BDF2 · RelTol=1e-3 · AbsTol=1e-6
# Matches ode15s default tolerances, BDF2 approximates ode15s BDF(1-5)
#
# Reference: F0AM ExampleSetup_Chamber.m with tutorial_5sp mechanism
# Gold file: vs_F0AM_tutorial5_box.csv

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
[]

[ICs]
  [A_ic]   type = ScalarConstantIC  variable = A   value = 100  []
  [B_ic]   type = ScalarConstantIC  variable = B   value = 0.05  []
  [C_ic]   type = ScalarConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ScalarConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ScalarConstantIC  variable = RO2 value = 0 []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 400
  dt = 0.1
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-6
  [TimeStepper]
    type = ConstantDT
    dt = 0.1
  []
[]

[Outputs]
  csv = true
  execute_on = 'initial timestep_end'
  file_base = 'vs_F0AM_tutorial5_box'
[]
