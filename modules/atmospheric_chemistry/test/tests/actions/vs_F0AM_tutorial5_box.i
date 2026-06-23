# vs F0AM ode15s — tutorial_5sp — box mode (ScalarVariable 0D ODE)
# Reproduces F0AM LearnF0AM_ODE.mlx: 5 species, 6 reactions
# Chamber conditions: T=298K, no photolysis, no transport
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

[Postprocessors]
  [A_val]  type = ScalarVariable  variable = A  []
  [B_val]  type = ScalarVariable  variable = B  []
  [C_val]  type = ScalarVariable  variable = C  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 0.1
  num_steps = 10
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  file_base = 'vs_F0AM_tutorial5_box'
[]
