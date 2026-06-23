# vs F0AM ode15s — tutorial_5sp — coupled mode (FEM 0D transport)
# Reproduces F0AM LearnF0AM_ODE.mlx: 5 species, 6 reactions
# Chamber conditions: T=298K, no photolysis, no transport
#
# Reference: F0AM ExampleSetup_Chamber.m with tutorial_5sp mechanism
# Gold file: vs_F0AM_tutorial5_coupled.csv

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 100  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.05  []
  [C_ic]   type = ConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ConstantIC  variable = RO2 value = 0 []
[]

[Postprocessors]
  [A_avg]  type = ElementAverageValue  variable = A  []
  [B_avg]  type = ElementAverageValue  variable = B  []
  [C_avg]  type = ElementAverageValue  variable = C  []
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
  file_base = 'vs_F0AM_tutorial5_coupled'
[]
