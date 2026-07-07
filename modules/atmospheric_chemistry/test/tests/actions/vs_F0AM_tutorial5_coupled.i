# vs F0AM ode15s — tutorial_5sp — coupled mode (FEM 0D transport)
# Reproduces F0AM LearnF0AM_ODE.mlx: 5 species, 6 reactions
# tspan=[0,10]s · BDF2 · RelTol=1e-3 · AbsTol=1e-6
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
  [Coupled]
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
  []
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 100  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.05  []
  [C_ic]   type = ConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ConstantIC  variable = RO2 value = 0 []
[]

[Postprocessors]
  [A]  type = ElementAverageValue  variable = A  []
  [B]  type = ElementAverageValue  variable = B  []
  [C]  type = ElementAverageValue  variable = C  []
  [RO2]  type = ElementAverageValue  variable = RO2  []
  [ONE]  type = ElementAverageValue  variable = ONE  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 400
  dt = 0.1
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = ConstantDT
    dt = 0.1
  []
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_tutorial5_coupled'
[]
