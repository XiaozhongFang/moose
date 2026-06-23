[Mesh]
  [gen]  type = GeneratedMeshGenerator  dim = 1  nx = 20  xmin = 0  xmax = 1  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
  include_transport = true
[]

[Kernels]
  [ONE_fix] type = MCMConstraintKernel  variable = ONE  function = '1.0'  []
  [RO2_fix] type = MCMConstraintKernel  variable = RO2  function = '0.0'  []
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 1  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.1  []
  [C_ic]   type = ConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ConstantIC  variable = RO2 value = 0 []
[]

[BCs]
  [A_left]  type = DirichletBC  variable = A  boundary = left   value = 1  []
  [B_left]  type = DirichletBC  variable = B  boundary = left   value = 0.1  []
  [C_left]  type = DirichletBC  variable = C  boundary = left   value = 0  []
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
  num_steps = 50
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  file_base = 'coupled_tutorial_out'
[]
