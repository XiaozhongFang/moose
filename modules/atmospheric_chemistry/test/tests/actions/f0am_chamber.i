# 纯0D化学对比测试: 5物种6反应 (vs F0AM ode15s gold)
# Chamber conditions: T=298K, P=1013 mbar, no transport
# Matching F0AM ExampleSetup_Chamber.m for 5-species mechanism

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[MCMFacsimileAction]
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 1  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.1  []
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
  file_base = 'f0am_chamber_out'
[]
