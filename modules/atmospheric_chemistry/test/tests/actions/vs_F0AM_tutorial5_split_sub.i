# Sub-app for operator-split test: chemistry box model
# Receives concentrations from parent (NS transport), solves ODE, returns results
# Uses mode=box (ScalarVariable + ChemistryODEKernel)

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
  dt = 0.001
  num_steps = 1
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  checkpoint = false
  csv = true
[]
