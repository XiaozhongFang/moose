# F0AM EKMA-style sensitivity scan example (tutorial 5sp)
# Simplified EKMA: scans VOC × NOx scaling factors for O3 peak
# Uses the small tutorial_5sp mechanism for fast turnaround

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = box
  mechanism_file = 'doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  mcm_photolysis_file = ''
  temperature = 298.15
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
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = ConstantDT
    dt = 0.1
  []
[]

[Outputs]
  csv = true
  execute_on = 'FINAL'
[]
