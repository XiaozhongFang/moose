# RO2 species detection validation test for a simplified mechanism.
#
# tutorial_5sp.fac does not contain an explicit "RO2 = ..." declaration. It
# declares RO2 directly as a mechanism species, so the parser must detect it via
# the fallback name-based path.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
    mcm_photolysis_file = ''
    temperature = 298
  []
[]

[ICs]
  [A_ic]   type = ScalarConstantIC  variable = A   value = 100  []
  [B_ic]   type = ScalarConstantIC  variable = B   value = 0.05 []
  [C_ic]   type = ScalarConstantIC  variable = C   value = 0    []
  [ONE_ic] type = ScalarConstantIC  variable = ONE value = 1    []
  [RO2_ic] type = ScalarConstantIC  variable = RO2 value = 0    []
[]

[VectorPostprocessors]
  [ro2_list]
    type = MCMRO2ListPostprocessor
    box_model = box_model
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 1
  dt = 1
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'test_ro2_detection_tutorial_5sp'
[]
