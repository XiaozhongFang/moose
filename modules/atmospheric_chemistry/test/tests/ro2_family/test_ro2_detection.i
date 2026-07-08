# RO2 species detection validation test — MCM Inorg + Isoprene — box mode
#
# Purpose: Verify that the MCMBoxModel / MCMFacsimileParser correctly
# detects RO2 (peroxy radical) species from the .fac mechanism file's
# explicit RO2 = ... declaration line.
#
# The MCMRO2ListPostprocessor outputs one column per detected RO2 species name.
# The gold CSV is generated directly from the mechanism file's RO2 = line and
# CSVDiff verifies both the count and species-name set (117 species).
#
# Additional validation: the parser's internal validation against the
# peroxy-radicals reference file prints warnings for any detected RO2
# species not found in the MCM v3.3.1 reference list.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac'
    temperature = 298.0
    air_density = 2.46e19
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = '../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
    jfac = 1.0
  []
[]

[ICs]
  [C5H8_ic]    type = ScalarConstantIC  variable = C5H8  value = 2.47e11   []  # 10 ppb
  [NO2_ic]     type = ScalarConstantIC  variable = NO2   value = 2.47e9    []  # 0.1 ppb
  [H2O2_ic]    type = ScalarConstantIC  variable = H2O2  value = 4.94e12   []  # 200 ppb
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
  end_time = 100
  dt = 100
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'test_ro2_detection'
[]
