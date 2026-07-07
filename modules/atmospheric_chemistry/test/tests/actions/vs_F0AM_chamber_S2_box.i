# vs F0AM Chamber — MCM Inorg + Isoprene — box mode
# Full reproduction of F0AM ExampleSetup_Chamber:
#   610 species, 1974 reactions (MCM v3.3.1 subset)
#   Isoprene + NO2 + H2O2 chamber oxidation
#   BottomUp photolysis with lab lamp spectrum
#   T=298K, P=1013mbar, RH=10%
#   3h simulation, C5H8=10ppb, NO2=0.1ppb, H2O2=200ppb

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistryBox]
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/MCMv331_Inorg_Isoprene.fac'
  temperature = 298.0
  air_density = 2.46e19
  water_vapor = 3.12e17
  press = 1013.0
  photolysis_scheme = BOTTOMUP
  lamp_flux_file = 'ExampleLightFlux.txt'
  bottomup_data_dir = '../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
  jfac = 1.0
  chem_solver = petsc_ts
  chem_solver_type = bdf
  chem_solver_rtol = 1e-2
  chem_solver_atol = 1e-6
[]

[ICs]
  [C5H8_ic]    type = ScalarConstantIC  variable = C5H8  value = 2.47e11   []  # 10 ppb
  [NO2_ic]     type = ScalarConstantIC  variable = NO2   value = 2.47e10   []  # 1 ppb (S2)
  [H2O2_ic]    type = ScalarConstantIC  variable = H2O2  value = 4.94e12   []  # 200 ppb
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 10800      # 3h
  dt = 100              # 100s steps for stiff chemistry
  l_max_its = 200
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  [TimeStepper]
    type = ConstantDT
    dt = 100
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_chamber_S2_box'
  time_step_interval = 10
  [checkpoint]
    type = Checkpoint
    execute_on = 'final'
  []
[]
