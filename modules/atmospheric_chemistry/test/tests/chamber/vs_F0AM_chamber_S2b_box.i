# vs F0AM Chamber — S2b (restart from S2, jcorr=10)
# Full reproduction of F0AM ExampleSetup_Chamber EVENTS section:
#   610 species, 1974 reactions (MCM v3.3.1 subset)
#   Restart from S2 (NO2=1ppb) end state
#   BottomUp photolysis with lab lamp spectrum ×10
#   T=298K, P=1013mbar, RH=10%
#   1h continuation, jcorr=10

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
    air_density = 2.4622e19
    water_vapor = 7.6114e16
    press = 1013.0
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = '../../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
    jfac = 10.0
    use_limiting_reagent = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  start_time = 0
  end_time = 3600       # 1h relative to restart
  dt = 100
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
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_chamber_S2b_box'
  time_step_interval = 10
[]

[Problem]
  restart_file_base = 'vs_F0AM_chamber_S2_box_cp/LATEST'
[]
