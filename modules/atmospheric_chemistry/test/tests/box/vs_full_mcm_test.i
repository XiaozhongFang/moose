# Full MCM v3.3.1 mechanism — quick scalability test
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]
[AtmosphericChemistry]
  [Box]
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_export_all.fac'
    temperature = 298.0
    air_density = 2.46e19
    press = 1013.0
    photolysis_scheme = MCM_SZA
    chem_solver = petsc_ts
    chem_solver_type = bdf
    chem_solver_rtol = 1e-1
    chem_solver_atol = 1e-4
    default_ic = 1.0e6
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 100
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
  file_base = 'vs_full_mcm_test'
  time_step_interval = 1
[]
