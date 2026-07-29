# Aerosol Module Smoke Test — box mode
# Quick validation of aerosol partitioning + wall loss C++ code

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
    mechanism_file = 'aerosol_smoke_mech.fac'
    temperature = 298.0
    air_density = 2.46e19
    water_vapor = 0.0
    press = 1013.0
    photolysis_scheme = MCM_SZA
    mcm_photolysis_file = 'doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
    jfac = 0.0
    use_limiting_reagent = true
    chem_solver = petsc_ts
    chem_solver_type = bdf
    chem_solver_rtol = 1e-6
    chem_solver_atol = 1e-12

    # Aerosol partitioning: A, B, C are gas species with particle counterparts
    aerosol_gas_species = 'A B C'
    aerosol_particle_species = 'A_p B_p C_p'
    aerosol_cstar = '0.001 10.0 0.1'
    aerosol_molecular_weights = '100 100 100'
    aerosol_cstar_cutoff = 100.0
    aerosol_alpha = 0.1
    aerosol_gas_diffusivity = 1.0e-5
    aerosol_particle_number = 1.0e10
    aerosol_seed_radius = 25.0e-9
    aerosol_organic_density = 1400.0
    aerosol_background_organic_mass = 0.01
    aerosol_vapor_wall_loss = 1.0e-4
    aerosol_particle_wall_loss = 0.0
  []
[]

[ICs]
  [ONE_ic] type = ScalarConstantIC  variable = ONE  value = 1.0e20  []
  [RO2_ic] type = ScalarConstantIC  variable = RO2  value = 1.0    []
  [A_ic]  type = ScalarConstantIC  variable = A  value = 1.0e11  []
  [B_ic]  type = ScalarConstantIC  variable = B  value = 1.0e9   []
  [C_ic]  type = ScalarConstantIC  variable = C  value = 1.0e8   []
  [A_p_ic]  type = ScalarConstantIC  variable = A_p  value = 0.0  []
  [B_p_ic]  type = ScalarConstantIC  variable = B_p  value = 0.0  []
  [C_p_ic]  type = ScalarConstantIC  variable = C_p  value = 0.0  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 100
  l_max_its = 200
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  dt = 1.0
[]

[Outputs]
  checkpoint = false
  console = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'aerosol_smoke'
  time_step_interval = 1
[]
