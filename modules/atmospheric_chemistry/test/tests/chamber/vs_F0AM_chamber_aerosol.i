# F0AM Chamber — MCM Inorg + Isoprene with Dynamic Aerosol Partitioning
#   - Gas-particle partitioning for isoprene hydroperoxides
#   - Wall loss for vapor and particle species
#   - Seed aerosol (50 nm, 10^4 cm^-3)
#   - Outputs gas + particle concentrations for paper figure comparison

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
    jfac = 1.0
    output_ro2_sum = true
    use_limiting_reagent = true
    chem_solver = petsc_ts
    chem_solver_type = bdf
    chem_solver_rtol = 5e-4
    chem_solver_atol = 5e-8

    # Aerosol partitioning for isoprene hydroperoxides (most likely to condense)
    aerosol_gas_species = 'ISOPAOOH ISOPBOOH ISOPCOOH ISOPDOOH IEPOXA'
    aerosol_particle_species = 'ISOPAOOH_p ISOPBOOH_p ISOPCOOH_p ISOPDOOH_p IEPOXA_p'
    aerosol_cstar = '1000 1000 1000 1000 500'
    aerosol_molecular_weights = '118 118 118 118 118'
    aerosol_cstar_cutoff = 1000.0
    aerosol_alpha = 0.1
    aerosol_gas_diffusivity = 1.0e-5
    aerosol_particle_number = 1.0e10
    aerosol_seed_radius = 25.0e-9
    aerosol_organic_density = 1400.0
    aerosol_background_organic_mass = 0.0
    aerosol_vapor_wall_loss = 1.0e-5
    aerosol_particle_wall_loss = 6.0e-5
  []
[]

[ICs]
  [C5H8_ic]    type = ScalarConstantIC  variable = C5H8  value = 4.92440663682e11   []  # 20 ppb
  [NO2_ic]     type = ScalarConstantIC  variable = NO2   value = 2.46220331841e9    []  # 0.1 ppb
  [H2O2_ic]    type = ScalarConstantIC  variable = H2O2  value = 1.23110165921e14   []  # 5 ppm
  [ISOPAOOH_p_ic]  type = ScalarConstantIC  variable = ISOPAOOH_p  value = 0.0  []
  [ISOPBOOH_p_ic]  type = ScalarConstantIC  variable = ISOPBOOH_p  value = 0.0  []
  [ISOPCOOH_p_ic]  type = ScalarConstantIC  variable = ISOPCOOH_p  value = 0.0  []
  [ISOPDOOH_p_ic]  type = ScalarConstantIC  variable = ISOPDOOH_p  value = 0.0  []
  [IEPOXA_p_ic]    type = ScalarConstantIC  variable = IEPOXA_p    value = 0.0  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 3600
  l_max_its = 200
  l_tol = 1e-4
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  dt = 300.0
[]

[Outputs]
  checkpoint = false
  console = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_F0AM_chamber_aerosol'
  time_step_interval = 1
[]
