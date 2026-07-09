# Chemistry sub-application for the chamber split-coupling smoke test.
#
# The parent application clones its street-canyon mesh into this sub-app. The
# KPP source kernels then evaluate the full F0AM chamber mechanism at each local
# grid degree of freedom; this is not a spatially averaged box model.

[Mesh]
  # Populated by the parent TransientMultiApp with clone_parent_mesh = true.
[]

[AtmosphericChemistry]
  [Coupled]
    mechanism_file = 'chamber/kpp_chamber/generated_mechanisms/chamber_mcm_rosenbrock/chamber_mcm_rosenbrock.kpp'
    temperature = 298.0
    air_density = 2.4622e19
    water_vapor = 7.6114e16
    press = 1013.0
    jfac = 1.0
    roof_open = true
    photolysis_scheme = BOTTOMUP
    lamp_flux_file = 'ExampleLightFlux.txt'
    bottomup_data_dir = '../../doc/content/modules/atmospheric_chemistry/database/photolysis/bottomup'
    chem_solver = kpp_rosenbrock
    units = molec_cm3
  []
[]

[ICs]
  [C5H8_ic]
    type = FunctionIC
    variable = C5H8
    function = c5h8_initial
  []
  [NO2_ic]
    type = FunctionIC
    variable = NO2
    function = no2_initial
  []
  [H2O2_ic]
    type = ConstantIC
    variable = H2O2
    value = 4.92440663682e12
  []
[]

[Functions]
  [c5h8_initial]
    type = ParsedFunction
    expression = '2.46220331841e11*(1.0 + 0.10*exp(-((x-10.0)^2+(y-1.0)^2)/18.0))'
  []
  [no2_initial]
    type = ParsedFunction
    expression = '2.46220331841e10*(1.0 + 0.75*exp(-((x-10.0)^2+(y-0.6)^2)/8.0))'
  []
[]

[Postprocessors]
  # Diagnostic averages only; chemistry itself is evaluated on the full mesh.
  [C5H8_avg]
    type = ElementAverageValue
    variable = C5H8
  []
  [NO2_avg]
    type = ElementAverageValue
    variable = NO2
  []
  [OH_avg]
    type = ElementAverageValue
    variable = OH
  []
  [HO2_avg]
    type = ElementAverageValue
    variable = HO2
  []
  [CH3O2_avg]
    type = ElementAverageValue
    variable = CH3O2
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = bdf2
  dt = 1.0
  num_steps = 3
  l_max_its = 200
  l_tol = 1e-8
  nl_max_its = 15
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
  # Chemistry is local to each distributed mesh partition, so use rank-local
  # block preconditioning instead of a global parallel direct solve.
  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'gmres bjacobi lu NONZERO'
[]

[Outputs]
  checkpoint = false
  csv = true
  console = false
  file_base = 'vs_F0AM_tutorial5_split_sub'
[]
