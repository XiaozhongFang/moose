# MCM Box Model Verification Test
# Zero-dimensional chemical kinetic ODE system (20 species, 48 reactions)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_export.fac'
  temperature = 298.15
  air_density = 2.46e19
  water_vapor = 2.46e17
  # SZA-based photolysis calculation uses defaults: lat=51.51, lon=0.13, date=2010-06-21
[]

[ICs]
  [O3_ic]
    type = ConstantIC
    variable = O3
    value = 7.4e11
  []
  [NO_ic]
    type = ConstantIC
    variable = NO
    value = 2.5e11
  []
  [NO2_ic]
    type = ConstantIC
    variable = NO2
    value = 5.0e11
  []
  [O_ic]
    type = ConstantIC
    variable = O
    value = 1.0e5
  []
  [O1D_ic]
    type = ConstantIC
    variable = O1D
    value = 1.0e3
  []
  [OH_ic]
    type = ConstantIC
    variable = OH
    value = 1.0e6
  []
  [HO2_ic]
    type = ConstantIC
    variable = HO2
    value = 1.0e8
  []
  [H2O2_ic]
    type = ConstantIC
    variable = H2O2
    value = 1.0e10
  []
  [NO3_ic]
    type = ConstantIC
    variable = NO3
    value = 1.0e8
  []
  [N2O5_ic]
    type = ConstantIC
    variable = N2O5
    value = 1.0e8
  []
  [HNO3_ic]
    type = ConstantIC
    variable = HNO3
    value = 1.0e10
  []
  [HO2NO2_ic]
    type = ConstantIC
    variable = HO2NO2
    value = 1.0e8
  []
  [HONO_ic]
    type = ConstantIC
    variable = HONO
    value = 1.0e8
  []
  [CO_ic]
    type = ConstantIC
    variable = CO
    value = 1.0e13
  []
  [SO2_ic]
    type = ConstantIC
    variable = SO2
    value = 1.0e10
  []
  [HSO3_ic]
    type = ConstantIC
    variable = HSO3
    value = 1.0e5
  []
  [SO3_ic]
    type = ConstantIC
    variable = SO3
    value = 1.0e5
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
  solve_type = PJFNK
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-5
  end_time = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0e-9
    optimal_iterations = 6
    iteration_window = 3
    growth_factor = 1.5
    cutback_factor_at_failure = 0.5
  []
[]

[Outputs]
  console = true
  print_linear_residuals = false
[]
