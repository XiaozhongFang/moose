# Smoke test for AtmosphericChemistryCoupled transport-kernel generation.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = 2
    ny = 2
  []
[]

[AuxVariables]
  [vel_x]
    family = MONOMIAL
    order = CONSTANT
  []
  [vel_y]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [vel_x_aux]
    type = ConstantAux
    variable = vel_x
    value = 0.02
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
  [vel_y_aux]
    type = ConstantAux
    variable = vel_y
    value = -0.01
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[AtmosphericChemistry]
  [Coupled]
    mechanism_file = '../chamber/aerosol_smoke_mech.fac'
    temperature = 298.15
    air_density = 2.46e19
    water_vapor = 2.46e17
    mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
    transport_velocity = 'vel_x vel_y'
    advection_upwinding_type = full
    diffusivity = 1e-4
  []
[]

[ICs]
  [ONE_ic]
    type = ConstantIC
    variable = ONE
    value = 1
  []
  [RO2_ic]
    type = ConstantIC
    variable = RO2
    value = 1
  []
  [A_ic]
    type = FunctionIC
    variable = A
    function = '1 + 0.1*x'
  []
  [B_ic]
    type = ConstantIC
    variable = B
    value = 2
  []
  [C_ic]
    type = ConstantIC
    variable = C
    value = 0
  []
[]

[Postprocessors]
  [A_avg]
    type = ElementAverageValue
    variable = A
  []
  [C_avg]
    type = ElementAverageValue
    variable = C
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    solve_type = NEWTON
  []
[]

[Executioner]
  type = Transient
  scheme = implicit-euler
  num_steps = 1
  dt = 0.1
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-10
  l_max_its = 50
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
[]

[Outputs]
  csv = true
  exodus = false
[]
