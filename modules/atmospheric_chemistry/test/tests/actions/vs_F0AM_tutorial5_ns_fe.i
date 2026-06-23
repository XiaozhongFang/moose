# vs F0AM — tutorial_5sp — FE NS lid-driven cavity + chemistry
# Fully-coupled: INS FE flow (manual kernels) + AtmosphericChemistry Action
# Chemistry species are advected by NS velocity and diffuse
#
# Requires: combined-opt
# If non-convergent: reduce nx/ny

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 1.0
    nx = 8
    ny = 8
    elem_type = QUAD9
  []
  [corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  []
[]

[Variables]
  [vel_x]  order = SECOND  family = LAGRANGE  []
  [vel_y]  order = SECOND  family = LAGRANGE  []
  [p]      order = FIRST   family = LAGRANGE  []
[]

[Kernels]
  # NS mass
  [mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    pressure = p
  []
  # x-momentum time
  [x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  []
  # x-momentum space
  [x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = p
    component = 0
  []
  # y-momentum time
  [y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  []
  # y-momentum space
  [y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    pressure = p
    component = 1
  []
[]

[BCs]
  [x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'bottom right left'
    value = 0.0
  []
  [lid]
    type = FunctionDirichletBC
    variable = vel_x
    boundary = 'top'
    function = 'lid_function'
  []
  [y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'bottom right top left'
    value = 0.0
  []
  [pressure_pin]
    type = DirichletBC
    variable = p
    boundary = 'pinned_node'
    value = 0
  []
[]

[Functions]
  [lid_function]
    type = ParsedFunction
    expression = '4*x*(1-x)'
  []
[]

[Materials]
  [const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '1  1'
  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
  temperature = 298
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 100  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.05  []
  [C_ic]   type = ConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ConstantIC  variable = RO2 value = 0 []
[]

[Kernels]
  # Species advection + diffusion
  [A_advect]  type = ConservativeAdvection  variable = A  velocity = 'vel_x vel_y'  []
  [B_advect]  type = ConservativeAdvection  variable = B  velocity = 'vel_x vel_y'  []
  [C_advect]  type = ConservativeAdvection  variable = C  velocity = 'vel_x vel_y'  []
  [ONE_advect]  type = ConservativeAdvection  variable = ONE  velocity = 'vel_x vel_y'  []
  [RO2_advect]  type = ConservativeAdvection  variable = RO2  velocity = 'vel_x vel_y'  []

  [A_diff]  type = Diffusion  variable = A  []
  [B_diff]  type = Diffusion  variable = B  []
  [C_diff]  type = Diffusion  variable = C  []
  [ONE_diff]  type = Diffusion  variable = ONE  []
  [RO2_diff]  type = Diffusion  variable = RO2  []
[]

[BCs]
  [A_top]  type = DirichletBC  variable = A  boundary = 'top'  value = 100  []
  [B_top]  type = DirichletBC  variable = B  boundary = 'top'  value = 0.05  []
  [C_top]  type = DirichletBC  variable = C  boundary = 'top'  value = 0  []
  [ONE_top]  type = DirichletBC  variable = ONE  boundary = 'top'  value = 1  []
  [RO2_top]  type = DirichletBC  variable = RO2  boundary = 'top'  value = 0  []
[]

[Postprocessors]
  [A_avg]  type = ElementAverageValue  variable = A  []
  [B_avg]  type = ElementAverageValue  variable = B  []
  [C_avg]  type = ElementAverageValue  variable = C  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  []
[]

[Executioner]
  type = Transient
  num_steps = 5
  dt = 0.01
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  file_base = 'vs_F0AM_tutorial5_ns_fe'
[]
