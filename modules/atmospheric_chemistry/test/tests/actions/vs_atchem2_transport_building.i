# vs AtChem2 — transport + chemistry with building emission — FE NS fully-coupled
# Fully-coupled: INS FE flow + AtmosphericChemistry Action (mode=coupled)
# SUPG+PSPG stabilization enables equal-order LAGRANGE elements
# 
# Domain: 20×10 m — lid-driven wind field (max 10 m/s at x=10)
# Building: no-slip bottom center region → recirculation zone
# Emission: Gaussian CO plume above building (x=8, y=2)
# Chemistry: atchem2_example.fac (29 species, 71 rx, MCM_SZA photolysis)
# Ref: vs_F0AM_tutorial5_ns_fe.i (flow template) + vs_AtChem2_inorg_coupled.i (chemistry)

[GlobalParams]
  gravity = '0 0 0'
  laplace = true
  integrate_p_by_parts = true
  family = LAGRANGE
  order = FIRST
  supg = true
  pspg = true
  alpha = 1e-1
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 20.0
    ymin = 0
    ymax = 10.0
    nx = 64
    ny = 32
  []
  [corner_node]
    type = ExtraNodesetGenerator
    new_boundary = 'pinned_node'
    nodes = '0'
    input = gen
  []
[]

[Variables]
  [vel_x]  []
  [vel_y]  []
  [p]      []
[]

[Kernels]
  [mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    pressure = p
  []
  [x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  []
  [x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = p
    component = 0
  []
  [y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  []
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
  # No-slip on bottom, left, right
  [x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'bottom right left'
    value = 0.0
  []
  # Lid-driven top: wind from left to right, max 10 m/s at center
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
    # Scaled lid function: 4*(x/20)*(1-x/20)*10 → max 10 m/s at x=10
    expression = '40*x*(20-x)/100'
  []
[]

[Materials]
  [const]
    type = GenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '1.2 0.15'
  []
[]

[AtmosphericChemistry]
  mode = coupled
  include_transport = true
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/atchem2_example.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  press = 1013.25
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  photolysis_scheme = MCM_SZA
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  # Species ICs from AtChem2 validation
  [CH4_ic]  type = ConstantIC  variable = CH4  value = 4.9e13 []
  [O3_ic]   type = ConstantIC  variable = O3   value = 5.2e11 []
  [NO2_ic]  type = ConstantIC  variable = NO2  value = 2.4e11 []
  # CO: background 3.6e12 + Gaussian plume near building center (x=8, y=2)
  [co_plume]
    type = FunctionIC
    variable = CO
    function = '3.6e12 + 2e12 * exp(-((x-8)^2+(y-2)^2)/4)'
  []
  # Flow ICs
  [vel_x_ic]  type = ConstantIC  variable = vel_x  value = 0  []
  [vel_y_ic]  type = ConstantIC  variable = vel_y  value = 0  []
[]

[Postprocessors]
  [O3_avg]   type = ElementAverageValue  variable = O3   []
  [NO2_avg]  type = ElementAverageValue  variable = NO2  []
  [CO_avg]   type = ElementAverageValue  variable = CO   []
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
  num_steps = 100
  dt = 0.5
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type -sub_pc_factor_levels'
  petsc_options_value = 'asm      2               ilu          4'
  line_search = 'none'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  nl_max_its = 8
  l_tol = 1e-6
  l_max_its = 500
[]

[Outputs]
  exodus = true
  csv = true
  file_base = 'vs_atchem2_transport_building'
[]
