# Chamber operator-splitting smoke test on a small urban street-canyon mesh.
#
# Parent app: 2D incompressible FV flow plus transported chamber species.
# Child app: full F0AM chamber KPP Rosenbrock chemistry on the cloned mesh.
#
# This mirrors CTM operator splitting: transport acts on distributed grid fields,
# then chemistry is advanced locally on the same grid. It deliberately uses one
# mesh-wide MultiApp, not one sub-application per grid point.

rho = 1.2
mu = 1.0
D_species = 0.05
velocity_interp_method = 'rc'
advected_interp_method = 'average'
species = 'C5H8 NO2 H2O2 OH HO2 CH3O2 O3 MVK MACR HCHO'

[GlobalParams]
  two_term_boundary_expansion = false
  rhie_chow_user_object = 'rc'
[]

[Mesh]
  parallel_type = distributed
  [street_canyon]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 20
    ymin = 0
    ymax = 8
    nx = 2
    ny = 2
    subdomain_name = 'street_air'
  []
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
    block = 'street_air'
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    block = 'street_air'
  []
  [vel_y]
    type = INSFVVelocityVariable
    block = 'street_air'
  []
  [pressure]
    type = INSFVPressureVariable
    block = 'street_air'
  []
  [lambda]
    family = SCALAR
    order = FIRST
  []

  [C5H8]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-11
  []
  [NO2]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-10
  []
  [H2O2]
    family = LAGRANGE
    order = FIRST
    scaling = 1e-12
  []
  [OH]
    family = LAGRANGE
    order = FIRST
  []
  [HO2]
    family = LAGRANGE
    order = FIRST
  []
  [CH3O2]
    family = LAGRANGE
    order = FIRST
  []
  [O3]
    family = LAGRANGE
    order = FIRST
  []
  [MVK]
    family = LAGRANGE
    order = FIRST
  []
  [MACR]
    family = LAGRANGE
    order = FIRST
  []
  [HCHO]
    family = LAGRANGE
    order = FIRST
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
  [OH_ic]
    type = ConstantIC
    variable = OH
    value = 0
  []
  [HO2_ic]
    type = ConstantIC
    variable = HO2
    value = 0
  []
  [CH3O2_ic]
    type = ConstantIC
    variable = CH3O2
    value = 0
  []
  [O3_ic]
    type = ConstantIC
    variable = O3
    value = 0
  []
  [MVK_ic]
    type = ConstantIC
    variable = MVK
    value = 0
  []
  [MACR_ic]
    type = ConstantIC
    variable = MACR
    value = 0
  []
  [HCHO_ic]
    type = ConstantIC
    variable = HCHO
    value = 0
  []
[]

[Functions]
  [lid_function]
    type = ParsedFunction
    expression = '0.5*4.0*(x/20.0)*(1.0-x/20.0)'
  []
  [c5h8_initial]
    type = ParsedFunction
    expression = '2.46220331841e11*(1.0 + 0.10*exp(-((x-10.0)^2+(y-1.0)^2)/18.0))'
  []
  [no2_initial]
    type = ParsedFunction
    expression = '2.46220331841e10*(1.0 + 0.75*exp(-((x-10.0)^2+(y-0.6)^2)/8.0))'
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    block = 'street_air'
  []
  [mean_zero_pressure]
    type = FVIntegralValueConstraint
    variable = pressure
    lambda = lambda
    block = 'street_air'
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'x'
    block = 'street_air'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    block = 'street_air'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
    block = 'street_air'
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'y'
    block = 'street_air'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
    block = 'street_air'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
    block = 'street_air'
  []
[]

[FVBCs]
  [top_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'top'
    function = 'lid_function'
  []
  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'left right bottom'
    function = 0
  []
  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'left right top bottom'
    function = 0
  []
[]

[Materials]
  [transport_velocity]
    type = VectorFromComponentVariablesMaterial
    vector_prop_name = 'transport_velocity'
    u = vel_x
    v = vel_y
    block = 'street_air'
  []
[]

[Kernels]
  [C5H8_time]
    type = TimeDerivative
    variable = C5H8
  []
  [C5H8_advect]
    type = ConservativeAdvection
    variable = C5H8
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [C5H8_diff]
    type = CoefDiffusion
    variable = C5H8
    coef = ${D_species}
  []
  [NO2_time]
    type = TimeDerivative
    variable = NO2
  []
  [NO2_advect]
    type = ConservativeAdvection
    variable = NO2
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [NO2_diff]
    type = CoefDiffusion
    variable = NO2
    coef = ${D_species}
  []
  [H2O2_time]
    type = TimeDerivative
    variable = H2O2
  []
  [H2O2_advect]
    type = ConservativeAdvection
    variable = H2O2
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [H2O2_diff]
    type = CoefDiffusion
    variable = H2O2
    coef = ${D_species}
  []
  [OH_time]
    type = TimeDerivative
    variable = OH
  []
  [OH_advect]
    type = ConservativeAdvection
    variable = OH
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [OH_diff]
    type = CoefDiffusion
    variable = OH
    coef = ${D_species}
  []
  [HO2_time]
    type = TimeDerivative
    variable = HO2
  []
  [HO2_advect]
    type = ConservativeAdvection
    variable = HO2
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [HO2_diff]
    type = CoefDiffusion
    variable = HO2
    coef = ${D_species}
  []
  [CH3O2_time]
    type = TimeDerivative
    variable = CH3O2
  []
  [CH3O2_advect]
    type = ConservativeAdvection
    variable = CH3O2
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [CH3O2_diff]
    type = CoefDiffusion
    variable = CH3O2
    coef = ${D_species}
  []
  [O3_time]
    type = TimeDerivative
    variable = O3
  []
  [O3_advect]
    type = ConservativeAdvection
    variable = O3
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [O3_diff]
    type = CoefDiffusion
    variable = O3
    coef = ${D_species}
  []
  [MVK_time]
    type = TimeDerivative
    variable = MVK
  []
  [MVK_advect]
    type = ConservativeAdvection
    variable = MVK
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [MVK_diff]
    type = CoefDiffusion
    variable = MVK
    coef = ${D_species}
  []
  [MACR_time]
    type = TimeDerivative
    variable = MACR
  []
  [MACR_advect]
    type = ConservativeAdvection
    variable = MACR
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [MACR_diff]
    type = CoefDiffusion
    variable = MACR
    coef = ${D_species}
  []
  [HCHO_time]
    type = TimeDerivative
    variable = HCHO
  []
  [HCHO_advect]
    type = ConservativeAdvection
    variable = HCHO
    velocity_material = 'transport_velocity'
    upwinding_type = full
  []
  [HCHO_diff]
    type = CoefDiffusion
    variable = HCHO
    coef = ${D_species}
  []
[]

[MultiApps]
  [chem]
    type = TransientMultiApp
    input_files = 'vs_F0AM_tutorial5_split_sub.i'
    clone_parent_mesh = true
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [species_to_chem]
    type = MultiAppCopyTransfer
    to_multi_app = chem
    source_variable = '${species}'
    variable = '${species}'
    execute_on = 'timestep_end'
  []
  [species_from_chem]
    type = MultiAppCopyTransfer
    from_multi_app = chem
    source_variable = '${species}'
    variable = '${species}'
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  # These are diagnostics only; CSV output is not used to drive the coupling.
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
  num_steps = 2
  l_max_its = 300
  l_tol = 1e-8
  nl_max_its = 20
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-8
  # Keep the split transport solve distributed; the chemistry sub-app also uses
  # rank-local block preconditioning for per-grid-cell reactions.
  petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'gmres bjacobi lu NONZERO'
[]

[Outputs]
  checkpoint = false
  csv = true
  exodus = true
  file_base = 'vs_F0AM_tutorial5_split_fv'
[]
