# vs F0AM — tutorial_5sp — FV NS + energy + chemistry
# Fully-coupled: NavierStokesFV + heat transfer + AtmosphericChemistry
# Temperature affects chemistry rate through MCMRatesMaterial
#
# Requires: combined-opt

mu = 0.01
rho = 1
k = 0.01
cp = 1
velocity_interp_method = 'rc'
advected_interp_method = 'average'

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 0.1
    ymin = 0
    ymax = 0.1
    nx = 10
    ny = 10
  []
[]

[Variables]
  [vel_x]    type = INSFVVelocityVariable  []
  [vel_y]    type = INSFVVelocityVariable  []
  [pressure] type = INSFVPressureVariable  []
  [T_fluid]  type = INSFVEnergyVariable  []
  [lambda]   family = SCALAR  order = FIRST  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []
  [mean_zero_pressure]
    type = FVIntegralValueConstraint
    variable = pressure
    lambda = lambda
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
  [temp_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T_fluid
  []
  [temp_advection]
    type = INSFVEnergyAdvection
    variable = T_fluid
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
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
  [T_hot]
    type = FVDirichletBC
    variable = T_fluid
    boundary = 'bottom'
    value = 300
  []
  [T_cold]
    type = FVDirichletBC
    variable = T_fluid
    boundary = 'top'
    value = 298
  []
[]

[FunctorMaterials]
  [functor_constants]
    type = ADGenericFunctorMaterial
    prop_names = 'cp k'
    prop_values = '${cp} ${k}'
  []
  [ins_fv]
    type = INSFVEnthalpyFunctorMaterial
    temperature = 'T_fluid'
    rho = ${rho}
  []
[]

[Functions]
  [lid_function]
    type = ParsedFunction
    expression = '4*x*(1-x)'
  []
[]

[AtmosphericChemistry]
  [Coupled]
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/tutorial_5sp.fac'
    temperature = 298
  []
[]

[ICs]
  [A_ic]   type = ConstantIC  variable = A   value = 100  []
  [B_ic]   type = ConstantIC  variable = B   value = 0.05  []
  [C_ic]   type = ConstantIC  variable = C   value = 0  []
  [ONE_ic] type = ConstantIC  variable = ONE value = 1 []
  [RO2_ic] type = ConstantIC  variable = RO2 value = 0 []
[]

[Kernels]
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

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  checkpoint = false
  csv = true
  file_base = 'vs_F0AM_tutorial5_ns_fv_heat'
[]
