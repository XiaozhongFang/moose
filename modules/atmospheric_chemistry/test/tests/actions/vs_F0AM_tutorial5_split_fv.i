# vs F0AM — tutorial_5sp — operator-split (MultiApp)
# Parent: FV NS flow + passive scalar transport (no chemistry)
# Child:  atmospheric_chemistry box model (mode=box, 0D ODE)
# Transfer: parent concentrations → child → child reacts → back to parent
#
# Requires: combined-opt (parent), atmospheric_chemistry-opt (child)
# If non-convergent: this is a conceptual demo of the split approach

mu = 0.01
rho = 1

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

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    density = ${rho}
    dynamic_viscosity = ${mu}
    initial_pressure = 0.0
    inlet_boundaries = 'top'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_functors = '1 0'
    wall_boundaries = 'left right bottom'
    momentum_wall_types = 'noslip noslip noslip'
    pin_pressure = true
    pinned_pressure_type = average
    pinned_pressure_value = 0
    mass_advection_interpolation = 'average'
    momentum_advection_interpolation = 'average'
  []
[]

[Variables]
  [A]  family = LAGRANGE  order = FIRST  []
  [B]  family = LAGRANGE  order = FIRST  []
  [C]  family = LAGRANGE  order = FIRST  []
[]

[ICs]
  [A_ic]  type = ConstantIC  variable = A  value = 100  []
  [B_ic]  type = ConstantIC  variable = B  value = 0.05  []
  [C_ic]  type = ConstantIC  variable = C  value = 0  []
[]

[Kernels]
  [A_time]  type = TimeDerivative  variable = A  []
  [B_time]  type = TimeDerivative  variable = B  []
  [C_time]  type = TimeDerivative  variable = C  []
  [A_advect]  type = ConservativeAdvection  variable = A  velocity = 'vel_x vel_y'  []
  [B_advect]  type = ConservativeAdvection  variable = B  velocity = 'vel_x vel_y'  []
  [C_advect]  type = ConservativeAdvection  variable = C  velocity = 'vel_x vel_y'  []
  [A_diff]  type = Diffusion  variable = A  []
  [B_diff]  type = Diffusion  variable = B  []
  [C_diff]  type = Diffusion  variable = C  []
[]

[BCs]
  [A_top]  type = DirichletBC  variable = A  boundary = 'top'  value = 100  []
  [B_top]  type = DirichletBC  variable = B  boundary = 'top'  value = 0.05  []
  [C_top]  type = DirichletBC  variable = C  boundary = 'top'  value = 0  []
[]

[MultiApps]
  [chem]
    type = TransientMultiApp
    input_files = 'vs_F0AM_tutorial5_split_sub.i'
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [to_chem_A]
    type = MultiAppCopyTransfer
    to_multi_app = chem
    source_variable = A
    variable = A
  []
  [to_chem_B]
    type = MultiAppCopyTransfer
    to_multi_app = chem
    source_variable = B
    variable = B
  []
  [to_chem_C]
    type = MultiAppCopyTransfer
    to_multi_app = chem
    source_variable = C
    variable = C
  []
  [from_chem_A]
    type = MultiAppCopyTransfer
    from_multi_app = chem
    source_variable = A
    variable = A
  []
  [from_chem_B]
    type = MultiAppCopyTransfer
    from_multi_app = chem
    source_variable = B
    variable = B
  []
  [from_chem_C]
    type = MultiAppCopyTransfer
    from_multi_app = chem
    source_variable = C
    variable = C
  []
[]

[Postprocessors]
  [A_avg]  type = ElementAverageValue  variable = A  []
  [B_avg]  type = ElementAverageValue  variable = B  []
  [C_avg]  type = ElementAverageValue  variable = C  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  dt = 0.001
  num_steps = 5
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
[]

[Outputs]
  csv = true
  file_base = 'vs_F0AM_tutorial5_split_fv'
[]
