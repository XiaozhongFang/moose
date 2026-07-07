# vs F0AM — tutorial_5sp — FV NS lid-driven cavity + chemistry
# Chemistry species solved on same mesh as NS (same nonlinear solve)
# No advection coupling (ConservativeAdvection velocity API changed)
# — demonstrates fully-coupled architecture: NS + chemistry in one system
# Requires: combined-opt

mu = 0.01
rho = 1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0  xmax = 0.1  ymin = 0  ymax = 0.1
    nx = 10  ny = 10
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

# ——— Chemistry (manual, combined-opt compatible) ———
# Species order: ONE(0) RO2(1) A(2) B(3) C(4)
# Rx1:A+B→C+B  Rx2:B→  Rx3:A→  Rx4:B+B→  Rx5:→A  Rx6:C→

[Variables]
  [ONE]  family = LAGRANGE  order = FIRST  []
  [RO2]  family = LAGRANGE  order = FIRST  []
  [A]    family = LAGRANGE  order = FIRST  []
  [B]    family = LAGRANGE  order = FIRST  []
  [C]    family = LAGRANGE  order = FIRST  []
[]

[ICs]
  [ONE_ic]  type = ConstantIC  variable = ONE  value = 1  []
  [RO2_ic]  type = ConstantIC  variable = RO2  value = 0  []
  [A_ic]    type = ConstantIC  variable = A    value = 100  []
  [B_ic]    type = ConstantIC  variable = B    value = 0.05  []
  [C_ic]    type = ConstantIC  variable = C    value = 0  []
[]

[Kernels]
  [ONE_td]  type = TimeDerivative  variable = ONE  []
  [RO2_td]  type = TimeDerivative  variable = RO2  []
  [A_td]    type = TimeDerivative  variable = A  []
  [B_td]    type = TimeDerivative  variable = B  []
  [C_td]    type = TimeDerivative  variable = C  []

  [ONE_src]  type = ChemicalSourceKernel  variable = ONE
    stoichiometric_row = '0 0 0 0 1 0'
    all_species = 'ONE RO2 A B C'
    species_reactants = '4 1'
  []
  [A_src]  type = ChemicalSourceKernel  variable = A
    stoichiometric_row = '-1 0 -1 0 1 0'
    all_species = 'ONE RO2 A B C'
    species_reactants = '0 1 2 1'
  []
  [B_src]  type = ChemicalSourceKernel  variable = B
    stoichiometric_row = '0 -1 0 -2 0 0'
    all_species = 'ONE RO2 A B C'
    species_reactants = '0 1 1 1 3 1'
  []
  [C_src]  type = ChemicalSourceKernel  variable = C
    stoichiometric_row = '1 0 0 0 0 -1'
    all_species = 'ONE RO2 A B C'
    species_reactants = '5 1'
  []
[]

[Materials]
  [mcm_rates]
    type = MCMRatesMaterial
    temperature = 298
    air_density = 2.46e19
    species_list = 'ONE RO2 A B C'
    species_variables = 'ONE RO2 A B C'
    reaction_rate_expressions = 'K1 K2 K3 K4 K5 K6'
    coefficient_names = 'K1 K2 K3 K4 K5 K6'
    coefficient_expressions = '1.0e-3 1.0e-2 1.0e-4 1.0e-1 5.0e-1 1.0e-4'
    reactant_matrix = '2 1 3 1; 3 1; 2 1; 3 1 3 1; 0 1; 4 1'
    j_cl_values = ''
    j_cmm_values = ''
    j_cnn_values = ''
  []
[]

[BCs]
  [ONE_bnd]  type = DirichletBC  variable = ONE  boundary = 'top'  value = 1  []
  [RO2_bnd]  type = DirichletBC  variable = RO2  boundary = 'top'  value = 0  []
  [A_bnd]    type = DirichletBC  variable = A    boundary = 'top'  value = 100  []
  [B_bnd]    type = DirichletBC  variable = B    boundary = 'top'  value = 0.05  []
  [C_bnd]    type = DirichletBC  variable = C    boundary = 'top'  value = 0  []
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
  checkpoint = false
  csv = true
  file_base = 'vs_F0AM_tutorial5_ns_fv'
[]
