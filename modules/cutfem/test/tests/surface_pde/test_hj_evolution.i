[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 40
    xmin = -1.0
    xmax = 1.0
    ymin = -1.0
    ymax = 1.0
  []
[]

[Variables]
  [phi]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [phi_advection]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  # Time-dependent LevelSet advection
  # ∂φ/∂t + V·∇φ = 0  (simplified H-J for constant V)
  [phi_time]
    type = TimeDerivative
    variable = phi
  []
  [phi_advection_kernel]
    type = ConservativeAdvection
    variable = phi
    velocity = '0.0 -0.1 0.0'
  []
[]

[BCs]
  [all]
    type = DirichletBC
    variable = phi
    boundary = 'left right top bottom'
    value = 1.0
  []
[]

[ICs]
  [phi_ic]
    type = FunctionIC
    variable = phi
    function = 'sqrt(x^2 + y^2) - 0.5'
  []
[]

[Executioner]
  type = Transient
  dt = 0.05
  num_steps = 10
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-10
[]

[Outputs]
  exodus = true
[]
