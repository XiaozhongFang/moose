[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 10
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
  []
  [subdomain1]
    type = SubdomainBoundingBoxGenerator
    input = gen
    block_id = 1
    bottom_left = '0.5 0 0'
    top_right = '1 1 0'
  []
  [interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = subdomain1
    primary_block = 0
    paired_block = 1
    new_boundary = 'interface'
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [diff_left]
    type = Diffusion
    variable = u
    block = 0
  []
  [diff_right]
    type = MatDiffusion
    variable = u
    block = 1
    diffusivity = 10.0
  []
[]

[Materials]
  [left]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'diffusivity'
    prop_values = '1.0'
  []
  [right]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'diffusivity'
    prop_values = '10.0'
  []
[]

[InterfaceKernels]
  [ghost_penalty]
    type = GhostPenaltyKernel
    variable = u
    neighbor_var = u
    boundary = 'interface'
    gamma = 1.0
    k = 1
    c_F = 1.0
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = 'left'
    value = 0.0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = 'right'
    value = 1.0
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  console = false
[]
