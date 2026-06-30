[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
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
    type = MatDiffusion
    variable = u
    block = 0
    diffusivity = diff_coef
  []
  [diff_right]
    type = MatDiffusion
    variable = u
    block = 1
    diffusivity = diff_coef
  []
[]

[Materials]
  [left]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'diff_coef'
    prop_values = '1.0'
  []
  [right]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'diff_coef'
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

[Functions]
  # Analytic solution for 1D bimaterial diffusion:
  # k1=1 on [0,0.5), k2=10 on (0.5,1]
  # u(0)=0, u(1)=1
  # u1(x) = (k2/(k1+k2)) * 2x   for x in [0,0.5]
  # u2(x) = (k2/(k1+k2)) + (k1/(k1+k2)) * (2x-1) for x in [0.5,1]
  # With k1=1, k2=10: u1 = 10/11 * x at x=0.5, u2 = 10/11*(0.5) + 1/11*(2x-1)
  [analytic]
    type = ParsedFunction
    # k1=1, k2=10. u1=20/11*x on [0,0.5]; u2=9/11+2/11*x on [0.5,1]
    expression = 'if(x<0.5, 20.0/11.0*x, 9.0/11.0 + 2.0/11.0*x)'
  []
[]

[Postprocessors]
  [L2_error]
    type = ElementL2Error
    variable = u
    function = analytic
  []
  [H1_error]
    type = ElementH1SemiError
    variable = u
    function = analytic
  []
  [h]
    type = AverageElementSize
  []
  [ndofs]
    type = NumDOFs
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_rtol'
  petsc_options_value = 'hypre boomeramg 1e-12'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
  console = false
[]
