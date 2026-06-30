[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
    xmin = -1.0
    xmax = 1.0
    ymin = -1.0
    ymax = 1.0
  []
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [cut_diff]
    type = CutFEMDiffusion
    variable = u
    cut_cell_quadrature = cut_cell_quad
  []
  # Note: Additional kernels on the same variable (e.g. BodyForce) require
  # the custom quadrature to be applied uniformly. In this single-kernel
  # setup, CutFEMDiffusion handles cut elements via _fe_problem.reinitElemPhys.
[]

[Functions]
  # Domain: circle of radius 0.8 centered at origin
  # φ(x) < 0 inside the circle (physical domain)
  [level_set_func]
    type = ParsedFunction
    expression = 'sqrt(x^2 + y^2) - 0.8'
  []
[]

[UserObjects]
  [cut_cell_quad]
    type = CutCellQuadratureUserObject
    level_set_function = level_set_func
    quadrature_order = 4
    level_set_tolerance = 1e-12
    execute_on = 'INITIAL'
  []
[]

[BCs]
  # Nitsche BC on the circle boundary is not yet implemented
  # For demonstration, use Dirichlet on the background mesh boundary
  [all]
    type = DirichletBC
    variable = u
    boundary = 'left right top bottom'
    value = 0.0
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-10
[]

[Outputs]
  exodus = true
  console = true
[]
