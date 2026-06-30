[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
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
  [diff]
    type = Diffusion
    variable = u
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
  [level_set_func]
    type = ParsedFunction
    # Circle of radius 0.3 centered at (0.5, 0.5)
    # φ < 0 inside the circle
    expression = 'sqrt((x-0.5)^2 + (y-0.5)^2) - 0.3'
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

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  console = true
[]
