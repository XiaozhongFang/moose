[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
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
  # Interface at x=0.5 which passes exactly through element midlines
  # Elements on either side: phi is constant sign → NO cut elements
  # Elements at x=0.5 line: phi=0 at all nodes (boundary-tangent)
  [level_set_func]
    type = ParsedFunction
    expression = 'x - 0.5'
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
