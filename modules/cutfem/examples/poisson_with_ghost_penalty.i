# CutFEM Phase 1: Ghost Penalty Stabilization
# Poisson equation on unfitted mesh with Nitsche boundary conditions
# 
# Problem:
#   -∇²u = 1  in Ω
#    u = 0   on ∂Ω
#
# Domain: Circle of radius 0.5 centered at (0.5, 0.5)
# Background mesh: Cartesian grid [0,1]²
#
# References:
# - Burman et al. (2015) "CutFEM for elliptic interface problems"
# - Wichrowski (2026) Section 3: Numerical example

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  elem_type = QUAD9           # 9-node quadrilateral for higher accuracy
[]

[AuxVariables]
  # Level set function: φ(x,y) = r - 0.5, where r = sqrt((x-0.5)² + (y-0.5)²)
  [phi]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [level_set]
    type = FunctionAux
    variable = phi
    function = level_set_func
    execute_on = 'INITIAL'
  []
[]

[Functions]
  # Implicit representation of circular domain
  # φ > 0: inside circle (domain)
  # φ < 0: outside circle (fictitious domain)
  # φ = 0: boundary
  [level_set_func]
    type = ParsedFunction
    value = 'sqrt((x-0.5)^2 + (y-0.5)^2) - 0.5'
  []
[]

[Variables]
  [u]
    order = SECOND          # P2 elements for better accuracy
    family = LAGRANGE
  []
[]

[Kernels]
  [diffusion]
    type = Diffusion
    variable = u
  []
  
  [source]
    type = BodyForce
    variable = u
    function = one
    block = 0
  []
[]

[Functions]
  [one]
    type = ConstantFunction
    value = 1.0
  []
[]

[BCs]
  # Nitsche's method for weak enforcement of Dirichlet BC
  # Boundary = entire background domain boundary (mesh boundary)
  [dirichlet_nitsche]
    type = DirichletBC
    variable = u
    boundary = 'left right bottom top'  # All external boundaries
    value = 0.0
  []
[]

[InterfaceKernels]
  # Ghost Penalty stabilization on internal faces
  # Acts on faces between elements that are both cut by the interface
  [ghost_penalty]
    type = GhostPenaltyKernel
    variable = u
    neighbor_var = u
    boundary = 'internal'              # Internal faces (detected automatically)
    
    # Stabilization parameters (from Wichrowski 2026)
    gamma = 1.0                         # Parameter γ ∈ [0,1]
                                        # 1.0 = weakest but sufficient
                                        # 0.0 = strongest (more diffusive)
    k = 1                               # Derivative order (k=1 for 2nd order PDE)
    c_F = 1.0                           # Face penalty constant (O(1))
  []
[]

[UserObjects]
  # Cut cell quadrature initialization (placeholder for phase 2)
  # In phase 1, standard quadrature is used
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
  
  # PETSc options for robust solving
  petsc_options_iname = '-pc_type -pc_sub_type -ksp_type -ksp_rtol -ksp_atol'
  petsc_options_value = 'asm ilu gmres 1e-8 1e-10'
  
  # Nonlinear solver
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-10
  nl_max_its = 25
  
  # Linear solver
  l_tol = 1e-5
  l_max_its = 100
[]

[Outputs]
  exodus = true
  csv = true
  
  # Monitor convergence
  [console]
    type = Console
    verbosity = VERBOSE
  []
[]

# Post-processing
[Postprocessors]
  [L2_error]
    type = ElementL2Error
    variable = u
    function = zero_func
    execute_on = 'TIMESTEP_END'
  []
  
  [max_u]
    type = ElementExtremeValue
    variable = u
    value_type = max
  []
  
  [min_u]
    type = ElementExtremeValue
    variable = u
    value_type = min
  []
[]

[Functions]
  [zero_func]
    type = ConstantFunction
    value = 0.0
  []
[]

###############################################################################
# Verification Information
###############################################################################
# 
# Expected behavior (from Burman et al. 2015):
# 1. Condition number of stiffness matrix: κ ≤ C h^{-2} (uniform in h)
#    Without GP: κ ≤ C h^{-2} (potentially much larger constant)
# 2. H1 error: ||u - u_h||_H1 ~ O(h) for P1 elements
# 3. L2 error: ||u - u_h||_L2 ~ O(h²) for P1 elements
#
# With 20×20 mesh (h ≈ 0.05):
#   - Expected H1 error: ~0.025
#   - Expected L2 error: ~0.0006
#
###############################################################################
