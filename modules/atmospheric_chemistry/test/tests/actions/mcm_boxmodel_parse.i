[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[Variables]
  [dummy]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = dummy
  []
[]

[UserObjects]
  [box_model]
    type = MCMBoxModel
    mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_export.fac'
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  console = true
[]