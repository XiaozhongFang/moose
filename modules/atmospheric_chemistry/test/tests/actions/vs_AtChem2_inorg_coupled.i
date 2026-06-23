# vs AtChem2 CVODE — MCM inorganic — coupled mode (FEM 0D transport)
# MCM inorganic mechanism subset (20 species, 48 reactions)
# Matching AtChem2 model_mcm_inorg conditions: 12h diurnal, SZA photolysis
#
# Note: Uses mcm_export.fac (simpler inorganic subset) because
# atchem2_example.fac has RO2-lumping expressions incompatible
# with fparser in coupled mode. Box mode (vs_AtChem2_inorg_box.i)
# uses the full atchem2_example.fac mechanism.
#
# Gold file: vs_AtChem2_inorg_coupled.csv

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  mode = coupled
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_export.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [O3_ic]
    type = ConstantIC
    variable = O3
    value = 5.2e11
  []
  [NO_ic]
    type = ConstantIC
    variable = NO
    value = 2.5e11
  []
  [NO2_ic]
    type = ConstantIC
    variable = NO2
    value = 2.4e11
  []
  [OH_ic]
    type = ConstantIC
    variable = OH
    value = 1.0e6
  []
  [HO2_ic]
    type = ConstantIC
    variable = HO2
    value = 1.0e8
  []
  [CO_ic]
    type = ConstantIC
    variable = CO
    value = 3.6e12
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 900
  end_time = 43200
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-3
  [TimeStepper]
    type = ConstantDT
    dt = 900
  []
[]

[Postprocessors]
  [O3_val]
    type = ElementAverageValue
    variable = O3
  []
  [NO_val]
    type = ElementAverageValue
    variable = NO
  []
  [NO2_val]
    type = ElementAverageValue
    variable = NO2
  []
  [OH_val]
    type = ElementAverageValue
    variable = OH
  []
  [HO2_val]
    type = ElementAverageValue
    variable = HO2
  []
  [CO_val]
    type = ElementAverageValue
    variable = CO
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_AtChem2_inorg_coupled'
  time_step_interval = 1
[]
