# vs AtChem2 CVODE — MCM inorganic — coupled mode (FEM 0D transport)
# Same mechanism, ICs, parameters as vs_AtChem2_inorg_box.i.
# Both shall produce identical species concentrations; only the
# implementation differs (ScalarVariable ODE vs FEM with 1-element mesh).
#
# Gold file: vs_AtChem2_inorg_box.csv (shared — CSVDiff compares columns
# present in both files, i.e. species + time)

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistryCoupled]
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/atchem2_example.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  press = 1013.25
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
[]

[ICs]
  [CH4_ic]  type = ConstantIC  variable = CH4  value = 4.9e13  []
  [CO_ic]   type = ConstantIC  variable = CO   value = 3.6e12  []
  [O3_ic]   type = ConstantIC  variable = O3   value = 5.2e11  []
  [NO2_ic]  type = ConstantIC  variable = NO2  value = 2.4e11  []
[]

[Postprocessors]
  # ── Photolysis rates (J1-J61, AtChem2 photolysisRates.output order) ──
  # Read from MCMRatesMaterial photolysis_rates property (coupled mode support)
  [J1]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 1   []
  [J2]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 2   []
  [J3]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 3   []
  [J4]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 4   []
  [J5]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 5   []
  [J6]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 6   []
  [J7]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 7   []
  [J8]   type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 8   []
  [J11]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 11  []
  [J12]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 12  []
  [J13]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 13  []
  [J14]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 14  []
  [J15]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 15  []
  [J16]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 16  []
  [J17]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 17  []
  [J18]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 18  []
  [J19]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 19  []
  [J20]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 20  []
  [J21]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 21  []
  [J22]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 22  []
  [J23]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 23  []
  [J24]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 24  []
  [J31]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 31  []
  [J32]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 32  []
  [J33]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 33  []
  [J34]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 34  []
  [J35]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 35  []
  [J41]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 41  []
  [J51]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 51  []
  [J52]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 52  []
  [J53]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 53  []
  [J54]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 54  []
  [J55]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 55  []
  [J56]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 56  []
  [J61]  type = MCMPhotolysisPostprocessor  material = mcm_rates_material  j_number = 61  []

  # ── Environment variables (AtChem2 environmentVariables.output order) ──
  [M_env]     type = Receiver  default = 2.546917e19 []
  [TEMP]      type = Receiver  default = 288.15     []
  [PRESS]     type = Receiver  default = 1013.25    []
  [RH]        type = Receiver  default = -1.0        []
  [H2O]       type = Receiver  default = 4.45e17    []
  [DEC]       type = Receiver  default = 0.40931542032971796  []
  [BLHEIGHT]  type = Receiver  default = -1.0        []
  [DILUTE]    type = Receiver  default = -1.0        []
  [JFAC]      type = Receiver  default = 1.0        []
  [ROOF]      type = Receiver  default = 1.0        []
  [ASA]       type = Receiver  default = -1.0        []
  [RO2_sum]   type = ElementAverageValue  variable = CH3O2  []

  # ── Solar parameters (time-varying via MCMSolarPostprocessor) ──
  [cosx]   type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = cosx    []
  [secx]   type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = secx    []
  [lha]    type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = lha     []
  [sinld]  type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = sinld   []
  [cosld]  type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = cosld   []
  [eqtime] type = MCMSolarPostprocessor  material = mcm_rates_material  solar_param = eqtime  []
  [lat]    type = Receiver  default = 51.51    []
  [lon]    type = Receiver  default = 0.13     []

  # ── Species concentrations — ElementAverageValue extracts the single-element
  # FEM value as a scalar for CSV output.  Column names match the box-mode
  # gold file so CSVDiff compares the intersection.
  [CH3NO3]    type = ElementAverageValue  variable = CH3NO3    []
  [CH3O]      type = ElementAverageValue  variable = CH3O      []
  [CH3O2]     type = ElementAverageValue  variable = CH3O2     []
  [CH3O2NO2]  type = ElementAverageValue  variable = CH3O2NO2  []
  [CH3OH]     type = ElementAverageValue  variable = CH3OH     []
  [CH3OOH]    type = ElementAverageValue  variable = CH3OOH    []
  [CH4]       type = ElementAverageValue  variable = CH4       []
  [CL]        type = ElementAverageValue  variable = CL        []
  [CO]        type = ElementAverageValue  variable = CO        []
  [H2]        type = ElementAverageValue  variable = H2        []
  [H2O2]      type = ElementAverageValue  variable = H2O2      []
  [HCHO]      type = ElementAverageValue  variable = HCHO      []
  [HNO3]      type = ElementAverageValue  variable = HNO3      []
  [HO2]       type = ElementAverageValue  variable = HO2       []
  [HO2NO2]    type = ElementAverageValue  variable = HO2NO2    []
  [HONO]      type = ElementAverageValue  variable = HONO      []
  [HSO3]      type = ElementAverageValue  variable = HSO3      []
  [N2O5]      type = ElementAverageValue  variable = N2O5      []
  [NA]        type = ElementAverageValue  variable = NA        []
  [NO]        type = ElementAverageValue  variable = NO        []
  [NO2]       type = ElementAverageValue  variable = NO2       []
  [NO3]       type = ElementAverageValue  variable = NO3       []
  [O]         type = ElementAverageValue  variable = O         []
  [O1D]       type = ElementAverageValue  variable = O1D       []
  [O3]        type = ElementAverageValue  variable = O3        []
  [OH]        type = ElementAverageValue  variable = OH        []
  [SA]        type = ElementAverageValue  variable = SA        []
  [SO2]       type = ElementAverageValue  variable = SO2       []
  [SO3]       type = ElementAverageValue  variable = SO3       []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  dt = 90
  end_time = 43200
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  [TimeStepper]
    type = ConstantDT
    dt = 90
  []
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_AtChem2_inorg_coupled'
  time_step_interval = 1
[]