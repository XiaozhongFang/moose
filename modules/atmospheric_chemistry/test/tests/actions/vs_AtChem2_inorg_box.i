# vs AtChem2 CVODE — MCM inorganic — box mode
# Full reproduction of AtChem2 model_mcm_inorg:
#   29 species, 71 reactions, 12h diurnal cycle (48 steps × 900s)
#   SZA photolysis, T=288.15K, P=1013.25mbar
#   lat=51.51, lon=0.13, 2010-06-21

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 1
  []
[]

[AtmosphericChemistry]
  [Box]
  mechanism_file = '../../../doc/content/modules/atmospheric_chemistry/database/atchem2_example.fac'
  temperature = 288.15
  air_density = 2.55e19
  water_vapor = 4.45e17
  mcm_photolysis_file = '../../../doc/content/modules/atmospheric_chemistry/database/mcm_photolysis_rates_v3.3.1.dat'
  latitude = 51.51
  longitude = 0.13
  day = 21
  month = 6
  year = 2010
  press = 1013.25
  []
[]

[ICs]
  [CH4_ic]  type = ScalarConstantIC  variable = CH4  value = 4.9e13  []
  [CO_ic]   type = ScalarConstantIC  variable = CO   value = 3.6e12  []
  [O3_ic]   type = ScalarConstantIC  variable = O3   value = 5.2e11  []
  [NO2_ic]  type = ScalarConstantIC  variable = NO2  value = 2.4e11  []
[]

[Postprocessors]
  # ── Photolysis rates — all J values from AtChem2 photolysisRates.output ──
  [J1]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 1   []
  [J2]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 2   []
  [J3]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 3   []
  [J4]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 4   []
  [J5]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 5   []
  [J6]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 6   []
  [J7]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 7   []
  [J8]   type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 8   []
  [J11]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 11  []
  [J12]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 12  []
  [J13]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 13  []
  [J14]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 14  []
  [J15]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 15  []
  [J16]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 16  []
  [J17]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 17  []
  [J18]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 18  []
  [J19]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 19  []
  [J20]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 20  []
  [J21]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 21  []
  [J22]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 22  []
  [J23]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 23  []
  [J24]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 24  []
  [J31]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 31  []
  [J32]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 32  []
  [J33]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 33  []
  [J34]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 34  []
  [J35]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 35  []
  [J41]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 41  []
  [J51]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 51  []
  [J52]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 52  []
  [J53]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 53  []
  [J54]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 54  []
  [J55]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 55  []
  [J56]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 56  []
  [J61]  type = MCMPhotolysisPostprocessor  box_model = box_model  j_number = 61  []

  # ── Environment variables (AtChem2 environmentVariables.output order) ──
  [M_env]     type = Receiver  default = 2.546917e19 []
  [TEMP]      type = Receiver  default = 288.15     []
  [PRESS]     type = Receiver  default = 1013.25    []
  [RH]        type = Receiver  default = -1.0        []
  [H2O]       type = Receiver  default = 4.45e17    []
  [DEC]       type = MCMSolarPostprocessor  box_model = box_model  solar_param = dec  []
  [BLHEIGHT]  type = Receiver  default = -1.0        []
  [DILUTE]    type = Receiver  default = -1.0        []
  [JFAC]      type = Receiver  default = 1.0        []
  [ROOF]      type = Receiver  default = 1.0        []
  [ASA]       type = Receiver  default = -1.0        []
  [RO2_sum]   type = MCMRO2Postprocessor  box_model = box_model  
    species_variables = 'HCHO CH3NO3 CH3OH O1D O3 HO2NO2 NO3 N2O5 H2O2 NO NA HO2 NO2 CH4 HSO3 CO CL O HNO3 SO3 SO2 CH3O OH H2 HONO CH3O2NO2 CH3OOH SA CH3O2'
  []

  # ── Solar parameters (AtChem2 zenith_data_mod) ──
  [cosx]   type = MCMSolarPostprocessor  box_model = box_model  solar_param = cosx    []
  [secx]   type = MCMSolarPostprocessor  box_model = box_model  solar_param = secx    []
  [lha]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lha     []
  [sinld]  type = MCMSolarPostprocessor  box_model = box_model  solar_param = sinld   []
  [cosld]  type = MCMSolarPostprocessor  box_model = box_model  solar_param = cosld   []
  [eqtime] type = MCMSolarPostprocessor  box_model = box_model  solar_param = eqtime  []
  [lat]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lat     []
  [lon]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lon     []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'bdf2'
  end_time = 43200
  dt = 90
  l_max_its = 50
  l_tol = 1e-5
  nl_max_its = 10
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  [TimeStepper]
    type = ConstantDT
    dt = 90
  []
[]

[Outputs]
  checkpoint = false
  csv = true
  execute_on = 'timestep_end'
  file_base = 'vs_AtChem2_inorg_box'
  time_step_interval = 1
[]
