# MCMSolarPostprocessor

!syntax description /Postprocessors/MCMSolarPostprocessor

`MCMSolarPostprocessor` outputs solar parameters computed by `MCMBoxModel`
during photolysis rate evaluation.  Mirrors AtChem2's `zenith_data_mod` variables
for validation and diagnostics.

## Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `box_model` | UserObjectName | Name of the MCMBoxModel UserObject |
| `solar_param` | MooseEnum | Parameter to output: `cosx` (cosine of SZA), `secx` (secant of SZA), `lha` (local hour angle), `sinld` ($\sin\lambda\sin\delta$), `cosld` ($\cos\lambda\cos\delta$), `eqtime` (equation of time), `lat` (latitude), `lon` (longitude) |

## Description

The solar parameters are computed using Madronich (1993) parameterization and
cached during `MCMBoxModel::evaluateCoefficients()`.  One `MCMSolarPostprocessor`
instance outputs a single parameter; use multiple instances to output all eight.

Values for `lat` and `lon` are the simulation parameters (degrees).  All other
parameters vary with simulation time and are updated at each timestep.

## Example Input File Syntax

```moose
[Postprocessors]
  [cosx]   type = MCMSolarPostprocessor  box_model = box_model  solar_param = cosx   []
  [secx]   type = MCMSolarPostprocessor  box_model = box_model  solar_param = secx   []
  [lha]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lha    []
  [sinld]  type = MCMSolarPostprocessor  box_model = box_model  solar_param = sinld  []
  [cosld]  type = MCMSolarPostprocessor  box_model = box_model  solar_param = cosld  []
  [eqtime] type = MCMSolarPostprocessor  box_model = box_model  solar_param = eqtime []
  [lat]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lat    []
  [lon]    type = MCMSolarPostprocessor  box_model = box_model  solar_param = lon    []
[]
```

!syntax parameters /Postprocessors/MCMSolarPostprocessor

!syntax inputs /Postprocessors/MCMSolarPostprocessor

!syntax children /Postprocessors/MCMSolarPostprocessor
