//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMSolarPostprocessor.h"

registerMooseObject("AtmosphericChemistryApp", MCMSolarPostprocessor);

InputParameters
MCMSolarPostprocessor::validParams()
{
  InputParameters params = GeneralPostprocessor::validParams();
  params.addParam<UserObjectName>("box_model", "",
      "Name of the MCMBoxModel UserObject (box mode)");
  params.addParam<MaterialPropertyName>("material", "",
      "Name of the MCMRatesMaterial block (coupled mode)");

  MooseEnum param_enum("cosx secx lha sinld cosld eqtime lat lon dec", "cosx");
  params.addRequiredParam<MooseEnum>("solar_param", param_enum,
      "Solar parameter to output: cosx, secx, lha, sinld, cosld, eqtime, lat, lon, dec");
  params.addClassDescription("Outputs solar parameters from MCMBoxModel (box mode) or MCMRatesMaterial (coupled mode).");
  return params;
}

MCMSolarPostprocessor::MCMSolarPostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(nullptr),
    _solar_prop(nullptr),
    _value(0.0)
{
  // Determine parameter type first
  MooseEnum p = getParam<MooseEnum>("solar_param");
  if (p == "cosx")   _param = COSX;
  else if (p == "secx")   _param = SECX;
  else if (p == "lha")    _param = LHA;
  else if (p == "sinld")  _param = SINLD;
  else if (p == "cosld")  _param = COSLD;
  else if (p == "eqtime") _param = EQTIME;
  else if (p == "lat")    _param = LAT;
  else if (p == "lon")    _param = LON;
  else if (p == "dec")    _param = DEC;

  // Load data source: box_model (box mode) or material (coupled mode)
  auto box_name = getParam<UserObjectName>("box_model");
  auto mat_name = getParam<MaterialPropertyName>("material");

  if (!box_name.empty())
    _box_model = &getUserObject<MCMBoxModel>("box_model");
  else if (!mat_name.empty())
  {
    // Map solar parameter → material property name
    switch (_param)
    {
      case COSX:   _solar_prop = &getMaterialProperty<Real>("solar_cosx");  break;
      case SECX:   _solar_prop = &getMaterialProperty<Real>("solar_secx");  break;
      case LHA:    _solar_prop = &getMaterialProperty<Real>("solar_lha");   break;
      case SINLD:  _solar_prop = &getMaterialProperty<Real>("solar_sinld"); break;
      case COSLD:  _solar_prop = &getMaterialProperty<Real>("solar_cosld"); break;
      case EQTIME: _solar_prop = &getMaterialProperty<Real>("solar_eqt");   break;
      case DEC:    _solar_prop = &getMaterialProperty<Real>("solar_dec");    break;
      case LAT:    _value = 51.51; break;
      case LON:    _value = 0.13;  break;
    }
  }
  else
    mooseError("MCMSolarPostprocessor: either 'box_model' or 'material' must be specified");
}

void
MCMSolarPostprocessor::execute()
{
  switch (_param)
  {
    case COSX:   _value = _box_model ? _box_model->getSolarCosX()   : (*_solar_prop)[0];   break;
    case SECX:   _value = _box_model ? _box_model->getSolarSecX()   : (*_solar_prop)[0];   break;
    case LHA:    _value = _box_model ? _box_model->getSolarLHA()    : (*_solar_prop)[0];   break;
    case SINLD:  _value = _box_model ? _box_model->getSolarSinLD()  : (*_solar_prop)[0];   break;
    case COSLD:  _value = _box_model ? _box_model->getSolarCosLD()  : (*_solar_prop)[0];   break;
    case EQTIME: _value = _box_model ? _box_model->getSolarEQT()    : (*_solar_prop)[0];   break;
    case DEC:    _value = _box_model ? _box_model->getDeclination()  : (*_solar_prop)[0];   break;
    case LAT:    _value = _box_model ? _box_model->getLatitude()    : _value; break;
    case LON:    _value = _box_model ? _box_model->getLongitude()   : _value; break;
  }
}

Real
MCMSolarPostprocessor::getValue() const
{
  return _value;
}
