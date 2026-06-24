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
  params.addRequiredParam<UserObjectName>("box_model", "Name of the MCMBoxModel UserObject");

  MooseEnum param_enum("cosx secx lha sinld cosld eqtime lat lon dec", "cosx");
  params.addRequiredParam<MooseEnum>("solar_param", param_enum,
      "Solar parameter to output: cosx, secx, lha, sinld, cosld, eqtime, lat, lon");
  params.addClassDescription("Outputs solar parameters (cosx, secx, lha, etc.) from MCMBoxModel.");
  return params;
}

MCMSolarPostprocessor::MCMSolarPostprocessor(const InputParameters & params)
  : GeneralPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _value(0.0)
{
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
}

void
MCMSolarPostprocessor::execute()
{
  switch (_param)
  {
    case COSX:   _value = _box_model.getSolarCosX();   break;
    case SECX:   _value = _box_model.getSolarSecX();   break;
    case LHA:    _value = _box_model.getSolarLHA();    break;
    case SINLD:  _value = _box_model.getSolarSinLD();  break;
    case COSLD:  _value = _box_model.getSolarCosLD();  break;
    case EQTIME: _value = _box_model.getSolarEQT();    break;
    case LAT:    _value = _box_model.getLatitude();     break;
    case LON:    _value = _box_model.getLongitude();    break;
    case DEC:    _value = _box_model.getDeclination();   break;
  }
}

Real
MCMSolarPostprocessor::getValue() const
{
  return _value;
}
