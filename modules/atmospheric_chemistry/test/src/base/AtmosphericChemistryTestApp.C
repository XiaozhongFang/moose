//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryTestApp.h"
#include "AtmosphericChemistryApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
AtmosphericChemistryTestApp::validParams()
{
  InputParameters params = AtmosphericChemistryApp::validParams();
  return params;
}

AtmosphericChemistryTestApp::AtmosphericChemistryTestApp(const InputParameters & parameters)
  : AtmosphericChemistryApp(parameters)
{
  AtmosphericChemistryTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

AtmosphericChemistryTestApp::~AtmosphericChemistryTestApp() {}

void
AtmosphericChemistryTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s,
                                          bool use_test_objects)
{
  AtmosphericChemistryApp::registerAll(f, af, s);
  if (use_test_objects)
  {
    Registry::registerObjectsTo(f, {"AtmosphericChemistryTestApp"});
    Registry::registerActionsTo(af, {"AtmosphericChemistryTestApp"});
  }
}

void
AtmosphericChemistryTestApp::registerApps()
{
  registerApp(AtmosphericChemistryApp);
  registerApp(AtmosphericChemistryTestApp);
}
