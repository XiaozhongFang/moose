//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AtmosphericChemistryApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "Capabilities.h"
#include "NavierStokesApp.h"

InputParameters
AtmosphericChemistryApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("automatic_automatic_scaling") = false;
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

registerKnownLabel("AtmosphericChemistryApp");

AtmosphericChemistryApp::AtmosphericChemistryApp(const InputParameters & parameters)
  : MooseApp(parameters)
{
  AtmosphericChemistryApp::registerAll(_factory, _action_factory, _syntax);
}

AtmosphericChemistryApp::~AtmosphericChemistryApp() {}

void
AtmosphericChemistryApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  // Register navier_stokes dependency (transitively: fluid_properties, heat_transfer)
  NavierStokesApp::registerAll(f, af, syntax);

  Registry::registerObjectsTo(f, {"AtmosphericChemistryApp"});
  Registry::registerActionsTo(af, {"AtmosphericChemistryApp"});

  registerSyntax("AtmosphericChemistryBoxAction", "AtmosphericChemistry/Box");
  registerSyntax("AtmosphericChemistryCoupledAction", "AtmosphericChemistry/Coupled");
}

void
AtmosphericChemistryApp::registerApps()
{
  registerApp(AtmosphericChemistryApp);
  NavierStokesApp::registerApps();
}
