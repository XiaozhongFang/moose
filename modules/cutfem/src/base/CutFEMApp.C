#include "CutFEMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
CutFEMApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  return params;
}

registerKnownLabel("CutFEMApp");

CutFEMApp::CutFEMApp(const InputParameters & parameters) : MooseApp(parameters)
{
  srand(processor_id());
  CutFEMApp::registerAll(_factory, _action_factory, _syntax);
}

CutFEMApp::~CutFEMApp() {}

void
CutFEMApp::registerApps()
{
  registerApp(CutFEMApp);
}

void
CutFEMApp::registerAll(Factory & f, ActionFactory & af, Syntax & /*s*/)
{
  Registry::registerObjectsTo(f, {"CutFEMApp"});
  Registry::registerActionsTo(af, {"CutFEMApp"});
}

extern "C" void
CutFEMApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  CutFEMApp::registerAll(f, af, s);
}

extern "C" void
CutFEMApp__registerApps()
{
  CutFEMApp::registerApps();
}
