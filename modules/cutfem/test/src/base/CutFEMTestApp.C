#include "CutFEMTestApp.h"
#include "CutFEMApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
CutFEMTestApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  return params;
}

registerKnownLabel("CutFEMTestApp");

CutFEMTestApp::CutFEMTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  srand(processor_id());
  CutFEMTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

CutFEMTestApp::~CutFEMTestApp() {}

void
CutFEMTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  CutFEMApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"CutFEMTestApp"});
    Registry::registerActionsTo(af, {"CutFEMTestApp"});
  }
}

void
CutFEMTestApp::registerApps()
{
  CutFEMApp::registerApps();
  registerApp(CutFEMTestApp);
}

extern "C" void
CutFEMTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  CutFEMTestApp::registerAll(f, af, s);
}

extern "C" void
CutFEMTestApp__registerApps()
{
  CutFEMTestApp::registerApps();
}
