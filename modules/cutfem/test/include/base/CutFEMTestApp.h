#pragma once

#include "MooseApp.h"

class CutFEMTestApp : public MooseApp
{
public:
  static InputParameters validParams();

  CutFEMTestApp(const InputParameters & parameters);
  virtual ~CutFEMTestApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s,
                          bool use_test_objs = false);
};
