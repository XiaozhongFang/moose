#pragma once

#include "MooseApp.h"

/**
 * CutFEM Module Application Class
 *
 * Cut Finite Element Methods (CutFEM) for unfitted mesh discretization
 * of PDEs with complex geometries and evolving interfaces.
 */
class CutFEMApp : public MooseApp
{
public:
  static InputParameters validParams();

  CutFEMApp(const InputParameters & parameters);
  virtual ~CutFEMApp();

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};
