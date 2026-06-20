//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "Function.h"

/**
 * Kernel that constrains a chemical species concentration to a
 * prescribed time-dependent function.  Used for species whose
 * concentrations are fixed by observations rather than solved by
 * the chemical ODE system (AtChem2 "constrained species" mode).
 *
 * The residual is R = u - f(t), which enforces u = f(t) at each
 * quadrature point/Newton iteration.
 */
class MCMConstraintKernel : public Kernel
{
public:
  static InputParameters validParams();
  MCMConstraintKernel(const InputParameters & params);

protected:
  Real computeQpResidual() override;
  Real computeQpJacobian() override;

  /// Time-dependent function providing the constrained value
  const Function & _func;
};
