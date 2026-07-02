//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BoxIntegrator.h"
#include "MCMBoxModel.h"

MooseImplicitIntegrator::MooseImplicitIntegrator(const MCMBoxModel & box_model)
  : _box(box_model)
{
}

Real
MooseImplicitIntegrator::computeResidual(unsigned int species_idx,
                                          const std::vector<Real> & C) const
{
  return _box.getDCdt(species_idx, C);
}

Real
MooseImplicitIntegrator::computeJacobianDiagonal(unsigned int species_idx,
                                                  const std::vector<Real> & C) const
{
  return _box.getJacobianDiagonal(species_idx, C);
}

Real
MooseImplicitIntegrator::computeJacobianOffDiagonal(unsigned int species_idx,
                                                     unsigned int jvar,
                                                     const std::vector<Real> & C) const
{
  return _box.getJacobianOffDiagonal(species_idx, jvar, C);
}

void
MooseImplicitIntegrator::reinit(Real time) const
{
  _box.markDirty();
  _box.setCurrentTime(time);
}

Real
MooseImplicitIntegrator::ppbToMolec() const
{
  return _box.ppbToMolec();
}
