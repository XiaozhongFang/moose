//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"
#include <string>
#include <vector>
#include <set>

/**
 * Manages AtChem2-style constrained species.
 * Constrained species are fixed to observed values and not solved.
 */
class ConstraintHandler
{
public:
  ConstraintHandler() = default;

  void setConstrained(const std::vector<std::string> & names,
                      const std::vector<std::string> & all_species);

  void updateValues(const std::vector<Real> & values);

  bool isConstrained(unsigned int idx) const { return _constrained_set.count(idx) > 0; }
  unsigned int nConstrained() const { return _constrained_set.size(); }
  const std::vector<Real> & values() const { return _constrained_values; }

private:
  std::set<unsigned int> _constrained_set;
  std::vector<Real> _constrained_values;
};
