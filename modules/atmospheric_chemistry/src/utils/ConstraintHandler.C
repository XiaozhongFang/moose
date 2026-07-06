//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConstraintHandler.h"
#include <algorithm>

void
ConstraintHandler::setConstrained(const std::vector<std::string> & names,
                                   const std::vector<std::string> & all_species)
{
  _constrained_set.clear();
  for (const auto & name : names)
  {
    auto it = std::find(all_species.begin(), all_species.end(), name);
    if (it != all_species.end())
      _constrained_set.insert(static_cast<unsigned int>(it - all_species.begin()));
  }
  _constrained_values.resize(_constrained_set.size(), 0.0);
}

void
ConstraintHandler::updateValues(const std::vector<Real> & values)
{
  if (values.size() == _constrained_set.size())
    _constrained_values = values;
}
