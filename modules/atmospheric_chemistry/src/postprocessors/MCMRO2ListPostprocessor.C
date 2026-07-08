//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMRO2ListPostprocessor.h"

#include <sstream>

registerMooseObject("AtmosphericChemistryApp", MCMRO2ListPostprocessor);

InputParameters
MCMRO2ListPostprocessor::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  params.addRequiredParam<UserObjectName>("box_model",
                                          "Name of the MCMBoxModel UserObject");
  params.addClassDescription(
      "Outputs the RO2 (peroxy radical) species list detected by MCMBoxModel "
      "from the parsed mechanism. Produces one vector per detected species name.");
  return params;
}

MCMRO2ListPostprocessor::MCMRO2ListPostprocessor(const InputParameters & params)
  : GeneralVectorPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _ro2_count(declareVector("ro2_count"))
{
}

void
MCMRO2ListPostprocessor::execute()
{
  const auto & names = _box_model.getRO2Species();

  _ro2_count.resize(1);
  _ro2_count[0] = (Real)names.size();

  if (_ro2_species_flags.empty())
    for (const auto & name : names)
      _ro2_species_flags.push_back(&declareVector(name));

  for (auto * flag : _ro2_species_flags)
  {
    flag->resize(1);
    (*flag)[0] = 1.0;
  }

  // Print RO2 species NAMES to console in a structured, parseable format.
  // The check_ro2.py --run-app mode parses this line to verify detection.
  std::ostringstream oss;
  for (std::size_t i = 0; i < names.size(); ++i)
  {
    if (i > 0) oss << ",";
    oss << names[i];
  }
  _console << "RO2_SPECIES(" << names.size() << "): " << oss.str() << std::endl;
}
