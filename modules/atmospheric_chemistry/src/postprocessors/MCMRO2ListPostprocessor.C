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
      "from the parsed mechanism. Produces vectors for both count and indices.");
  return params;
}

MCMRO2ListPostprocessor::MCMRO2ListPostprocessor(const InputParameters & params)
  : GeneralVectorPostprocessor(params),
    _box_model(getUserObject<MCMBoxModel>("box_model")),
    _ro2_count(declareVector("ro2_count")),
    _ro2_species(declareVector("ro2_species"))
{
}

void
MCMRO2ListPostprocessor::execute()
{
  const auto & indices = _box_model.getRO2Indices();
  const auto & names = _box_model.getRO2Species();

  // RO2 count (single-element vector)
  _ro2_count.resize(1);
  _ro2_count[0] = (Real)indices.size();

  // RO2 species indices
  _ro2_species.resize(indices.size());
  for (std::size_t i = 0; i < indices.size(); ++i)
    _ro2_species[i] = (Real)indices[i];

  // Print RO2 species NAMES to console in a structured, parseable format.
  // The check_ro2.py --run mode parses this line to verify detection.
  std::ostringstream oss;
  for (std::size_t i = 0; i < names.size(); ++i)
  {
    if (i > 0) oss << ",";
    oss << names[i];
  }
  _console << "RO2_SPECIES(" << names.size() << "): " << oss.str() << std::endl;
}
