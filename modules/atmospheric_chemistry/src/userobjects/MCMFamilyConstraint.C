//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MCMFamilyConstraint.h"

registerMooseObject("AtmosphericChemistryApp", MCMFamilyConstraint);

InputParameters
MCMFamilyConstraint::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Defines chemical family conservation groups using the DAE method. "
      "Each family's first member is treated as an algebraic slack variable "
      "to enforce conservation of the total family concentration. "
      "Analogous to F0AM Mass_eval.m / dydt_eval.m family section.");

  params.addParam<std::vector<std::string>>("family_names", {},
      "Names of chemical families (e.g. 'NOx', 'Ox')");
  params.addParam<std::vector<std::vector<std::string>>>(
      "family_members", {},
      "Member species for each family. First member is the DAE slack variable.");
  params.addParam<std::vector<std::vector<Real>>>(
      "family_scaling", {},
      "Scaling/weighting factors for each family member.");
  params.addParam<std::vector<std::string>>("species_list", {},
      "Full list of species names in mechanism order for index resolution.");

  return params;
}

MCMFamilyConstraint::MCMFamilyConstraint(const InputParameters & params)
  : GeneralUserObject(params)
{
  auto names = getParam<std::vector<std::string>>("family_names");
  auto members = getParam<std::vector<std::vector<std::string>>>("family_members");
  auto scaling = getParam<std::vector<std::vector<Real>>>("family_scaling");
  _all_species = getParam<std::vector<std::string>>("species_list");

  if (names.size() != members.size())
    mooseError("MCMFamilyConstraint: family_names (", names.size(),
               ") and family_members (", members.size(), ") must have same length");
  if (names.size() != scaling.size())
    mooseError("MCMFamilyConstraint: family_names (", names.size(),
               ") and family_scaling (", scaling.size(), ") must have same length");

  for (unsigned int i = 0; i < names.size(); ++i)
  {
    if (members[i].empty())
      mooseError("MCMFamilyConstraint: Family '", names[i], "' has no members");
    if (members[i].size() != scaling[i].size())
      mooseError("MCMFamilyConstraint: Family '", names[i],
                 "' has ", members[i].size(), " members but ",
                 scaling[i].size(), " scaling factors");
    addFamily(names[i], members[i], scaling[i]);
  }
}

void
MCMFamilyConstraint::resolveIndices()
{
  _slack_set.clear();

  for (auto & pair : _families)
  {
    Family & fam = pair.second;
    fam.member_indices.clear();

    for (const auto & name : fam.member_names)
    {
      auto it = std::find(_all_species.begin(), _all_species.end(), name);
      if (it != _all_species.end())
        fam.member_indices.push_back((unsigned int)(it - _all_species.begin()));
      else
        mooseError("MCMFamilyConstraint: Species '", name, "' not found in mechanism");
    }

    if (fam.member_indices.empty())
      mooseError("MCMFamilyConstraint: Family '", fam.name, "' has no valid members");

    fam.slack_idx = fam.member_indices[0];
    _slack_set.insert(fam.slack_idx);
  }
}

void
MCMFamilyConstraint::addFamily(const std::string & name,
                                const std::vector<std::string> & member_names,
                                const std::vector<Real> & scaling_factors)
{
  Family fam;
  fam.name = name;
  fam.member_names = member_names;
  fam.scaling = scaling_factors;
  if (!member_names.empty())
    fam.slack_idx = 0; // first member is slack, resolved later
  else
    fam.slack_idx = 0;
  _families[name] = fam;
  _family_names.push_back(name);
}

int
MCMFamilyConstraint::slackIndex(const std::string & family_name) const
{
  auto it = _families.find(family_name);
  if (it == _families.end())
    return -1;
  return (int)it->second.slack_idx;
}

const std::vector<unsigned int> &
MCMFamilyConstraint::memberIndices(const std::string & family_name) const
{
  static const std::vector<unsigned int> empty;
  auto it = _families.find(family_name);
  if (it != _families.end())
    return it->second.member_indices;
  return empty;
}

const std::vector<Real> &
MCMFamilyConstraint::scalingFactors(const std::string & family_name) const
{
  static const std::vector<Real> empty;
  auto it = _families.find(family_name);
  if (it != _families.end())
    return it->second.scaling;
  return empty;
}
