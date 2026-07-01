//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include <string>
#include <vector>
#include <map>
#include <set>

/**
 * Defines chemical family conservation groups using the DAE method.
 *
 * F0AM/AtChem2 analog: family conservation (NOx, Ox, NOy, etc.)
 * Reference: F0AM Mass_eval.m + dydt_eval.m family section
 *
 * Each family is defined as a weighted sum of member species:
 *   F_total = sum(member_i * scaling_i)
 *
 * The first member of each family is the "DAE slack variable":
 * its dC/dt is algebraically determined to enforce d(F_total)/dt = 0.
 * Other members evolve normally via the chemical ODE.
 *
 * Usage:
 *   NOx = NO2 + NO          → slack=NO2, members={NO2, NO}, scaling={1,1}
 *   Ox  = O3 + NO2 + 2*NO3  → slack=O3, members={O3, NO2, NO3}, scaling={1,1,2}
 */
class MCMFamilyConstraint : public GeneralUserObject
{
public:
  static InputParameters validParams();
  MCMFamilyConstraint(const InputParameters & params);

  /** Define a chemical family. Called during initialization. */
  void addFamily(const std::string & name, const std::vector<std::string> & member_names,
                 const std::vector<Real> & scaling_factors);

  /** Resolve member names to species indices using stored species list. */
  void resolveIndices();

  /** Set the list of all species names (in mechanism order) for index resolution. */
  void setSpeciesList(const std::vector<std::string> & all_species) { _all_species = all_species; }

  /** Get family names. */
  const std::vector<std::string> & familyNames() const { return _family_names; }

  /** Get DAE slack variable index for a family. Returns -1 if not found. */
  int slackIndex(const std::string & family_name) const;

  /** Get all family member indices. */
  const std::vector<unsigned int> & memberIndices(const std::string & family_name) const;

  /** Get scaling factors for each member. */
  const std::vector<Real> & scalingFactors(const std::string & family_name) const;

  /** Whether species idx is a DAE slack variable. */
  bool isSlack(unsigned int idx) const { return _slack_set.count(idx); }

  /** Species indices that are DAE slack variables. */
  const std::set<unsigned int> & slackIndices() const { return _slack_set; }

  /** Total number of family constraints. */
  unsigned int nFamilies() const { return _families.size(); }

  // GeneralUserObject interface
  void initialize() override { resolveIndices(); }
  void execute() override {}
  void finalize() override {}

private:
  struct Family
  {
    std::string name;
    std::vector<std::string> member_names;
    std::vector<unsigned int> member_indices;
    std::vector<Real> scaling;
    unsigned int slack_idx;
  };

  std::vector<std::string> _family_names;
  std::map<std::string, Family> _families;
  std::set<unsigned int> _slack_set;
  std::vector<std::string> _all_species;
};
