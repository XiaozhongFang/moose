//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MechanismLoader.h"
#include "MCMFacsimileParser.h"

#include <fstream>
#include <sstream>
#include <algorithm>

MechanismData
MechanismLoader::load(const std::string & mechanism_file,
                       const std::string & photo_path,
                       const std::string & mcm_version,
                       const std::string & peroxy_path,
                       const std::vector<std::string> & input_file_dirs)
{
  MechanismData data;

  // ── Parse the .fac mechanism file ────────────────────────────────────
  MCMFacsimileParser parser;
  parser.setMCMVersion(mcm_version);
  ParsedMechanism mech = parser.parse(mechanism_file, photo_path, peroxy_path);

  // ── Copy species ─────────────────────────────────────────────────────
  data.species = mech.species;
  data.ro2_species = mech.ro2_species;

  // ── Copy reactions (convert ParsedReaction → MechanismData::Reaction) ─
  for (auto & r : mech.reactions)
  {
    Reaction rx;
    rx.rate_expression = r.rate_expression;
    rx.reactants = r.reactants;
    rx.products = r.products;
    data.reactions.push_back(rx);
  }
  data.stoichiometric_matrix = mech.stoichiometry;

  // ── Copy rate coefficients ───────────────────────────────────────────
  for (unsigned int i = 0; i < mech.coefficient_names.size(); ++i)
  {
    data.rate_coefficients[mech.coefficient_names[i]] = mech.coefficient_expressions[i];
    data.converted_coefficients[mech.coefficient_names[i]] = mech.coefficient_expressions[i];
    data.coefficient_names.insert(mech.coefficient_names[i]);
  }
  data.eval_order = mech.coefficient_names;

  // ── Copy reaction rate expressions ───────────────────────────────────
  data.reaction_rate_expressions.resize(mech.reactions.size());
  for (unsigned int i = 0; i < mech.reactions.size(); ++i)
    data.reaction_rate_expressions[i] = mech.reactions[i].rate_expression;

  // ── Copy photolysis parameters (mechanism-referenced only) ───────────
  for (unsigned int i = 0; i < mech.j_numbers.size(); ++i)
  {
    std::string jkey = "J<" + std::to_string(mech.j_numbers[i]) + ">";
    data.j_numbers.push_back(mech.j_numbers[i]);
    data.j_CL.push_back(mech.j_CL[i]);
    data.j_CMM.push_back(mech.j_CMM[i]);
    data.j_CNN.push_back(mech.j_CNN[i]);
  }

  // ── Build base variable set ──────────────────────────────────────────
  data.base_variables = {"TEMP", "M", "O2", "N2", "H2O"};
  for (auto & [jname, _] : data.rate_coefficients)
  {
    // jname starts with "J<" for photolysis rate coefficients
    if (jname.size() > 2 && jname[0] == 'J' && jname[1] == '<')
      data.base_variables.insert(jname);
  }

  // ── Resolve photolysis file path ─────────────────────────────────────
  {
    // Determine the mechanism file directory for fallback resolution
    std::string mech_dir;
    auto pos = mechanism_file.find_last_of("/\\");
    if (pos != std::string::npos)
      mech_dir = mechanism_file.substr(0, pos);

    data.resolved_photo_path = resolvePhotoPath(photo_path, input_file_dirs, mech_dir);
  }

  // ── Load full photolysis parameter set ───────────────────────────────
  // The parser only transfers J<N> entries referenced by the mechanism.
  // We need the full set for MCMPhotolysisPostprocessor / MCMRatesMaterial.
  if (!data.resolved_photo_path.empty())
  {
    loadFullPhotolysisSet(data.resolved_photo_path,
                           data.j_numbers_all,
                           data.j_cl_values,
                           data.j_cmm_values,
                           data.j_cnn_values);
  }

  return data;
}

std::string
MechanismLoader::resolvePhotoPath(const std::string & photo_path,
                                   const std::vector<std::string> & input_file_dirs,
                                   const std::string & mech_file_dir)
{
  // If empty, nothing to resolve
  if (photo_path.empty())
    return photo_path;

  // Try as-is (relative to working directory)
  {
    std::ifstream test_file(photo_path);
    if (test_file.good())
      return photo_path;
  }

  // Try relative to each input file directory
  for (const auto & input_file : input_file_dirs)
  {
    auto pos = input_file.find_last_of("/\\");
    if (pos != std::string::npos)
    {
      std::string resolved = input_file.substr(0, pos) + "/" + photo_path;
      std::ifstream test_file(resolved);
      if (test_file.good())
        return resolved;
    }
  }

  // Try relative to the mechanism file's directory
  if (!mech_file_dir.empty())
  {
    auto bname_pos = photo_path.find_last_of("/\\");
    std::string base = (bname_pos != std::string::npos)
                           ? photo_path.substr(bname_pos + 1)
                           : photo_path;
    std::string resolved = mech_file_dir + "/" + base;
    std::ifstream test_file(resolved);
    if (test_file.good())
      return resolved;
  }

  // Fall back to the original path
  return photo_path;
}

void
MechanismLoader::loadFullPhotolysisSet(const std::string & resolved_path,
                                        std::vector<unsigned int> & j_numbers,
                                        std::vector<Real> & j_cl,
                                        std::vector<Real> & j_cmm,
                                        std::vector<Real> & j_cnn)
{
  std::ifstream pfile(resolved_path);
  if (!pfile.good())
    return;

  std::string line;
  std::getline(pfile, line); // skip header
  while (std::getline(pfile, line))
  {
    if (line.empty() || line[0] == '#')
      continue;

    // Convert Fortran D-notation (6.073D-05) → E-notation (6.073E-05)
    std::replace(line.begin(), line.end(), 'D', 'E');
    std::replace(line.begin(), line.end(), 'd', 'e');

    // Parse "j l m n name tau" columns (we only need first 4)
    std::istringstream iss(line);
    unsigned int jn;
    double cl, cmm, cnn;
    std::string unused1, unused2;
    if (iss >> jn >> cl >> cmm >> cnn)
    {
      j_numbers.push_back(jn);
      j_cl.push_back(cl);
      j_cmm.push_back(cmm);
      j_cnn.push_back(cnn);
    }
  }
}
