//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemistryMechanismSpec.h"
#include "MechanismLoader.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <set>

ChemistryMechanismSpec::ChemistryMechanismSpec(const std::string & mechanism_file,
                                               const std::string & chem_solver,
                                               const std::string & mcm_version,
                                               const std::string & photo_path,
                                               const std::string & peroxy_path,
                                               const std::vector<std::string> & input_files)
{
  _is_kpp = (chem_solver == "kpp_rosenbrock" || chem_solver == "kpp_sdirk" ||
             chem_solver == "kpp_runge_kutta");
  if (!_is_kpp && mechanism_file.size() > 4 &&
      mechanism_file.substr(mechanism_file.size() - 4) == ".kpp")
    _is_kpp = true;

  if (_is_kpp)
  {
    _species = parseKPPSpecies(mechanism_file);
  }
  else
  {
    MechanismData data = MechanismLoader::load(
        mechanism_file, photo_path, mcm_version, peroxy_path, input_files);
    _mech_data = data;
    _species = data.species;
    _ro2_species = data.ro2_species;

    std::set<std::string> base_vars = {"TEMP", "M", "O2", "N2", "H2O"};
    for (auto jn : data.j_numbers)
      base_vars.insert("J<" + std::to_string(jn) + ">");
    _base_vars = base_vars;
  }
}

std::vector<std::string>
ChemistryMechanismSpec::parseKPPSpecies(const std::string & kpp_file) const
{
  std::vector<std::string> species;
  std::set<std::string> seen;

  std::function<void(const std::string &)> resolve = [&](const std::string & fpath) {
    // Determine mechanism file directory for relative path resolution
    auto slash = kpp_file.find_last_of('/');
    std::string mech_dir = (slash != std::string::npos) ? kpp_file.substr(0, slash + 1) : "";

    // Search order: 1) relative to mechanism dir, 2) KPP_HOME/models/, 3) cwd
    std::string abs_path = fpath;
    if (!mech_dir.empty() && fpath.find('/') == std::string::npos)
      abs_path = mech_dir + fpath;

    std::ifstream file(abs_path);
    if (!file.is_open())
    {
      // Second attempt: search relative to KPP_HOME/models/ (for #MODEL directives)
      const char * kpp_home = std::getenv("KPP_HOME");
      if (kpp_home)
      {
        abs_path = std::string(kpp_home) + "/models/" + fpath;
        file.open(abs_path);
      }
    }
    if (!file.is_open())
    {
      // Third attempt: search current working directory
      file.open(fpath);
    }
    if (!file.is_open())
      mooseError("ChemistryMechanismSpec: cannot open KPP include file '", fpath,
                 "' (searched relative to mechanism directory, KPP_HOME/models/, and cwd).");

    std::string line;
    bool in_defvar = false;
    while (std::getline(file, line))
    {
      if (line.find("#INCLUDE") == 0 || line.find("#include") == 0)
      {
        auto pos = line.find_first_of(" \t");
        if (pos != std::string::npos)
        {
          std::string incl = line.substr(pos + 1);
          incl.erase(0, incl.find_first_not_of(" \t"));
          incl.erase(incl.find_last_not_of(" \t") + 1);
          resolve(incl);
        }
        continue;
      }
      if (line.find("#MODEL") == 0 || line.find("#model") == 0)
      {
        auto pos = line.find_first_of(" \t");
        if (pos != std::string::npos)
        {
          std::string model = line.substr(pos + 1);
          model.erase(0, model.find_first_not_of(" \t"));
          model.erase(model.find_last_not_of(" \t") + 1);
          resolve(model + ".spc");
          resolve(model + ".eqn");
        }
        continue;
      }
      if (line.find("#DEFVAR") == 0) { in_defvar = true; continue; }
      if (line.find("#DEFFIX") == 0) { in_defvar = false; continue; }
      if (line.find('#') == 0) continue;
      if (!in_defvar) continue;

      auto eq = line.find('=');
      if (eq == std::string::npos) continue;
      std::string name = line.substr(0, eq);
      name.erase(0, name.find_first_not_of(" \t"));
      name.erase(name.find_last_not_of(" \t") + 1);
      if (!name.empty() && seen.insert(name).second)
        species.push_back(name);
    }
  };

  resolve(kpp_file);
  if (species.empty())
    mooseError("ChemistryMechanismSpec: no species found in ", kpp_file,
               ". Ensure the .kpp mechanism file has a #DEFVAR section.");
  return species;
}
