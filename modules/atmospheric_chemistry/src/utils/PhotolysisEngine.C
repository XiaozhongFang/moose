//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PhotolysisEngine.h"

#include "MCMPhoto.h"
#include "HybridPhoto.h"
#include "BottomUpPhoto.h"

std::unique_ptr<PhotolysisEngine>
PhotolysisEngine::create(Method method, const PhotolysisParams & pparams)
{
  switch (method)
  {
    case MCM_SZA:
      return std::make_unique<MCMPhoto>(pparams.j_CL,
                                         pparams.j_CMM,
                                         pparams.j_CNN,
                                         pparams.j_numbers);

    case HYBRID:
      return std::make_unique<HybridPhoto>(pparams.hybrid_table_dir,
                                            pparams.hybrid_j_numbers);

    case BOTTOMUP:
      return std::make_unique<BottomUpPhoto>(pparams.bottomup_data_dir,
                                              pparams.bottomup_flux_file,
                                              pparams.j_numbers);

    default:
      return nullptr;
  }
}
