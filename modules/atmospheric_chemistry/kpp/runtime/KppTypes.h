// kpp/runtime/KppTypes.h
#pragma once

#include "Moose.h"

/// Type mapping between MOOSE Real and KPP_REAL (double).
/// KPP generates code using KPP_REAL = double.
using KppReal = double;

/// KPP integer type for indexing.
using KppInt = int;
