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

#include <vector>
#include <cstdlib>

// Forward declaration
struct ParsedMechanism;

/**
 * Lightweight sparse stoichiometric matrix with pluggable storage format.
 *
 * Follows PETSc's parameter-driven format selection pattern (cf. -mat_type aij|dense):
 * a Format enum selects the storage backend at construction time.  All compute
 * methods use a single template iteration interface (forEachInRow) that the
 * compiler fully inlines — zero virtual-call overhead regardless of format.
 *
 * Supported formats:
 *   CSR   — Compressed Sparse Row, net stoichiometric coefficients.
 *           PETSc AIJ-compatible, compact, optimal for HPC with large mechanisms.
 *   COO   — AtChem2-style split reactant / product vectors (clhs/crhs).
 *           Enables separate loss-rate / production-rate diagnostics.
 *   DENSE — Dense 2D array.  Simple, cache-friendly, zero indirect-addressing
 *           overhead.  Best for tiny mechanisms (< ~50 species) where memory is
 *           irrelevant and instruction-level efficiency matters.
 *   CSC   — Compressed Sparse Column: species-major storage.  Each column
 *           lists (reaction, coeff) for all reactions involving that species.
 *           Enables column queries ("which reactions involve species X?").
 *           Row iteration builds a lightweight forward index at build() time.
 */
struct StoichMatrix
{
  enum Format { CSR, COO, DENSE, CSC };

  Format format;

  unsigned int nSpecies = 0, nReactions = 0;

  // ---- CSR data ----
  std::vector<int>    csr_cols;
  std::vector<Real>   csr_vals;
  std::vector<size_t> csr_row_ptr;

  // ---- COO data (AtChem2-style clhs / crhs) ----
  std::vector<int>    lhs_species, rhs_species;
  std::vector<Real>   lhs_coeff,   rhs_coeff;
  std::vector<size_t> lhs_row_ptr, rhs_row_ptr;

  // ---- DENSE data ----
  /// dense[r][s] = net stoichiometric coefficient (0.0 for non-participating)
  std::vector<std::vector<Real>> dense;

  // ---- CSC data (species-major, column queries) ----
  /// Column s spans csc_col_ptr[s] .. csc_col_ptr[s+1]-1.
  /// Element k: reaction = csc_rows[k], coefficient = csc_c_vals[k].
  std::vector<int>    csc_rows;
  std::vector<Real>   csc_c_vals;
  std::vector<size_t> csc_col_ptr;

  /// Build from parsed mechanism.
  void build(const ParsedMechanism & mech, Format fmt);

  /**
   * Iterate non-zero stoichiometric entries for reaction r.
   * Calls fn(species_index, net_coefficient) for each participating species.
   */
  template <typename F>
  void forEachInRow(unsigned int r, F && fn) const
  {
    switch (format)
    {
      case CSR:
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          fn(csr_cols[k], csr_vals[k]);
        break;
      case COO:
        for (size_t k = lhs_row_ptr[r]; k < lhs_row_ptr[r + 1]; ++k)
          fn(lhs_species[k], -lhs_coeff[k]);
        for (size_t k = rhs_row_ptr[r]; k < rhs_row_ptr[r + 1]; ++k)
          fn(rhs_species[k], rhs_coeff[k]);
        break;
      case DENSE:
        for (unsigned int s = 0; s < nSpecies; ++s)
          if (std::abs(dense[r][s]) > 1e-30)
            fn((int)s, dense[r][s]);
        break;
      case CSC:
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          fn(csr_cols[k], csr_vals[k]);
        break;
    }
  }

  /// O(k) lookup — k = entries per row (typically 2-10).  Diagnostic use only.
  Real get(unsigned int r, unsigned int s) const
  {
    if (r >= nReactions || s >= nSpecies)
      return 0.0;

    switch (format)
    {
      case CSR:
        for (size_t k = csr_row_ptr[r]; k < csr_row_ptr[r + 1]; ++k)
          if ((unsigned int)csr_cols[k] == s) return csr_vals[k];
        return 0.0;
      case COO:
        for (size_t k = lhs_row_ptr[r]; k < lhs_row_ptr[r + 1]; ++k)
          if ((unsigned int)lhs_species[k] == s) return -lhs_coeff[k];
        for (size_t k = rhs_row_ptr[r]; k < rhs_row_ptr[r + 1]; ++k)
          if ((unsigned int)rhs_species[k] == s) return rhs_coeff[k];
        return 0.0;
      case DENSE:
        return dense[r][s];
      case CSC:
        for (size_t k = csc_col_ptr[s]; k < csc_col_ptr[s + 1]; ++k)
          if ((unsigned int)csc_rows[k] == r) return csc_c_vals[k];
        return 0.0;
    }
    return 0.0;
  }
};
