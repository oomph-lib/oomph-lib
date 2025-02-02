// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================


#include "matrices.h"
#include "double_vector.h"
#include "sum_of_matrices.h"


namespace oomph
{
  // =================================================================
  /// Matrix-vector multiplication for a sumofmatrices class. Just
  /// delegate each multiplication to the appropriate class then add up the
  /// results.
  // =================================================================
  void SumOfMatrices::multiply(const DoubleVector& x, DoubleVector& soln) const
  {
    // We assume that appropriate checks and initialisation on x and soln are
    // carried out within the individual matrix multiplys.

    // Multiply for the main matrix
    Main_matrix_pt->multiply(x, soln);

    // Now add contribution for the added matrices
    for (unsigned i_matrix = 0; i_matrix < Added_matrix_pt.size(); i_matrix++)
    {
      // If possible copy the matrix distribution, otherwise it isn't
      // distributed so make a serial LinearAlgebraDistribution object.
      LinearAlgebraDistribution col_dist, row_dist;
      OomphCommunicator serial_comm; // Serial communcator (does nothing)
      col_dist.build(&serial_comm, added_matrix_pt(i_matrix)->ncol(), false);
      row_dist.build(&serial_comm, added_matrix_pt(i_matrix)->nrow(), false);

      // Create temporary output DoubleVector
      DoubleVector temp_soln(row_dist);

      // Create a const iterator for the map (faster than .find() or []
      // access, const means can't change the map via the iterator).
      std::map<unsigned, unsigned>::const_iterator it;

      // Pull out the appropriate values into a temp vector
      //??ds not parallel
      DoubleVector temp_x(col_dist);
      for (it = Col_map_pt[i_matrix]->main_to_added_mapping_pt()->begin();
           it != Col_map_pt[i_matrix]->main_to_added_mapping_pt()->end();
           it++)
      {
        temp_x[it->second] = x[it->first];
      }

      // Perform the multiplication
      Added_matrix_pt[i_matrix]->multiply(temp_x, temp_soln);

      // Add result to solution vector
      //??ds not parallel
      for (it = Row_map_pt[i_matrix]->main_to_added_mapping_pt()->begin();
           it != Row_map_pt[i_matrix]->main_to_added_mapping_pt()->end();
           it++)
      {
        soln[it->first] += temp_soln[it->second];
      }
    }
  }

} // namespace oomph
