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
// Functions for OomphCommunicator

// MPI headers
#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// Oomph-lib error handler
#include "communicator.h"
#include "matrices.h"

namespace oomph
{
#ifdef OOMPH_HAS_MPI

  //=============================================================================
  /// A broadcast function for DenseMatrix<double>
  //=============================================================================
  void OomphCommunicator::broadcast(const int& source, DenseMatrix<double>& x)
  {
    // Get number of entries on processor source (where the matrix exists)
    unsigned nrow, ncol;
    if (this->my_rank() == source)
    {
      nrow = x.nrow();
      if (nrow > 0)
      {
        ncol = x.ncol();
      }
      else
      {
        ncol = 0;
      }
    }

    // Broadcast to everybody how many entries to expect
    MPI_Bcast(&nrow, 1, MPI_UNSIGNED_LONG, source, this->mpi_comm());
    MPI_Bcast(&ncol, 1, MPI_UNSIGNED_LONG, source, this->mpi_comm());

    if (ncol != 0 && nrow != 0)
    {
      // convert into a C-style array
      double* x_bcast = new double[nrow * ncol];

      if (this->my_rank() == source)
        for (unsigned long i = 0; i < nrow; i++)
        {
          for (unsigned long j = 0; j < ncol; j++)
          {
            x_bcast[i * ncol + j] = x(i, j);
          }
        }

      // broadcast the array
      MPI_Bcast(x_bcast, ncol * nrow, MPI_DOUBLE, source, this->mpi_comm());

      // Now convert back into matrix (everywhere apart from source)
      if (this->my_rank() != source)
      {
        x.resize(nrow, ncol);
        for (unsigned long i = 0; i < nrow; i++)
          for (unsigned long j = 0; j < ncol; j++)
          {
            x(i, j) = x_bcast[i * ncol + j];
          }
      }
      // delete C-style array
      delete[] x_bcast;
    }
    else
    {
      x.resize(0, 0);
    }
  }

#endif

} // namespace oomph
