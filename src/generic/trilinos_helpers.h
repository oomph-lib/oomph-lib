//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
#ifndef OOMPH_TRILINOS_HELPERS_HEADER
#define OOMPH_TRILINOS_HELPERS_HEADER

// trilinos headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_MatrixMatrix.h"
#include "ml_epetra_utils.h"
#ifdef OOMPH_HAS_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// oomph-lib headers
#include "../generic/matrices.h"

// mpi includes
#ifdef OOMPH_HAS_MPI
#include "double_vector.h"
#include "linear_algebra_distribution.h"
#endif

namespace oomph
{

//// forward declaration of oomph-lib matrices
class CRDoubleMatrix;
class DoubleMatrixBase;

//=============================================================================
// Helper functions for use with the Trilinos library - primarily functions to
// convert oomph-lib matrices and vectors to Trilinos (Epetra) matrices and 
// vectors.
//=============================================================================
namespace TrilinosHelpers
{

 /// \short creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt. \n
 /// NOTE. the Epetra_Map and the LinearAlgebraDistribution object in the
 /// vector must descibe the same distribution \n
 void create_epetra_vector(const DoubleVector& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt,
                           bool view = false);

 /// \short creates an empty trilinos vector with Epetra_Map row_map_pt. \n
 /// NOTE. elements are NOT set to zero
 void create_epetra_vector(const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt);

 /// \short Helper function to copy the contents of a Trilinos vector to an
 /// oomph-lib distributed vector. The distribution of the two vectors must
 /// be identical \n
 void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                              DoubleVector& oomph_lib_v);

 /// \short Helper function to create a distributed Epetra_CrsMatrix from a 
 /// from a CRDoubleMatrix \n
 /// NOTE 1. This function constructs an Epetra_CrsMatrix using new,
 /// "delete trilinos_matrix_pt;" is NOT called. \n
 void create_epetra_matrix(DoubleMatrixBase* oomph_matrix,
                           const Epetra_Map* epetra_range_map_pt,
                           const Epetra_Map* epetra_domain_map_pt,
                           const Epetra_Map* epetra_col_map_pt,
                           Epetra_CrsMatrix* &epetra_matrix_pt,
                           bool view = false);
 
 /// \short Helper function to create a distributed Epetra_CrsMatrix from a 
 /// from a CRDoubleMatrix \n
 /// NOTE 1. This function constructs an Epetra_CrsMatrix using new,
 /// "delete trilinos_matrix_pt;" is NOT called. \n
 /// NOTE 2. Specialisation for SQUARE matrices. \n
 void create_epetra_matrix(DoubleMatrixBase* oomph_matrix,
                           const Epetra_Map* epetra_range_pt,
                           const Epetra_Map* epetra_col_map_pt,
                           Epetra_CrsMatrix* &epetra_matrix_pt);

 /// \short  Function to perform a matrix-vector multiplication on a 
 /// oomph-lib matrix and vector using Trilinos functionality.\n
 /// NOTE 1. the matrix (matrix) and the vectors (x and soln) must have the 
 /// same communicator.\n
 /// NOTE 2. The vector (soln) will be returned with the same distribution
 /// as the matrix, unless a distribution is predefined in the solution
 /// vector in which case the vector will be returned with that distribution.\n
 void multiply(CRDoubleMatrix &matrix,
               const DoubleVector& x,
               DoubleVector &soln);

 /// \short Function to perform a matrix-matrix multiplication on oomph-lib
 /// matrices by using Trilinos functionality.\n
 /// \b NOTE 1. There are two Trilinos matrix-matrix multiplication methods 
 /// available, using either the EpetraExt::MatrixMatrix class (if use_ml == 
 /// false) or using ML (Epetra_MatrixMult method)
 /// \b NOTE 2. the solution matrix (matrix_soln) will be returned with the 
 /// same distribution as matrix1
 /// \b NOTE 3. All matrices must share the same communicator. 
 void multiply(CRDoubleMatrix &matrix_1,
               CRDoubleMatrix &matrix_2,
               CRDoubleMatrix &matrix_soln,
               const bool& use_ml = false);

#ifdef OOMPH_HAS_MPI
 /// \short takes a oomph-lib distribution and returns a trilinos mpi map 
 /// (Epetra_Map)
 void create_epetra_map(const LinearAlgebraDistribution* distribution_pt,
                        const Epetra_MpiComm* epetra_comm_pt,
                        Epetra_Map* &epetra_map_pt, int* &my_global_rows);
#else
  /// \short creates a serial epetra_map of size nrow
 void create_epetra_map(const LinearAlgebraDistribution* distribution_pt,
                        const Epetra_SerialComm* epetra_comm_pt,
                        Epetra_Map* &epetra_map_pt);
#endif
   
#ifdef PARANOID
 /// \short Function to compare a distribution described by an oomph-lib 
 /// DistributionInfo object to a distribution described by a Trilinos
 /// Epetra_Map object. Returns true if they describe the same distribution.
 /// Used as a PARANOID check in other TrilinosHelpers methods
 bool compare_maps(const LinearAlgebraDistribution* oomph_distribution_pt, 
                   const Epetra_Map& epetra_map);
#endif
} // end of trilinos helpers namespace
} // end of namspace oomph
#endif
