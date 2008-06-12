//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include "distributed_vector.h"
#include "distribution_info.h"
#endif

namespace oomph
{

//// forward declaration of oomph-lib matrices
class CRDoubleMatrix;
class DistributedCRDoubleMatrix;
class DoubleMatrixBase;

//=============================================================================
// Helper functions for use with the Trilinos library - primarily functions to
// convert oomph-lib matrices and vectors to Trilinos (Epetra) matrices and 
// vectors.
//=============================================================================
namespace TrilinosHelpers
{



#ifdef OOMPH_HAS_MPI
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // DISTRIBUTED Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////




 /// \short creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt. \n
 /// \b NOTE 1. the Epetra_Map and the DistributionInfo objects must descibe
 /// the same distribution
 /// \b NOTE 2. if the bool copy is true (default) then the values in the 
 /// oomph-lib vector (oomph_v) will be copied into the epetra vector, or if 
 /// it is false then the epetra vector will 'view' (i.e. point to) the 
 /// contents of the oomph-lib vector
 void create_epetra_vector(DistributedVector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt,
                           const bool& copy = true);


 /// \short creates a distributed Epetra_Vector that has the same contents as
 /// oomph_v with distribution row_map_pt, the contents of the trilinos vector 
 /// are copied from the oomph-lib vector. \n
 /// \b NOTE 1. duplication of helper due to the const oomph-lib vector
 /// \b NOTE 2. the Epetra_Map and the DistributionInfo objects must descibe
 /// the same distribution
 void create_epetra_vector(const DistributedVector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt);

 /// \short Helper function to copy the contents of a Trilinos vector to an
 /// oomph-lib distributed vector. The my_distribution DistributionInfo should 
 /// be the same as the distribution of epetra_v_pt
 void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                              const DistributionInfo my_distribution, 
                              DistributedVector<double>& oomph_lib_v);


 /// \short  Function to perform a matrix-vector multiplication on a 
 /// distributed matrix and a distributed vector using Trilinos functionality.
 /// \n
 /// \b NOTE 1. the matrix (matrix) and the vectors (x and soln) must have the 
 /// same communicator.
 /// \b NOTE 2. The vector (soln) will be returned with the same distribution
 /// as the matrix, unless a distribution is predefined in the solution
 /// vector in which case the vector will be returned with that distribution
 void multiply(DistributedCRDoubleMatrix &matrix,
               const DistributedVector<double>& x,
               DistributedVector<double> &soln);

 /// \short Function to perform a matrix-matrix multiplication on distributed 
 /// matrices by using Trilinos functionality.\n
 /// \b NOTE 1. There are two Trilinos matrix-matrix multiplication methods 
 /// available, using either the EpetraExt::MatrixMatrix class (if use_ml == 
 /// false) or using ML (Epetra_MatrixMult method)
 /// \b NOTE 2. the solution matrix (matrix_soln) will be returned with the 
 /// same distribution as matrix1
 /// \b NOTE 3. All matrices must share the same communicator. 
 void multiply(DistributedCRDoubleMatrix &matrix_1,
               DistributedCRDoubleMatrix &matrix_2,
               DistributedCRDoubleMatrix &matrix_soln,
               const bool& use_ml = false);

 /// \short takes a oomph-lib distribution (DistributionInfo) and returns
 /// a trilinos map (Epetra_Map)
 void create_epetra_map(const DistributionInfo& my_distribution,
                        const Epetra_MpiComm* epetra_comm_pt,
                        Epetra_Map* &epetra_map_pt, int* &my_global_rows);
                       
 
 /// \short Function to compare a distribution described by an oomph-lib 
 /// DistributionInfo object to a distribution described by a Trilinos
 /// Epetra_Map object. Returns true if they describe the same distribution.
 /// Used as a PARANOID check in other TrilinosHelpers methods
 bool compare_maps(const DistributionInfo& oomph_distribution, 
                   const Epetra_Map& epetra_map);


#else




 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // SERIAL Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////




 /// \short creates a serial epetra_map of size nrow
 void create_epetra_map(const unsigned& nrow, 
                        const Epetra_SerialComm* epetra_comm_pt,
                        Epetra_Map* &epetra_map_pt);
                       



#endif 




 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////
 // SERIAL and DISTRIBUTED Helpers
 //////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////




 /// \short creates an Epetra_Vector from a serial oomph-lib vector oomph_v
 /// with map row_map_pt.
 void create_epetra_vector(const Vector<double>& oomph_v,
                           const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt);

 /// creates an empty trilinos vector with Epetra_Map row_map_pt
 void create_epetra_vector(const Epetra_Map* row_map_pt,
                           Epetra_Vector*& epetra_v_pt);

 /// \short Helper function to copy the contents of a epetra vector to an
 /// oomph-lib distributed vector.
 void copy_to_oomphlib_vector(const Epetra_Vector* epetra_v_pt,
                              Vector<double>& oomph_lib_v);

 /// \short Helper function to create a distributed Epetra_CrsMatrix from a 
 /// from a DistributedCRDoubleMatrix or a CRDoubleMatrix
 /// \b NOTE 1. This function constructs an Epetra_CrsMatrix using new,
 /// "delete trilinos_matrix_pt;" is NOT called.
 void create_epetra_matrix(DoubleMatrixBase* oomph_matrix,
                           const Epetra_Map* epetra_row_map_pt,
                           const Epetra_Map* epetra_col_map_pt,
                           Epetra_CrsMatrix* &epetra_matrix_pt);
 
 /// \short Helper function to create a distributed Epetra_CrsMatrix from a 
 /// from a DistributedCRDoubleMatrix or a CRDoubleMatrix
 /// \b NOTE 1. This function constructs an Epetra_CrsMatrix using new,
 /// "delete trilinos_matrix_pt;" is NOT called.
 /// \b NOTE 2. Specialisation for SQUARE matrices.
 void create_epetra_matrix(DoubleMatrixBase* oomph_matrix,
                           const Epetra_Map* epetra_row_pt,
                           Epetra_CrsMatrix* &epetra_matrix_pt);
}
}
#endif
