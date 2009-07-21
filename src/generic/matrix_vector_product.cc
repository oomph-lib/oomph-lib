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
#include "matrix_vector_product.h"

namespace oomph
{


 //============================================================================
 /// \short Setup the matrix vector product operator.\n
 /// WARNING: This class is wrapper to Trilinos Epetra matrix vector
 /// multiply methods, if Trilinos is not installed then this class will
 /// function as expected, but there will be no computational speed gain.\n
 /// By default the Epetra_CrsMatrix::multiply(...) are employed.\n
 //============================================================================
 void MatrixVectorProduct::setup(CRDoubleMatrix* matrix_pt)
 {
  // clean memory
  this->clean_up_memory();

  //
  Distribution_pt->rebuild(matrix_pt->distribution_pt());

#ifdef HAVE_TRILINOS

  // create the communicator
#ifdef OOMPH_HAS_MPI
  Epetra_comm_pt = new Epetra_MpiComm(
   matrix_pt->distribution_pt()->communicator_pt()->mpi_comm());
#else
  Epetra_comm_pt = new Epetra_SerialComm();
#endif

  // create the rows map
#ifdef OOMPH_HAS_MPI
  TrilinosHelpers::create_epetra_map(matrix_pt->distribution_pt(),
                                     Epetra_comm_pt,
                                     Epetra_range_map_pt,Global_rows);
#else
  TrilinosHelpers::create_epetra_map(matrix_pt->distribution_pt(),
                                     Epetra_comm_pt,Epetra_range_map_pt);
#endif
     
  // create the cols map
  Column_distribution_pt = new LinearAlgebraDistribution
   (matrix_pt->distribution_pt()->communicator_pt(),
    matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());
#ifdef OOMPH_HAS_MPI
  TrilinosHelpers::create_epetra_map(Column_distribution_pt,Epetra_comm_pt,
                                     Epetra_domain_map_pt,Global_cols); 
#else
  TrilinosHelpers::create_epetra_map(Column_distribution_pt,Epetra_comm_pt,
                                     Epetra_domain_map_pt); 
#endif

  // convert epetra matrix
  Epetra_col_map_pt = new Epetra_Map(matrix_pt->ncol(),
                                     matrix_pt->ncol(),
                                     0,*Epetra_comm_pt);
  double t_start = TimingHelpers::timer();
  TrilinosHelpers::create_epetra_matrix(matrix_pt,Epetra_range_map_pt,
                                        Epetra_domain_map_pt,
                                        Epetra_col_map_pt,
                                        Epetra_matrix_pt,false);
  double t_end = TimingHelpers::timer();
  oomph_info << "Time to build epetra matrix [sec] : "
             << t_end - t_start << std::endl;

  // store the number of columns
  Ncol = matrix_pt->ncol();


#else
  // create the cols map
  Column_distribution_pt = new LinearAlgebraDistribution
   (matrix_pt->distribution_pt()->communicator_pt(),
    matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());

  // if no trilinos then copy the oomph-lib matrix
  Oomph_matrix_pt = new CRDoubleMatrix(matrix_pt->distribution_pt());
  double* values_pt = matrix_pt->value();
  int* column_indices = matrix_pt->column_index();
  int* row_start = matrix_pt->row_start();
  unsigned nnz = matrix_pt->nnz();
  unsigned nrow = matrix_pt->nrow();
  double* my_values_pt = new double[nnz];
  int* my_column_indices = new int[nnz];
  int* my_row_start = new int[nrow+1];
  for (unsigned i = 0; i < nnz; i++)
   {
    my_values_pt[i] = values_pt[i];
   }
  for (unsigned i = 0; i < nnz; i++)
   {
    my_column_indices[i] = column_indices[i];
   }
  for (unsigned i = 0; i <= nrow; i++)
   {
    my_row_start[i] = row_start[i];
   }
  Ncol = matrix_pt->ncol();
  Oomph_matrix_pt->rebuild_matrix_without_copy(Ncol,nnz,my_values_pt,
                                               my_column_indices,my_row_start);
  Ncol = matrix_pt->ncol();
#endif
 }

 //============================================================================
 /// \short Apply the operator to the vector x and return the result in 
 /// the vector y
 //============================================================================
 void MatrixVectorProduct::multiply(const DoubleVector& x, 
                                    DoubleVector& y)
 {
#ifdef PARANOID
  // check that the distribution of x is setup
  if (!x.distribution_setup())
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "The distribution of the vector x must be setup";
    throw OomphLibError(error_message_stream.str(),
                        "MatrixVectorProduct::multiply()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  // Check to see if x.size() = ncol().
  if (*this->Column_distribution_pt != *x.distribution_pt())
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "This class assumes that the x vector has a uniform "
     << "distributed distribution.";
    throw OomphLibError(error_message_stream.str(),
                        "MatrixVectorProduct::multiply()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  // if y is setup then it should have the same distribution as x
  if (y.distribution_setup())
   {
    if (!(*y.distribution_pt() == *this->Distribution_pt))
     {
      std::ostringstream error_message_stream;
      error_message_stream 
       << "The y vector is setup and therefore must have the same "
       << "distribution as the vector x";
      throw OomphLibError(error_message_stream.str(),
                          "MatrixVectorProduct::multiply()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // if y is not setup then setup the distribution
  if (!y.distribution_setup())
   {
    // Resize and initialize the solution vector
    y.build(this->distribution_pt(),0.0);
   }

#ifdef HAVE_TRILINOS

  // convert x to Trilinos vector
  Epetra_Vector* epetra_x_pt;
  TrilinosHelpers::create_epetra_vector(x,Epetra_domain_map_pt,epetra_x_pt,true);
  // create Trilinos vector for soln ( 'viewing' the contents of the oomph-lib
  // matrix)
  Epetra_Vector* epetra_soln_pt;
  TrilinosHelpers::create_epetra_vector(y,Epetra_range_map_pt,epetra_soln_pt,true);

  // do the multiply
  int epetra_error_flag = Epetra_matrix_pt->Multiply(false,*epetra_x_pt,
                                                     *epetra_soln_pt);

  // throw error if there is an epetra error
#if PARANOID
  if (epetra_error_flag != 0)
   {
    std::ostringstream error_message;
    error_message 
     << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
     << epetra_error_flag;
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpersMPI::multiply(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // do something with the error flag to keep the compiler from
  // complaining in non-paranoid mode
  epetra_error_flag=0;

  // clean up
  delete epetra_x_pt;
  delete epetra_soln_pt;

#else
  // just multiply using the copied oomph-lib matrix
  Oomph_matrix_pt->multiply(x,y);
#endif
 }

 //============================================================================
 /// \short Apply the transpose of the operator to the vector x and return 
 /// the result in the vector y
 //============================================================================
 void MatrixVectorProduct::multiply_transpose(const DoubleVector& x, 
                                              DoubleVector& y)
 {

#ifdef PARANOID
  // check that the distribution of x is setup
  if (!x.distribution_setup())
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "The distribution of the vector x must be setup";
    throw OomphLibError(error_message_stream.str(),
                        "MatrixVectorProduct::multiply()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  // Check to see if x.size() = ncol()
  if (*this->Distribution_pt != *x.distribution_pt())
   {
    std::ostringstream error_message_stream;
    error_message_stream 
     << "This class assumes that the y vector has a uniform "
     << "distributed distribution.";
    throw OomphLibError(error_message_stream.str(),
                        "MatrixVectorProduct::multiply()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  // if y is setup then it should have the same distribution as x
  if (y.distribution_setup())
   {
    if (!(*y.distribution_pt() == *this->Column_distribution_pt))
     {
      std::ostringstream error_message_stream;
      error_message_stream 
       << "The y vector is setup and therefore must have the same "
       << "distribution as the vector x";
      throw OomphLibError(error_message_stream.str(),
                          "MatrixVectorProduct::multiply()",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // if y is not setup then setup the distribution
  if (!y.distribution_setup())
   {
    // Resize and initialize the solution vector
    y.build(this->Column_distribution_pt,0.0);
   }

#ifdef HAVE_TRILINOS

  // convert x to Trilinos vector
  Epetra_Vector* epetra_x_pt;
  TrilinosHelpers::create_epetra_vector(x,Epetra_range_map_pt,epetra_x_pt,true);

  // create Trilinos vector for soln ( 'viewing' the contents of the oomph-lib
  // matrix)
  Epetra_Vector* epetra_soln_pt;
  TrilinosHelpers::create_epetra_vector(y,Epetra_domain_map_pt,epetra_soln_pt,
                                        true);

  // do the multiply
  int epetra_error_flag = Epetra_matrix_pt->Multiply(true,*epetra_x_pt,
                                                     *epetra_soln_pt);

  // throw error if there is an epetra error
#if PARANOID
  if (epetra_error_flag != 0)
   {
    std::ostringstream error_message;
    error_message 
     << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
     << epetra_error_flag;
    throw OomphLibError(error_message.str(),
                        "TrilinosHelpersMPI::multiply(...)",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // do something with the error flag to keep the compiler from
  // complaining in non-paranoid mode
  epetra_error_flag=0;

  // clean up
  delete epetra_x_pt;
  delete epetra_soln_pt;

#else
  // just multiply using the copied oomph-lib matrix
  Oomph_matrix_pt->multiply_transpose(x,y);
#endif
 }
}
