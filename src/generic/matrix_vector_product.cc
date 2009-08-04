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

  // (re)build distribution_pt
  Distribution_pt->rebuild(matrix_pt->distribution_pt());

#ifdef HAVE_TRILINOS

  // Has MPI been initialised?
  if (MPI_Helpers::MPI_has_been_initialised)
   {
#ifdef OOMPH_HAS_MPI // MPI defined and initialised - can use TRILINOS
    // create the communicator
    Epetra_comm_pt = new Epetra_MpiComm(
     matrix_pt->distribution_pt()->communicator_pt()->mpi_comm());

    // create the rows map
    TrilinosHelpers::create_epetra_map(matrix_pt->distribution_pt(),
                                       Epetra_comm_pt,
                                       Epetra_range_map_pt,Global_rows);
     
    // create the cols map
    Column_distribution_pt = new LinearAlgebraDistribution
     (matrix_pt->distribution_pt()->communicator_pt(),
      matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());
    TrilinosHelpers::create_epetra_map(Column_distribution_pt,Epetra_comm_pt,
                                       Epetra_domain_map_pt,Global_cols);

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
#endif
   }
  else
   {
#ifdef OOMPH_HAS_MPI // MPI is defined but not initialised - can't use TRILINOS

    setup_oomph_method_helper(matrix_pt);

//     // create the cols map
//     Column_distribution_pt = new LinearAlgebraDistribution
//      (matrix_pt->distribution_pt()->communicator_pt(),
//       matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());

//     // No trilinos, so copy the oomph-lib matrix
//     Oomph_matrix_pt = new CRDoubleMatrix(matrix_pt->distribution_pt());
//     double* values_pt = matrix_pt->value();
//     int* column_indices = matrix_pt->column_index();
//     int* row_start = matrix_pt->row_start();
//     unsigned nnz = matrix_pt->nnz();
//     unsigned nrow = matrix_pt->nrow();
//     double* my_values_pt = new double[nnz];
//     int* my_column_indices = new int[nnz];
//     int* my_row_start = new int[nrow+1];
//     for (unsigned i = 0; i < nnz; i++)
//      {
//       my_values_pt[i] = values_pt[i];
//      }
//     for (unsigned i = 0; i < nnz; i++)
//      {
//       my_column_indices[i] = column_indices[i];
//      }
//     for (unsigned i = 0; i <= nrow; i++)
//      {
//       my_row_start[i] = row_start[i];
//      }
//     Ncol = matrix_pt->ncol();
//     Oomph_matrix_pt->build_matrix_without_copy(Ncol,nnz,my_values_pt,
//                                                my_column_indices,my_row_start);
//     Ncol = matrix_pt->ncol();

#else // MPI has not been defined or initialised - use serial TRILINOS

    // create the communicator
    Epetra_comm_pt = new Epetra_SerialComm();

    // create the rows map
    TrilinosHelpers::create_epetra_map(matrix_pt->distribution_pt(),
                                       Epetra_comm_pt,Epetra_range_map_pt);
     
    // create the cols map
    Column_distribution_pt = new LinearAlgebraDistribution
     (matrix_pt->distribution_pt()->communicator_pt(),
      matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());
    TrilinosHelpers::create_epetra_map(Column_distribution_pt,Epetra_comm_pt,
                                       Epetra_domain_map_pt); 

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
#endif
   }

#else // TRILINOS not being used, so use oomph methods

  setup_oomph_method_helper(matrix_pt);
//   // create the cols map
//   Column_distribution_pt = new LinearAlgebraDistribution
//    (matrix_pt->distribution_pt()->communicator_pt(),
//     matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());

//   // if no trilinos then copy the oomph-lib matrix
//   Oomph_matrix_pt = new CRDoubleMatrix(matrix_pt->distribution_pt());
//   double* values_pt = matrix_pt->value();
//   int* column_indices = matrix_pt->column_index();
//   int* row_start = matrix_pt->row_start();
//   unsigned nnz = matrix_pt->nnz();
//   unsigned nrow = matrix_pt->nrow();
//   double* my_values_pt = new double[nnz];
//   int* my_column_indices = new int[nnz];
//   int* my_row_start = new int[nrow+1];
//   for (unsigned i = 0; i < nnz; i++)
//    {
//     my_values_pt[i] = values_pt[i];
//    }
//   for (unsigned i = 0; i < nnz; i++)
//    {
//     my_column_indices[i] = column_indices[i];
//    }
//   for (unsigned i = 0; i <= nrow; i++)
//    {
//     my_row_start[i] = row_start[i];
//    }
//   Ncol = matrix_pt->ncol();
//   Oomph_matrix_pt->build_matrix_without_copy(Ncol,nnz,my_values_pt,
//                                              my_column_indices,my_row_start);
//   Ncol = matrix_pt->ncol();
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

  // Only use (parallel) Trilinos if MPI has been defined and initialised
  if (MPI_Helpers::MPI_has_been_initialised)
   {
#ifdef OOMPH_HAS_MPI
    multiply_helper(x,y);

//     // convert x to Trilinos vector
//     Epetra_Vector* epetra_x_pt;
//     TrilinosHelpers::create_epetra_vector(x,Epetra_domain_map_pt,
//                                           epetra_x_pt,true);

//     // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
//     // matrix)
//     Epetra_Vector* epetra_soln_pt;
//     TrilinosHelpers::create_epetra_vector(y,Epetra_range_map_pt,
//                                           epetra_soln_pt,true);

//     // do the multiply
//     int epetra_error_flag = Epetra_matrix_pt->Multiply(false,*epetra_x_pt,
//                                                        *epetra_soln_pt);

//     // throw error if there is an epetra error
// #if PARANOID
//     if (epetra_error_flag != 0)
//      {
//       std::ostringstream error_message;
//       error_message 
//        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
//        << epetra_error_flag;
//       throw OomphLibError(error_message.str(),
//                           "TrilinosHelpersMPI::multiply(...)",
//                           OOMPH_EXCEPTION_LOCATION);
//      }
// #endif

//     // do something with the error flag to keep the compiler from
//     // complaining in non-paranoid mode
//     epetra_error_flag=0;

//     // clean up
//     delete epetra_x_pt;
//     delete epetra_soln_pt;
#endif
   }
  else // MPI has not been initialised
   {
#ifdef OOMPH_HAS_MPI // MPI has been defined - can't use Trilinos
    Oomph_matrix_pt->multiply(x,y);
#else // MPI not defined - can use Trilinos (serial)

    multiply_helper(x,y);

//     // convert x to Trilinos vector
//     Epetra_Vector* epetra_x_pt;
//     TrilinosHelpers::create_epetra_vector(x,Epetra_domain_map_pt,
//                                           epetra_x_pt,true);
//     // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
//     // matrix)
//     Epetra_Vector* epetra_soln_pt;
//     TrilinosHelpers::create_epetra_vector(y,Epetra_range_map_pt,
//                                           epetra_soln_pt,true);

//     // do the multiply
//     int epetra_error_flag = Epetra_matrix_pt->Multiply(false,*epetra_x_pt,
//                                                        *epetra_soln_pt);

//     // throw error if there is an epetra error
// #if PARANOID
//     if (epetra_error_flag != 0)
//      {
//       std::ostringstream error_message;
//       error_message 
//        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
//        << epetra_error_flag;
//       throw OomphLibError(error_message.str(),
//                           "TrilinosHelpersMPI::multiply(...)",
//                           OOMPH_EXCEPTION_LOCATION);
//      }
// #endif

//     // do something with the error flag to keep the compiler from
//     // complaining in non-paranoid mode
//     epetra_error_flag=0;

//     // clean up
//     delete epetra_x_pt;
//     delete epetra_soln_pt;
#endif
   }

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

  // Only use (parallel) Trilinos if MPI has been defined and initialised
  if (MPI_Helpers::MPI_has_been_initialised)
   {
#ifdef OOMPH_HAS_MPI

    multiply_transpose_helper(x,y);

//     // convert x to Trilinos vector
//     Epetra_Vector* epetra_x_pt;
//     TrilinosHelpers::create_epetra_vector(x,Epetra_range_map_pt,
//                                           epetra_x_pt,true);

//     // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
//     // matrix)
//     Epetra_Vector* epetra_soln_pt;
//     TrilinosHelpers::create_epetra_vector(y,Epetra_domain_map_pt,
//                                           epetra_soln_pt,true);

//     // do the multiply
//     int epetra_error_flag = Epetra_matrix_pt->Multiply(true,*epetra_x_pt,
//                                                        *epetra_soln_pt);

//     // throw error if there is an epetra error
// #if PARANOID
//     if (epetra_error_flag != 0)
//      {
//       std::ostringstream error_message;
//       error_message 
//        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
//        << epetra_error_flag;
//       throw OomphLibError(error_message.str(),
//                           "TrilinosHelpersMPI::multiply(...)",
//                           OOMPH_EXCEPTION_LOCATION);
//      }
// #endif

//     // do something with the error flag to keep the compiler from
//     // complaining in non-paranoid mode
//     epetra_error_flag=0;

//     // clean up
//     delete epetra_x_pt;
//     delete epetra_soln_pt;
#endif
   }
  else // MPI has not been initialised
   {
#ifdef OOMPH_HAS_MPI // MPI is defined - can't use Trilinos
    Oomph_matrix_pt->multiply_transpose(x,y);
#else // MPI is not defined - can use (serial) Trilinos

    multiply_transpose_helper(x,y);

//     // convert x to Trilinos vector
//     Epetra_Vector* epetra_x_pt;
//     TrilinosHelpers::create_epetra_vector(x,Epetra_range_map_pt,
//                                           epetra_x_pt,true);

//     // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
//     // matrix)
//     Epetra_Vector* epetra_soln_pt;
//     TrilinosHelpers::create_epetra_vector(y,Epetra_domain_map_pt,
//                                           epetra_soln_pt,true);

//     // do the multiply
//     int epetra_error_flag = Epetra_matrix_pt->Multiply(true,*epetra_x_pt,
//                                                        *epetra_soln_pt);

//     // throw error if there is an epetra error
// #if PARANOID
//     if (epetra_error_flag != 0)
//      {
//       std::ostringstream error_message;
//       error_message 
//        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = " 
//        << epetra_error_flag;
//       throw OomphLibError(error_message.str(),
//                           "TrilinosHelpersMPI::multiply(...)",
//                           OOMPH_EXCEPTION_LOCATION);
//      }
// #endif

//     // do something with the error flag to keep the compiler from
//     // complaining in non-paranoid mode
//     epetra_error_flag=0;

//     // clean up
//     delete epetra_x_pt;
//     delete epetra_soln_pt;
#endif
   }

#else
  // just multiply using the copied oomph-lib matrix
  Oomph_matrix_pt->multiply_transpose(x,y);
#endif
 }


 //============================================================================
 /// \short Setup the matrix vector product operator for oomph_method.\n
 //============================================================================
 void MatrixVectorProduct::setup_oomph_method_helper(CRDoubleMatrix* matrix_pt)
 {
  // create the cols map
  Column_distribution_pt = new LinearAlgebraDistribution
   (matrix_pt->distribution_pt()->communicator_pt(),
    matrix_pt->ncol(),matrix_pt->distribution_pt()->distributed());

  // No trilinos, so copy the oomph-lib matrix
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
  Oomph_matrix_pt->build_matrix_without_copy(Ncol,nnz,my_values_pt,
                                             my_column_indices,my_row_start);
  Ncol = matrix_pt->ncol();
 }

#ifdef HAVE_TRILINOS
 //============================================================================
 /// \short Apply the operator to the vector x and return
 /// the result in the vector y (helper function)
 //============================================================================
 void MatrixVectorProduct::multiply_helper(const DoubleVector& x, 
                                           DoubleVector& y)
 {
  // convert x to Trilinos vector
  Epetra_Vector* epetra_x_pt;
  TrilinosHelpers::create_epetra_vector(x,Epetra_domain_map_pt,
                                        epetra_x_pt,true);

  // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
  // matrix)
  Epetra_Vector* epetra_soln_pt;
  TrilinosHelpers::create_epetra_vector(y,Epetra_range_map_pt,
                                        epetra_soln_pt,true);

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
 }

 //============================================================================
 /// \short Apply the transpose of the operator to the vector x and return
 /// the result in the vector y (helper function)
 //============================================================================
 void MatrixVectorProduct::multiply_transpose_helper(const DoubleVector& x, 
                                                     DoubleVector& y)
 {

  // convert x to Trilinos vector
  Epetra_Vector* epetra_x_pt;
  TrilinosHelpers::create_epetra_vector(x,Epetra_range_map_pt,
                                        epetra_x_pt,true);

  // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
  // matrix)
  Epetra_Vector* epetra_soln_pt;
  TrilinosHelpers::create_epetra_vector(y,Epetra_domain_map_pt,
                                        epetra_soln_pt,true);

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
 }
#endif

}
