// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#include "matrix_vector_product.h"

namespace oomph
{
  //============================================================================
  /// Setup the matrix vector product operator.
  /// WARNING: This class is wrapper to Trilinos Epetra matrix vector
  /// multiply methods, if Trilinos is not installed then this class will
  /// function as expected, but there will be no computational speed gain.
  /// By default the Epetra_CrsMatrix::multiply(...) are employed.
  /// The optional argument col_dist_pt is the distribution of:
  /// x if using multiply(...) or y if using multiply_transpose(...)
  /// where this is A x = y. By default, this is assumed to the uniformly
  /// distributed based on matrix_pt->ncol().
  //============================================================================
  void MatrixVectorProduct::setup(CRDoubleMatrix* matrix_pt,
                                  const LinearAlgebraDistribution* col_dist_pt)
  {
    // clean memory
    this->clean_up_memory();

    // (re)build distribution_pt
    this->build_distribution(matrix_pt->distribution_pt());

    // store the number of columns
    Ncol = matrix_pt->ncol();

    // determine whether we are using trilinos
    Using_trilinos = false;
#ifdef OOMPH_HAS_TRILINOS
#ifdef OOMPH_HAS_MPI
    if (MPI_Helpers::mpi_has_been_initialised())
    {
      Using_trilinos = true;
    }
#else
    Using_trilinos = true;
#endif
#endif

    // create the column distribution map, if a distribution is not provided,
    // create a uniformly distributed based on matrix_pt->ncol().
    // Otherwise, use the provided distribution.
    if (col_dist_pt == 0)
    {
      Column_distribution_pt = new LinearAlgebraDistribution(
        matrix_pt->distribution_pt()->communicator_pt(),
        matrix_pt->ncol(),
        matrix_pt->distribution_pt()->distributed());
    }
    else
    {
      Column_distribution_pt = new LinearAlgebraDistribution(col_dist_pt);
    }

    // setup the operator
    if (Using_trilinos)
    {
#ifdef OOMPH_HAS_TRILINOS
      double t_start = TimingHelpers::timer();
      Epetra_matrix_pt =
        TrilinosEpetraHelpers::create_distributed_epetra_matrix(
          matrix_pt, Column_distribution_pt);
      double t_end = TimingHelpers::timer();
      oomph_info << "Time to build epetra matrix [sec] : " << t_end - t_start
                 << std::endl;
#endif
    }
    else
    {
      double t_start = TimingHelpers::timer();
      Oomph_matrix_pt = new CRDoubleMatrix(*matrix_pt);
      double t_end = TimingHelpers::timer();
      oomph_info << "Time to copy CRDoubleMatrix [sec] : " << t_end - t_start
                 << std::endl;
    }
  }

  //============================================================================
  /// Apply the operator to the vector x and return the result in
  /// the vector y
  //============================================================================
  void MatrixVectorProduct::multiply(const DoubleVector& x,
                                     DoubleVector& y) const
  {
#ifdef PARANOID
    // check that the distribution of x is setup
    if (!x.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The distribution of the vector x must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Check to see if the distribution of the matrix column is the same as the
    // distribution of the vector to be operated on.
    if (*this->Column_distribution_pt != *x.distribution_pt())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The distribution of the x Vector is not the same as"
        << " the column distribution.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // if y is setup then it should have the same distribution as
    // the matrix used to set this up.
    if (y.built())
    {
      if (!(*y.distribution_pt() == *this->distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The y vector is setup and therefore must have the same "
          << "distribution as the matrix used to set up the "
             "MatrixVectorProduct";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if y is not setup then setup the distribution
    if (!y.built())
    {
      // Resize and initialize the solution vector
      y.build(this->distribution_pt(), 0.0);
    }

    // apply the operator
    if (Using_trilinos)
    {
#ifdef OOMPH_HAS_TRILINOS
      trilinos_multiply_helper(x, y);
#endif
    }
    else
    {
      Oomph_matrix_pt->multiply(x, y);
    }
  }

  //============================================================================
  /// Apply the transpose of the operator to the vector x and return
  /// the result in the vector y
  //============================================================================
  void MatrixVectorProduct::multiply_transpose(const DoubleVector& x,
                                               DoubleVector& y) const
  {
#ifdef PARANOID
    // check that the distribution of x is setup
    if (!x.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The distribution of the vector x must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // Check to see if x.size() = ncol()
    if (*this->distribution_pt() != *x.distribution_pt())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "This class assumes that the y vector has a uniform "
        << "distributed distribution.";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // if y is setup then it should have the same distribution as x
    if (y.built())
    {
      if (!(*y.distribution_pt() == *this->Column_distribution_pt))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The y vector is setup and therefore must have the same "
          << "distribution as the vector x";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if y is not setup then setup the distribution
    if (!y.built())
    {
      // Resize and initialize the solution vector
      y.build(this->Column_distribution_pt, 0.0);
    }

    // apply the transpose operator
    if (Using_trilinos)
    {
#ifdef OOMPH_HAS_TRILINOS
      trilinos_multiply_transpose_helper(x, y);
#endif
    }
    else
    {
      Oomph_matrix_pt->multiply_transpose(x, y);
    }
  }

#ifdef OOMPH_HAS_TRILINOS
  //============================================================================
  /// Apply the operator to the vector x and return
  /// the result in the vector y (helper function)
  //============================================================================
  void MatrixVectorProduct::trilinos_multiply_helper(const DoubleVector& x,
                                                     DoubleVector& y) const
  {
    // convert x to a Trilinos Epetra vector.
    // x is const so it much be copied.
    Epetra_Vector* epetra_x_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(x);

    // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
    // matrix)
    Epetra_Vector* epetra_soln_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(y);

    // do the multiply
#ifdef PARANOID
    int epetra_error_flag = 0;
    epetra_error_flag =
#endif
      Epetra_matrix_pt->Multiply(false, *epetra_x_pt, *epetra_soln_pt);

    // throw error if there is an epetra error
#ifdef PARANOID
    if (epetra_error_flag != 0)
    {
      std::ostringstream error_message;
      error_message
        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = "
        << epetra_error_flag;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // return solution
    TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_soln_pt, y);

    // clean up
    delete epetra_x_pt;
    delete epetra_soln_pt;
  }

  //============================================================================
  /// Apply the transpose of the operator to the vector x and return
  /// the result in the vector y (helper function)
  //============================================================================
  void MatrixVectorProduct::trilinos_multiply_transpose_helper(
    const DoubleVector& x, DoubleVector& y) const
  {
    // convert x to a Trilinos Epetra vector.
    // x is const so it much be copied.
    Epetra_Vector* epetra_x_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(x);

    // create Trilinos vector for soln ('viewing' the contents of the oomph-lib
    // matrix)
    Epetra_Vector* epetra_soln_pt =
      TrilinosEpetraHelpers::create_distributed_epetra_vector(y);

    // do the multiply
#ifdef PARANOID
    int epetra_error_flag = 0;
    epetra_error_flag =
#endif
      Epetra_matrix_pt->Multiply(true, *epetra_x_pt, *epetra_soln_pt);

    // throw error if there is an epetra error
#ifdef PARANOID
    if (epetra_error_flag != 0)
    {
      std::ostringstream error_message;
      error_message
        << "Epetra Matrix Vector Multiply Error : epetra_error_flag = "
        << epetra_error_flag;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // copy to solution vector
    TrilinosEpetraHelpers::copy_to_oomphlib_vector(epetra_soln_pt, y);

    // clean up
    delete epetra_x_pt;
    delete epetra_soln_pt;
  }
#endif
} // namespace oomph
