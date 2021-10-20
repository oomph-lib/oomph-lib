// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
#include "trilinos_helpers.h"

namespace oomph
{
  // VECTOR METHODS
  // =============================================================


  //=============================================================================
  ///  create an Epetra_Vector from an oomph-lib DoubleVector.
  /// If oomph_vec is NOT distributed (i.e. locally replicated) and
  /// on more than one processor, then the returned Epetra_Vector will be
  /// uniformly distributed. If the oomph_vec is distributed then the
  /// Epetra_Vector returned will have the same distribution as oomph_vec.
  //=============================================================================
  Epetra_Vector* TrilinosEpetraHelpers::create_distributed_epetra_vector(
    const DoubleVector& oomph_vec)
  {
#ifdef PARANOID
    // check the the oomph lib vector is setup
    if (!oomph_vec.built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib vector (oomph_v) must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // create the corresponding Epetra_Map
    LinearAlgebraDistribution* dist_pt = 0;
    if (oomph_vec.distributed())
    {
      dist_pt = new LinearAlgebraDistribution(oomph_vec.distribution_pt());
    }
    else
    {
      dist_pt = new LinearAlgebraDistribution(
        oomph_vec.distribution_pt()->communicator_pt(), oomph_vec.nrow(), true);
    }
    Epetra_Map* epetra_map_pt = create_epetra_map(dist_pt);

    // first first coefficient of the oomph vector to be inserted into the
    // Epetra_Vector
    unsigned offset = 0;
    if (!oomph_vec.distributed())
    {
      offset = dist_pt->first_row();
    }

    // copy the values into the oomph-lib vector
    // const_cast OK because Epetra_Vector construction is Copying values and
    // therefore does not modify data.
    double* v_pt = const_cast<double*>(oomph_vec.values_pt());
    Epetra_Vector* epetra_vec_pt =
      new Epetra_Vector(Copy, *epetra_map_pt, v_pt + offset);

    // clean up
    delete epetra_map_pt;
    delete dist_pt;

    // return
    return epetra_vec_pt;
  }

  //=============================================================================
  ///  create an Epetra_Vector based on the argument oomph-lib
  /// LinearAlgebraDistribution
  /// If dist is NOT distributed and
  /// on more than one processor, then the returned Epetra_Vector will be
  /// uniformly distributed. If dist is distributed then the Epetra_Vector
  /// returned will have the same distribution as dist.
  /// The coefficient values are not set.
  //=============================================================================
  Epetra_Vector* TrilinosEpetraHelpers::create_distributed_epetra_vector(
    const LinearAlgebraDistribution* dist_pt)
  {
#ifdef PARANOID
    // check the the oomph lib vector is setup
    if (!dist_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The LinearAlgebraDistribution dist_pt must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // create the corresponding Epetra_Map
    LinearAlgebraDistribution* target_dist_pt = 0;
    if (dist_pt->distributed())
    {
      target_dist_pt = new LinearAlgebraDistribution(dist_pt);
    }
    else
    {
      target_dist_pt = new LinearAlgebraDistribution(
        dist_pt->communicator_pt(), dist_pt->nrow(), true);
    }
    Epetra_Map* epetra_map_pt = create_epetra_map(target_dist_pt);

    // create epetra_vector
    Epetra_Vector* epetra_vec_pt = new Epetra_Vector(*epetra_map_pt, false);

    // clean up
    delete epetra_map_pt;
    delete target_dist_pt;

    // return
    return epetra_vec_pt;
  }

  //=============================================================================
  /// Create an Epetra_Vector equivalent of DoubleVector
  /// The argument DoubleVector must be built.
  /// The Epetra_Vector will point to, and NOT COPY the underlying data in the
  /// DoubleVector.
  /// The oomph-lib DoubleVector and the returned Epetra_Vector will have the
  /// the same distribution.
  //=============================================================================
  Epetra_Vector* TrilinosEpetraHelpers::create_epetra_vector_view_data(
    DoubleVector& oomph_vec)
  {
#ifdef PARANOID
    // check the the oomph lib vector is setup
    if (!oomph_vec.built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib vector (oomph_v) must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // create the corresponding Epetra_Map
    Epetra_Map* epetra_map_pt = create_epetra_map(oomph_vec.distribution_pt());

    // copy the values into the oomph-lib vector
    double* v_pt = oomph_vec.values_pt();
    Epetra_Vector* epetra_vec_pt =
      new Epetra_Vector(View, *epetra_map_pt, v_pt);

    // clean up
    delete epetra_map_pt;

    // return
    return epetra_vec_pt;
  }

  //=============================================================================
  ///  Helper function to copy the contents of a Trilinos vector to an
  /// oomph-lib distributed vector. The distribution of the two vectors must
  /// be identical
  //=============================================================================
  void TrilinosEpetraHelpers::copy_to_oomphlib_vector(
    const Epetra_Vector* epetra_vec_pt, DoubleVector& oomph_vec)
  {
#ifdef PARANOID
    // check the the oomph lib vector is setup
    if (!oomph_vec.built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib vector (oomph_v) must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // if the oomph-lib vector is distributed
    if (oomph_vec.distributed())
    {
      // extract values from epetra_v_pt
      double* v_values;
      epetra_vec_pt->ExtractView(&v_values);

      // copy the values
      unsigned nrow_local = oomph_vec.nrow_local();
      for (unsigned i = 0; i < nrow_local; i++)
      {
        oomph_vec[i] = v_values[i];
      }
    }

    // else teh oomph-lib vector is not distributed
    else
    {
      // get the values vector
#ifdef OOMPH_HAS_MPI
      int nproc = epetra_vec_pt->Map().Comm().NumProc();
      if (nproc == 1)
      {
        epetra_vec_pt->ExtractCopy(oomph_vec.values_pt());
      }
      else
      {
        // get the local values
        double* local_values;
        epetra_vec_pt->ExtractView(&local_values);

        // my rank
        int my_rank = epetra_vec_pt->Map().Comm().MyPID();

        // number of local rows
        Vector<int> nrow_local(nproc);
        nrow_local[my_rank] = epetra_vec_pt->MyLength();

        // gather the First_row vector
        int my_nrow_local_copy = nrow_local[my_rank];
        MPI_Allgather(
          &my_nrow_local_copy,
          1,
          MPI_INT,
          &nrow_local[0],
          1,
          MPI_INT,
          oomph_vec.distribution_pt()->communicator_pt()->mpi_comm());

        // number of local rows
        Vector<int> first_row(nproc);
        first_row[my_rank] = epetra_vec_pt->Map().MyGlobalElements()[0];

        // gather the First_row vector
        int my_first_row = first_row[my_rank];
        MPI_Allgather(
          &my_first_row,
          1,
          MPI_INT,
          &first_row[0],
          1,
          MPI_INT,
          oomph_vec.distribution_pt()->communicator_pt()->mpi_comm());

        // gather the local solution values
        MPI_Allgatherv(
          local_values,
          nrow_local[my_rank],
          MPI_DOUBLE,
          oomph_vec.values_pt(),
          &nrow_local[0],
          &first_row[0],
          MPI_DOUBLE,
          oomph_vec.distribution_pt()->communicator_pt()->mpi_comm());
      }
#else
      epetra_vec_pt->ExtractCopy(oomph_vec.values_pt());
#endif
    }
  }

  // MATRIX METHODS
  // =============================================================


  //=============================================================================
  ///  create an Epetra_CrsMatrix from an oomph-lib CRDoubleMatrix.
  /// If oomph_matrix_pt is NOT distributed (i.e. locally replicated) and
  /// on more than one processor, then the returned Epetra_Vector will be
  /// uniformly distributed. If the oomph_matrix_pt is distributed then the
  /// Epetra_CrsMatrix returned will have the same distribution as
  /// oomph_matrix_pt.
  /// The LinearAlgebraDistribution argument dist_pt should specify the
  /// distribution of the object this matrix will operate on.
  //=============================================================================
  Epetra_CrsMatrix* TrilinosEpetraHelpers::create_distributed_epetra_matrix(
    const CRDoubleMatrix* oomph_matrix_pt,
    const LinearAlgebraDistribution* dist_pt)
  {
#ifdef PARANOID
    if (!oomph_matrix_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib matrix must be built.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!oomph_matrix_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib matrix must be built.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!oomph_matrix_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib matrix must be built.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // get pointers to the matrix values, column indices etc
    // const_cast is safe because we use the Epetra_Vector "Copy" construction
    // method
    int* column = const_cast<int*>(oomph_matrix_pt->column_index());
    double* value = const_cast<double*>(oomph_matrix_pt->value());
    int* row_start = const_cast<int*>(oomph_matrix_pt->row_start());

    // create the corresponding Epetra_Map
    LinearAlgebraDistribution* target_dist_pt = 0;
    if (oomph_matrix_pt->distributed())
    {
      target_dist_pt =
        new LinearAlgebraDistribution(oomph_matrix_pt->distribution_pt());
    }
    else
    {
      target_dist_pt = new LinearAlgebraDistribution(
        oomph_matrix_pt->distribution_pt()->communicator_pt(),
        oomph_matrix_pt->nrow(),
        true);
    }
    Epetra_Map* epetra_map_pt = create_epetra_map(target_dist_pt);

    // first first coefficient of the oomph vector to be inserted into the
    // Epetra_Vector
    unsigned offset = 0;
    if (!oomph_matrix_pt->distributed())
    {
      offset = target_dist_pt->first_row();
    }

    // get my nrow_local and first_row
    unsigned nrow_local = target_dist_pt->nrow_local();
    unsigned first_row = target_dist_pt->first_row();

    // store the number of non zero entries per row
    int* nnz_per_row = new int[nrow_local];
    for (unsigned row = 0; row < nrow_local; row++)
    {
      nnz_per_row[row] = row_start[row + offset + 1] - row_start[offset + row];
    }

    // create the matrix
    Epetra_CrsMatrix* epetra_matrix_pt =
      new Epetra_CrsMatrix(Copy, *epetra_map_pt, nnz_per_row, true);

    // insert the values
    for (unsigned row = 0; row < nrow_local; row++)
    {
      // get pointer to this row in values/columns
      int ptr = row_start[row + offset];
#ifdef PARANOID
      int err = 0;
      err =
#endif
        epetra_matrix_pt->InsertGlobalValues(
          first_row + row, nnz_per_row[row], value + ptr, column + ptr);
#ifdef PARANOID
      if (err != 0)
      {
        std::ostringstream error_message;
        error_message
          << "Epetra Matrix Insert Global Values : epetra_error_flag = " << err;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // complete the build of the trilinos matrix
    LinearAlgebraDistribution* target_col_dist_pt = 0;
    if (dist_pt->distributed())
    {
      target_col_dist_pt = new LinearAlgebraDistribution(dist_pt);
    }
    else
    {
      target_col_dist_pt = new LinearAlgebraDistribution(
        dist_pt->communicator_pt(), dist_pt->nrow(), true);
    }
    Epetra_Map* epetra_domain_map_pt = create_epetra_map(target_col_dist_pt);
#ifdef PARANOID
    int err = 0;
    err =
#endif
      epetra_matrix_pt->FillComplete(*epetra_domain_map_pt, *epetra_map_pt);
#ifdef PARANOID
    if (err != 0)
    {
      std::ostringstream error_message;
      error_message
        << "Epetra Matrix Fill Complete Error : epetra_error_flag = " << err;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // tidy up memory
    delete[] nnz_per_row;
    delete epetra_map_pt;
    delete epetra_domain_map_pt;
    delete target_dist_pt;
    delete target_col_dist_pt;

    // return
    return epetra_matrix_pt;
  }


  //=============================================================================
  /// Class to allow sorting of column indices in conversion to epetra matrix
  //=============================================================================
  class DistributionPredicate
  {
  public:
    /// Constructor: Pass number of first column and the number of local columns
    DistributionPredicate(const int& first_col, const int& ncol_local)
      : First_col(first_col), Last_col(first_col + ncol_local - 1)
    {
    }

    ///  Comparison operator: is column col in the range
    /// between (including) First_col and Last_col
    bool operator()(const int& col)
    {
      if (col >= First_col && col <= Last_col)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    /// First column held locally
    int First_col;

    /// Last colum held locally
    int Last_col;
  };


  //=============================================================================
  ///  create and Epetra_CrsMatrix from an oomph-lib CRDoubleMatrix.
  /// Specialisation for Trilinos AztecOO.
  /// If oomph_matrix_pt is NOT distributed (i.e. locally replicated) and
  /// on more than one processor, then the returned Epetra_Vector will be
  /// uniformly distributed. If the oomph_matrix_pt is distributed then the
  /// Epetra_CrsMatrix returned will have the same distribution as
  /// oomph_matrix_pt.
  /// For AztecOO, the column map is ordered such that the local rows are
  /// first.
  //=============================================================================
  Epetra_CrsMatrix* TrilinosEpetraHelpers::
    create_distributed_epetra_matrix_for_aztecoo(
      CRDoubleMatrix* oomph_matrix_pt)
  {
#ifdef PARANOID
    if (!oomph_matrix_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The oomph-lib matrix must be built.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // get pointers to the matrix values, column indices etc
    // const_cast is safe because we use the Epetra_Vector "Copy" construction
    // method
    int* column = const_cast<int*>(oomph_matrix_pt->column_index());
    double* value = const_cast<double*>(oomph_matrix_pt->value());
    int* row_start = const_cast<int*>(oomph_matrix_pt->row_start());

    // create the corresponding Epetra_Map
    LinearAlgebraDistribution* target_dist_pt = 0;
    if (oomph_matrix_pt->distributed())
    {
      target_dist_pt =
        new LinearAlgebraDistribution(oomph_matrix_pt->distribution_pt());
    }
    else
    {
      target_dist_pt = new LinearAlgebraDistribution(
        oomph_matrix_pt->distribution_pt()->communicator_pt(),
        oomph_matrix_pt->nrow(),
        true);
    }
    Epetra_Map* epetra_map_pt = create_epetra_map(target_dist_pt);

    // create the epetra column map

#ifdef OOMPH_HAS_MPI
    int first_col = oomph_matrix_pt->first_row();
    int ncol_local = oomph_matrix_pt->nrow_local();

    // Build colum map
    Epetra_Map* epetra_col_map_pt = 0;
    {
      // Vector of column indices; on processor goes first
      std::vector<int> col_index_vector;
      col_index_vector.reserve(oomph_matrix_pt->nnz() + ncol_local);
      col_index_vector.resize(ncol_local);

      // Global column indices corresponding to on-processor rows
      for (int c = 0; c < ncol_local; ++c)
      {
        col_index_vector[c] = c + first_col;
      }

      // Remember where the on-processor rows (columns) end
      std::vector<int>::iterator mid = col_index_vector.end();

      // Now insert ALL column indices of ALL entries
      col_index_vector.insert(mid, column, column + oomph_matrix_pt->nnz());

      // Loop over the newly added entries and remove them if they
      // refer to on-processor columns
      std::vector<int>::iterator end =
        std::remove_if(mid,
                       col_index_vector.end(),
                       DistributionPredicate(first_col, ncol_local));

      // Now sort the newly added entries
      std::sort(mid, end);

      //...and remove duplicates
      end = std::unique(mid, end);

      // Make the map
      epetra_col_map_pt = new Epetra_Map(
        -1,
        end - col_index_vector.begin(),
        &col_index_vector[0],
        0,
        Epetra_MpiComm(
          oomph_matrix_pt->distribution_pt()->communicator_pt()->mpi_comm()));

      // Hack to clear memory
      std::vector<int>().swap(col_index_vector);
    }

#else

    int ncol = oomph_matrix_pt->ncol();
    Epetra_Map* epetra_col_map_pt =
      new Epetra_LocalMap(ncol, 0, Epetra_SerialComm());

#endif

    // first first coefficient of the oomph vector to be inserted into the
    // Epetra_Vector
    unsigned offset = 0;
    if (!oomph_matrix_pt->distributed())
    {
      offset = target_dist_pt->first_row();
    }

    // get my nrow_local and first_row
    unsigned nrow_local = target_dist_pt->nrow_local();
    unsigned first_row = target_dist_pt->first_row();

    // store the number of non zero entries per row
    int* nnz_per_row = new int[nrow_local];
    for (unsigned row = 0; row < nrow_local; ++row)
    {
      nnz_per_row[row] = row_start[row + offset + 1] - row_start[offset + row];
    }

    // create the matrix
    Epetra_CrsMatrix* epetra_matrix_pt = new Epetra_CrsMatrix(
      Copy, *epetra_map_pt, *epetra_col_map_pt, nnz_per_row, true);

    // insert the values
    for (unsigned row = 0; row < nrow_local; row++)
    {
      // get pointer to this row in values/columns
      int ptr = row_start[row + offset];
#ifdef PARANOID
      int err = 0;
      err =
#endif
        epetra_matrix_pt->InsertGlobalValues(
          first_row + row, nnz_per_row[row], value + ptr, column + ptr);
#ifdef PARANOID
      if (err != 0)
      {
        std::ostringstream error_message;
        error_message
          << "Epetra Matrix Insert Global Values : epetra_error_flag = " << err;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // complete the build of the trilinos matrix
#ifdef PARANOID
    int err = 0;
    err =
#endif
      epetra_matrix_pt->FillComplete();

#ifdef PARANOID
    if (err != 0)
    {
      std::ostringstream error_message;
      error_message
        << "Epetra Matrix Fill Complete Error : epetra_error_flag = " << err;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // tidy up memory
    delete[] nnz_per_row;
    delete epetra_map_pt;
    delete epetra_col_map_pt;
    delete target_dist_pt;

    // return
    return epetra_matrix_pt;
  }


  // MATRIX OPERATION METHODS ==================================================


  //============================================================================
  ///   Function to perform a matrix-vector multiplication on a
  /// oomph-lib matrix and vector using Trilinos functionality.
  /// NOTE 1. the matrix and the vectors must have the same communicator.
  /// NOTE 2. The vector will be returned with the same distribution
  /// as the matrix, unless a distribution is predefined in the solution
  /// vector in which case the vector will be returned with that distribution.
  //============================================================================
  void TrilinosEpetraHelpers::multiply(const CRDoubleMatrix* oomph_matrix_pt,
                                       const DoubleVector& oomph_x,
                                       DoubleVector& oomph_y)
  {
#ifdef PARANOID
    // check that this matrix is built
    if (!oomph_matrix_pt->built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // check that the distribution of the matrix and the soln are the same
    if (oomph_y.built())
    {
      if (!(*oomph_matrix_pt->distribution_pt() == *oomph_y.distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The soln vector and this matrix must have the same distribution.";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    // check that the distribution of the oomph-lib vector x is setup
    if (!oomph_x.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "The x vector must be setup";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // setup the distribution
    if (!oomph_y.distribution_pt()->built())
    {
      oomph_y.build(oomph_matrix_pt->distribution_pt(), 0.0);
    }

    // convert matrix1 to epetra matrix
    Epetra_CrsMatrix* epetra_matrix_pt = create_distributed_epetra_matrix(
      oomph_matrix_pt, oomph_x.distribution_pt());

    // convert x to Trilinos vector
    Epetra_Vector* epetra_x_pt = create_distributed_epetra_vector(oomph_x);

    // create Trilinos vector for soln ( 'viewing' the contents of the oomph-lib
    // matrix)
    Epetra_Vector* epetra_y_pt = create_distributed_epetra_vector(oomph_y);

    // do the multiply
#ifdef PARANOID
    int epetra_error_flag = 0;
    epetra_error_flag =
#endif
      epetra_matrix_pt->Multiply(false, *epetra_x_pt, *epetra_y_pt);

    // return the solution
    copy_to_oomphlib_vector(epetra_y_pt, oomph_y);

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

    // clean up
    delete epetra_matrix_pt;
    delete epetra_x_pt;
    delete epetra_y_pt;
  }

  //=============================================================================
  ///  Function to perform a matrix-matrix multiplication on oomph-lib
  /// matrices by using Trilinos functionality.
  /// \b NOTE 1. There are two Trilinos matrix-matrix multiplication methods
  /// available, using either the EpetraExt::MatrixMatrix class (if use_ml ==
  /// false) or using ML (Epetra_MatrixMult method)
  /// \b NOTE 2. the solution matrix (matrix_soln) will be returned with the
  /// same distribution as matrix1
  /// \b NOTE 3. All matrices must share the same communicator.
  //=============================================================================
  void TrilinosEpetraHelpers::multiply(const CRDoubleMatrix& matrix_1,
                                       const CRDoubleMatrix& matrix_2,
                                       CRDoubleMatrix& matrix_soln,
                                       const bool& use_ml)
  {
#ifdef PARANOID
    // check that matrix 1 is built
    if (!matrix_1.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix matrix_1 has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check that matrix 2 is built
    if (!matrix_2.built())
    {
      std::ostringstream error_message_stream;
      error_message_stream << "This matrix matrix_2 has not been built";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // check matrix dimensions are compatable
    if (matrix_1.ncol() != matrix_2.nrow())
    {
      std::ostringstream error_message;
      error_message
        << "Matrix dimensions incompatible for matrix-matrix multiplication"
        << "ncol() for first matrix: " << matrix_1.ncol()
        << "nrow() for second matrix: " << matrix_2.nrow();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // check that the have the same communicator
    OomphCommunicator temp_comm(matrix_1.distribution_pt()->communicator_pt());
    if (temp_comm != *matrix_2.distribution_pt()->communicator_pt())
    {
      std::ostringstream error_message;
      error_message
        << "Matrix dimensions incompatible for matrix-matrix multiplication"
        << "ncol() for first matrix: " << matrix_1.ncol()
        << "nrow() for second matrix: " << matrix_2.nrow();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // check that the distribution of the matrix and the soln are the same
    if (matrix_soln.distribution_pt()->built())
    {
      if (!(*matrix_soln.distribution_pt() == *matrix_1.distribution_pt()))
      {
        std::ostringstream error_message_stream;
        error_message_stream << "The solution matrix and matrix_1 must have "
                                "the same distribution.";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // setup the distribution
    if (!matrix_soln.distribution_pt()->built())
    {
      matrix_soln.build(matrix_1.distribution_pt());
    }

    // temporary fix
    // ML MM method only appears to work for square matrices
    // Should be investigated further.
    bool temp_use_ml = false;
    if ((*matrix_1.distribution_pt() == *matrix_2.distribution_pt()) &&
        (matrix_1.ncol() == matrix_2.ncol()))
    {
      temp_use_ml = use_ml;
    }

    // create matrix 1
    Epetra_CrsMatrix* epetra_matrix_1_pt =
      create_distributed_epetra_matrix(&matrix_1, matrix_2.distribution_pt());

    // create matrix 2
    LinearAlgebraDistribution matrix_2_column_dist(
      matrix_2.distribution_pt()->communicator_pt(), matrix_2.ncol(), true);
    Epetra_CrsMatrix* epetra_matrix_2_pt =
      create_distributed_epetra_matrix(&matrix_2, &matrix_2_column_dist);

    // create the Trilinos epetra matrix to hold solution - will have same map
    // (and number of rows) as matrix_1
    Epetra_CrsMatrix* solution_pt;

    // do the multiplication
    // ---------------------
    if (temp_use_ml)
    {
      // there is a problem using this function, many pages of
      // warning messages are issued....
      // "tmpresult->InsertGlobalValues returned 3"
      // and
      // "Result_epet->InsertGlobalValues returned 3"
      // from function ML_back_to_epetraCrs(...) in
      // file ml_epetra_utils.cpp unless the
      // relevant lines are commented out. However this function
      // is much faster (at least on Biowulf) than the alternative
      // below.
      solution_pt = Epetra_MatrixMult(epetra_matrix_1_pt, epetra_matrix_2_pt);
    }
    else
    {
      // this method requires us to pass in the solution matrix
      solution_pt = new Epetra_CrsMatrix(Copy, epetra_matrix_1_pt->RowMap(), 0);
#ifdef PARANOID
      int epetra_error_flag = 0;
      epetra_error_flag =
#endif
        EpetraExt::MatrixMatrix::Multiply(
          *epetra_matrix_1_pt, false, *epetra_matrix_2_pt, false, *solution_pt);
#ifdef PARANOID
      if (epetra_error_flag != 0)
      {
        std::ostringstream error_message;
        error_message << "error flag from Multiply(): " << epetra_error_flag
                      << " from TrilinosHelpers::multiply" << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    // extract values and put into solution
    // ------------------------------------

    // find
    int nnz_local = solution_pt->NumMyNonzeros();
    int nrow_local = matrix_1.nrow_local();

    // do some checks
#ifdef PARANOID
    // check number of global rows in soluton matches that in matrix_1
    if ((int)matrix_1.nrow() != solution_pt->NumGlobalRows())
    {
      std::ostringstream error_message;
      error_message << "Incorrect number of global rows in solution matrix. "
                    << "nrow() for first input matrix: " << matrix_1.nrow()
                    << " nrow() for solution: " << solution_pt->NumGlobalRows();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // check number of local rows in soluton matches that in matrix_1
    if (static_cast<int>(matrix_1.nrow_local()) != solution_pt->NumMyRows())
    {
      std::ostringstream error_message;
      error_message << "Incorrect number of local rows in solution matrix. "
                    << "nrow_local() for first input matrix: "
                    << matrix_1.nrow_local() << " nrow_local() for solution: "
                    << solution_pt->NumMyRows();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // check number of global columns in soluton matches that in matrix_2
    if ((int)matrix_2.ncol() != solution_pt->NumGlobalCols())
    {
      std::ostringstream error_message;
      error_message << "Incorrect number of global columns in solution matrix. "
                    << "ncol() for second input matrix: " << matrix_2.ncol()
                    << " ncol() for solution: " << solution_pt->NumGlobalCols();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // check global index of the first row matches
    if (static_cast<int>(matrix_1.first_row()) != solution_pt->GRID(0))
    {
      std::ostringstream error_message;
      error_message
        << "Incorrect global index for first row of solution matrix. "
        << "first_row() for first input matrix : " << matrix_1.first_row()
        << " first_row() for solution: " << solution_pt->GRID(0);
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // extract values from Epetra matrix row by row
    double* value = new double[nnz_local];
    int* column_index = new int[nnz_local];
    int* row_start = new int[nrow_local + 1];
    int ptr = 0;
    int num_entries = 0;
    int first = matrix_soln.first_row();
    int last = first + matrix_soln.nrow_local();
    for (int row = first; row < last; row++)
    {
      row_start[row - first] = ptr;
      solution_pt->ExtractGlobalRowCopy(
        row, nnz_local, num_entries, value + ptr, column_index + ptr);
      ptr += num_entries;
    }
    row_start[nrow_local] = ptr;

    // delete Trilinos objects
    delete epetra_matrix_1_pt;
    delete epetra_matrix_2_pt;
    delete solution_pt;

    // Build the Oomph-lib solution matrix using build function
    matrix_soln.build(matrix_1.distribution_pt());
    matrix_soln.build_without_copy(
      matrix_2.ncol(), nnz_local, value, column_index, row_start);
  }


  // HELPER METHODS
  // =============================================================


  //=============================================================================
  /// create an Epetra_Map corresponding to the LinearAlgebraDistribution
  //=============================================================================
  Epetra_Map* TrilinosEpetraHelpers::create_epetra_map(
    const LinearAlgebraDistribution* const dist_pt)
  {
#ifdef OOMPH_HAS_MPI
    if (dist_pt->distributed())
    {
      unsigned first_row = dist_pt->first_row();
      unsigned nrow_local = dist_pt->nrow_local();
      int* my_global_rows = new int[nrow_local];
      for (unsigned i = 0; i < nrow_local; ++i)
      {
        my_global_rows[i] = first_row + i;
      }
      Epetra_Map* epetra_map_pt =
        new Epetra_Map(dist_pt->nrow(),
                       nrow_local,
                       my_global_rows,
                       0,
                       Epetra_MpiComm(dist_pt->communicator_pt()->mpi_comm()));
      delete[] my_global_rows;
      return epetra_map_pt;
    }
    else
    {
      return new Epetra_LocalMap(
        int(dist_pt->nrow()),
        int(0),
        Epetra_MpiComm(dist_pt->communicator_pt()->mpi_comm()));
    }
#else
    return new Epetra_LocalMap(
      int(dist_pt->nrow()), int(0), Epetra_SerialComm());
#endif
  }

} // namespace oomph
