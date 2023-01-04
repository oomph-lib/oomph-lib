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
// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib includes
#include "general_purpose_preconditioners.h"


namespace oomph
{
  //=============================================================
  /// Setup diagonal preconditioner: Store the inverse of the
  /// diagonal entries from the fully
  /// assembled matrix.
  //=============================================================
  void MatrixBasedDiagPreconditioner::setup()
  {
    // first attempt to cast to DistributableLinearAlgebraObject
    DistributableLinearAlgebraObject* dist_matrix_pt =
      dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt());

    // if it is a distributable matrix
    if (dist_matrix_pt != 0)
    {
      // cache the number of first_rows and nrow_local
      unsigned nrow_local = dist_matrix_pt->nrow_local();
      unsigned first_row = dist_matrix_pt->first_row();

      // resize the inverse diagonal storage
      Inv_diag.resize(nrow_local);

      // Extract the diagonal entries
      for (unsigned i = 0; i < nrow_local; i++)
      {
        unsigned index = i + first_row;
#ifdef PARANOID
        if ((*matrix_pt())(i, index) == 0.0)
        {
          throw OomphLibError(
            "Zero diagonal in matrix --> Cannot use diagonal preconditioner.",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
#endif
        Inv_diag[i] = 1.0 / (*matrix_pt())(i, index);
      }

      // store the distribution
      this->build_distribution(dist_matrix_pt->distribution_pt());
    }

    // else it is not a distributable matrix
    else
    {
      // # of rows in the matrix
      unsigned n_row = matrix_pt()->nrow();

      // Resize the Inv_diag vector to accommodate the # of
      // diagonal entries
      Inv_diag.resize(n_row);

      // Extract the diagonal entries
      for (unsigned i = 0; i < n_row; i++)
      {
#ifdef PARANOID
        if ((*matrix_pt())(i, i) == 0.0)
        {
          throw OomphLibError(
            "Zero diagonal in matrix --> Cannot use diagonal preconditioner.",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
        }
        else
#endif
        {
          Inv_diag[i] = 1.0 / (*matrix_pt())(i, i);
        }
      }

      // create the distribution
      LinearAlgebraDistribution dist(comm_pt(), n_row, false);
      this->build_distribution(dist);
    }
  }


  //=============================================================================
  /// Apply preconditioner: Multiply r by the inverse of the diagonal.
  //=============================================================================
  void MatrixBasedDiagPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
#ifdef PARANOID
    if (*r.distribution_pt() != *this->distribution_pt())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The r vector must have the same distribution as the "
           "preconditioner. "
        << "(this is the same as the matrix passed to setup())";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (z.built())
    {
      if (*z.distribution_pt() != *this->distribution_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The z vector distribution has been setup; it must have the "
          << "same distribution as the r vector (and preconditioner).";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // if z has not been setup then rebuild it
    if (!z.built())
    {
      z.build(this->distribution_pt(), 0.0);
    }

    // apply the preconditioner
    const double* r_values = r.values_pt();
    double* z_values = z.values_pt();
    unsigned nrow_local = this->nrow_local();
    for (unsigned i = 0; i < nrow_local; i++)
    {
      z_values[i] = Inv_diag[i] * r_values[i];
    }
  }

  //=============================================================================
  /// Setup the lumped preconditioner. Specialisation for CCDoubleMatrix.
  //=============================================================================
  template<>
  void MatrixBasedLumpedPreconditioner<CCDoubleMatrix>::setup()
  {
    // # of rows in the matrix
    Nrow = matrix_pt()->nrow();

    // Create the vector for the inverse lumped
    if (Inv_lumped_diag_pt != 0)
    {
      delete[] this->Inv_lumped_diag_pt;
    }
    Inv_lumped_diag_pt = new double[this->Nrow];

    // zero the vector
    for (unsigned i = 0; i < Nrow; i++)
    {
      Inv_lumped_diag_pt[i] = 0.0;
    }

    // cast the Double Base Matrix to Compressed Column Double Matrix
    CCDoubleMatrix* cc_matrix_pt = dynamic_cast<CCDoubleMatrix*>(matrix_pt());

    // get the matrix
    int* m_row_index = cc_matrix_pt->row_index();
    double* m_value = cc_matrix_pt->value();
    unsigned m_nnz = cc_matrix_pt->nnz();

    // intially set positive matrix to true
    Positive_matrix = true;

    // lump the matrix
    for (unsigned i = 0; i < m_nnz; i++)
    {
      // if the matrix contains negative coefficient the matrix not positive
      if (m_value[i] < 0.0)
      {
        Positive_matrix = false;
      }

      // computed lumped matrix - temporarily stored in Inv_lumped_diag_pt
      Inv_lumped_diag_pt[m_row_index[i]] += m_value[i];
    }

    // invert the lumped matrix
    for (unsigned i = 0; i < Nrow; i++)
    {
      Inv_lumped_diag_pt[i] = 1.0 / Inv_lumped_diag_pt[i];
    }

    // create the distribution
    LinearAlgebraDistribution dist(comm_pt(), Nrow, false);
    this->build_distribution(dist);
  }

  //=============================================================================
  /// Setup the lumped preconditioner. Specialisation for CRDoubleMatrix.
  //=============================================================================
  template<>
  void MatrixBasedLumpedPreconditioner<CRDoubleMatrix>::setup()
  {
    // first attempt to cast to CRDoubleMatrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

    // # of rows in the matrix
    Nrow = cr_matrix_pt->nrow_local();

    // Create the vector for the inverse lumped
    if (Inv_lumped_diag_pt != 0)
    {
      delete[] this->Inv_lumped_diag_pt;
    }
    Inv_lumped_diag_pt = new double[this->Nrow];

    // zero the vector
    for (unsigned i = 0; i < Nrow; i++)
    {
      Inv_lumped_diag_pt[i] = 0.0;
    }

    // get the matrix
    int* m_row_start = cr_matrix_pt->row_start();
    double* m_value = cr_matrix_pt->value();

    // intially set positive matrix to true
    Positive_matrix = true;

    // lump and invert matrix
    for (unsigned i = 0; i < Nrow; i++)
    {
      Inv_lumped_diag_pt[i] = 0.0;
      for (int j = m_row_start[i]; j < m_row_start[i + 1]; j++)
      {
        // if the matrix contains negative coefficient the matrix not positive
        if (m_value[j] < 0.0)
        {
          Positive_matrix = false;
        }

        // computed lumped coef (i,i)- temporarily stored in Inv_lumped_diag_pt
        Inv_lumped_diag_pt[i] += m_value[j];
      }
      // invert coef (i,i)
      Inv_lumped_diag_pt[i] = 1.0 / Inv_lumped_diag_pt[i];
    }

    // store the distribution
    this->build_distribution(cr_matrix_pt->distribution_pt());
  }


  //=============================================================================
  /// Apply preconditioner: Multiply r by the inverse of the lumped matrix
  //=============================================================================
  template<typename MATRIX>
  void MatrixBasedLumpedPreconditioner<MATRIX>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
#ifdef PARANOID
    if (*r.distribution_pt() != *this->distribution_pt())
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "The r vector must have teh same distribution as the "
           "preconditioner. "
        << "(this is the same as the matrix passed to setup())";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    if (z.built())
    {
      if (*z.distribution_pt() != *this->distribution_pt())
      {
        std::ostringstream error_message_stream;
        error_message_stream
          << "The z vector distribution has been setup; it must have the "
          << "same distribution as the r vector (and preconditioner).";
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    z.build(r.distribution_pt(), 0.0);
    for (unsigned i = 0; i < Nrow; i++)
    {
      z[i] = Inv_lumped_diag_pt[i] * r[i];
    }
  }


  // ensure the lumped preconditioner get built
  template class MatrixBasedLumpedPreconditioner<CCDoubleMatrix>;
  template class MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;


  //=============================================================================
  /// setup ILU(0) preconditioner for Matrices of CCDoubleMatrix type
  //=============================================================================
  void ILUZeroPreconditioner<CCDoubleMatrix>::setup()
  {
    // cast the Double Base Matrix to Compressed Column Double Matrix
    CCDoubleMatrix* cc_matrix_pt = dynamic_cast<CCDoubleMatrix*>(matrix_pt());

#ifdef PARANOID
    if (cc_matrix_pt == 0)
    {
      std::ostringstream error_msg;
      error_msg << "Failed to conver matrix_pt to CCDoubleMatrix*.";
      throw OomphLibError(
        error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // number of rows in matrix
    int n_row = cc_matrix_pt->nrow();

    // set the distribution
    LinearAlgebraDistribution dist(comm_pt(), n_row, false);
    this->build_distribution(dist);

    // declares variables to store number of non zero entires in L and U
    int l_nz = 0;
    int u_nz = 0;

    // create space for m matrix
    int* m_column_start;
    int* m_row_index;
    double* m_value;

    // get the m matrix
    m_column_start = cc_matrix_pt->column_start();
    m_row_index = cc_matrix_pt->row_index();
    m_value = cc_matrix_pt->value();

    // find number non zero entries in L and U
    for (int i = 0; i < n_row; i++)
    {
      for (int j = m_column_start[i]; j < m_column_start[i + 1]; j++)
      {
        if (m_row_index[j] > i)
        {
          l_nz++;
        }
        else
        {
          u_nz++;
        }
      }
    }

    // resize vectors to store the data for the lower prior to building the
    // matrices
    L_column_start.resize(n_row + 1);
    L_row_entry.resize(l_nz);

    // and the upper matrix
    U_column_start.resize(n_row + 1);
    U_row_entry.resize(u_nz);

    // set first column pointers to zero
    L_column_start[0] = 0;
    U_column_start[0] = 0;

    // split the matrix into L and U
    for (int i = 0; i < n_row; i++)
    {
      L_column_start[i + 1] = L_column_start[i];
      U_column_start[i + 1] = U_column_start[i];
      for (int j = m_column_start[i]; j < m_column_start[i + 1]; j++)
      {
        if (m_row_index[j] > i)
        {
          int k = L_column_start[i + 1]++;
          L_row_entry[k].index() = m_row_index[j];
          L_row_entry[k].value() = m_value[j];
        }
        else
        {
          int k = U_column_start[i + 1]++;
          U_row_entry[k].index() = m_row_index[j];
          U_row_entry[k].value() = m_value[j];
        }
      }
    }

    // sort each row entry vector into row index order for each column

    // loop over the columns
    for (unsigned i = 0; i < unsigned(n_row); i++)
    {
      // sort the columns of the L matrix
      std::sort(L_row_entry.begin() + L_column_start[i],
                L_row_entry.begin() + L_column_start[i + 1]);

      // sort the columns of the U matrix
      std::sort(U_row_entry.begin() + U_column_start[i],
                U_row_entry.begin() + U_column_start[i + 1]);
    }


    // factorise matrix
    int i;
    unsigned j, pn, qn, rn;
    pn = 0;
    qn = 0;
    rn = 0;
    double multiplier;
    for (i = 0; i < n_row - 1; i++)
    {
      multiplier = U_row_entry[U_column_start[i + 1] - 1].value();
      for (j = L_column_start[i]; j < L_column_start[i + 1]; j++)
        L_row_entry[j].value() /= multiplier;
      for (j = U_column_start[i + 1]; j < U_column_start[i + 2] - 1; j++)
      {
        multiplier = U_row_entry[j].value();
        qn = j + 1;
        rn = L_column_start[i + 1];
        for (pn = L_column_start[U_row_entry[j].index()];
             pn < L_column_start[U_row_entry[j].index() + 1] &&
             static_cast<int>(L_row_entry[pn].index()) <= i + 1;
             pn++)
        {
          while (qn < U_column_start[i + 2] &&
                 U_row_entry[qn].index() < L_row_entry[pn].index())
            qn++;
          if (qn < U_column_start[i + 2] &&
              L_row_entry[pn].index() == U_row_entry[qn].index())
            U_row_entry[qn].value() -= multiplier * L_row_entry[pn].value();
        }
        for (; pn < L_column_start[U_row_entry[j].index() + 1]; pn++)
        {
          while (rn < L_column_start[i + 2] &&
                 L_row_entry[rn].index() < L_row_entry[pn].index())
            rn++;
          if (rn < L_column_start[i + 2] &&
              L_row_entry[pn].index() == L_row_entry[rn].index())
            L_row_entry[rn].value() -= multiplier * L_row_entry[pn].value();
        }
      }
    }
  }


  //=============================================================================
  /// setup ILU(0) preconditioner for Matrices of CRDoubleMatrix Type
  //=============================================================================
  void ILUZeroPreconditioner<CRDoubleMatrix>::setup()
  {
    // cast the Double Base Matrix to Compressed Column Double Matrix
    CRDoubleMatrix* cr_matrix_pt = dynamic_cast<CRDoubleMatrix*>(matrix_pt());

#ifdef PARANOID
    if (cr_matrix_pt == 0)
    {
      std::ostringstream error_msg;
      error_msg << "Failed to conver matrix_pt to CRDoubleMatrix*.";
      throw OomphLibError(
        error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // if the matrix is distributed then build global version
    bool built_global = false;
    if (cr_matrix_pt->distributed())
    {
      // get the global matrix
      CRDoubleMatrix* global_matrix_pt = cr_matrix_pt->global_matrix();

      // store at cr_matrix pointer
      cr_matrix_pt = global_matrix_pt;

      // set the flag so we can delete later
      built_global = true;
    }

    // store the Distribution
    this->build_distribution(cr_matrix_pt->distribution_pt());

    // number of rows in matrix
    int n_row = cr_matrix_pt->nrow();

    // declares variables to store number of non zero entires in L and U
    int l_nz = 0;
    int u_nz = 0;

    // create space for m matrix
    int* m_row_start;
    int* m_column_index;
    double* m_value;

    // get the m matrix
    m_row_start = cr_matrix_pt->row_start();
    m_column_index = cr_matrix_pt->column_index();
    m_value = cr_matrix_pt->value();

    // find number non zero entries in L and U
    for (int i = 0; i < n_row; i++)
    {
      for (int j = m_row_start[i]; j < m_row_start[i + 1]; j++)
      {
        if (m_column_index[j] < i)
        {
          l_nz++;
        }
        else
        {
          u_nz++;
        }
      }
    }

    // resize vectors to store the data for the lower prior to building the
    // matrices
    L_row_start.resize(n_row + 1);
    L_row_entry.resize(l_nz);

    // and the upper matrix
    U_row_start.resize(n_row + 1);
    U_row_entry.resize(u_nz);

    // set first column pointers to zero
    L_row_start[0] = 0;
    U_row_start[0] = 0;

    // split the matrix into L and U
    for (int i = 0; i < n_row; i++)
    {
      L_row_start[i + 1] = L_row_start[i];
      U_row_start[i + 1] = U_row_start[i];
      for (int j = m_row_start[i]; j < m_row_start[i + 1]; j++)
      {
        if (m_column_index[j] < i)
        {
          int k = L_row_start[i + 1]++;
          L_row_entry[k].value() = m_value[j];
          L_row_entry[k].index() = m_column_index[j];
        }
        else
        {
          int k = U_row_start[i + 1]++;
          U_row_entry[k].value() = m_value[j];
          U_row_entry[k].index() = m_column_index[j];
        }
      }
    }


    // factorise matrix
    unsigned i, j, pn, qn, rn;
    pn = 0;
    qn = 0;
    rn = 0;
    double multiplier;
    for (i = 1; i < static_cast<unsigned>(n_row); i++)
    {
      for (j = L_row_start[i]; j < L_row_start[i + 1]; j++)
      {
        pn = U_row_start[L_row_entry[j].index()];
        multiplier = (L_row_entry[j].value() /= U_row_entry[pn].value());
        qn = j + 1;
        rn = U_row_start[i];
        for (pn++; pn < U_row_start[L_row_entry[j].index() + 1] &&
                   U_row_entry[pn].index() < i;
             pn++)
        {
          while (qn < L_row_start[i + 1] &&
                 L_row_entry[qn].index() < U_row_entry[pn].index())
            qn++;
          if (qn < L_row_start[i + 1] &&
              U_row_entry[pn].index() == L_row_entry[qn].index())
            L_row_entry[qn].value() -= multiplier * U_row_entry[pn].value();
        }
        for (; pn < U_row_start[L_row_entry[j].index() + 1]; pn++)
        {
          while (rn < U_row_start[i + 1] &&
                 U_row_entry[rn].index() < U_row_entry[pn].index())
            rn++;
          if (rn < U_row_start[i + 1] &&
              U_row_entry[pn].index() == U_row_entry[rn].index())
            U_row_entry[rn].value() -= multiplier * U_row_entry[pn].value();
        }
      }
    }

    // if we built the global matrix then delete it
    if (built_global)
    {
      delete cr_matrix_pt;
    }
  }


  //=============================================================================
  /// Apply ILU(0) preconditioner for CCDoubleMatrix: Solve Ly=r then
  /// Uz=y and return z
  //=============================================================================
  void ILUZeroPreconditioner<CCDoubleMatrix>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // # of rows in the matrix
    int n_row = r.nrow();

    // store the distribution of z
    LinearAlgebraDistribution* z_dist = 0;
    if (z.built())
    {
      z_dist = new LinearAlgebraDistribution(z.distribution_pt());
    }

    // copy r to z
    z = r;

    // if z is distributed then change to global
    if (z.distributed())
    {
      z.redistribute(this->distribution_pt());
    }

    // solve Ly=r (note L matrix is unit and diagonal is not stored)
    for (unsigned i = 0; i < static_cast<unsigned>(n_row); i++)
    {
      for (unsigned j = L_column_start[i]; j < L_column_start[i + 1]; j++)
      {
        z[L_row_entry[j].index()] =
          z[L_row_entry[j].index()] - z[i] * L_row_entry[j].value();
      }
    }

    // solve Uz=y
    double x;
    for (int i = n_row - 1; i >= 0; i--)
    {
      x = z[i] / U_row_entry[U_column_start[i + 1] - 1].value();
      z[i] = x;
      for (unsigned j = U_column_start[i]; j < U_column_start[i + 1] - 1; j++)
      {
        z[U_row_entry[j].index()] =
          z[U_row_entry[j].index()] - x * U_row_entry[j].value();
      }
    }

    // if the distribution of z was preset the redistribute to original
    if (z_dist != 0)
    {
      z.redistribute(z_dist);
      delete z_dist;
    }
  }

  //=============================================================================
  /// Apply ILU(0) preconditioner for CRDoubleMatrix: Solve Ly=r then
  /// Uz=y
  ///  and return z
  //=============================================================================
  void ILUZeroPreconditioner<CRDoubleMatrix>::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // # of rows in the matrix
    int n_row = r.nrow();

    // store the distribution of z
    LinearAlgebraDistribution* z_dist = 0;
    if (z.built())
    {
      z_dist = new LinearAlgebraDistribution(z.distribution_pt());
    }

    // copy r to z
    z = r;

    // if z is distributed then change to global
    if (z.distributed())
    {
      z.redistribute(this->distribution_pt());
    }

    // solve Ly=r (note L matrix is unit and diagonal is not stored)
    double t;
    for (int i = 0; i < n_row; i++)
    {
      t = 0;
      for (unsigned j = L_row_start[i]; j < L_row_start[i + 1]; j++)
      {
        t = t + L_row_entry[j].value() * z[L_row_entry[j].index()];
      }
      z[i] = z[i] - t;
    }

    // solve Uz=y
    for (int i = n_row - 1; i >= 0; i--)
    {
      t = 0;
      for (unsigned j = U_row_start[i] + 1; j < U_row_start[i + 1]; j++)
      {
        t = t + U_row_entry[j].value() * z[U_row_entry[j].index()];
      }
      z[i] = z[i] - t;
      z[i] = z[i] / U_row_entry[U_row_start[i]].value();
    }

    // if the distribution of z was preset the redistribute to original
    if (z_dist != 0)
    {
      z.redistribute(z_dist);
      delete z_dist;
    }
  }
} // namespace oomph
