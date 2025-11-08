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
#include "double_vector.h"
#include "matrices.h"


namespace oomph
{
  //============================================================================
  /// Just copys the argument DoubleVector
  //============================================================================
  void DoubleVector::build(const DoubleVector& old_vector)
  {
    if (!(*this == old_vector))
    {
      // the vector owns the internal data
      Internal_values = true;

      // reset the distribution and resize the data
      this->build(old_vector.distribution_pt(), 0.0);

      // copy the data
      if (this->distribution_built())
      {
        unsigned nrow_local = this->nrow_local();
        const double* old_vector_values = old_vector.values_pt();
        std::copy(old_vector_values, old_vector_values + nrow_local, Values_pt);
      }
    }
  }

  //============================================================================
  /// Assembles a DoubleVector with distribution dist, if v is specified
  /// each row is set to v
  //============================================================================
  void DoubleVector::build(const LinearAlgebraDistribution* const& dist_pt,
                           const double& v)
  {
    // clean the memory
    this->clear();

    // the vector owns the internal data
    Internal_values = true;

    // Set the distribution
    this->build_distribution(dist_pt);

    // update the values
    if (dist_pt->built())
    {
      unsigned nrow_local = this->nrow_local();
      Values_pt = new double[nrow_local];

      std::fill_n(Values_pt, nrow_local, v);
      Built = true;
    }
    else
    {
      Built = false;
    }
  }

  //============================================================================
  /// Assembles a DoubleVector with a distribution dist and coefficients
  /// taken from the vector v.
  /// Note. The vector v MUST be of length nrow()
  //============================================================================
  void DoubleVector::build(const LinearAlgebraDistribution* const& dist_pt,
                           const Vector<double>& v)
  {
    // clean the memory
    this->clear();

    // the vector owns the internal data
    Internal_values = true;

    // Set the distribution
    this->build_distribution(dist_pt);

    // update the values
    if (dist_pt->built())
    {
      // re-allocate memory which was deleted by clear()
      unsigned nrow_local = this->nrow_local();
      Values_pt = new double[nrow_local];

      // use the initialise method to populate the vector
      this->initialise(v);
      Built = true;
    }
    else
    {
      Built = false;
    }
  }

  //============================================================================
  /// initialise the whole vector with value v
  //============================================================================
  void DoubleVector::initialise(const double& v)
  {
    if (Built)
    {
      // cache nrow local
      unsigned nrow_local = this->nrow_local();

      std::fill_n(Values_pt, nrow_local, v);
    }
  }

  //============================================================================
  /// initialise the vector with coefficient from the vector v.
  /// Note: The vector v must be of length
  //============================================================================
  void DoubleVector::initialise(const Vector<double> v)
  {
#ifdef PARANOID
    if (v.size() != this->nrow())
    {
      std::ostringstream error_message;
      error_message << "The vector passed to initialise(...) must be of length "
                    << "nrow()";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    unsigned begin_first_row = this->first_row();
    unsigned end = begin_first_row + this->nrow_local();

    std::copy(v.begin() + begin_first_row, v.begin() + end, Values_pt);
  }

  //============================================================================
  /// The contents of the vector are redistributed to match the new
  /// distribution. In a non-MPI build this method works, but does nothing.
  /// \b NOTE 1: The current distribution and the new distribution must have
  /// the same number of global rows.
  /// \b NOTE 2: The current distribution and the new distribution must have
  /// the same Communicator.
  //============================================================================
  void DoubleVector::redistribute(
    const LinearAlgebraDistribution* const& dist_pt)
  {
#ifdef OOMPH_HAS_MPI
#ifdef PARANOID
    if (!Internal_values)
    {
      // if this vector does not own the double* values then it cannot be
      // distributed.
      // note: this is not stictly necessary - would just need to be careful
      // with delete[] below.
      std::ostringstream error_message;
      error_message << "This vector does not own its data (i.e. it has been "
                    << "passed in via set_external_values() and therefore "
                    << "cannot be redistributed";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the nrows for both distributions is the
    // same
    if (dist_pt->nrow() != this->nrow())
    {
      std::ostringstream error_message;
      error_message << "The number of global rows in the new distribution ("
                    << dist_pt->nrow() << ") is not equal to the number"
                    << " of global rows in the current distribution ("
                    << this->nrow() << ").\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // paranoid check that the current distribution and the new distribution
    // have the same Communicator
    OomphCommunicator temp_comm(*dist_pt->communicator_pt());
    if (!(temp_comm == *this->distribution_pt()->communicator_pt()))
    {
      std::ostringstream error_message;
      error_message << "The new distribution and the current distribution must "
                    << "have the same communicator.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // check the distributions are not the same
    if (!((*this->distribution_pt()) == *dist_pt))
    {
      // get the rank and the number of processors
      int my_rank = this->distribution_pt()->communicator_pt()->my_rank();
      int nproc = this->distribution_pt()->communicator_pt()->nproc();

      // if both vectors are distributed
      if (this->distributed() && dist_pt->distributed())
      {
        // new nrow_local and first_row data
        Vector<unsigned> new_first_row_data(nproc);
        Vector<unsigned> new_nrow_local_data(nproc);
        Vector<unsigned> current_first_row_data(nproc);
        Vector<unsigned> current_nrow_local_data(nproc);
        for (int i = 0; i < nproc; i++)
        {
          new_first_row_data[i] = dist_pt->first_row(i);
          new_nrow_local_data[i] = dist_pt->nrow_local(i);
          current_first_row_data[i] = this->first_row(i);
          current_nrow_local_data[i] = this->nrow_local(i);
        }

        // compute which local rows are expected to be received from each
        // processor / sent to each processor
        Vector<unsigned> new_first_row_for_proc(nproc);
        Vector<unsigned> new_nrow_local_for_proc(nproc);
        Vector<unsigned> new_first_row_from_proc(nproc);
        Vector<unsigned> new_nrow_local_from_proc(nproc);

        // for every processor compute first_row and nrow_local that will
        // will sent and received by this processor
        for (int p = 0; p < nproc; p++)
        {
          // start with data to be sent
          if ((new_first_row_data[p] < (current_first_row_data[my_rank] +
                                        current_nrow_local_data[my_rank])) &&
              (current_first_row_data[my_rank] <
               (new_first_row_data[p] + new_nrow_local_data[p])))
          {
            new_first_row_for_proc[p] =
              std::max(current_first_row_data[my_rank], new_first_row_data[p]);
            new_nrow_local_for_proc[p] =
              std::min((current_first_row_data[my_rank] +
                        current_nrow_local_data[my_rank]),
                       (new_first_row_data[p] + new_nrow_local_data[p])) -
              new_first_row_for_proc[p];
          }

          // and data to be received
          if ((new_first_row_data[my_rank] <
               (current_first_row_data[p] + current_nrow_local_data[p])) &&
              (current_first_row_data[p] <
               (new_first_row_data[my_rank] + new_nrow_local_data[my_rank])))
          {
            new_first_row_from_proc[p] =
              std::max(current_first_row_data[p], new_first_row_data[my_rank]);
            new_nrow_local_from_proc[p] =
              std::min(
                (current_first_row_data[p] + current_nrow_local_data[p]),
                (new_first_row_data[my_rank] + new_nrow_local_data[my_rank])) -
              new_first_row_from_proc[p];
          }
        }

        // temporary storage for the new data
        double* temp_data = new double[new_nrow_local_data[my_rank]];

        // "send to self" or copy Data that does not need to be sent else where
        // to temp_data
        if (new_nrow_local_for_proc[my_rank] != 0)
        {
          unsigned j =
            new_first_row_for_proc[my_rank] - current_first_row_data[my_rank];
          unsigned k =
            new_first_row_for_proc[my_rank] - new_first_row_data[my_rank];
          for (unsigned i = 0; i < new_nrow_local_for_proc[my_rank]; i++)
          {
            temp_data[k + i] = Values_pt[j + i];
          }
        }

        // send and receive circularly
        for (int p = 1; p < nproc; p++)
        {
          // next processor to send to
          unsigned dest_p = (my_rank + p) % nproc;

          // next processor to receive from
          unsigned source_p = (nproc + my_rank - p) % nproc;

          // send and receive the value
          MPI_Status status;
          MPI_Sendrecv(Values_pt + new_first_row_for_proc[dest_p] -
                         current_first_row_data[my_rank],
                       new_nrow_local_for_proc[dest_p],
                       MPI_DOUBLE,
                       dest_p,
                       1,
                       temp_data + new_first_row_from_proc[source_p] -
                         new_first_row_data[my_rank],
                       new_nrow_local_from_proc[source_p],
                       MPI_DOUBLE,
                       source_p,
                       1,
                       this->distribution_pt()->communicator_pt()->mpi_comm(),
                       &status);
        }

        // copy from temp data to Values_pt
        delete[] Values_pt;
        unsigned nrow_local = dist_pt->nrow_local();
        Values_pt = new double[nrow_local];
        for (unsigned i = 0; i < nrow_local; i++)
        {
          Values_pt[i] = temp_data[i];
        }
        delete[] temp_data;
      }

      // if this vector is distributed but the new distributed is global
      else if (this->distributed() && !dist_pt->distributed())
      {
        // copy existing Values_pt to temp_data
        unsigned nrow_local = this->nrow_local();
        double* temp_data = new double[nrow_local];
        for (unsigned i = 0; i < nrow_local; i++)
        {
          temp_data[i] = Values_pt[i];
        }

        // clear and resize Values_pt
        delete[] Values_pt;
        Values_pt = new double[this->nrow()];

        // create a int vector of first rows
        int* dist_first_row = new int[nproc];
        int* dist_nrow_local = new int[nproc];
        for (int p = 0; p < nproc; p++)
        {
          dist_first_row[p] = this->first_row(p);
          dist_nrow_local[p] = this->nrow_local(p);
        }

        // gather the local vectors from all processors on all processors
        int my_nrow_local(this->nrow_local());
        MPI_Allgatherv(temp_data,
                       my_nrow_local,
                       MPI_DOUBLE,
                       Values_pt,
                       dist_nrow_local,
                       dist_first_row,
                       MPI_DOUBLE,
                       this->distribution_pt()->communicator_pt()->mpi_comm());

        // update the distribution
        this->build_distribution(dist_pt);

        // delete the temp_data
        delete[] temp_data;

        // clean up
        delete[] dist_first_row;
        delete[] dist_nrow_local;
      }

      // if this vector is not distrubted but the target vector is
      else if (!this->distributed() && dist_pt->distributed())
      {
        // cache the new nrow_local
        unsigned nrow_local = dist_pt->nrow_local();

        // and first_row
        unsigned first_row = dist_pt->first_row();

        // temp storage for the new data
        double* temp_data = new double[nrow_local];

        // copy the data
        for (unsigned i = 0; i < nrow_local; i++)
        {
          temp_data[i] = Values_pt[first_row + i];
        }

        // copy to Values_pt
        delete[] Values_pt;
        Values_pt = temp_data;

        // update the distribution
        this->build_distribution(dist_pt);
      }

      // copy the Distribution
      this->build_distribution(dist_pt);
    }
#endif
  }

  //============================================================================
  /// [] access function to the (local) values of this vector
  //============================================================================
  double& DoubleVector::operator[](int i)
  {
#ifdef RANGE_CHECKING
    if (i >= int(this->nrow_local()))
    {
      std::ostringstream error_message;
      error_message << "Range Error: " << i << " is not in the range (0,"
                    << this->nrow_local() - 1 << ")";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return Values_pt[i];
  }

  //============================================================================
  /// == operator
  //============================================================================
  bool DoubleVector::operator==(const DoubleVector& v)
  {
    // if v is not setup return false
    if (v.built() && !this->built())
    {
      return false;
    }
    else if (!v.built() && this->built())
    {
      return false;
    }
    else if (!v.built() && !this->built())
    {
      return true;
    }
    else
    {
      const double* v_values_pt = v.values_pt();
      unsigned nrow_local = this->nrow_local();
      for (unsigned i = 0; i < nrow_local; i++)
      {
        if (Values_pt[i] != v_values_pt[i])
        {
          return false;
        }
      }
      return true;
    }
  }

  //============================================================================
  /// += operator
  //============================================================================
  void DoubleVector::operator+=(const DoubleVector& v)
  {
#ifdef PARANOID
    // PARANOID check that this vector is setup
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "This vector must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // PARANOID check that the vector v is setup
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // PARANOID check that the vectors have the same distribution
    if (!(*v.distribution_pt() == *this->distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The vector v and this vector must have the same "
                    << "distribution.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif //

    // cache nrow_local
    double* v_values_pt = v.values_pt();
    unsigned nrow_local = this->nrow_local();

    // Decided to keep this as a loop rather than use std::transform, because
    // this is a very simple loop and should compile to the same code.
    for (unsigned i = 0; i < nrow_local; i++)
    {
      Values_pt[i] += v_values_pt[i];
    }
  }

  //============================================================================
  /// -= operator
  //============================================================================
  void DoubleVector::operator-=(const DoubleVector& v)
  {
#ifdef PARANOID
    // PARANOID check that this vector is setup
    if (!this->distribution_built())
    {
      std::ostringstream error_message;
      error_message << "This vector must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // PARANOID check that the vector v is setup
    if (!v.built())
    {
      std::ostringstream error_message;
      error_message << "The vector v must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    // PARANOID check that the vectors have the same distribution
    if (!(*v.distribution_pt() == *this->distribution_pt()))
    {
      std::ostringstream error_message;
      error_message << "The vector v and this vector must have the same "
                    << "distribution.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // cache nrow_local
    double* v_values_pt = v.values_pt();
    unsigned nrow_local = this->nrow_local();

    // Decided to keep this as a loop rather than use std::transform, because
    // this is a very simple loop and should compile to the same code.
    for (unsigned i = 0; i < nrow_local; i++)
    {
      Values_pt[i] -= v_values_pt[i];
    }
  }


  //============================================================================
  /// Multiply by double
  //============================================================================
  void DoubleVector::operator*=(const double& d)
  {
#ifdef PARANOID
    if (!this->distribution_built())
    {
      std::ostringstream error_msg;
      error_msg << "DoubleVector must be set up.";
      throw OomphLibError(
        error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Decided to keep this as a loop rather than use std::transform, because
    // this is a very simple loop and should compile to the same code.
    for (unsigned i = 0, ni = this->nrow_local(); i < ni; i++)
    {
      Values_pt[i] *= d;
    }
  }

  //============================================================================
  /// Divide by double
  //============================================================================
  void DoubleVector::operator/=(const double& d)
  {
    // PARANOID checks are done inside operator *=

    // Decided to keep this as a loop rather than use std::transform, because
    // this is a very simple loop and should compile to the same code.
    double divisor = (1.0 / d);
    this->operator*=(divisor);
  }

  //============================================================================
  /// [] access function to the (local) values of this vector
  //============================================================================
  const double& DoubleVector::operator[](int i) const
  {
#ifdef RANGE_CHECKING
    if (i >= int(this->nrow_local()))
    {
      std::ostringstream error_message;
      error_message << "Range Error: " << i << " is not in the range (0,"
                    << this->nrow_local() - 1 << ")";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return Values_pt[i];
  }

  //============================================================================
  /// returns the maximum coefficient
  //============================================================================
  double DoubleVector::max() const
  {
    // the number of local rows
    unsigned nrow = this->nrow_local();

    // get the local maximum
    double max = 0.0;
    for (unsigned i = 0; i < nrow; i++)
    {
      if (std::fabs(Values_pt[i]) > std::fabs(max))
      {
        max = std::fabs(Values_pt[i]);
      }
    }

    // now return the maximum
#ifdef OOMPH_HAS_MPI
    // if this vector is not distributed then the local maximum is the global
    // maximum
    if (!this->distributed())
    {
      return max;
    }
    // else if the vector is distributed but only on a single processor
    // then the local maximum is the global maximum
    else if (this->distribution_pt()->communicator_pt()->nproc() == 1)
    {
      return max;
    }
    // otherwise use MPI_Allreduce to find the global maximum
    else
    {
      double local_max = max;
      MPI_Allreduce(&local_max,
                    &max,
                    1,
                    MPI_DOUBLE,
                    MPI_MAX,
                    this->distribution_pt()->communicator_pt()->mpi_comm());
      return max;
    }
#else
    return max;
#endif
  }

  //============================================================================
  /// output the contents of the vector
  //============================================================================
  void DoubleVector::output(std::ostream& outfile,
                            const int& output_precision) const
  {
    // temp pointer to values
    double* temp;

    // number of global row
    unsigned nrow = this->nrow();

#ifdef OOMPH_HAS_MPI

    // number of local rows
    int nrow_local = this->nrow_local();

    // gather from all processors
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
      // number of processors
      int nproc = this->distribution_pt()->communicator_pt()->nproc();

      // number of gobal row
      unsigned nrow = this->nrow();

      // get the vector of first_row s and nrow_local s
      int* dist_first_row = new int[nproc];
      int* dist_nrow_local = new int[nproc];
      for (int p = 0; p < nproc; p++)
      {
        dist_first_row[p] = this->first_row(p);
        dist_nrow_local[p] = this->nrow_local(p);
      }

      // gather
      temp = new double[nrow];
      MPI_Allgatherv(Values_pt,
                     nrow_local,
                     MPI_DOUBLE,
                     temp,
                     dist_nrow_local,
                     dist_first_row,
                     MPI_DOUBLE,
                     this->distribution_pt()->communicator_pt()->mpi_comm());

      // clean up
      delete[] dist_first_row;
      delete[] dist_nrow_local;
    }
    else
    {
      temp = Values_pt;
    }
#else
    temp = Values_pt;
#endif

    // output
    // Store the precision so we can revert it.
    std::streamsize old_precision = 0;
    if (output_precision > 0)
    {
      old_precision = outfile.precision();
      outfile << std::setprecision(output_precision);
    }

    for (unsigned i = 0; i < nrow; i++)
    {
      outfile << i << " " << temp[i] << std::endl;
    }

    // Revert the precision.
    if (output_precision > 0)
    {
      outfile << std::setprecision(old_precision);
    }

    // clean up if requires
#ifdef OOMPH_HAS_MPI
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
      delete[] temp;
    }
#endif
  }

  //============================================================================
  /// output the local contents of the vector
  //============================================================================
  void DoubleVector::output_local_values(std::ostream& outfile,
                                         const int& output_precision) const
  {
    // Number of local rows.
    unsigned nrow_local = this->nrow_local();

    // output
    // Store the precision so we can revert it.
    std::streamsize old_precision = 0;
    if (output_precision > 0)
    {
      old_precision = outfile.precision();
      outfile << std::setprecision(output_precision);
    }

    for (unsigned i = 0; i < nrow_local; i++)
    {
      outfile << i << " " << Values_pt[i] << std::endl;
    }

    // Revert the precision.
    if (output_precision > 0)
    {
      outfile << std::setprecision(old_precision);
    }
  }

  //============================================================================
  /// output the local contents of the vector with the first row offset.
  //============================================================================
  void DoubleVector::output_local_values_with_offset(
    std::ostream& outfile, const int& output_precision) const
  {
    // Number of local rows.
    unsigned nrow_local = this->nrow_local();

    // First row on this processor.
    unsigned first_row = this->first_row();

    // output
    // Store the precision so we can revert it.
    std::streamsize old_precision = 0;
    if (output_precision > 0)
    {
      old_precision = outfile.precision();
      outfile << std::setprecision(output_precision);
    }

    for (unsigned i = 0; i < nrow_local; i++)
    {
      outfile << (i + first_row) << " " << Values_pt[i] << std::endl;
    }

    // Revert the precision.
    if (output_precision > 0)
    {
      outfile << std::setprecision(old_precision);
    }
  }

  //============================================================================
  /// compute the dot product of this vector with the vector vec
  //============================================================================
  double DoubleVector::dot(const DoubleVector& vec) const
  {
#ifdef PARANOID
    // paranoid check that the vector is setup
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "This vector must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!vec.built())
    {
      std::ostringstream error_message;
      error_message << "The input vector be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*this->distribution_pt() != *vec.distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of this vector and the vector vec "
                    << "must be the same."
                    << "\n\n  this: " << *this->distribution_pt()
                    << "\n  vec:  " << *vec.distribution_pt();
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // compute the local norm
    unsigned nrow_local = this->nrow_local();
    double n = 0.0;
    const double* vec_values_pt = vec.values_pt();
    for (unsigned i = 0; i < nrow_local; i++)
    {
      n += Values_pt[i] * vec_values_pt[i];
    }

    // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
    double n2 = n;
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
      MPI_Allreduce(&n,
                    &n2,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    this->distribution_pt()->communicator_pt()->mpi_comm());
    }
    n = n2;
#endif

    // and return;
    return n;
  }

  //============================================================================
  /// compute the 2 norm of this vector
  //============================================================================
  double DoubleVector::norm() const
  {
#ifdef PARANOID
    // paranoid check that the vector is setup
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "This vector must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // compute the local norm
    unsigned nrow_local = this->nrow_local();
    double n = 0;
    for (unsigned i = 0; i < nrow_local; i++)
    {
      n += Values_pt[i] * Values_pt[i];
    }

    // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
    double n2 = n;
    if (this->distributed() &&
        this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
      MPI_Allreduce(&n,
                    &n2,
                    1,
                    MPI_DOUBLE,
                    MPI_SUM,
                    this->distribution_pt()->communicator_pt()->mpi_comm());
    }
    n = n2;
#endif

    // sqrt the norm
    n = sqrt(n);

    // and return
    return n;
  }

  //============================================================================
  /// compute the A-norm using the matrix at matrix_pt
  //============================================================================
  double DoubleVector::norm(const CRDoubleMatrix* matrix_pt) const
  {
#ifdef PARANOID
    // paranoid check that the vector is setup
    if (!this->built())
    {
      std::ostringstream error_message;
      error_message << "This vector must be setup.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!matrix_pt->built())
    {
      std::ostringstream error_message;
      error_message << "The input matrix be built.";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (*this->distribution_pt() != *matrix_pt->distribution_pt())
    {
      std::ostringstream error_message;
      error_message << "The distribution of this vector and the matrix at "
                    << "matrix_pt must be the same";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // compute the matrix norm
    DoubleVector x(this->distribution_pt(), 0.0);
    matrix_pt->multiply(*this, x);
    return sqrt(this->dot(x));
  }

  /// output operator
  std::ostream& operator<<(std::ostream& out, const DoubleVector& v)
  {
    // Do the first value outside the loop to get the ", "s right.
    out << "[" << v[0];

    for (unsigned i = 1, ni = v.nrow_local(); i < ni; i++)
    {
      out << ", " << v[i];
    }
    out << "]";

    return out;
  }

  //=================================================================
  /// Namespace for helper functions for DoubleVectors
  //=================================================================
  namespace DoubleVectorHelpers
  {
    //===========================================================================
    /// Concatenate DoubleVectors.
    /// Takes a Vector of DoubleVectors. If the out vector is built, we will not
    /// build a new distribution. Otherwise we build a uniform distribution.
    ///
    /// The rows of the out vector is seen "as it is" in the in vectors.
    /// For example, if we have DoubleVectors with distributions A and B,
    /// distributed across two processors (p0 and p1),
    ///
    /// A: [a0] (on p0)    B: [b0] (on p0)
    ///    [a1] (on p1)       [b1] (on P1),
    ///
    /// then the out_vector is
    ///
    /// [a0  (on p0)
    ///  a1] (on p0)
    /// [b0]  (on p1)
    ///  b1] (on p1),
    ///
    /// Communication is required between processors. The sum of the global
    /// number of rows in the in vectors must equal to the global number of rows
    /// in the out vector. This condition must be met if one is to supply an out
    /// vector with a distribution, otherwise we can let the function generate
    /// the out vector distribution itself.
    //===========================================================================
    void concatenate(const Vector<DoubleVector*>& in_vector_pt,
                     DoubleVector& out_vector)
    {
      // How many in vectors to concatenate?
      unsigned nvectors = in_vector_pt.size();

      // PARANIOD checks which involves the in vectors only
#ifdef PARANOID
      // Check that there is at least one vector.
      if (nvectors == 0)
      {
        std::ostringstream error_message;
        error_message << "There is no vector to concatenate...\n"
                      << "Perhaps you forgot to fill in_vector_pt?\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Does this vector need concatenating?
      if (nvectors == 1)
      {
        std::ostringstream warning_message;
        warning_message << "There is only one vector to concatenate...\n"
                        << "This does not require concatenating...\n";
        OomphLibWarning(warning_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all the DoubleVectors in in_vector_pt are built
      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        if (!in_vector_pt[vec_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The vector in position " << vec_i
                        << " is not built.\n"
                        << "I cannot concatenate an unbuilt vector.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // The communicator pointer for the first in vector.
      const OomphCommunicator* const comm_pt =
        in_vector_pt[0]->distribution_pt()->communicator_pt();

      // Check if the first in vector is distributed.
      bool distributed = in_vector_pt[0]->distributed();

      // If the out vector is not built, build it with a uniform distribution.
      if (!out_vector.built())
      {
        // Nrow for the out vector is the sum of the nrow of the in vectors.
        unsigned tmp_nrow = 0;
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          tmp_nrow += in_vector_pt[vec_i]->nrow();
        }

        // Build the out vector with uniform distribution.
        out_vector.build(
          LinearAlgebraDistribution(comm_pt, tmp_nrow, distributed), 0.0);
      }
      else
      {
#ifdef PARANOID
        // Check that the sum of nrow of in vectors match the nrow in the out
        // vectors.
        unsigned in_nrow = 0;
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          in_nrow += in_vector_pt[vec_i]->nrow();
        }

        if (in_nrow != out_vector.nrow())
        {
          std::ostringstream error_message;
          error_message << "The sum of nrow of the in vectors does not match\n"
                        << "the nrow of the out vector.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }

#ifdef PARANOID
      // Check that all communicators of the vectors to concatenate are the same
      // by comparing all communicators against the out vector.
      const OomphCommunicator out_comm =
        *(out_vector.distribution_pt()->communicator_pt());

      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        // Get the Communicator for the current vector.
        const OomphCommunicator in_comm =
          *(in_vector_pt[vec_i]->distribution_pt()->communicator_pt());

        if (out_comm != in_comm)
        {
          std::ostringstream error_message;
          error_message << "The vector in position " << vec_i << " has a\n"
                        << "different communicator from the out vector.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that the distributed boolean is the same for all vectors.
      if (out_comm.nproc() != 1)
      {
        const bool out_distributed = out_vector.distributed();
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          if (out_distributed != in_vector_pt[vec_i]->distributed())
          {
            std::ostringstream error_message;
            error_message << "The vector in position " << vec_i << " has a\n"
                          << "different distributed boolean from "
                          << "the out vector.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif


      // Now we do the concatenation.
      if ((comm_pt->nproc() == 1) || !distributed)
      {
        // Serial version of the code.
        // This is trivial, we simply loop through the in vectors and
        // fill in the out vector.

        // Out vector index.
        unsigned out_i = 0;

        // Out vector values.
        double* out_value_pt = out_vector.values_pt();

        // Loop through the in vectors.
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          // Nrow of current in vector.
          unsigned in_nrow = in_vector_pt[vec_i]->nrow();

          // In vector values.
          double* in_value_pt = in_vector_pt[vec_i]->values_pt();

          // Loop through the entries of this in vector.
          for (unsigned i = 0; i < in_nrow; i++)
          {
            out_value_pt[out_i++] = in_value_pt[i];
          }
        }
      }
      // Otherwise we are dealing with a distributed vector.
      else
      {
#ifdef OOMPH_HAS_MPI
        // Get the number of processors
        unsigned nproc = comm_pt->nproc();

        // My rank
        unsigned my_rank = comm_pt->my_rank();

        // Storage for the data (per processor) to send
        Vector<Vector<double>> values_to_send(nproc);

        // The sum of the nrow for the in vectors (so far). This is used as an
        // offset to calculate the global equation number in the out vector
        unsigned long sum_of_vec_nrow = 0;

        // Loop over the in vectors and work out:
        // out_p: the rank the of receiving processor
        // out_local_eqn: the local equation number of the receiving processor
        //
        // Then put the value and out_local_eqn at out_p in values_to_send

        LinearAlgebraDistribution* out_distribution_pt =
          out_vector.distribution_pt();
        for (unsigned in_vec_i = 0; in_vec_i < nvectors; in_vec_i++)
        {
          // Loop through the local equations
          unsigned in_vec_nrow_local = in_vector_pt[in_vec_i]->nrow_local();
          unsigned in_vec_first_row = in_vector_pt[in_vec_i]->first_row();

          for (unsigned in_row_i = 0; in_row_i < in_vec_nrow_local; in_row_i++)
          {
            // Calculate the global equation number for this in_row_i
            unsigned out_global_eqn =
              in_row_i + in_vec_first_row + sum_of_vec_nrow;

            // Get the processor that this global row belongs to.
            // The rank_of_global_row(...) function loops through all the
            // processors and does two unsigned comparisons. Since we have to do
            // this for every row, it may be better to store a list mapping for
            // very large number of processors.
            unsigned out_p =
              out_distribution_pt->rank_of_global_row(out_global_eqn);
            //         unsigned out_p = out_distribution_pt
            //           ->rank_of_global_row_map(out_global_eqn);

            // Knowing out_p enables us to work out the out_first_row and
            // out_local_eqn.
            unsigned out_first_row = out_distribution_pt->first_row(out_p);
            unsigned out_local_eqn = out_global_eqn - out_first_row;

            // Now push back the out_local_eqn and the value
            values_to_send[out_p].push_back(out_local_eqn);
            values_to_send[out_p].push_back(
              (*in_vector_pt[in_vec_i])[in_row_i]);
          }

          // Update the offset.
          sum_of_vec_nrow += in_vector_pt[in_vec_i]->nrow();
        }

        // Prepare to send the data!

        // Storage for the number of data to be sent to each processor.
        Vector<int> send_n(nproc, 0);

        // Storage for all the values to be send to each processor.
        Vector<double> send_values_data;

        // Storage location within send_values_data
        Vector<int> send_displacement(nproc, 0);

        // Get the total amount of data which needs to be sent, so we can
        // reserve space for it.
        unsigned total_ndata = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          if (rank != my_rank)
          {
            total_ndata += values_to_send[rank].size();
          }
        }

        // Now we don't have to re-allocate data/memory when push_back is
        // called. Nb. Using push_back without reserving memory may cause
        // multiple re-allocation behind the scenes, this is expensive.
        send_values_data.reserve(total_ndata);

        // Loop over all the processors to "flat pack" the data for sending.
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          // Set the offset for the current processor
          send_displacement[rank] = send_values_data.size();

          // Don't bother to do anything if
          // the processor in the loop is the current processor.
          if (rank != my_rank)
          {
            // Put the values into the send data vector.
            unsigned n_data = values_to_send[rank].size();
            for (unsigned j = 0; j < n_data; j++)
            {
              send_values_data.push_back(values_to_send[rank][j]);
            } // Loop over the data
          } // if rank != my_rank

          // Find the number of data to be added to the vector.
          send_n[rank] = send_values_data.size() - send_displacement[rank];
        } // Loop over processors

        // Storage for the number of data to be received from each processor.
        Vector<int> receive_n(nproc, 0);
        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     comm_pt->mpi_comm());

        // Prepare the data to be received
        // by working out the displacement from the received data.
        Vector<int> receive_displacement(nproc, 0);
        int receive_data_count = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors.
        // Make sure that it has size of at least one.
        if (receive_data_count == 0)
        {
          receive_data_count++;
        }
        Vector<double> receive_values_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault.
        if (send_values_data.size() == 0)
        {
          send_values_data.resize(1);
        }

        // Now send the data between all processors
        MPI_Alltoallv(&send_values_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_DOUBLE,
                      &receive_values_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_DOUBLE,
                      comm_pt->mpi_comm());

        // Data from all other processors are stored in:
        // receive_values_data
        // Data already on this processor is stored in:
        // values_to_send[my_rank]

        // Loop through the data on this processor.
        unsigned location_i = 0;
        unsigned my_values_to_send_size = values_to_send[my_rank].size();
        while (location_i < my_values_to_send_size)
        {
          out_vector[unsigned(values_to_send[my_rank][location_i])] =
            values_to_send[my_rank][location_i + 1];

          location_i += 2;
        }

        // Before we loop through the data on other processors, we need to check
        // if any data has been received.
        bool data_has_been_received = false;
        unsigned send_rank = 0;
        while (send_rank < nproc)
        {
          if (receive_n[send_rank] > 0)
          {
            data_has_been_received = true;
            break;
          }
          send_rank++;
        }

        location_i = 0;
        if (data_has_been_received)
        {
          unsigned receive_values_data_size = receive_values_data.size();
          while (location_i < receive_values_data_size)
          {
            out_vector[unsigned(receive_values_data[location_i])] =
              receive_values_data[location_i + 1];
            location_i += 2;
          }
        }
#else
        {
          std::ostringstream error_message;
          error_message << "I don't know what to do with distributed vectors\n"
                        << "without MPI... :(";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    } // function concatenate

    //===========================================================================
    /// Wrapper around the other concatenate(...) function.
    /// Be careful with Vector of vectors. If the DoubleVectors are resized,
    /// there could be reallocation of memory. If we wanted to use the function
    /// which takes a Vector of pointers to DoubleVectors, we would either have
    /// to invoke new and remember to delete, or create a temporary Vector to
    /// store pointers to the DoubleVector objects.
    /// This wrapper is meant to make life easier for the user by avoiding calls
    /// to new/delete AND without creating a temporary vector of pointers to
    /// DoubleVectors.
    /// If we had C++ 11, this would be so much nicer since we can use smart
    /// pointers which will delete themselves, so we do not have to remember
    /// to delete!
    //===========================================================================
    void concatenate(Vector<DoubleVector>& in_vector, DoubleVector& out_vector)
    {
      const unsigned n_in_vector = in_vector.size();

      Vector<DoubleVector*> in_vector_pt(n_in_vector, 0);

      for (unsigned i = 0; i < n_in_vector; i++)
      {
        in_vector_pt[i] = &in_vector[i];
      }

      DoubleVectorHelpers::concatenate(in_vector_pt, out_vector);
    } // function concatenate

    //===========================================================================
    /// Split a DoubleVector into the out DoubleVectors.
    /// Let vec_A be the in Vector, and let vec_B and vec_C be the out vectors.
    /// Then the splitting of vec_A is depicted below:
    /// vec_A: [a0  (on p0)
    ///         a1] (on p0)
    ///        [a2  (on p1)
    ///         a3] (on p1)
    ///
    /// vec_B: [a0] (on p0)    vec_C: [a2] (on p0)
    ///        [a1] (on p1)           [a3] (on p1)
    ///
    /// Communication is required between processors.
    /// The out_vector_pt must contain pointers to DoubleVector which has
    /// already been built with the correct distribution; the sum of the number
    /// of global row of the out vectors must be the same the number of global
    /// rows of the in vector.
    //===========================================================================
    void split(const DoubleVector& in_vector,
               Vector<DoubleVector*>& out_vector_pt)
    {
      // How many out vectors do we have?
      unsigned nvec = out_vector_pt.size();
#ifdef PARANOID

      // Check that the in vector is built.
      if (!in_vector.built())
      {
        std::ostringstream error_message;
        error_message << "The in_vector is not built.\n"
                      << "Please build it!.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all the out vectors are built.
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        if (!out_vector_pt[vec_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The vector at position " << vec_i
                        << "  is not built.\n"
                        << "Please build it!.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that the sum of the nrow from out vectors is the same as the
      // nrow from in_vector.
      unsigned out_nrow_sum = 0;
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        out_nrow_sum += out_vector_pt[vec_i]->nrow();
      }

      if (in_vector.nrow() != out_nrow_sum)
      {
        std::ostringstream error_message;
        error_message << "The global number of rows in the in_vector\n"
                      << "is not equal to the sum of the global nrows\n"
                      << "of the in vectors.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all communicators are the same. We use a communicator to
      // get the number of processors and my_rank. So we would like them to be
      // the same for in_vector and all out vectors.
      const OomphCommunicator in_vector_comm =
        *(in_vector.distribution_pt()->communicator_pt());
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        const OomphCommunicator dist_i_comm =
          *(out_vector_pt[vec_i]->distribution_pt()->communicator_pt());

        if (in_vector_comm != dist_i_comm)
        {
          std::ostringstream error_message;
          error_message << "The communicator for the distribution in the \n"
                        << "position " << vec_i
                        << " is not the same as the in_vector\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that the distributed boolean is the same for all vectors.
      bool para_distributed = in_vector.distributed();

      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        if (para_distributed != out_vector_pt[vec_i]->distributed())
        {
          std::ostringstream error_message;
          error_message
            << "The vector in position " << vec_i << " does not \n"
            << " have the same distributed boolean as the in_vector\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

#endif

      // The communicator.
      const OomphCommunicator* const comm_pt =
        in_vector.distribution_pt()->communicator_pt();

      // Is this distributed?
      bool distributed = in_vector.distributed();

      // The serial code.
      if ((comm_pt->nproc() == 1) || !distributed)
      {
        // Serial version of the code: loop through all the out vectors and
        // insert the elements of in_vector.

        // index for in vector, and in vector values.
        unsigned in_vec_i = 0;
        double* in_value_pt = in_vector.values_pt();

        // Fill in the out vectors.
        for (unsigned out_vec_i = 0; out_vec_i < nvec; out_vec_i++)
        {
          // out vector nrow and values.
          unsigned out_nrow = out_vector_pt[out_vec_i]->nrow();
          double* out_value_pt = out_vector_pt[out_vec_i]->values_pt();

          // Fill in the current out vector.
          for (unsigned out_val_i = 0; out_val_i < out_nrow; out_val_i++)
          {
            out_value_pt[out_val_i] = in_value_pt[in_vec_i++];
          }
        }
      }
      // Otherwise we are dealing with a distributed vector.
      else
      {
#ifdef OOMPH_HAS_MPI
        // For each entry in the in_vector, we need to work out:
        // 1) Which out vector this entry belongs to,
        // 2) which processor to send the data to and
        // 3) the local equation number in the out vector.
        //
        // We know the in_local_eqn, we can work out the in_global_eqn.
        //
        // From in_global_eqn we can work out the out vector and
        // the out_global_eqn.
        //
        // The out_global_eqn allows us to determine which processor to send to.
        // With the out_p (processor to send data to) and out vector, we get the
        // out_first_row which then allows us to work out the out_local_eqn.


        // Get the number of processors
        unsigned nproc = comm_pt->nproc();

        // My rank
        unsigned my_rank = comm_pt->my_rank();

        // Storage for the data (per processor) to send.
        Vector<Vector<double>> values_to_send(nproc);

        // Sum of the nrow of the out vectors so far. This is used to work out
        // which out_vector a in_global_eqn belongs to.
        Vector<unsigned> sum_of_out_nrow(nvec + 1);
        for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
        {
          sum_of_out_nrow[vec_i + 1] =
            sum_of_out_nrow[vec_i] + out_vector_pt[vec_i]->nrow();
        }

        // Loop through the in_vector local values.
        unsigned in_nrow_local = in_vector.nrow_local();
        for (unsigned in_local_eqn = 0; in_local_eqn < in_nrow_local;
             in_local_eqn++)
        {
          // The global equation number of this row.
          unsigned in_global_eqn = in_local_eqn + in_vector.first_row();

          // Which out_vector does this in_global_eqn belong to?
          unsigned out_vector_i = 0;
          while (in_global_eqn < sum_of_out_nrow[out_vector_i] ||
                 in_global_eqn >= sum_of_out_nrow[out_vector_i + 1])
          {
            out_vector_i++;
          }

          // The out_global_eqn
          // (this is the global equation in the current out vector)
          unsigned out_global_eqn =
            in_global_eqn - sum_of_out_nrow[out_vector_i];

          // The processor to send this row to.
          unsigned out_p =
            out_vector_pt[out_vector_i]->distribution_pt()->rank_of_global_row(
              out_global_eqn);

          // The local_eqn in the out_vector_i
          unsigned out_local_eqn =
            out_global_eqn -
            out_vector_pt[out_vector_i]->distribution_pt()->first_row(out_p);


          // Fill in the data to send

          // Which out vector to put this data in.
          values_to_send[out_p].push_back(out_vector_i);

          // The local equation of the data.
          values_to_send[out_p].push_back(out_local_eqn);

          // The actual data.
          values_to_send[out_p].push_back(in_vector[in_local_eqn]);
        }

        // Prepare to send the data!

        // Storage for the number of data to be sent to each processor.
        Vector<int> send_n(nproc, 0);

        // Storage for all the values to be send to each processor.
        Vector<double> send_values_data;

        // Storage location within send_values_data
        Vector<int> send_displacement(nproc, 0);

        // Get the total amount of data which needs to be sent, so we can
        // reserve space for it.
        unsigned total_ndata = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          if (rank != my_rank)
          {
            total_ndata += values_to_send[rank].size();
          }
        }

        // Now we don't have to re-allocate data/memory when push_back is
        // called. Nb. Using push_back without reserving memory may cause
        // multiple re-allocation behind the scenes, this is expensive.
        send_values_data.reserve(total_ndata);

        // Loop over all the processors to "flat pack" the data for sending.
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          // Set the offset for the current processor
          send_displacement[rank] = send_values_data.size();

          // Don't bother to do anything if
          // the processor in the loop is the current processor.
          if (rank != my_rank)
          {
            // Put the values into the send data vector.
            unsigned n_data = values_to_send[rank].size();
            for (unsigned j = 0; j < n_data; j++)
            {
              send_values_data.push_back(values_to_send[rank][j]);
            } // Loop over the data
          } // if rank != my_rank

          // Find the number of data to be added to the vector.
          send_n[rank] = send_values_data.size() - send_displacement[rank];
        } // Loop over processors

        // Storage for the number of data to be received from each processor.
        Vector<int> receive_n(nproc, 0);
        MPI_Alltoall(&send_n[0],
                     1,
                     MPI_INT,
                     &receive_n[0],
                     1,
                     MPI_INT,
                     comm_pt->mpi_comm());

        // Prepare the data to be received
        // by working out the displacement from the received data.
        Vector<int> receive_displacement(nproc, 0);
        int receive_data_count = 0;
        for (unsigned rank = 0; rank < nproc; rank++)
        {
          receive_displacement[rank] = receive_data_count;
          receive_data_count += receive_n[rank];
        }

        // Now resize the receive buffer for all data from all processors.
        // Make sure that it has size of at least one.
        if (receive_data_count == 0)
        {
          receive_data_count++;
        }
        Vector<double> receive_values_data(receive_data_count);

        // Make sure that the send buffer has size at least one
        // so that we don't get a segmentation fault.
        if (send_values_data.size() == 0)
        {
          send_values_data.resize(1);
        }

        // Now send the data between all processors
        MPI_Alltoallv(&send_values_data[0],
                      &send_n[0],
                      &send_displacement[0],
                      MPI_DOUBLE,
                      &receive_values_data[0],
                      &receive_n[0],
                      &receive_displacement[0],
                      MPI_DOUBLE,
                      comm_pt->mpi_comm());

        // Data from all other processors are stored in:
        // receive_values_data
        // Data already on this processor is stored in:
        // values_to_send[my_rank]
        //

        // Index for values_to_send Vector.
        unsigned location_i = 0;
        // Loop through the data on this processor
        unsigned my_values_to_send_size = values_to_send[my_rank].size();
        while (location_i < my_values_to_send_size)
        {
          // The vector to put the values in.
          unsigned out_vector_i =
            unsigned(values_to_send[my_rank][location_i++]);

          // Where to put the value.
          unsigned out_local_eqn =
            unsigned(values_to_send[my_rank][location_i++]);

          // The actual value!
          double out_value = values_to_send[my_rank][location_i++];

          // Insert the value in the out vector.
          (*out_vector_pt[out_vector_i])[out_local_eqn] = out_value;
        }

        // Before we loop through the data on other processors, we need to check
        // if any data has been received. This is because the
        // receive_values_data has been resized to at least one, even if no data
        // is sent.
        bool data_has_been_received = false;
        unsigned send_rank = 0;
        while (send_rank < nproc)
        {
          if (receive_n[send_rank] > 0)
          {
            data_has_been_received = true;
            break;
          }
          send_rank++;
        }

        // Reset the index, it is now being used to index the
        // receive_values_data vector.
        location_i = 0;
        if (data_has_been_received)
        {
          // Extract the data and put it into the out vector.
          unsigned receive_values_data_size = receive_values_data.size();
          while (location_i < receive_values_data_size)
          {
            // Which out vector to put the value in?
            unsigned out_vector_i = unsigned(receive_values_data[location_i++]);

            // Where in the out vector to put the value?
            unsigned out_local_eqn =
              unsigned(receive_values_data[location_i++]);

            // The value to put in.
            double out_value = receive_values_data[location_i++];

            // Insert the value in the out vector.
            (*out_vector_pt[out_vector_i])[out_local_eqn] = out_value;
          }
        }
#else
        {
          std::ostringstream error_message;
          error_message << "You have a distributed vector but with no mpi...\n"
                        << "I don't know what to do :( \n";
          throw OomphLibError(
            error_message.str(), "RYARAYERR", OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    } // function split(...)

    //===========================================================================
    /// Wrapper around the other split(...) function.
    /// Be careful with Vector of vectors. If the DoubleVectors are resized,
    /// there could be reallocation of memory. If we wanted to use the function
    /// which takes a Vector of pointers to DoubleVectors, we would either have
    /// to invoke new and remember to delete, or create a temporary Vector to
    /// store pointers to the DoubleVector objects.
    /// This wrapper is meant to make life easier for the user by avoiding calls
    /// to new/delete AND without creating a temporary vector of pointers to
    /// DoubleVectors.
    /// If we had C++ 11, this would be so much nicer since we can use smart
    /// pointers which will delete themselves, so we do not have to remember
    /// to delete!
    //===========================================================================
    void split(const DoubleVector& in_vector, Vector<DoubleVector>& out_vector)
    {
      const unsigned n_out_vector = out_vector.size();
      Vector<DoubleVector*> out_vector_pt(n_out_vector, 0);

      for (unsigned i = 0; i < n_out_vector; i++)
      {
        out_vector_pt[i] = &out_vector[i];
      }

      DoubleVectorHelpers::split(in_vector, out_vector_pt);
    } // function split(...)

    //===========================================================================
    /// Concatenate DoubleVectors.
    /// Takes a Vector of DoubleVectors. If the out vector is built, we will not
    /// build a new distribution. Otherwise a new distribution will be built
    /// using LinearAlgebraDistribution::concatenate(...).
    ///
    /// The out vector has its rows permuted according to the individual
    /// distributions of the in vectors. For example, if we have DoubleVectors
    /// with distributions A and B, distributed across two processors
    /// (p0 and p1),
    ///
    /// A: [a0] (on p0)    B: [b0] (on p0)
    ///    [a1] (on p1)       [b1] (on P1),
    ///
    /// then the out_vector is
    ///
    /// [a0  (on p0)
    ///  b0] (on p0)
    /// [a1  (on p1)
    ///  b1] (on p1),
    ///
    /// as opposed to
    ///
    /// [a0  (on p0)
    ///  a1] (on p0)
    /// [b0  (on p1)
    ///  b1] (on p1).
    ///
    /// Note (1): The out vector may not be uniformly distributed even
    /// if the in vectors have uniform distributions. The nrow_local of the
    /// out vector will be the sum of the nrow_local of the in vectors.
    /// Try this out with two distributions of global rows 3 and 5, uniformly
    /// distributed across two processors. Compare this against a distribution
    /// of global row 8 distributed across two processors.
    ///
    /// There are no MPI send and receive, the data stays on the processor
    /// as defined by the distributions from the in vectors.
    //===========================================================================
    void concatenate_without_communication(
      const Vector<DoubleVector*>& in_vector_pt, DoubleVector& out_vector)
    {
      // How many in vectors do we want to concatenate?
      unsigned nvectors = in_vector_pt.size();

      // PARANOID checks which involves the in vectors only.
#ifdef PARANOID
      // Check that there is at least one vector.
      if (nvectors == 0)
      {
        std::ostringstream error_message;
        error_message << "There is no vector to concatenate...\n"
                      << "Perhaps you forgot to fill in_vector_pt?\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Does this vector need concatenating?
      if (nvectors == 1)
      {
        std::ostringstream warning_message;
        warning_message << "There is only one vector to concatenate...\n"
                        << "This does not require concatenating...\n";
        OomphLibWarning(warning_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all the DoubleVectors in in_vector_pt are built
      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        if (!in_vector_pt[vec_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The vector in position " << vec_i
                        << " is not built.\n"
                        << "I cannot concatenate an unbuilt vector.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // If the out vector is not built, build it with the correct distribution.
      if (!out_vector.built())
      {
        Vector<LinearAlgebraDistribution*> in_distribution_pt(nvectors, 0);
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          in_distribution_pt[vec_i] = in_vector_pt[vec_i]->distribution_pt();
        }

        LinearAlgebraDistribution tmp_distribution;
        LinearAlgebraDistributionHelpers::concatenate(in_distribution_pt,
                                                      tmp_distribution);
        out_vector.build(tmp_distribution, 0.0);
      }

      // PARANOID checks which involves all in vectors and out vectors.
#ifdef PARANOID

      // Check that all communicators of the vectors to concatenate are the same
      // by comparing all communicators against the out vector.
      const OomphCommunicator out_comm =
        *(out_vector.distribution_pt()->communicator_pt());

      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        // Get the Communicator for the current vector.
        const OomphCommunicator in_comm =
          *(in_vector_pt[vec_i]->distribution_pt()->communicator_pt());

        if (out_comm != in_comm)
        {
          std::ostringstream error_message;
          error_message << "The vector in position " << vec_i << " has a\n"
                        << "different communicator from the out vector.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that the distributed boolean is the same for all vectors.
      if (out_comm.nproc() > 1)
      {
        const bool out_distributed = out_vector.distributed();
        for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
        {
          if (out_distributed != in_vector_pt[vec_i]->distributed())
          {
            std::ostringstream error_message;
            error_message << "The vector in position " << vec_i << " has a\n"
                          << "different distributed boolean from the "
                          << "out vector.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }

      // Check that the distribution from the out vector is indeed the
      // same as the one created by
      // LinearAlgebraDistributionHelpers::concatenate(...). This test is
      // redundant if the out_vector is not built to begin with.

      // Create tmp_distribution, a concatenation of all distributions from
      // the in vectors.
      Vector<LinearAlgebraDistribution*> in_distribution_pt(nvectors, 0);
      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        in_distribution_pt[vec_i] = in_vector_pt[vec_i]->distribution_pt();
      }

      LinearAlgebraDistribution tmp_distribution;
      LinearAlgebraDistributionHelpers::concatenate(in_distribution_pt,
                                                    tmp_distribution);
      // The the distribution from the out vector.
      LinearAlgebraDistribution out_distribution =
        *(out_vector.distribution_pt());

      // Compare them!
      if (tmp_distribution != out_distribution)
      {
        std::ostringstream error_message;
        error_message << "The distribution of the out vector is not correct.\n"
                      << "Please call the function with a cleared out vector,\n"
                      << "or compare the distribution of the out vector with\n"
                      << "the distribution created by\n"
                      << "LinearAlgebraDistributionHelpers::concatenate(...)\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Do not need these distributions.
      tmp_distribution.clear();
      out_distribution.clear();
#endif


      unsigned out_value_offset = 0;

      double* out_value_pt = out_vector.values_pt();

      // Loop through the vectors.
      for (unsigned vec_i = 0; vec_i < nvectors; vec_i++)
      {
        // Get the nrow_local and
        // pointer to the values for the current in vector.
        unsigned in_vector_nrow_local = in_vector_pt[vec_i]->nrow_local();
        double* in_vector_value_pt = in_vector_pt[vec_i]->values_pt();

        // Loop through the local values and inset them into the out_vector.
        for (unsigned val_i = 0; val_i < in_vector_nrow_local; val_i++)
        {
          out_value_pt[out_value_offset + val_i] = in_vector_value_pt[val_i];
        }

        // Update the offset.
        out_value_offset += in_vector_nrow_local;
      }
    } // function concatenate_without_communication

    //===========================================================================
    /// Wrapper around the other concatenate_without_communication(...)
    /// function.
    /// Be careful with Vector of vectors. If the DoubleVectors are resized,
    /// there could be reallocation of memory. If we wanted to use the function
    /// which takes a Vector of pointers to DoubleVectors, we would either have
    /// to invoke new and remember to delete, or create a temporary Vector to
    /// store pointers to the DoubleVector objects.
    /// This wrapper is meant to make life easier for the user by avoiding calls
    /// to new/delete AND without creating a temporary vector of pointers to
    /// DoubleVectors.
    /// If we had C++ 11, this would be so much nicer since we can use smart
    /// pointers which will delete themselves, so we do not have to remember
    /// to delete!
    //===========================================================================
    void concatenate_without_communication(Vector<DoubleVector>& in_vector,
                                           DoubleVector& out_vector)
    {
      const unsigned n_in_vector = in_vector.size();

      Vector<DoubleVector*> in_vector_pt(n_in_vector, 0);

      for (unsigned i = 0; i < n_in_vector; i++)
      {
        in_vector_pt[i] = &in_vector[i];
      }

      DoubleVectorHelpers::concatenate_without_communication(in_vector_pt,
                                                             out_vector);
    } // function concatenate_without_communication

    //===========================================================================
    /// Split a DoubleVector into the out DoubleVectors.
    /// Data stays on its current processor, no data is sent between processors.
    /// This results in our vectors which are a permutation of the in vector.
    ///
    /// Let vec_A be the in Vector, and let vec_B and vec_C be the out vectors.
    /// Then the splitting of vec_A is depicted below:
    /// vec_A: [a0  (on p0)
    ///         a1] (on p0)
    ///        [a2  (on p1)
    ///         a3] (on p1)
    ///
    /// vec_B: [a0] (on p0)    vec_C: [a1] (on p0)
    ///        [a2] (on p1)           [a3] (on p1).
    ///
    /// This means that the distribution of the in vector MUST be a
    /// concatenation of the out vector distributions, refer to
    /// LinearAlgebraDistributionHelpers::concatenate(...) to concatenate
    /// distributions.
    //===========================================================================
    void split_without_communication(const DoubleVector& in_vector,
                                     Vector<DoubleVector*>& out_vector_pt)
    {
      // How many out vectors do we need?
      unsigned nvec = out_vector_pt.size();

#ifdef PARANOID
      // Check that in_vector is built
      if (!in_vector.built())
      {
        std::ostringstream error_message;
        error_message << "The in_vector is not built.\n"
                      << "Please build it!.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Check that all out vectors are built.
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        if (!out_vector_pt[vec_i]->built())
        {
          std::ostringstream error_message;
          error_message << "The vector at position " << vec_i
                        << " is not built.\n"
                        << "Please build it!.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that the concatenation of distributions from the out vectors is
      // the same as the distribution from in_vector.

      // Create the distribution from out_distribution.
      Vector<LinearAlgebraDistribution*> tmp_out_distribution_pt(nvec, 0);
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        tmp_out_distribution_pt[vec_i] =
          out_vector_pt[vec_i]->distribution_pt();
      }

      LinearAlgebraDistribution tmp_distribution;
      LinearAlgebraDistributionHelpers::concatenate(tmp_out_distribution_pt,
                                                    tmp_distribution);
      // Compare the distributions
      if (tmp_distribution != *(in_vector.distribution_pt()))
      {
        std::ostringstream error_message;
        error_message << "The distribution from the in vector is incorrect.\n"
                      << "It must be a concatenation of all the distributions\n"
                      << "from the out vectors.\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      // Clear the distribution.
      tmp_distribution.clear();

      // Check that all communicators are the same. We use a communicator to
      // get the number of processors and my_rank. So we would like them to be
      // the same for the in vector and all the out vectors.
      const OomphCommunicator in_vector_comm =
        *(in_vector.distribution_pt()->communicator_pt());
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        const OomphCommunicator vec_i_comm =
          *(out_vector_pt[vec_i]->distribution_pt()->communicator_pt());

        if (in_vector_comm != vec_i_comm)
        {
          std::ostringstream error_message;
          error_message << "The communicator for the vector in position\n"
                        << vec_i << " is not the same as the in_vector\n"
                        << "communicator.";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Check that if the in vector is distributed, then all the out vectors
      // are also distributed.
      if (in_vector_comm.nproc() > 1)
      {
        bool in_distributed = in_vector.distributed();
        for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
        {
          if (in_distributed != out_vector_pt[vec_i]->distributed())
          {
            std::ostringstream error_message;
            error_message << "The vector in position " << vec_i
                          << " does not have\n"
                          << "the same distributed boolean as the in vector";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
      }
#endif

      // Loop through the sub vectors and insert the values from the
      // in vector.
      double* in_value_pt = in_vector.values_pt();
      unsigned in_value_offset = 0;
      for (unsigned vec_i = 0; vec_i < nvec; vec_i++)
      {
        // The nrow_local and values for the current out vector.
        unsigned out_nrow_local = out_vector_pt[vec_i]->nrow_local();
        double* out_value_pt = out_vector_pt[vec_i]->values_pt();

        // Loop through the local values of out vector.
        for (unsigned val_i = 0; val_i < out_nrow_local; val_i++)
        {
          out_value_pt[val_i] = in_value_pt[in_value_offset + val_i];
        }

        // Update the offset.
        in_value_offset += out_nrow_local;
      }
    } // function split_distribution_vector

    //===========================================================================
    /// Wrapper around the other split_without_communication(...)
    /// function.
    /// Be careful with Vector of vectors. If the DoubleVectors are resized,
    /// there could be reallocation of memory. If we wanted to use the function
    /// which takes a Vector of pointers to DoubleVectors, we would either have
    /// to invoke new and remember to delete, or create a temporary Vector to
    /// store pointers to the DoubleVector objects.
    /// This wrapper is meant to make life easier for the user by avoiding calls
    /// to new/delete AND without creating a temporary vector of pointers to
    /// DoubleVectors.
    /// If we had C++ 11, this would be so much nicer since we can use smart
    /// pointers which will delete themselves, so we do not have to remember
    /// to delete!
    //===========================================================================
    void split_without_communication(const DoubleVector& in_vector,
                                     Vector<DoubleVector>& out_vector)
    {
      const unsigned n_out_vector = out_vector.size();

      Vector<DoubleVector*> out_vector_pt(n_out_vector, 0);

      for (unsigned i = 0; i < n_out_vector; i++)
      {
        out_vector_pt[i] = &out_vector[i];
      }

      DoubleVectorHelpers::split_without_communication(in_vector,
                                                       out_vector_pt);

    } // function split_distribution_vector

  } // namespace DoubleVectorHelpers

} // namespace oomph
