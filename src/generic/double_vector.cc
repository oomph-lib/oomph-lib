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
     this->build(old_vector.distribution_pt(),0.0);

     // copy the data
     if (this->distribution_built())
      {
       unsigned nrow_local = this->nrow_local();
       const double* old_vector_values = old_vector.values_pt();
       for (unsigned i = 0; i < nrow_local; i++)
        {
         Values_pt[i] = old_vector_values[i];
        }
      }
    }
  }

 //============================================================================
 /// Assembles a DoubleVector with distribution dist, if v is specified
 /// each row is set to v
 //============================================================================
 void DoubleVector::build(const LinearAlgebraDistribution* const &dist_pt,
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
     unsigned nrow = this->nrow_local();
     Values_pt = new double[nrow];

     // set the data
     for (unsigned i = 0; i < nrow; i++)
      {
       Values_pt[i] = v;
      }
     Built=true;
    }
   else
    {
     Built=false;
    }
  }

 //============================================================================
 /// \short Assembles a DoubleVector with a distribution dist and coefficients
 /// taken from the vector v.\n
 /// Note. The vector v MUST be of length nrow()
 //============================================================================
 void DoubleVector::build(const LinearAlgebraDistribution* const &dist_pt,
                          const Vector<double>& v)
 {
  // clean the memory
   this->clear();

   // the vector owns the internal data
   Internal_values = true;

   // Set the distribution
   this->build_distribution(dist_pt);

   // use the initialise method to populate the vector
   this->initialise(v);

   // indicate that its built
   Built=true;
 }

 //============================================================================
 /// \short initialise the whole vector with value v
 //============================================================================
 void DoubleVector::initialise(const double& v)
  {
   if (Built)
    {
     // cache nrow local
     unsigned nrow_local = this->nrow_local();

     // set the residuals
     for (unsigned i = 0; i < nrow_local; i++)
      {
       Values_pt[i] = v;
      }
    }
  }

 //============================================================================
 /// \short initialise the vector with coefficient from the vector v.\n
 /// Note: The vector v must be of length
 //============================================================================
 void DoubleVector::initialise(const Vector<double> v)
 {
#ifdef PARANOID
  if (v.size()!=this->nrow())
   {
    std::ostringstream error_message;
    error_message << "The vector passed to initialise(...) must be of length "
                  << "nrow()";
    throw OomphLibError(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  unsigned begin = this->first_row();
  unsigned end = begin + this->nrow_local();
  for (unsigned i = begin; i < end; i++)
   {
    Values_pt[i-begin] = v[i];
   }
 }

 //============================================================================
 /// The contents of the vector are redistributed to match the new
 /// distribution. In a non-MPI build this method works, but does nothing. \n
 /// \b NOTE 1: The current distribution and the new distribution must have
 /// the same number of global rows.\n
 /// \b NOTE 2: The current distribution and the new distribution must have
 /// the same Communicator.
 //============================================================================
 void DoubleVector::redistribute(const LinearAlgebraDistribution*
                                 const& dist_pt)
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // paranoid check that the current distribution and the new distribution
   // have the same Communicator
   OomphCommunicator temp_comm(*dist_pt->communicator_pt());
   if (!(temp_comm == *this->distribution_pt()->communicator_pt()))
    {
     std::ostringstream error_message;
     error_message << "The new distribution and the current distribution must "
                   << "have the same communicator.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
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
             (current_first_row_data[my_rank] < (new_first_row_data[p] +
                                                 new_nrow_local_data[p])))
         {
          new_first_row_for_proc[p] =
           std::max(current_first_row_data[my_rank],
                    new_first_row_data[p]);
          new_nrow_local_for_proc[p] =
           std::min((current_first_row_data[my_rank] +
                     current_nrow_local_data[my_rank]),
                    (new_first_row_data[p] +
                     new_nrow_local_data[p])) - new_first_row_for_proc[p];
         }

         // and data to be received
         if ((new_first_row_data[my_rank] < (current_first_row_data[p] +
                                                 current_nrow_local_data[p]))
             && (current_first_row_data[p] < (new_first_row_data[my_rank] +
                                              new_nrow_local_data[my_rank])))
         {
          new_first_row_from_proc[p] =
           std::max(current_first_row_data[p],
                    new_first_row_data[my_rank]);
          new_nrow_local_from_proc[p] =
           std::min((current_first_row_data[p] +
                     current_nrow_local_data[p]),
                    (new_first_row_data[my_rank] +
                     new_nrow_local_data[my_rank]))-new_first_row_from_proc[p];
         }
        }

       // temporary storage for the new data
       double* temp_data = new double[new_nrow_local_data[my_rank]];

       // "send to self" or copy Data that does not need to be sent else where
       // to temp_data
       if (new_nrow_local_for_proc[my_rank] != 0)
        {
         unsigned j = new_first_row_for_proc[my_rank] -
          current_first_row_data[my_rank];
         unsigned k = new_first_row_for_proc[my_rank] -
          new_first_row_data[my_rank];
         for (unsigned i = 0; i < new_nrow_local_for_proc[my_rank]; i++)
          {
           temp_data[k+i] = Values_pt[j+i];
          }
        }

       // send and receive circularly
       for (int p = 1; p < nproc; p++)
        {
         // next processor to send to
         unsigned dest_p = (my_rank + p)%nproc;

         // next processor to receive from
         unsigned source_p = (nproc + my_rank - p)%nproc;

         // send and receive the value
         MPI_Status status;
         MPI_Sendrecv(Values_pt + new_first_row_for_proc[dest_p] -
                      current_first_row_data[my_rank],
                      new_nrow_local_for_proc[dest_p],MPI_DOUBLE,dest_p,1,
                      temp_data + new_first_row_from_proc[source_p] -
                      new_first_row_data[my_rank],
                      new_nrow_local_from_proc[source_p],MPI_DOUBLE,source_p,1,
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
       int* dist_nrow_local =  new int[nproc];
       for (int p = 0; p < nproc; p++)
        {
         dist_first_row[p] = this->first_row(p);
         dist_nrow_local[p] = this->nrow_local(p);
        }

       // gather the local vectors from all processors on all processors
       int my_nrow_local(this->nrow_local());
       MPI_Allgatherv(temp_data,my_nrow_local,MPI_DOUBLE,
                      Values_pt,dist_nrow_local,dist_first_row,MPI_DOUBLE,
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

       //copy to Values_pt
       delete[] Values_pt;
       Values_pt= temp_data;

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
     error_message << "Range Error: " << i
                   << " is not in the range (0,"
                   << this->nrow_local()-1 << ")";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   return Values_pt[i];
  }

 //============================================================================
 /// \short == operator
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
 /// \short += operator
 //============================================================================
 void DoubleVector::operator+=(const DoubleVector& v)
  {
#ifdef PARANOID
   // PARANOID check that this vector is setup
   if (!this->built())
    {
     std::ostringstream error_message;
     error_message << "This vector must be setup.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // PARANOID check that the vector v is setup
   if (!v.built())
    {
     std::ostringstream error_message;
     error_message << "The vector v must be setup.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // PARANOID check that the vectors have the same distribution
   if (!(*v.distribution_pt() == *this->distribution_pt()))
    {
     std::ostringstream error_message;
     error_message << "The vector v and this vector must have the same "
                   << "distribution.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif//

   // cache nrow_local
   double* v_values_pt = v.values_pt();
   unsigned nrow_local = this->nrow_local();
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // PARANOID check that the vector v is setup
   if (!v.built())
    {
     std::ostringstream error_message;
     error_message << "The vector v must be setup.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   // PARANOID check that the vectors have the same distribution
   if (!(*v.distribution_pt() == *this->distribution_pt()))
    {
     std::ostringstream error_message;
     error_message << "The vector v and this vector must have the same "
                   << "distribution.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // cache nrow_local
   double* v_values_pt = v.values_pt();
   unsigned nrow_local = this->nrow_local();
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
  for(unsigned i=0, ni=this->nrow_local(); i<ni; i++)
   {
    Values_pt[i] *= d;
   }
 }

 //============================================================================
 /// Divide by double
 //============================================================================
 void DoubleVector::operator/=(const double& d)
 {

  for(unsigned i=0, ni=this->nrow_local(); i<ni; i++)
   {
    Values_pt[i] /= d;
   }

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
     error_message << "Range Error: " << i
                   << " is not in the range (0,"
                   << this->nrow_local()-1 << ")";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
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
     MPI_Allreduce(&local_max,&max,1,MPI_DOUBLE,MPI_MAX,
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
 void DoubleVector::output(std::ostream &outfile)
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
     int* dist_nrow_local =  new int[nproc];
     for (int p = 0; p < nproc; p++)
      {
       dist_first_row[p] = this->first_row(p);
       dist_nrow_local[p] = this->nrow_local(p);
      }

     // gather
     temp = new double[nrow];
     MPI_Allgatherv(Values_pt,nrow_local,MPI_DOUBLE,
                    temp,dist_nrow_local,dist_first_row,MPI_DOUBLE,
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
   for (unsigned i = 0; i < nrow; i++)
    {
     outfile << i << " " << temp[i] << std::endl;
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (!vec.built())
    {
     std::ostringstream error_message;
     error_message << "The input vector be setup.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (*this->distribution_pt() != *vec.distribution_pt())
    {
     std::ostringstream error_message;
     error_message << "The distribution of this vector and the vector vec "
                   << "must be the same."
                   << "\n\n  this: " << *this->distribution_pt()
                   << "\n  vec:  " << *vec.distribution_pt();
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // compute the local norm
   unsigned nrow_local = this->nrow_local();
   double n = 0.0;
   const double* vec_values_pt = vec.values_pt();
   for (unsigned i = 0; i < nrow_local; i++)
    {
     n += Values_pt[i]*vec_values_pt[i];
    }

   // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
   double n2 = n;
   if (this->distributed() &&
       this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
     MPI_Allreduce(&n,&n2,1,MPI_DOUBLE,MPI_SUM,
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // compute the local norm
   unsigned nrow_local = this->nrow_local();
   double n = 0;
   for (unsigned i = 0; i < nrow_local; i++)
    {
     n += Values_pt[i]*Values_pt[i];
    }

   // if this vector is distributed and on multiple processors then gather
#ifdef OOMPH_HAS_MPI
   double n2 = n;
   if (this->distributed() &&
       this->distribution_pt()->communicator_pt()->nproc() > 1)
    {
     MPI_Allreduce(&n,&n2,1,MPI_DOUBLE,MPI_SUM,
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
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (!matrix_pt->built())
    {
     std::ostringstream error_message;
     error_message << "The input matrix be built.";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (*this->distribution_pt() != *matrix_pt->distribution_pt())
    {
     std::ostringstream error_message;
     error_message << "The distribution of this vector and the matrix at "
                   << "matrix_pt must be the same";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // compute the matrix norm
   DoubleVector x(this->distribution_pt(),0.0);
   matrix_pt->multiply(*this,x);
   return sqrt(this->dot(x));
  }


 /// \short output operator
 std::ostream& operator<< (std::ostream &out, const DoubleVector& v)
 {
  // Do the first value outside the loop to get the ", "s right.
  out << "[" << v[0];

  for(unsigned i=1, ni=v.nrow_local(); i<ni; i++)
   {
    out << ", " << v[i];
   }
  out << "]";

  return out;
 }



}//end of oomph namespace
