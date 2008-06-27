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
#ifndef OOMPH_SuperLU_Preconditioner_HEADER
#define OOMPH_SuperLU_Preconditioner_HEADER

//oomph-lib headers

#include "linear_solver.h"
#include "preconditioner.h"

namespace oomph
{

//====================================================================
/// An interface to allow SuperLU to be used as an (exact) Preconditioner
//====================================================================
class SuperLU_Preconditioner : public Preconditioner, public SuperLU
{
 public:

 /// Constructor.
 SuperLU_Preconditioner()
  {}
 
 /// Destructor.
 ~SuperLU_Preconditioner()
  {}
 
  /// Broken copy constructor.
  SuperLU_Preconditioner(const SuperLU_Preconditioner&)
  {
   BrokenCopy::broken_copy("SuperLU_Preconditioner");
  }


  /// Broken assignment operator.
  void operator=(const SuperLU_Preconditioner&)
  {
   BrokenCopy::broken_assign("SuperLU_Preconditioner");
  }

  /// \short Function to set up a preconditioner for the linear
  /// system defined by matrix_pt. This function must be called
  /// before using preconditioner_solve.
  /// Note: matrix_pt must point to an object of class
  /// CRDoubleMatrix or CCDoubleMatrix
  void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
   {
    oomph_info << "Setting up SuperLU (exact) preconditioner" 
               << std::endl;
    factorise(matrix_pt);
   }
  
  /// \short Function applies SuperLU to vector r for (exact) preconditioning,
  /// this requires a call to setup(...) first.
  void preconditioner_solve(const Vector<double>&r,
                            Vector<double> &z)
   {
    resolve(r, z);
   }
  

  /// \short Clean up memory -- forward the call to the version in
  /// SuperLU in its LinearSolver incarnation.
  virtual void clean_up_memory()
   {
    SuperLU::clean_up_memory();
   }

};

#ifdef OOMPH_HAS_MPI
//====================================================================
/// \short An interface to allow SuperLU_dist to be used as
/// an (exact) Preconditioner for CRDoubleMatrix and
/// DistributedCRDoubleMatrix classes.
//====================================================================
class SuperLU_dist_Preconditioner : public Preconditioner, public SuperLU_dist
{
  public:
 
 /// Constructor.
 SuperLU_dist_Preconditioner()
  {
   Preconditioning_distributed_system = true;
  }
 
 /// Destructor.
 ~SuperLU_dist_Preconditioner()
  {}
 
 /// Broken copy constructor.
 SuperLU_dist_Preconditioner(const SuperLU_dist_Preconditioner&)
  {
   BrokenCopy::broken_copy("SuperLU_dist_Preconditioner");
  }
 
 
 /// Broken assignment operator.
 void operator=(const SuperLU_dist_Preconditioner&)
  {
   BrokenCopy::broken_assign("SuperLU_dist_Preconditioner");
  }
 
 /// \short Function to set up a preconditioner for the linear
 /// system defined by matrix_pt. This function must be called
 /// before using preconditioner_solve.
 /// Note: matrix_pt must point to an object of class
 /// CRDoubleMatrix or CCDoubleMatrix
 void setup(Problem* problem_pt, DoubleMatrixBase* matrix_pt)
  {
   oomph_info << "Setting up SuperLU_dist (exact) preconditioner"
              << std::endl;
   
   // Test the matrix type and record whether we are dealing with a 
   // distributed system or not
   if(dynamic_cast<DistributedCRDoubleMatrix*>(matrix_pt))
    {
     Preconditioning_distributed_system = true;
    }
   else if (dynamic_cast<CRDoubleMatrix*>(matrix_pt))
    {
     Preconditioning_distributed_system = false;
    }
   else
    {
     std::ostringstream error_message_stream;
     error_message_stream << "SuperLU_dist_Preconditioner implemented only "
                          << "for CRDoubleMatrix\n"
                          << "and DistributedCRDoubleMatrix matrices\n";
     throw OomphLibError(error_message_stream.str(),
                         "SuperLU_dist_Preconditioner::setup()",
                         OOMPH_EXCEPTION_LOCATION);
    }
   factorise(matrix_pt);
  }
  
 /// \short Function applies SuperLU_dist to vector r for (exact)
 /// preconditioning, this requires a call to setup(...) first.
 void preconditioner_solve(const Vector<double>&r, Vector<double> &z)
  {
   // Switch off output for preconditioner solve
   bool copy_of_Doc_time = Doc_time;
   Doc_time = false;
   
   // Solve the system
   if (Preconditioning_distributed_system)
    {
     // SuperLU_dist returns the global solution vector but we must
     // return only the 'local' values 
     Vector<double> global_z(ndof());
     resolve(r, global_z);

     // Copy local values from global_z
     const unsigned nrow_local = ndof_local();
     const unsigned first_local_row = first_local_dof();
     z.resize(nrow_local);

     for (unsigned i=0; i<nrow_local; i++)
      {
       z[i] = global_z[i+first_local_row];
      }
    }
   else
    {
      // Call resolve, this returns the global Vector
     resolve(r, z);
    }
   
   // Reset Doc_time
   Doc_time = copy_of_Doc_time;
  }

  /// \short Clean up memory -- forward the call to the version in
  /// SuperLU_dist in its LinearSolver incarnation.
  virtual void clean_up_memory()
   {
    SuperLU_dist::clean_up_memory();
   }

  private:

  /// Flag to record whether we are dealing with a distributed system
  bool Preconditioning_distributed_system;

};
#endif

}
#endif
