//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
class SuperLUPreconditioner : public Preconditioner
{
 public:

 /// Constructor.
 SuperLUPreconditioner()
  {
   Solver.disable_doc_stats();
   Solver.disable_doc_time();
  }
 
 /// Destructor.
 ~SuperLUPreconditioner()
  {}
 
  /// Broken copy constructor.
  SuperLUPreconditioner(const SuperLUPreconditioner&)
  {
   BrokenCopy::broken_copy("SuperLUPreconditioner");
  }


  /// Broken assignment operator.
  void operator=(const SuperLUPreconditioner&)
  {
   BrokenCopy::broken_assign("SuperLUPreconditioner");
  }

  /// \short Function to set up a preconditioner for the linear
  /// system defined by matrix_pt. This function must be called
  /// before using preconditioner_solve.
  /// Note: matrix_pt must point to an object of class
  /// CRDoubleMatrix or CCDoubleMatrix
 void setup()
   {
    oomph_info << "Setting up SuperLU (exact) preconditioner" 
               << std::endl;
    if (dynamic_cast<DistributableLinearAlgebraObject*>(matrix_pt()) != 0)
     {
      LinearAlgebraDistribution dist
       ((dynamic_cast<DistributableLinearAlgebraObject*>
         (matrix_pt()))->distribution_pt());
      this->build_distribution(dist);
      Solver.factorise(matrix_pt());
     }
    else
     {
      std::ostringstream error_message_stream;                         
      error_message_stream                                        
       << "SuperLUPreconditioner can only be applied to matrices derived \n"
       << "DistributableLinearAlgebraObject.\n"
       << "You are most likely to be here because you are using the\n "
       << "soon to be obsolete CCDoubleMatrix\n";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,             
                         OOMPH_EXCEPTION_LOCATION);        
     }
   }
  
  /// \short Function applies SuperLU to vector r for (exact) preconditioning,
  /// this requires a call to setup(...) first.
  void preconditioner_solve(const DoubleVector&r,
                            DoubleVector &z)
   {
    Solver.resolve(r, z);
   }
  

  /// \short Clean up memory -- forward the call to the version in
  /// SuperLU in its LinearSolver incarnation.
  virtual void clean_up_memory()
   {
    Solver.clean_up_memory();
   }

  private:

  /// \short the SuperLU solver emplyed by this preconditioner
  SuperLUSolver Solver;
};

}
#endif
