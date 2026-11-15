//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

// Oomph-lib includes
#include "generic.h"

using namespace oomph;
using namespace std;



//======================================================================
/// Dummy class containing mutable workspaces
//======================================================================
class MyClass
{
public:

 /// Constructor.
 /// The size of Vector_workspace's content is established here and only here
 MyClass() : Vector_workspace(2),
             Matrix_workspace(2,2),
             Rank_three_tensor_workspace(2,2,2), 
             Rank_four_tensor_workspace(2,2,2,2)
  {}


 /// Function that uses the workspace
 void do_stuff() const
  {
   // Get a workspace.
   //
   // In a PARANOID build this records that associated GuardedWorkspaces
   // are now in active use.
   Workspace<Vector<double>> v_workspace = Vector_workspace.workspace();
   v_workspace.content()[0] = 1.0;
   v_workspace.content()[1] = 2.0;

   Workspace<DenseMatrix<double>> m_workspace = Matrix_workspace.workspace();
   m_workspace.content()(0,0) = 1.0;
   
   Workspace<RankThreeTensor<double>> r3_workspace =
    Rank_three_tensor_workspace.workspace();
   r3_workspace.content()(0,0,0) = 1.0;
   
   Workspace<RankFourTensor<double>> r4_workspace =
    Rank_four_tensor_workspace.workspace();
   r4_workspace.content()(0,0,0,0) = 1.0;


   // Perform whatever operations are required with the content.
   // [...]

   
   // There is no explicit release operation. When workspaces go out
   // of scope, their destructor records that the relevant GuardedWorkspaces
   // are no longer in use.
  }

 
 /// Function that uses the workspace and then calls another function
 /// that attempts to use the same workspace
 void do_stuff_again() const
  {
   Workspace<Vector<double>> v_workspace = Vector_workspace.workspace();
   v_workspace.content()[0] = 3.0;
   v_workspace.content()[1] = 4.0;

   Workspace<DenseMatrix<double>> m_workspace = Matrix_workspace.workspace();
   m_workspace.content()(0,0) = 5.0;
   
   Workspace<RankThreeTensor<double>> r3_workspace =
    Rank_three_tensor_workspace.workspace();
   r3_workspace.content()(0,0,0) = 6.0;
   
   Workspace<RankFourTensor<double>> r4_workspace =
    Rank_four_tensor_workspace.workspace();
   r4_workspace.content()(0,0,0,0) = 7.0;

   
   // In a PARANOID build, do_stuff() attempts to create a second
   // Workspace referring to Vector_workspace while the Workspace above
   // is still alive. The overlapping use is therefore detected.
   do_stuff();
  }
 

protected:
 
 /// Mutable workspace containing a dynamically sized vector whose
 /// size is set (and then remains) equal to two when the
 /// GuardedWorkspace is constructed in the constructor of the
 /// current class
 mutable GuardedWorkspace<Vector<double>> Vector_workspace;
 
 /// Mutable workspace containing a dynamically sized matrix whose
 /// size is set (and then remains) equal to 2x2 when the
 /// GuardedWorkspace is constructed in the constructor of the
 /// current class
 mutable GuardedWorkspace<DenseMatrix<double>> Matrix_workspace;
 
 /// Mutable workspace containing a dynamically sized rank 3 tensor whose
 /// size is set (and then remains) equal to 2x2x2 when the
 /// GuardedWorkspace is constructed in the constructor of the
 /// current class
 mutable GuardedWorkspace<RankThreeTensor<double>> Rank_three_tensor_workspace;
 
 /// Mutable workspace containing a dynamically sized rank 4 tensor whose
 /// size is set (and then remains) equal to 2x2x2x2 when the
 /// GuardedWorkspace is constructed in the constructor of the
 /// current class
 mutable GuardedWorkspace<RankFourTensor<double>> Rank_four_tensor_workspace;
 
};




//======================================================================
/// Driver
//======================================================================
int main()
{
  MyClass my_object;

  oomph_info << "About to do stuff" << std::endl;
  my_object.do_stuff();
  oomph_info << "Done stuff" << std::endl;


  // Clever checking of overlapping use of workspace only
  // happens in paranoid mode!
#ifdef PARANOID
  
  oomph_info << "About to do stuff again" << std::endl;
  bool expected_exception_was_caught = false;
  try
  {
    my_object.do_stuff_again();
  }
  catch (const OomphLibError& exception)
  {
   expected_exception_was_caught = true;   
   oomph_info << "Caught expected exception: "
              << exception.what() << std::endl;
  }
  if (!expected_exception_was_caught)
   {
    throw OomphLibError(
     "Overlapping workspace use was not detected.",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
  oomph_info << "Done stuff again" << std::endl;

  // Demonstrate that exception handling correctly released the workspace.
  //
  // If stack unwinding had failed to clear the In_use state, this call
  // would incorrectly report that the workspace was still in use.
  oomph_info << "About to use workspace after exception" << std::endl;

  my_object.do_stuff();
  oomph_info << "Successfully used workspace after exception" << std::endl;
  
#else

  // In a non-PARANOID build the diagnostic machinery is absent.
  my_object.do_stuff_again();
 

#endif



  // Output celebratory message if we get here. (Don't edit -- this is our
  // validata!)
  oomph_info << "Reached the end! That all worked well." << std::endl;
}
