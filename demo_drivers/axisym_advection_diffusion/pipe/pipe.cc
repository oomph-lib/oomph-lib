//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
//Driver for a simple axisymmetric pipe-flow problem 
//with adaptive mesh refinement

//Generic routines
#include "generic.h"

// The axisymmetric advection-diffusion  equations
#include "axisym_advection_diffusion.h"

// The mesh 
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//======start_of_namespace============================================
/// Namespace for physical variables
//====================================================================
namespace Global_Physical_Variables
{

 /// Peclet number
 double Peclet=200.0;

 /// Length of the pipe
 double Length=10.0;

 /// Radius of the pipe
 double Radius=1.0;

 /// Wind
 void wind_function(const Vector<double>& x, Vector<double>& wind)
 {
  //Poiseulle flow no radial flow, but quadratic axial flow
  //plus a possible swirl component which will do nothing
  //because the concentration cannot have gradients in the swirl 
  //direction
  wind[0] = 0.0;
  wind[1] = (Radius*Radius - x[0]*x[0]);
  wind[2] = 0.0;
 }
 
} // end of namespace

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// AdvectionDiffusion problem on rectangular domain, discretised 
/// with refineable Axisymmetric QAdvectionDiffusion elements. 
/// The specific type of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableAdvectionDiffusionPipeProblem : public Problem
{

public:

 /// Constructor: 
 RefineableAdvectionDiffusionPipeProblem();

 /// Destructor. Empty
 ~RefineableAdvectionDiffusionPipeProblem(){}

 ///  Set the inlet concentration condition
 void set_inlet_concentration();

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the tanh solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt: Document the solution
 void actions_before_adapt()
  {
   // Doc the solution
   doc_solution();
   
   // Increment label for output files
   Doc_info.number()++;
  }

 ///  Doc the solution.
 void doc_solution();

 ///  Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

private:

 /// DocInfo object
 DocInfo Doc_info;

}; // end of problem class



//=====start_of_constructor===============================================
///  Constructor for AdvectionDiffusion problem
//========================================================================
template<class ELEMENT>
RefineableAdvectionDiffusionPipeProblem<ELEMENT>::
RefineableAdvectionDiffusionPipeProblem()
{ 

 // Set output directory
 Doc_info.set_directory("RESLT");

 // Setup mesh

 // # of elements in r-direction
 unsigned n_r=4;

 // # of elements in z-direction
 unsigned n_z=4;

 // Domain length in r-direction
 double l_r=Global_Physical_Variables::Radius;

 // Domain length in z-direction
 double l_z=Global_Physical_Variables::Length;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_r,n_z,l_r,l_z);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here
 // Assume specified inlet concentration (boundary 0)
 {
  unsigned b=0;
  unsigned num_nod= mesh_pt()->nboundary_node(b);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    mesh_pt()->boundary_node_pt(b,inod)->pin(0); 
   }
 } 
 
 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the wind function pointer
   el_pt->wind_fct_pt() = Global_Physical_Variables::wind_function;

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//=============================start_of_actions_before_newton_solve=======
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the tanh solution.
//========================================================================
template<class ELEMENT>
void RefineableAdvectionDiffusionPipeProblem<ELEMENT>::
actions_before_newton_solve()
{
 this->set_inlet_concentration();
}  // end of actions before solve


//=============================start_of_set_concentration_profile=======
///Set a specified inlet concentration profile
//========================================================================
template<class ELEMENT>
void RefineableAdvectionDiffusionPipeProblem<ELEMENT>::
set_inlet_concentration()
{
 unsigned b=0;
 unsigned n_node = mesh_pt()->nboundary_node(b);
 for(unsigned n=0;n<n_node;n++)
  {
   Node* nod_pt = this->mesh_pt()->boundary_node_pt(b,n);
   //Get the radial value
   double r = nod_pt->x(0);
   
   //Now chose an exponetially decaying profile
   double c = exp(-10.0*r*r);

   nod_pt->set_value(0,c);
  }
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableAdvectionDiffusionPipeProblem<ELEMENT>::doc_solution()
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc


//===== start_of_main=====================================================
/// Driver code for 2D AdvectionDiffusion problem
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 // Create the problem with axisymmetric nine-node refineable elements from the
 // RefineableQAdvectionDiffusionElement family.
 RefineableAdvectionDiffusionPipeProblem<
 RefineableQAxisymAdvectionDiffusionElement<3> > problem;

 // Solve the problem, performing up to 4 adaptive refinements
 problem.newton_solve(4);

 //Output the solution
 problem.doc_solution();
 
} // end of main









