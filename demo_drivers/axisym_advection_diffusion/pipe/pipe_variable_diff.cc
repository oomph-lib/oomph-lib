//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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

 /// Peclet number multiplied by Strouhal number
 double Peclet_St = 1.0;

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

 //Conserved bit (the swimming velocity)
 //You can include three velocity components but only two will
 //matter in the axisymmetric advection-diffusion equations
 void swimming(const Vector<double> &x, Vector<double> &swim)
 {
  //Radial direction
  swim[0] = 0.0;
  //Vertical direction
  swim[1] = -5*x[1];
  //Azimuthal direction (ignored)
  swim[2] = 0.0;
 }

 /// Diffusivity tensor (again 3x3 but only the upper 2x2 block affects
 /// the axisymmetric advection diffusion equations)
 void diff_function(const Vector<double> &x, DenseMatrix<double> &D)
 {
  //Radial component
  D(0,0) = 0.5;
  //Radial-Axial component
  D(0,1) = 0.0;
  //Axial-Radial component
  D(1,0) = 0.0;
  //Axial component
  D(1,1) = 1.0;

  //These will not affect the equation
  D(0,2) = 0.0;
  D(1,2) = 0.0;
  D(2,0) = 0.0;
  D(2,1) = 0.0;
  D(2,2) = 0.0;
 }

 
} // end of namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////

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

 /// Set the initial state of the system
 void set_initial_condition();

 /// Doc the solution.
 void doc_solution();

 /// Overloaded version of the problem's access function to 
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
/// Constructor for AdvectionDiffusion problem
//========================================================================
template<class ELEMENT>
RefineableAdvectionDiffusionPipeProblem<ELEMENT>::
RefineableAdvectionDiffusionPipeProblem()
{ 
 //Add the timestepper
 this->add_time_stepper_pt(new BDF<2>);

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
  new RefineableRectangularQuadMesh<ELEMENT>(n_r,n_z,l_r,l_z,
                                             this->time_stepper_pt());

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  // Make the mesh periodic in the z-direction by setting the nodes on
 // right boundary (boundary 0) to be the periodic counterparts of
 // those on the left one (boundary 2).
 unsigned n_node = mesh_pt()->nboundary_node(0);
 for(unsigned n=0;n<n_node;n++)
  {
   mesh_pt()->boundary_node_pt(0,n)
    ->make_periodic(mesh_pt()->boundary_node_pt(2,n));
  }
 
 // Now establish the new neighbours (created by "wrapping around"
 // the domain) in the TreeForst representation of the mesh

 // Get pointers to tree roots associated with elements on the 
 // top and bottom  boundaries
 Vector<TreeRoot*> top_root_pt(n_r);
 Vector<TreeRoot*> bottom_root_pt(n_r);
 for(unsigned i=0;i<n_r;i++) 
  {
   top_root_pt[i] = 
    dynamic_cast<ELEMENT*>(mesh_pt()->element_pt((n_z-1)*n_r + i))->
    tree_pt()->root_pt();
   bottom_root_pt[i] = 
    dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i))->
    tree_pt()->root_pt();
  }

  // Switch on QuadTreeNames for enumeration of directions
   using namespace QuadTreeNames;

  //Set the neighbour and periodicity
  for(unsigned i=0;i<n_r;i++) 
   {
    // The northern neighbours of the elements on the top
    // boundary are those on the bottom
    top_root_pt[i]->neighbour_pt(N) = bottom_root_pt[i];
    top_root_pt[i]->set_neighbour_periodic(N); 
    
    // The southern neighbours of the elements on the bottom
    // boundary are those on the top
    bottom_root_pt[i]->neighbour_pt(S) = top_root_pt[i];
    bottom_root_pt[i]->set_neighbour_periodic(S);     
   } // done
  
 //No boundary conditions because we have periodicity in the axial direction
 //and no-flux in the radial direction
  
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

   //Set the swimming (conserved wind) 
   //el_pt->conserved_wind_fct_pt() = Global_Physical_Variables::swimming;

   //Set the diffusivity
   el_pt->diff_fct_pt() = Global_Physical_Variables::diff_function;

   //Set the Peclet Strouhal number
   el_pt->pe_st_pt() = &Global_Physical_Variables::Peclet_St;

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=============================start_of_set_concentration_profile=======
/// Set a specified inlet concentration profile
//========================================================================
template<class ELEMENT>
void RefineableAdvectionDiffusionPipeProblem<ELEMENT>::
set_initial_condition()
{
 //Lets have a little blob of swimmers along the axis
 unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   Node* nod_pt = mesh_pt()->node_pt(n);
   double r = nod_pt->x(0);
   double z = nod_pt->x(1);
   
   //Exponentially decaying in two directions
   double c = exp(-10*r*r
                  -2.0*(z-0.5*Global_Physical_Variables::Length)*
                  (z-0.5*Global_Physical_Variables::Length));
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

 //Increase the number
 ++Doc_info.number();

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
 RefineableQGeneralisedAxisymAdvectionDiffusionElement<3> > problem;

 //Refine uniformly once
 problem.refine_uniformly();
 //Set the initial distribution
 problem.set_initial_condition();
 problem.doc_solution();

 //Set the timestep
 double dt = 0.001;
 problem.assign_initial_values_impulsive(dt);

 //Now take a few timesteps
 bool first = true;
 unsigned max_adapt=2;
 for(unsigned t=0;t<3;t++)
  {
   problem.unsteady_newton_solve(dt,max_adapt,first);
   problem.doc_solution();
   if(t==0) {first=false;}
  }
 
} // end of main









