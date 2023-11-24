//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver for a simple 1D reaction-diffusion problem 
// with adaptive mesh refinement

// Generic oomph-lib routines
#include "generic.h"
//Advection diffusion reaction
#include "advection_diffusion_reaction.h"
// The mesh
#include "meshes/one_d_mesh.h"

//Standard C++ namespace
using namespace std;
//The oomph-lib namespace
using namespace oomph;


//Define Global variables in a namespace to keep things neat
namespace GlobalVariables
{
 //The two timescale parameters (both set to one initially)
 Vector<double> Tau(2,1.0);
 
 //The two diffusion parameters (both set to one initially)
 Vector<double> D(2,1.0);

 //The bifurcation parameter (source of inhibitor)
 double A = -0.1;

 //Simple reaction kinetics
 void activator_inhibitor_reaction(const Vector<double> &C, Vector<double> &R)
 {
  //Inhibitor loss is linearly proportional to concentrations of activator
  //and inhibitor
  R[0] = C[0] + C[1] - A;
  //Activator growth is linearly proportional to activator and inhibitor
  //and loss is self-catalysed at a cubic rate
  R[1] = C[1]*C[1]*C[1] - C[1] - C[0];
 }

 /// Derivative of simple reaction kinetics above
 void activator_inhibitor_reaction_derivative(const Vector<double> &C,
                                 DenseMatrix<double> &dRdC)
 {
  dRdC(0,0) = 1.0; dRdC(0,1) = 1.0;
  dRdC(1,0) = -1.0; dRdC(1,1) = 3.0*C[1]*C[1] - 1.0;
 }

}


//======start_of_problem_class============================================
/// 1D AdvectionDiffusionReaction problem discretised with refineable 
/// 1D QAdvectionDiffusionReaction elements.
/// The specific type of element is specified via the template parameter.
/// (The bit in the angle brackets)
//========================================================================
template<class ELEMENT> 
class RefineableOneDAdvectionDiffusionReactionProblem : public Problem
{

public:

 /// Constructor. No arguments
 RefineableOneDAdvectionDiffusionReactionProblem();
 
 /// Destructor (empty)
 ~RefineableOneDAdvectionDiffusionReactionProblem() {}

 /// Set the initial condition 
 void set_initial_condition();

 /// Perform nstep timesteps of size dt
 void timestep(const double &dt, const unsigned &nstep);

 /// Overloaded Problem's access function to the mesh.
 /// Recasts the pointer to the base Mesh object to the actual mesh type.
 /// This is required so that we can call specific RefineableMesh functions
 RefineableOneDMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableOneDMesh<ELEMENT>*>(Problem::mesh_pt());
  }

private:

 /// Internal storage for the timestep
 double Dt;
 
}; // End of problem class




//=====start_of_constructor===============================================
/// Constructor for AdvectionDiffusionReaction problem:
//========================================================================
template<class ELEMENT>
RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
RefineableOneDAdvectionDiffusionReactionProblem()
{ 
 //Allocate the timestepper (second order implicit)
 add_time_stepper_pt(new BDF<2>); 
 
 // Set up the mesh

 // Number of elements initially
 const unsigned n = 2;

 // Domain length
 const double length = 1.0;

 // Build and assign the refineable mesh, need to pass in number of
 // elements, length and the timestepper
 Problem::mesh_pt() = 
  new RefineableOneDMesh<ELEMENT>(n,length,Problem::time_stepper_pt());

 // Create/set error estimator (default)
 mesh_pt()->spatial_error_estimator_pt() = new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem. 
 // Make the domain periodic by setting the values at the left-hand boundary
 // equal to those on the right
 mesh_pt()->boundary_node_pt(0,0)
  ->make_periodic(mesh_pt()->boundary_node_pt(1,0));

 // Loop over the elements to set up element-specific things that cannot
 // be handled by the (argument-free!) ELEMENT constructor: Pass pointer
 // to source function
 const unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the timescales
   elem_pt->tau_pt() = &GlobalVariables::Tau;

   //Set the diffusion coefficients
   elem_pt->diff_pt() = &GlobalVariables::D;

   //Set the reaction terms
   elem_pt->reaction_fct_pt() = &GlobalVariables::activator_inhibitor_reaction;

    //And their derivatives
   elem_pt->reaction_deriv_fct_pt() = 
    &GlobalVariables::activator_inhibitor_reaction_derivative;
  }

 // Set up equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

} // End of constructor


//=====================================================================
/// Set the initial conditions to be a single "hot" spot
//=====================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::
set_initial_condition()
{
 //Set the number of hot spots
 const unsigned n_spot = 1;
 //Set the centre of the hot spot
 const double centre_x[n_spot] = {0.5};
 
 //Set the initial concentrations of the reagents --- the homogeneous state
 //Plus an exponential hot spot
 unsigned n_node = mesh_pt()->nnode();
 //Loop over the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Local pointer to the node
   Node* nod_pt = mesh_pt()->node_pt(n);
   //This is |A|^{1/3}
   double a13 = pow(std::abs(GlobalVariables::A),(1.0/3.0));

   //Set the value of the inhibitor from the homogeneous solution
   //This comes from an analytic solution of the equations
   nod_pt->set_value(0,a13*(1.0-a13*a13));
   
   //Set a localised hot spot by making an exponential
   //hump in the activator
   double x = nod_pt->x(0);    
   double spot = 0.0;
   //Loop over the number of spots
   for(unsigned s=0; s<n_spot; s++)
    {
     //Find the square distance from the centre of the 
     //hot spot
     double r2 = (x - centre_x[s])*(x - centre_x[s]);
     
     //Add an exponential hump to the value of the spot variable
     //The "10" in the exponential sets the steepness of the spot
     spot += 2.0*exp(-10.0*r2);
    }
   
   //Set the value of the activator
   nod_pt->set_value(1,-a13 + spot);
  }

 //Document the initial solution
 ofstream filename("initial.dat");
 //Plot the solution with 5 points per element
 mesh_pt()->output(filename,5);
 filename.close();
 
 //Set the initial values impulsive 
 //i.e. assume that the solution has been at the initial condition for all
 //previous times
 assign_initial_values_impulsive(Dt);
} 


//====================================================================
/// Timestep the problem for nstep timesteps of length dt
//===================================================================
template<class ELEMENT>
void RefineableOneDAdvectionDiffusionReactionProblem<ELEMENT>::timestep(
 const double &dt, const unsigned &nstep)
{
 //Set the problem's Dt to use in the inital condition
 Dt = dt;

 //Maximum adaptation for the first timestep
 unsigned max_adapt = 2;

 //Take the first timestep
 bool first = true;
 
 //Set the initial condition
 set_initial_condition();

 //Solve the first step (you need the additional flag first so that
 //the initial conditions are reset when you adapt in space)
 unsteady_newton_solve(dt,max_adapt,first);
 //Output the result
 {
  unsigned i=0;
  char file1[100];
  sprintf(file1,"step%i.dat",i+1);
  ofstream out1(file1);
  mesh_pt()->output(out1,5);
  out1.close();
 }

 //Now set so that only one round of adaptation is performed each timestep
 max_adapt = 1;
 //This is not the first timestep, so we shouldn't use the initial conditions
 first = false;

 //Loop over timesteps
 for(unsigned i=1;i<nstep;i++)
  {
   //Take a timestep
   unsteady_newton_solve(dt,max_adapt,first);
   //Output the result
   char file1[100];
   sprintf(file1,"step%i.dat",i+1);
   ofstream out1(file1);
   mesh_pt()->output(out1,5);
   out1.close();
  }
}




//======start_of_main=====================================================
/// Driver code for 1D AdvectionDiffusionReaction problem
//========================================================================
int main()
{
 //Diffusive length-scale
 double epsilon = 0.05;
 //Set the control parameters
 GlobalVariables::Tau[1] = 0.2;
 GlobalVariables::D[1] = epsilon*epsilon;
 GlobalVariables::A = -0.4;
 
 //Set the (largeish) timestep
 double dt = 0.1;

 //Set up the problem
 //------------------
 // DREIGIAU: There's an inherent problem with the 1D elements. Doesn't
 // seem to be a problem with 2D elements but we can't even initialise
 // a new object of the type used to template the problem below.
 
 // Create the problem with 1D three-node refineable elements from the
 // RefineableLineAdvectionDiffusionReactionElement family. 
 RefineableOneDAdvectionDiffusionReactionProblem<
  RefineableQAdvectionDiffusionReactionElement<2,1,3> > problem;

 //Take four levels of uniform refinement to start things off
 for(unsigned i=0;i<4;i++) { problem.refine_uniformly(); }

 //Now timestep the problem
 problem.timestep(dt,10);
 
} // End of main









