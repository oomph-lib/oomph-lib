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
//Driver for a simple 2D reaction diffusioni problem with
//adaptive mesh refinement
//The problem is the "hot spot" expansion given in Muratov and Osipov
//(see two_d_act_inhibit.cc for full reference)
//Starting with a radially symmetric spot gives slightly slower growth
//that with an eccentric one, but the general dynamics are the same.

//Generic routines
#include "generic.h"
#include "advection_diffusion_reaction.h"
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


namespace GlobalVariables
{
 ///The vector of timescales
 Vector<double> Tau(2,1.0);
 
 ///The vector of diffusion coefficients
 Vector<double> D(2,1.0);

 ///Bifurcation (loss) parameter
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

 ///Derivative of simple reaction kinetics
 void activator_inhibitor_reaction_derivative(const Vector<double> &C,
                                 DenseMatrix<double> &dRdC)
 {
  dRdC(0,0) = 1.0; dRdC(0,1) = 1.0;
  dRdC(1,0) = -1.0; dRdC(1,1) = 3.0*C[1]*C[1] - 1.0;
 }

}


//====== start_of_problem_class=======================================
/// 2D ActivatorInhibitor problem on rectangular domain, discretised 
/// with refineable 2D QActivatorInhibitor elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableActivatorInhibitorProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source and wind functions
 RefineableActivatorInhibitorProblem();

 /// Destructor. Empty
 ~RefineableActivatorInhibitorProblem() {}

 /// Empty
 void actions_before_newton_solve() {} 

 // Empty
 void actions_after_newton_solve(){}

 //Set the initial condition
 void set_initial_condition();

 //Set the timestep
 void timestep(const double &dt, const unsigned &nstep);

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

private:

 /// Storage for the timestep
 double Dt;

}; // end of problem class



//=====start_of_constructor===============================================
/// Constructor for ActivatorInhibitor problem
//========================================================================
template<class ELEMENT>
RefineableActivatorInhibitorProblem<ELEMENT>::
RefineableActivatorInhibitorProblem() : Dt(0.1)
{ 
 //Allocate the timestepper 
 add_time_stepper_pt(new BDF<2>);

 // Setup initial mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=10.0;

 // Domain length in y-direction
 double l_y=10.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(
   n_x,n_y,l_x,l_y,Problem::time_stepper_pt());
 
 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 //Totally free (zero-flux) boundary conditions, so we do nothing
 
 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
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

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=====================================================================
/// Set the initial conditions to be a single "hot" spot
//=====================================================================
template<class ELEMENT>
void RefineableActivatorInhibitorProblem<ELEMENT>::set_initial_condition()
{
 //Set the number of hot spots
 const unsigned n_spot = 1;
 //Set the centre of the hot spot
 const double centre_x[n_spot] = {6.0};
 const double centre_y[n_spot] = {7.0};
 
 //Set the initial concentrations of the reagents --- the homogeneous state
 //Plus an exponential hot spot
 unsigned n_node = mesh_pt()->nnode();
 //Loop over the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Local pointer to the node
   Node* nod_pt = mesh_pt()->node_pt(n);
   //Get the absolute value of A
   double a13 = pow(std::abs(GlobalVariables::A),(1.0/3.0));

   //Set the value of the inhibitor from the homogeneous solution
   nod_pt->set_value(0,a13*(1.0-a13*a13));
   
   //Set a localised hot spot by making an exponential
   //hump in the activator
   double x = nod_pt->x(0);    double y = nod_pt->x(1);
   double spot = 0.0;

   for(unsigned s=0; s<n_spot; s++)
    {
     //Find the square distance from the centre of the 
     //hot spot
     double r2 = (x - centre_x[s])*(x - centre_x[s]) + 
                 (y - centre_y[s])*(y - centre_y[s]);
     
     //Add an exponential hump to the value of the spot variable
     spot += 2.0*exp(-2.0*r2);
    }
   
   //Set the value of the activator
   nod_pt->set_value(1,-a13 + spot);
  }

 //Document the initial solution
 ofstream filename("input.dat");
 mesh_pt()->output(filename,5);
 filename.close();

 //Set the initial values impulsive
 assign_initial_values_impulsive(Dt);
} 


//====================================================================
/// Timestep the problem
//===================================================================
template<class ELEMENT>
void RefineableActivatorInhibitorProblem<ELEMENT>::timestep(
 const double &dt, const unsigned &nstep)
{
 //Set the problem's Dt for the inital condition bit
 Dt = dt;

 //One uniform refinement so that we stand a chance of getting 
 //something
 this->refine_uniformly();
 //Maximum adaptation for the first timestep
 unsigned max_adapt = 2;

 //Take the first timestep
 bool first = true;
 
 //Set the initial condition
 set_initial_condition();

 //Solve the first step
 unsteady_newton_solve(dt,max_adapt,first);
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
 first = false;

 //Loop over possibilities
 for(unsigned i=1;i<nstep;i++)
  {
   unsteady_newton_solve(dt,max_adapt,first);
   char file1[100];
   sprintf(file1,"step%i.dat",i+1);
   ofstream out1(file1);
   mesh_pt()->output(out1,5);
   out1.close();
  }
}

//===== start_of_main=====================================================
/// Driver code for 2D ActivatorInhibitor problem
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
 double dt = 0.25;
 
 //Construct the problem
 RefineableActivatorInhibitorProblem 
  <RefineableQAdvectionDiffusionReactionElement<2,2,3> >problem;

 //Timestep it
 problem.timestep(dt,2);
}









