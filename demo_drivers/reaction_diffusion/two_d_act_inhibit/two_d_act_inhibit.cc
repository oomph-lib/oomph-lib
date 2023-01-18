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
//A driver code that solves the activator-inhibitor reaction-diffusion
//problem specified in Muratov & Osipov (Phys Rev E 54, 5, pp 4860--4879)
//Periodic boundary conditions are applied in this case
#include <iostream>
#include <fstream>
#include <cmath>

#include "generic.h"
#include "navier_stokes.h"
#include "advection_diffusion_reaction.h"
#include "meshes/rectangular_quadmesh.template.cc"

using namespace std;

using namespace oomph;

namespace GlobalVariables
{
 //The constant pi
 const double pi = MathematicalConstants::Pi;

 //The two timescale parameters (both set to one)
 Vector<double> Tau(2,1.0);
 
 //The two diffusion parameters (both set to one)
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

//=======================================================================
/// The ActivatorInhibitor problem class solved on a 2x2 square
/// with periodic boundary conditions
//======================================================================
template<class ELEMENT>
class ActivatorInhibitorProblem : public Problem
{
public:

 /// Constructor, number of elements in the x and y directions
 ActivatorInhibitorProblem(const unsigned &nx, const unsigned &ny);

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve: 
 void actions_before_newton_solve() {}

 /// Timestep the problem with constant timestep dt for nstep steps
 void timestep(const double &dt, const unsigned &nstep);

};


//========================================================================
/// Constructor pass the number of elements in the x and y directions
//========================================================================
 template<class ELEMENT>
 ActivatorInhibitorProblem<ELEMENT>::ActivatorInhibitorProblem(
  const unsigned &nx, const unsigned &ny)
{ 
 //Allocate the timestepper 
 add_time_stepper_pt(new BDF<2>); //SEG FAULT for Steady<0>

 // Build simple rectangular mesh on 2x2 domain
 Problem::mesh_pt()= 
  new RectangularQuadMesh<ELEMENT>(nx,ny,2.0,2.0,
                                   Problem::time_stepper_pt());
 
 //Sort out the periodic boundaries (miss out the ends)
 //Storage for the corners
 Vector<Node*> corners_pt(4);
  //Find the number of nodes on the vertical boundary
 unsigned n_node = mesh_pt()->nboundary_node(1);
 if(n_node > 0)
  {
   corners_pt[0] = mesh_pt()->boundary_node_pt(1,0);
   corners_pt[1] = mesh_pt()->boundary_node_pt(1,n_node-1);
   corners_pt[2] = mesh_pt()->boundary_node_pt(3,0);
   corners_pt[3] = mesh_pt()->boundary_node_pt(3,n_node-1);
  }
 
 //Make the side boundaries periodic (miss out the ends)
 for(unsigned n=1;n<n_node-1;n++)
  {
   mesh_pt()->boundary_node_pt(1,n)
    ->make_periodic(mesh_pt()->boundary_node_pt(3,n));
  }
 
 //Now make the lower boundaries periodic (miss out the ends)
 n_node = mesh_pt()->nboundary_node(0);
 for(unsigned n=1;n<n_node-1;n++)
  {
   mesh_pt()->boundary_node_pt(0,n)
    ->make_periodic(mesh_pt()->boundary_node_pt(2,n));
  }
 
 //Now sort out the corners
 corners_pt[0]->make_periodic_nodes(corners_pt); 

//Add the element specific information
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   
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

 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}

//=======================================================================
/// Function to timestep the problem
//======================================================================
template<class ELEMENT>
void ActivatorInhibitorProblem<ELEMENT>::timestep(const double &dt,
                                        const unsigned &nstep)
 
{ 
 //Set the initial concentrations of the reagents --- the homogeneous state
 //plus a perturbation of amplitude amp

 //Find the Maximum unsigned value from the RNG
 unsigned long Max = RAND_MAX;
 //Seed the random number generator with a constant value
 srand(10273);
 
 //Add some small noise of amplitude amp
 const double amp = 1.0e-4;

 //Store the global value of pi
 const double pi = GlobalVariables::pi;

 //Loop over all nodes in the mesh
 unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   //Local pointer to the node
   Node* nod_pt = mesh_pt()->node_pt(n);
   //Get the absolute value of the parameter A
   double a13 = pow(std::abs(GlobalVariables::A),(1.0/3.0));
   
   //Get two random values from a uniform distribution
   double u1 = rand()/(double)(Max);
   double u2 = rand()/(double)(Max);
   
   //Reapeat until u1 is not too small
   while(u1 < 1.0e-10) {u1 = rand()/(double)(Max);}
   //Now do box muller to get to normal variables
   double Amp = amp*sqrt(-2.0*log(u1));
   double n1 = Amp*cos(2.0*pi*u2);
   double n2 = Amp*sin(2.0*pi*u2);
   
   //Set the values with Gaussian noise
   nod_pt->set_value(0, a13*(1.0-a13*a13) + n1);
   nod_pt->set_value(1, -a13 + n2);
  }
 
 //Document the solution
 ofstream filename("input.dat");
 mesh_pt()->output(filename,5);
 filename.close();

 //Set the initial values
 assign_initial_values_impulsive(dt);

 //Loop over possibilities
 for(unsigned i=0;i<nstep;i++)
  {
   unsteady_newton_solve(dt);
   char file1[100];
   // sprintf(file1,"step%g.dat",time_pt()->time());
   sprintf(file1,"step%i.dat",i+1);
   ofstream out1(file1);
   mesh_pt()->output(out1,5);
   out1.close();
  }
}


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


 int main()
  {
   double epsilon = 0.05;
   //Set the control parameters (should produce a spot)
   GlobalVariables::Tau[1] = 0.5;
   GlobalVariables::D[1] = epsilon*epsilon;
   GlobalVariables::A = -0.1;

   //Set the timestep
   double dt = 0.1;

   ActivatorInhibitorProblem 
    <QAdvectionDiffusionReactionElement<2,2,3> >problem(20,20);
   problem.timestep(dt,5);
  }
