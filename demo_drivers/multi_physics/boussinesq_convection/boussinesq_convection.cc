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
//Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations, 
//giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"
#include "multi_physics.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/rectangular_quadmesh.h"

// Use the oomph and std namespaces 
using namespace oomph;
using namespace std;


//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{
 /// Peclet number (identically one from our non-dimensionalisation)
 double Peclet=1.0;

 /// 1/Prandtl number
 double Inverse_Prandtl=1.0;

 /// Rayleigh number, set to be greater than 
 /// the threshold for linear instability
 double Rayleigh = 1800.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(2);
  
} // end_of_namespace

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//====== start_of_problem_class=======================================
/// 2D Convection  problem on rectangular domain, discretised 
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class ConvectionProblem : public Problem
{

public:

 ///Constructor
 ConvectionProblem();

 /// Destructor. Empty
 ~ConvectionProblem() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt:(empty)
 void actions_before_adapt(){}

 /// Actions before the timestep (update the the time-dependent 
 /// boundary conditions)
 void actions_before_implicit_timestep() 
  {
   set_boundary_conditions(time_pt()->time());
  }

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to specific element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// Doc the solution.
 void doc_solution();

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class ELEMENT>
ConvectionProblem<ELEMENT>::ConvectionProblem()
{
 //Allocate a timestepper
 add_time_stepper_pt(new BDF<2>);

 // Set output directory
 Doc_info.set_directory("RESLT");
 
 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=3.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build a standard rectangular quadmesh
 Problem::mesh_pt() = 
  new RectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=3;
   //If we are on the side-walls, the v-velocity and temperature
   //satisfy natural boundary conditions, so we only pin the
   //first value
   if((ibound==1) || (ibound==3)) {val_max=1;}

   //Loop over the number of nodes on the boundry
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor.
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class ELEMENT>
void ConvectionProblem<ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //Set the pinned velocities to zero
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
      {
       //Set the temperature to -0.5 (cooled)
       nod_pt->set_value(2,-0.5);

       //Add small velocity imperfection if desired
       double epsilon = 0.01;

       //Read out the x position
       double x = nod_pt->x(0);

       //Set a sinusoidal perturbation in the vertical velocity
       //This perturbation is mass conserving
       double value = sin(2.0*MathematicalConstants::Pi*x/3.0)*
        epsilon*time*exp(-time);
       nod_pt->set_value(1,value);
      }

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(2,0.5);}
    }
  }
} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ConvectionProblem<ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
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

 Doc_info.number()++;
} // end of doc


//=======start_of_main================================================
/// Driver code for 2D Boussinesq convection problem
//====================================================================
int main(int argc, char **argv)
{

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;

 //Construct our problem
 ConvectionProblem<BuoyantQCrouzeixRaviartElement<2> > problem;

 // Apply the boundary condition at time zero
 problem.set_boundary_conditions(0.0);

 //Perform a single steady Newton solve
 problem.steady_newton_solve();

 //Document the solution
 problem.doc_solution();

 //Set the timestep
 double dt = 0.1;

 //Initialise the value of the timestep and set an impulsive start
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 200;

 //If we have a command line argument, perform fewer steps 
 //(used for self-test runs)
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution();
  }

} // end of main









