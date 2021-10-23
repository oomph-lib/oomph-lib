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
//Driver for a multi-physics problem that couples a Navier--Stokes
//mesh to an advection diffusion mesh, giving Boussinesq convection


//Oomph-lib headers, we require the generic, advection-diffusion,
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"
#include "multi_physics.h"

// Both meshes are the standard rectangular quadmesh
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
/// 2D Convection  problem on two rectangular domains, discretised 
/// with Navier-Stokes and Advection-Diffusion elements. The specific type
/// of elements is specified via the template parameters.
//====================================================================
template<class NST_ELEMENT,class AD_ELEMENT> 
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
   dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e))->
    fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// Doc the solution.
 void doc_solution();

 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);

 /// Access function to the Navier-Stokes mesh
 RectangularQuadMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<NST_ELEMENT>*>(Nst_mesh_pt);
  }

 /// Access function to the Advection-Diffusion mesh
 RectangularQuadMesh<AD_ELEMENT>* adv_diff_mesh_pt() 
  {
   return dynamic_cast<RectangularQuadMesh<AD_ELEMENT>*>(Adv_diff_mesh_pt);
  }
 
private:
 
 /// DocInfo object
 DocInfo Doc_info;

protected:

 /// Mesh of Navier Stokes elements
 RectangularQuadMesh<NST_ELEMENT>* Nst_mesh_pt;

 /// Mesh of advection diffusion elements
 RectangularQuadMesh<AD_ELEMENT>* Adv_diff_mesh_pt;

}; // end of problem class

//===========start_of_constructor=========================================
/// Constructor for convection problem
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::ConvectionProblem()
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

 // Build two standard rectangular quadmesh
 Nst_mesh_pt = 
  new RectangularQuadMesh<NST_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());
 Adv_diff_mesh_pt = 
  new RectangularQuadMesh<AD_ELEMENT>(n_x,n_y,l_x,l_y,time_stepper_pt());

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (both velocity values by default)
   unsigned val_max=2;

   //Loop over the number of nodes on the boundry
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the v-velocity satisfies natural
     //boundary conditions, so we only pin the u-velocity
     if ((ibound==1) || (ibound==3)) 
      {
       val_max=1;
      }

     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 //Pin the zero-th pressure dof in element 0 and set its value to
 //zero:
 fix_pressure(0,0,0.0);

 //Loop over the boundaries of the AD mesh
 num_bound = adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (the temperature value by default)
   unsigned val_max=1;

   //Loop over the number of nodes on the boundry
   unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //If we are on the side-walls, the temperature
     //satisfies natural boundary conditions, so we don't pin anything
     // in this mesh
     if ((ibound==1) || (ibound==3)) 
      {
       val_max=0;
      }
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }


 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructors. 
 unsigned n_nst_element = nst_mesh_pt()->nelement();
 for(unsigned i=0;i<n_nst_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   NST_ELEMENT *el_pt = dynamic_cast<NST_ELEMENT*>
    (nst_mesh_pt()->element_pt(i));

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Inverse_Prandtl;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;

   // We can ignore the external geometric data in the "external"
   // advection diffusion element when computing the Jacobian matrix
   // because the interaction does not involve spatial gradients of 
   // the temperature (and also because the mesh isn't moving!)
   el_pt->ignore_external_geometric_data();

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();

  }

 unsigned n_ad_element = adv_diff_mesh_pt()->nelement();
 for(unsigned i=0;i<n_ad_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   AD_ELEMENT *el_pt = dynamic_cast<AD_ELEMENT*>
    (adv_diff_mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Peclet number multiplied by the Strouhal number
   el_pt->pe_st_pt() =&Global_Physical_Variables::Peclet;

   //The mesh is fixed, so we can disable ALE
   el_pt->disable_ALE();

   // We can ignore the external geometric data in the "external"
   // advection diffusion element when computing the Jacobian matrix
   // because the interaction does not involve spatial gradients of 
   // the temperature (and also because the mesh isn't moving!)
   el_pt->ignore_external_geometric_data();
  }

 // combine the submeshes
 add_sub_mesh(Nst_mesh_pt);
 add_sub_mesh(Adv_diff_mesh_pt);
 build_global_mesh();

 // Set sources
 Multi_domain_functions::
  setup_multi_domain_interactions<NST_ELEMENT,AD_ELEMENT>
  (this,nst_mesh_pt(),adv_diff_mesh_pt());

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl; 

} // end of constructor



//===========start_of_set_boundary_conditions================
/// Set the boundary conditions as a function of continuous 
/// time
//===========================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::set_boundary_conditions(
 const double &time)
{
 // Loop over all the boundaries on the NST mesh
 unsigned num_bound=nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=nst_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=nst_mesh_pt()->boundary_node_pt(ibound,inod);

     //Set the number of velocity components to be pinned
     //(both by default)
     unsigned vel_max=2;

     //If we are on the side walls we only set the x-velocity.
     if((ibound==1) || (ibound==3)) {vel_max = 1;}

     //Set the pinned velocities to zero on NST mesh
     for(unsigned j=0;j<vel_max;j++) {nod_pt->set_value(j,0.0);}

     //If we are on the top boundary
     if(ibound==2) 
      {
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

    }
  }

 // Loop over all the boundaries on the AD mesh
 num_bound=adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary 
   unsigned num_nod=adv_diff_mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=adv_diff_mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the top boundary, set the temperature 
     //to -0.5 (cooled)
     if(ibound==2) {nod_pt->set_value(0,-0.5);}

     //If we are on the bottom boundary, set the temperature
     //to 0.5 (heated)
     if(ibound==0) {nod_pt->set_value(0,0.5);}
    }
  }


} // end_of_set_boundary_conditions

//===============start_doc_solution=======================================
/// Doc the solution
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
void ConvectionProblem<NST_ELEMENT,AD_ELEMENT>::doc_solution()
{ 
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output Navier-Stokes solution
 sprintf(filename,"%s/fluid_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nst_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output advection diffusion solution
 sprintf(filename,"%s/temperature_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 adv_diff_mesh_pt()->output(some_file,npts);
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

#define NEW
//#undef NEW
#ifdef NEW

 //Construct our problem
 ConvectionProblem<
 NavierStokesBoussinesqElement<QCrouzeixRaviartElement<2>,
                               QAdvectionDiffusionElement<2,3> > ,
  AdvectionDiffusionBoussinesqElement<QAdvectionDiffusionElement<2,3>,
                                      QCrouzeixRaviartElement<2> >
  > 
  problem;

#else

//Construct our problem
ConvectionProblem<QCrouzeixRaviartBoussinesqElement<2>,
 QAdvectionDiffusionBoussinesqElement<2> > 
 problem;

#endif

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

 problem.refine_uniformly();

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









