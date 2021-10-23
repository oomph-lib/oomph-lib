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
//Refineable version of a problem of convection in a spherical shell with
//an electromagnetic (rather than gravitational) central body force. 
//The problem is described in Feudel et al. (2011) Phys Rev E 83, 046304
//There is a loss of spherical symmetry at Ra=2491, which leads to an
//axisymmetric branch of solutions at a transcritical bifurcation.


//Generic includes
#include "generic.h"
#include "spherical_advection_diffusion.h"
#include "spherical_buoyant_navier_stokes.h"

//Standard rectangular mesh
#include "meshes/rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
 //Reynolds number 
 double Re  = 1.0;

 /// Peclet number
 double Pe = 1.0;

 /// Rayleigh number
 double Ra = 10.0;
 
 /// Gravity
 Vector<double> G;

} // end of GPV namespace

//==start_of_problem_class===============================================
/// Convection in a spherical shell
//=======================================================================
template<class ELEMENT>
class RefineableSphereConvectionProblem : public Problem
{

public:


 /// Constructor
 RefineableSphereConvectionProblem();

 /// Destructor to clean up memory
 ~RefineableSphereConvectionProblem();

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Set the boundary conditions
 void set_boundary_conditions(const double &time);
 
 void actions_before_implicit_timestep()
  {
   this->set_boundary_conditions(this->time());
  }

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableSphericalNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableSphericalNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   //Reset the boundary conditions
   set_boundary_conditions(this->time());
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);
  }


 // Access function for the specific mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Document the solution
 void doc_solution(DocInfo& doc_info);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RefineableSphereConvection problem 
//========================================================================
template<class ELEMENT>
RefineableSphereConvectionProblem<ELEMENT>::RefineableSphereConvectionProblem()
{
 
 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------
 this->add_time_stepper_pt(new BDF<2>);

 // pi definition
 const double pi = MathematicalConstants::Pi;
     
 // # of elements in r-direction
 unsigned n_r=4;

 // # of elements in theta-direction
 unsigned n_theta=4;

 // Radius of inner sphere
 double R_inner = 0.5;

 // Radius of outer sphere
 double R_outer = 1.0;
 
  // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_r,n_theta,R_inner,R_outer,
                                             0.0,pi,this->time_stepper_pt());

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 
 // Pin the temperature and velocities on boundaries 1 and 3
 for(unsigned ibound=1;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<4;++i)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries 1 and 3
 

 // Now pin the theta and phi velocities on boundaries 0 and 2
 for(unsigned ibound=0;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over the theta- and phi-velocities
     for (unsigned i=1; i<3; i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
      }
    }
  } // end loop over boundaries 0 and 2
  
  // end of set boundary conditions
  
   
  
 // Complete the build of all elements so they are fully functional
 //================================================================

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));


    // Set the Peclet number
    el_pt->pe_pt() = &Global_Physical_Variables::Pe;
    
    // Set the Peclet number multiplied by the Strouhal number
    el_pt->pe_st_pt() =&Global_Physical_Variables::Pe;
    
    // Set the Reynolds number (1/Pr in our non-dimensionalisation)
    el_pt->re_pt() = &Global_Physical_Variables::Re;
    
   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Re;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Ra;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::G;

   //Disable ALE
   el_pt->disable_ALE();

  } // end loop over elements


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
} // end_of_constructor


//=========start of set_boundary_conditions===============================
/// Set the boundary conditions so that the inner sphere has
/// a constant angular rotation of angular velocity one.
//========================================================================
template<class ELEMENT>
void RefineableSphereConvectionProblem<ELEMENT>::set_boundary_conditions(
 const double &time)
{
 //Set velocity and temperature for boundary 1 - outer wall (zero)
  unsigned ibound=1;
 
   // Loop over the nodes on boundary
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     for(unsigned i=0;i<4;++i)
      {
       nod_pt->set_value(i,0.0);
      }

     //Set a wobble on the azimuthal velocity
     double theta = nod_pt->x(1);
     double value = sin(theta)*0.001*time*exp(-time);
     nod_pt->set_value(2,value);
    }
  
   //Setting for boundary 3 (inner sphere that is heated)
   ibound=3;
  
   // Loop over the nodes on boundary
   num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     for(unsigned i=0;i<3;++i)
      {
       nod_pt->set_value(i,0.0);
      }
     nod_pt->set_value(3,1.0);
    }
  
  
} // end of actions_before_implicit_timestep


//==start_of_destructor===================================================
/// Destructor for RefineableSphereConvection problem 
//========================================================================
template<class ELEMENT>
RefineableSphereConvectionProblem<ELEMENT>::~RefineableSphereConvectionProblem()
{ 
 //Delete the error estimator
 delete mesh_pt()->spatial_error_estimator_pt();

 //Clean up the memory allocated for the mesh
 delete Problem::mesh_pt();

} // end_of_destructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableSphereConvectionProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 
 
 
 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for Ra =" << Global_Physical_Variables::Ra 
      << std::endl;
 cout << "=================================================" << std::endl;


 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%g_%g.dat",doc_info.directory().c_str(),
         Global_Physical_Variables::Ra,this->time());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 
 some_file.close();
} // end_of_doc_solution


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for RefineableSphereConvection test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main(int argc, char **argv)
{
 //Set gravity
 Global_Physical_Variables::G.resize(3,0.0);
 Global_Physical_Variables::G[0] = -1.0;

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT_adapt");

 // ---------------
 // end of Set up doc info
 
 // Build the problem with QCrouzeixRaviartElements
 RefineableSphereConvectionProblem<
 RefineableBuoyantQSphericalCrouzeixRaviartElement> 
  problem;
 
 // Over-ride the maximum and minimum permitted errors
 problem.mesh_pt()->max_permitted_error() = 1.0e-3; //Default = 1.0e-3
 problem.mesh_pt()->min_permitted_error() = 1.0e-5; //Default = 1.0e-5
 
 // Over-ride the maximum and minimum permitted refinement levels
 problem.mesh_pt()->max_refinement_level() = 4;//maximum_ref_level;
 problem.mesh_pt()->min_refinement_level() = 1;//minimum_ref_level;
 
  //Set the boundary conditions
 problem.set_boundary_conditions(0.0);
 
 //Set the maximum adaptation
 unsigned max_adapt = 2;
  
 // Solve the problem
 problem.steady_newton_solve(max_adapt);
 problem.doc_solution(doc_info);

 //Crank up the Rayleigh number and kick it
 Global_Physical_Variables::Ra = 2500.0;
 
 //Set the timestep

 //Set the timestep
 double dt = 1.0;

 //Initialise the value of the timestep and set an impulsive start
 problem.assign_initial_values_impulsive(dt);

 //Set the number of timesteps to our default value
 unsigned n_steps = 100; //This is enough to get to steady state

 //If we have a command line argument, perform fewer steps
 if(argc > 1) {n_steps = 5;}

 //Perform n_steps timesteps
 for(unsigned i=0;i<n_steps;++i)
  {
   problem.unsteady_newton_solve(dt);
   problem.doc_solution(doc_info);
  }



} // end_of_main
