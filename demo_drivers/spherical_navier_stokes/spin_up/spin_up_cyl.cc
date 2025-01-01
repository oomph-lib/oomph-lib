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
//Driver for the spin up of a viscous fluid within a sphere, 
//calculated using axisymmetric cylindrical polar coordinates
//and imposing up-down symmetry
//Note that although this is for a fixed mesh, the mesh
//is obtained by uniformly refining a refineable mesh,
//hence we must use refineable elements.

//Generic includes
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "meshes/quarter_circle_sector_mesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 // Reynolds number
 double Re=10;
 
 // Womersley number
 double ReSt=10;
 
 // Eccentricity parameter
 double lambda=1.0;

} // end of GPV namespace



//===================================================
// Namespace for solution vectors
//===================================================
namespace Boundary_Items
{    	
 // Provide a scalar value for the velocity at boundary 1
 void boundary_velocity(const double& time, const Vector<double>& x, double& u)
 {    	  
  // u assignment - spin-up problem
  u = x[0]*(1.0 - std::exp(-100*time));
 }
	
} // End of Boundary_Items namespace




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Spin-up problem in  a sphere
//====================================================================
template<class ELEMENT>
class SphericalSpinUpProblem : public Problem
{

public:


 /// Constructor
 SphericalSpinUpProblem();

 /// Destructor to clean up memory
 ~SphericalSpinUpProblem();


 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure



 /// Set the boundary conditions
 void set_boundary_conditions();

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve
 void actions_before_newton_solve(){}
  
 /// Update the problem specs before next timestep: 
 /// Set Dirichlet boundary conditions from exact solution.
 void actions_before_implicit_timestep() {set_boundary_conditions();}

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   //Reset the boundary conditions
   set_boundary_conditions();
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);
  }
 
 
 // Access function for the specific mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, std::ofstream&);

 /// Perform a timestepping study
 void timestep(const double &dt, const unsigned &nstep, 
               const string &output_dir);

private:

 /// Geometric object that defines the boundary of the domain
 Ellipse* Curved_boundary_pt;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for SphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
SphericalSpinUpProblem<ELEMENT>::SphericalSpinUpProblem()
{ 
 /// Build the geometric object that describes the outer wall
 Curved_boundary_pt = new Ellipse(1.0,1.0);
 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timesteps. 
 add_time_stepper_pt(new BDF<2>);
 
 // The curved boundary of the mesh is defined by the geometric object
 // What follows are the start and end coordinates on the geometric object:
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);
 
 // Fraction along geometric object at which the radial dividing line
 // is placed
 double fract_mid=0.5;
 
 //Now create the mesh using the geometric object
 Problem::mesh_pt() = new RefineableQuarterCircleSectorMesh<ELEMENT>(
  Curved_boundary_pt,xi_lo,fract_mid,xi_hi,time_stepper_pt());
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 
 //The base is a symmetry boundary
 {
  const unsigned ibound = 0;
  const unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Pin the z-velocity only
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
   }
 }


 //The axis (boundary 2) is a symmetry boundary
 {
  const unsigned ibound = 2;
  const unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Pin the r-velocity and phi velocity
    mesh_pt()->boundary_node_pt(2,inod)->pin(0); 
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(2); 
   }
 }
  
 //The outside has pure Dirichlet conditions
 {
  const unsigned ibound = 1;
  const unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Pin the all velocity components
    for(unsigned i=0;i<3;i++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
     }
   }
 }
  
  
 // Complete the build of all elements so they are fully functional
 //================================================================

 //Find number of elements in mesh
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   //Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   
   //Disable ALE
   el_pt->disable_ALE();

  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=========start of actions_before_implicit_timestep======================
/// Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void SphericalSpinUpProblem<ELEMENT>::set_boundary_conditions()
{
 // Get current time
 const double time=time_pt()->time();
 
 //Setting for boundary 0 - zero z velocity
 unsigned ibound=0;
 
 // Loop over the nodes on boundary
 unsigned num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   //Set the z value to zero
   nod_pt->set_value(1,0.0);
  }
  
  
 //Set velocity for boundary 1 - driving wall, varies with theta
 ibound=1;
 
 // Loop over the nodes on boundary
 num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   double u;
   
   Vector<double> x(2);
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   // Get current values of the boundary conditions from the
   // exact solution
   nod_pt->set_value(0,0.0);
   nod_pt->set_value(1,0.0);
   Boundary_Items::boundary_velocity(time,x,u);
   nod_pt->set_value(2,u);
  }
    
    
 //Setting for boundary 2 - zero r and phi velocities
 ibound=2;
 
 // Loop over the nodes on boundary
 num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   
   nod_pt->set_value(0,0.0);
   nod_pt->set_value(2,0.0);
  }
    
    
} // end of actions_before_implicit_timestep


//==start_of_destructor===================================================
/// Destructor for SphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
SphericalSpinUpProblem<ELEMENT>::~SphericalSpinUpProblem()
{ 
 
 //Delete the mesh
 delete Problem::mesh_pt();
 //Delete the geometric object
 delete Curved_boundary_pt;

} // end_of_destructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void SphericalSpinUpProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream&)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 const unsigned npts=5;
 
 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for t=" << time_pt()->time() << std::endl;
 cout << "=================================================" << std::endl;


 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);

 // Print out the time details for data
 some_file << "time = " << time_pt()->time() << "\"";
  
 some_file.close();
 } // end_of_doc_solution



/// Timestep the problem with a given (fixed) timestep dt for nstep steps
template<class ELEMENT>
void SphericalSpinUpProblem<ELEMENT>::timestep(const double &dt, 
                                               const unsigned &nstep,
                                               const string &output_dir)
{
 // Set up doc info
 // ---------------
 
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory(output_dir);
 
 // Step number
 doc_info.number()=0;

 // ---------------
 // end of Set up doc info

 // Open a trace file
 std::ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/time_trace.dat",doc_info.directory().c_str());
 trace_file.open(filename);
 trace_file << "time " << "u " << "v " << "w "
            << std::endl;


 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 this->initialise_dt(dt);
 
 // Set IC
 this->assign_initial_values_impulsive();
 
 //Output initial condition
 this->doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
 doc_info.number()++;

 //Store nodes for output
 Node* nod_pt =  this->mesh_pt()->boundary_node_pt(0,1);
 Node* nod2_pt = this->mesh_pt()->boundary_node_pt(0,2);
 Node* nod3_pt = this->mesh_pt()->boundary_node_pt(0,3);
 
 // Timestepping loop
 for(unsigned istep=0;istep<nstep;istep++)
  {
   std::cout << " Timestep " << istep << std::endl;
   
   // Take timestep
   this->unsteady_newton_solve(dt);
   
   //Output solution
   this->doc_solution(doc_info,trace_file);
   
   //Output azimuthal velocities on the equator
   trace_file << this->time() << " " 
              << nod_pt->x(0) << " " << nod_pt->x(1)
              << " " << nod_pt->value(2) << " " 
              << nod2_pt->x(0) << " " << nod2_pt->x(1)
              << " " << nod2_pt->value(2) << " " 
              << nod3_pt->x(0) << " " << nod3_pt->x(1)
              << " " << nod3_pt->value(2) << std::endl; 
   
   
   //Increment counter for solutions 
   doc_info.number()++;
  }
 
 // Close trace file
 trace_file.close();
}



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for SphericalSpinUp test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{
  
 // Choose simulation interval and timestep
 double dt=0.05; 
 unsigned nstep=3;

 // Doing QTaylorHoodElements
 {
  
  // Build the problem with QTaylorHoodElements
  SphericalSpinUpProblem<RefineableAxisymmetricQTaylorHoodElement > 
   problem;
  cout << "Doing QTaylorHoodElement" << std::endl;
  
  //Refine uniformly a couple of times
  for(unsigned i=0;i<3;i++) {problem.refine_uniformly();}

  problem.timestep(dt,nstep,"RESLT_TH");

 } // end of QTaylorHoodElements

 
 // Doing QCrouzeixRaviartElements
 {
  
  // Build the problem with QCrouzeixRaviartElements
  SphericalSpinUpProblem<RefineableAxisymmetricQCrouzeixRaviartElement > 
   problem;
  cout << "Doing QCrouzeixRaviartElement" << std::endl;
  
  //Refine uniformly a couple of times
  for(unsigned i=0;i<3;i++) {problem.refine_uniformly();}

  problem.timestep(dt,nstep,"RESLT_CR");

 } // end of QCrouzeixRaviartElements


} // end_of_main
