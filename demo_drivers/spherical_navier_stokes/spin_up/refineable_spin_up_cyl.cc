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
//Driver for the spin-up inside a sphere of radius 1, but in cylindrical
//polar coordinates

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
 void get_exact_u_b1(double& time, const Vector<double>& x, double& u)
 {    	  
  // u assignment - spin-up problem
  u = x[0]*(1.0 - std::exp(-time));
  
  // u assignment - DSDJ problem
  //u = x[0]*sin(time)*sin(x[1]);
 }
	
} // End of Boundary_Items namespace




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain - time dependent version
//====================================================================
template<class ELEMENT>
class RefineableSphericalSpinUpProblem : public Problem
{

public:


 /// Constructor
 RefineableSphericalSpinUpProblem();

 /// Destructor to clean up memory
 ~RefineableSphericalSpinUpProblem();


 ///Fix pressure in element e at pressure dof pdof and set to pvalue
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
 
 
  /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 
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
 
 /// Set initial conditions: Set all nodal velocities to zero and
 /// initialise the previous velocities to correspond to an impulsive
 /// start
 void set_initial_condition()
  {
   // Determine number of nodes in mesh
   const unsigned n_node = mesh_pt()->nnode();

   // Loop over all nodes in mesh
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       // Set velocity component i of node n to zero
       mesh_pt()->node_pt(n)->set_value(i,0.0);
      }
    }

   // Initialise the previous velocity values for timestepping
   // corresponding to an impulsive start
   assign_initial_values_impulsive();

  } // End of set_initial_condition
 
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

private:

 ///Geometric object that defines the boundary of the domain
 Ellipse* Curved_boundary_pt;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for RefineableSphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
RefineableSphericalSpinUpProblem<ELEMENT>::RefineableSphericalSpinUpProblem()
{ 
 ///Build the geometric object that describes the outer wall
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
 

 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------
 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
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
 unsigned n_element = mesh_pt()->nelement();

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
void RefineableSphericalSpinUpProblem<ELEMENT>::set_boundary_conditions()
{
 // Get current time
 double time=time_pt()->time();
 
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
   Boundary_Items::get_exact_u_b1(time,x,u);
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
/// Destructor for RefineableSphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
RefineableSphericalSpinUpProblem<ELEMENT>::~RefineableSphericalSpinUpProblem()
{ 

 // Mesh gets killed in general problem destructor

} // end_of_destructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableSphericalSpinUpProblem<ELEMENT>::doc_solution(DocInfo& doc_info, std::ofstream&)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 
 
 
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



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for RefineableSphericalSpinUp test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{
  
 // Choose simulation interval and timestep
 unsigned nstep =1;
 double dt=0.25;
	
 // Set up doc info
 // ---------------

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;

 // ---------------
 // end of Set up doc info
 

 // Doing QCrouzeixRaviartElements
 {
  
  // Build the problem with QCrouzeixRaviartElements
  RefineableSphericalSpinUpProblem<RefineableAxisymmetricQCrouzeixRaviartElement > 
   problem;
  cout << "Doing QCrouzeixRaviartElement" << std::endl;
    
    // Open a trace file
    ofstream trace_file;
    char filename[100];   
    sprintf(filename,"%s/trace0.dat",doc_info.directory().c_str());
    trace_file.open(filename);
    trace_file << "time " << "u " << "v " << "w "
               << std::endl;

    //Refine uniformly once
    problem.refine_uniformly();


    // Initialise timestep -- also sets the weights for all timesteppers
    // in the problem.
    problem.initialise_dt(dt);
 
    // Set initial conditions
    problem.set_initial_condition();
    
    // Over-ride the maximum and minimum permitted errors
    problem.mesh_pt()->max_permitted_error() = 1.0e-2; //Default = 1.0e-3
    problem.mesh_pt()->min_permitted_error() = 1.0e-3; //Default = 1.0e-5
    
    // Over-ride the maximum and minimum permitted refinement levels
    problem.mesh_pt()->max_refinement_level() = 4;//maximum_ref_level;
    problem.mesh_pt()->min_refinement_level() = 1;//minimum_ref_level;
 
    //Specify the normalising factor explicitly
    Z2ErrorEstimator* error_pt = dynamic_cast<Z2ErrorEstimator*>(
     problem.mesh_pt()->spatial_error_estimator_pt());
    error_pt->reference_flux_norm() = 0.01;

    //Output initial condition
    problem.doc_solution(doc_info,trace_file);
 
    //Increment counter for solutions 
    doc_info.number()++;

    unsigned max_adapt = 3;

    bool first = true;

    // Timestepping loop
    for (unsigned istep=0;istep<nstep;istep++)
     {
      cout << " Timestep " << istep << std::endl;
   
      // Take timestep
      problem.unsteady_newton_solve(dt,max_adapt,first);
   
      //Set first to false
      first=false;
      max_adapt=1;

      //Output solution
      problem.doc_solution(doc_info,trace_file);

      Node* nod_pt = //problem.mesh_pt()->finite_element_pt(435)->node_pt(0);
       problem.mesh_pt()->node_pt(50);

      trace_file << problem.time() << " " 
                 << nod_pt->x(0) << " " << nod_pt->x(1)
                 << " " << nod_pt->value(2) << std::endl;

      //Increment counter for solutions 
      doc_info.number()++;
  }
 
 // Close trace file
 trace_file.close();

 } // end of QCrouzeixRaviartElements


} // end_of_main
