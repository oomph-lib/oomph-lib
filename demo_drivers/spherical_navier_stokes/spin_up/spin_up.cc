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
//Driver for the spin up of a viscous fluid within a sphere, calculated
//using axisymmetric spherical polar coordinates

//Generic includes
#include "generic.h"
#include "spherical_navier_stokes.h"
#include "meshes/simple_rectangular_quadmesh.h"


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
 // Provide a scalar value for the azimuthal velocity
 void boundary_velocity(double& time, const Vector<double>& x, double& u)
 {    	  
  //Spin-up with a "fast" exponential forcing
  u = x[0]*sin(x[1])*(1.0 - std::exp(-100*time));
 }
 
} // End of Boundary_Items namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// SpinUp problem in a sphere
//====================================================================
template<class ELEMENT>
class SphericalSpinUpProblem : public Problem
{

public:

 /// Constructor
 SphericalSpinUpProblem();

 /// Destructor to clean up memory
 ~SphericalSpinUpProblem();

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}
 
 ///  Update the problem specs before solve
 void actions_before_newton_solve(){}
 
 ///  Update the problem specs before next timestep: 
 /// Set Dirichlet boundary conditions from exact solution.
 void actions_before_implicit_timestep();
 
 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, std::ofstream&);

 /// Perform a timestepping study
 void timestep(const double &dt, const unsigned &nstep, 
               const string &output_dir);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for SphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
SphericalSpinUpProblem<ELEMENT>::SphericalSpinUpProblem()
{ 
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timesteps. 
 add_time_stepper_pt(new BDF<2>);
 
 
 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------
 // pi definition
 double pi = MathematicalConstants::Pi;
 
 // # of elements in r-direction
 unsigned n_x=20;
 
 // # of elements in theta-direction
 unsigned n_y=20;
 
 // Domain length in r-direction
 double l_x=1.0;
 
 // Domain length in theta-direction
 double l_y=pi;
 
 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y, time_stepper_pt());
 
 //-------------------------------------------------------------------------
 // Move the nodes to an elliptical shape
 
 // Loop over all the nodes in the mesh
 unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   // Eccentricity variable definition
   double L = Global_Physical_Variables::lambda;
   
   // Redefine the radial coordinate for the elliptical case
   //Get the current values
   double r = mesh_pt()->node_pt(n)->x(0);
   double theta = mesh_pt()->node_pt(n)->x(1);
   //Tell me what you are doing
   //std::cout << "I'm chaging r from " << r << " ";
   r *= sqrt(cos(theta)*cos(theta) + L*L*sin(theta)*sin(theta));
   //std::cout << " to " << r << "\n";

   //Now set the r-position of the node to the new position
   mesh_pt()->node_pt(n)->x(0) = r;
  }
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 const unsigned num_bound = mesh_pt()->nboundary();
 
 // Pin all three velocities on boundaries 1 and 3
 for(unsigned ibound=1;ibound<num_bound;ibound = ibound + 2)
  {
   const unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u/v/w velocities)
     for (unsigned i=0;i<3;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries 1 and 3
  
  // Now pin the theta and phi velocities on boundaries 0 and 2
 for(unsigned ibound=0;ibound<num_bound;ibound = ibound + 2)
  {
   const unsigned num_nod= mesh_pt()->nboundary_node(ibound);
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

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=========start of actions_before_implicit_timestep======================
///  Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void SphericalSpinUpProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Get current time
 double time=time_pt()->time();
 
 
 //Setting for boundary 0 - zero theta and phi velocities
 unsigned ibound=0;
 
 // Loop over the nodes on boundary
 unsigned num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   nod_pt->set_value(1,0);
   nod_pt->set_value(2,0);
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
   nod_pt->set_value(0,0);
   nod_pt->set_value(1,0);
   Boundary_Items::boundary_velocity(time,x,u);
   nod_pt->set_value(2,u);
  }
 
 //Setting for boundary 2 - zero theta and phi velocities
 ibound=2;
 
 // Loop over the nodes on boundary
 num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

   nod_pt->set_value(1,0);
   nod_pt->set_value(2,0);
  }
 
 
 //Setting for boundary 3
 ibound=3;
 
 // Loop over the nodes on boundary
 num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   
   nod_pt->set_value(0,0);
   nod_pt->set_value(1,0);
   nod_pt->set_value(2,0);
  }
  
  
} // end of actions_before_implicit_timestep


//==start_of_destructor===================================================
/// Destructor for SphericalSpinUp problem 
//========================================================================
template<class ELEMENT>
SphericalSpinUpProblem<ELEMENT>::~SphericalSpinUpProblem()
{ 
 //Delete the allocated mesh
 delete Problem::mesh_pt();
} // end_of_destructor


//==start_of_doc_solution=================================================
/// Document the solution
//========================================================================
template<class ELEMENT>
void SphericalSpinUpProblem<ELEMENT>::doc_solution(DocInfo& doc_info, 
                                                   std::ofstream&)
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
 Node* nod_pt = this->mesh_pt()->finite_element_pt(205)->node_pt(0);
 Node* nod2_pt = this->mesh_pt()->finite_element_pt(210)->node_pt(0);
 Node* nod3_pt = this->mesh_pt()->finite_element_pt(215)->node_pt(0);
 
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



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for SphericalSpinUp test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{


 //Set the timestep
 double dt=0.05;
 //Set the number of steps
 unsigned nstep = 3;

 // Doing QTaylorHoodElements
 {
  // Build the problem with QTaylorHoodElements
  SphericalSpinUpProblem<QSphericalTaylorHoodElement > problem;
  cout << "Doing QTaylorHoodElement" << std::endl;

  problem.timestep(dt,nstep,"RESLT_TH");
 } // end of QTaylorHoodElements

 
 // Doing QCrouzeixRaviartElements
 {
  // Build the problem with QCrouzeixRaviartElements
  SphericalSpinUpProblem<QSphericalCrouzeixRaviartElement > problem;
  cout << "Doing QCrouzeixRaviartElement" << std::endl;

  problem.timestep(dt,nstep,"RESLT_CR");
 } // end of QCrouzeixRaviartElements


} // end_of_main
