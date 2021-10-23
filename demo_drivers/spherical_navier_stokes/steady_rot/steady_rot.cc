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
//Driver for the steady rigid body rotation
//of a viscous fluid within a sphere, calculated
//using axisymmetric spherical polar coordinates
//The main aim of this driver is to check convergence rates of
//the spherical navier stokes elements against an exact solution.

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
 double Re=0.0;

} // end of GPV namespace


//=======================================================
///Namespace for the exact rigid-body-rotation solution
//=======================================================
namespace ExactSolution
{
 ///The rigid-body rotation solution
 void rigid_body_rotation(const Vector<double> &x,
                          Vector<double> &u)
 {
  //store r sin(theta)
  double r_sin_theta = x[0]*sin(x[1]);
  u[0] = 0.0;
  u[1] = 0.0;
  //Only azimuthal velocity
  u[2] = r_sin_theta;
  //The pressure is the only part the involves the Reynolds number
  u[3] = 0.5*Global_Physical_Variables::Re*r_sin_theta*r_sin_theta;
 }

 ///The r-derivative of the rigid-body rotation solution
 void rigid_body_rotation_dr(const Vector<double> &x,
                              Vector<double> &u)
 {
  //store  sin(theta)
  double sin_theta = sin(x[1]);
  u[0] = 0.0;
  u[1] = 0.0;
  //Only azimuthal velocity
  u[2] = sin_theta;
  //The pressure is the only part the involves the Reynolds number
  u[3] = Global_Physical_Variables::Re*x[0]*sin_theta*sin_theta;
 }


 ///The theta-derivative of the rigid-body rotation solution
 void rigid_body_rotation_dtheta(const Vector<double> &x,
                              Vector<double> &u)
 {
  //store  sin(theta)
  u[0] = 0.0;
  u[1] = 0.0;
  //Only azimuthal velocity
  u[2] = x[0]*cos(x[1]);
  //The pressure is the only part the involves the Reynolds number
  u[3] = Global_Physical_Variables::Re*x[0]*x[0]*sin(x[1])*cos(x[1]);
 }
     
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// SteadyRotation problem in a sphere
//====================================================================
template<class ELEMENT>
class SphericalSteadyRotationProblem : public Problem
{

public:

 /// Constructor
 SphericalSteadyRotationProblem(const unsigned &nr, const unsigned &ntheta);

 /// Destructor to clean up memory
 ~SphericalSteadyRotationProblem();

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
 
 /// Update the problem specs before solve
 void actions_before_newton_solve(){}
 
 // Access function for the specific mesh
 SimpleRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<SimpleRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Perform a timestepping study
 void parameter_study(std::ofstream &trace, const string &output_dir);

private:

 //Storage for number of elements in the radial direction
 unsigned Nr;
 
 ///Storage for number of elements in the theta direction
 unsigned Ntheta;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for SphericalSteadyRotation problem 
//========================================================================
template<class ELEMENT>
SphericalSteadyRotationProblem<ELEMENT>::SphericalSteadyRotationProblem(
 const unsigned &nr, const unsigned &ntheta) : Nr(nr), Ntheta(ntheta)
{ 
 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------
 // pi definition
 double pi = MathematicalConstants::Pi;
 
 // # of elements in r-direction
 unsigned n_x=nr;
 
 // # of elements in theta-direction
 unsigned n_y=ntheta;
 
 // Domain length in r-direction
 double l_x=1.0;
 
 // Domain length in theta-direction
 double l_y=pi;
 
 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);
 
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
     Node* nod_pt = mesh_pt()->boundary_node_pt(ibound,inod);
     // Loop over values (u/v/w velocities)
     for (unsigned i=0;i<3;i++) {nod_pt->pin(i);}

     //Set the velocity on the outer boundary
     if(ibound==1)
      {
       double r = nod_pt->x(0); double theta = nod_pt->x(1);
       nod_pt->set_value(2,r*sin(theta));
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
   
   //Disable ALE
   el_pt->disable_ALE();

  } // end loop over elements

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//==start_of_destructor===================================================
/// Destructor for SphericalSteadyRotation problem 
//========================================================================
template<class ELEMENT>
SphericalSteadyRotationProblem<ELEMENT>::~SphericalSteadyRotationProblem()
{ 
 //Delete the allocated mesh
 delete Problem::mesh_pt();
} // end_of_destructor


/// Timestep the problem with a given (fixed) timestep dt for nstep steps
template<class ELEMENT>
void SphericalSteadyRotationProblem<ELEMENT>::parameter_study(
 std::ofstream &trace, const string &output_dir)
{
 Global_Physical_Variables::Re = 0.0;

 //Store the number of elements
 const unsigned n_element = this->mesh_pt()->nelement();

 ofstream junk;

 for(unsigned i=0;i<2;i++)
  {
   //Solve the problem
   this->newton_solve();
   
   
   //Output solution
   ofstream some_file;
   char filename[100];
   
   // Number of plot points
   const unsigned npts=5;
   
   // Output solution 
   //-----------------
   sprintf(filename,"soln_%s_%ix%i_%g.dat",
           output_dir.c_str(),Nr,Ntheta,Global_Physical_Variables::Re);
   some_file.open(filename);
   mesh_pt()->output(some_file,npts);
   some_file.close();
   
 
   //Storage for the global velocity and pressure errors
   double u_error = 0.0, u_norm = 0.0;
   double p_error = 0.0, p_norm = 0.0;

   //Loop over the elements and compute the appropriate errors 
   for(unsigned e=0;e<n_element;e++)
    {
     //Local storage for the velocity and pressure errors and norms
     double el_u_error=0.0, el_u_norm=0.0;
     double el_p_error=0.0, el_p_norm=0.0;
     dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(e))->
      compute_error_e(junk,ExactSolution::rigid_body_rotation,
                      ExactSolution::rigid_body_rotation_dr,
                      ExactSolution::rigid_body_rotation_dtheta,
                      el_u_error,el_u_norm,el_p_error,el_p_norm);

     //Add the elemental contribution to the global errors
     u_error += el_u_error;
     u_norm += el_u_norm;
     p_error += el_p_error;
     p_norm += el_p_norm;
    }
   
   
   //Output the velocity and pressure errors
   trace << Global_Physical_Variables::Re << " " 
         << Nr << " " << Ntheta << " "
         << std::sqrt(u_error) << " " << std::sqrt(u_norm) << " " 
         << std::sqrt(p_error) << " " << std::sqrt(p_norm) << std::endl;
   
   
   //Increment the Reynolds number
   Global_Physical_Variables::Re += 10;
 
  }
 
}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for SphericalSteadyRotation test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{

 //Number of elements in each direction
 unsigned N = 2;

 ofstream trace("trace_TH.dat");

 // Doing QTaylorHoodElements
 for(unsigned i=0;i<2;i++)
 {
  N *= 2;
  // Build the problem with QTaylorHoodElements
  SphericalSteadyRotationProblem<QSphericalTaylorHoodElement > problem(N,N);
  cout << "Doing QTaylorHoodElement" << std::endl;

  problem.parameter_study(trace,"TH");
 } // end of QTaylorHoodElements

 trace.close();
 trace.open("trace_CR.dat");

 N = 2;
 
 // Doing QCrouzeixRaviartElements
 for(unsigned i=0;i<2;i++)
 {
  N *= 2;
  // Build the problem with QCrouzeixRaviartElements
  SphericalSteadyRotationProblem<QSphericalCrouzeixRaviartElement > problem(N,N);
  cout << "Doing QCrouzeixRaviartElement" << std::endl;

  problem.parameter_study(trace,"CR");
 } // end of QCrouzeixRaviartElements

 trace.close();

} // end_of_main
