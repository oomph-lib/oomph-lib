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
//Driver for concentration field from a sphere falling in still fluid 
//secreting mass at a constant concentration (C=1).
//Not properly validated against anything analytic, but passes the eyeball
//test

//Generic includes
#include "generic.h"
#include "spherical_advection_diffusion.h"

//Standard rectangular mesh
#include "meshes/rectangular_quadmesh.h"


using namespace std;

using namespace oomph;
 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 // Peclet number
 double Pe=10;
} // end of GPV namespace


namespace Wind_Function
{
 /// Wind that represents a constantly translating sphere
 /// in spherical polar coordinates
 void get_wind(const Vector<double>& x, Vector<double>& wind)
 {
  wind[0]=(-1.0 + 1.5/x[0] - 0.5/(x[0]*x[0]*x[0]))*cos(x[1]);
  wind[1]=(1.0 - 0.75/x[0] - 0.25/(x[0]*x[0]*x[0]))*sin(x[1]);
 }
}



//==start_of_problem_class============================================
//A sphere that falls and elutes a chemical
//====================================================================
template<class ELEMENT>
class EultingSphereProblem : public Problem
{

public:


 /// Constructor
 EultingSphereProblem();

 /// Destructor to clean up memory
 ~EultingSphereProblem();

 /// Set the boundary conditions
 void set_boundary_conditions();
 
 // Access function for the specific mesh
 RectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<RectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Document the solution
 void doc_solution(DocInfo& doc_info);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for EultingSphere problem 
//========================================================================
template<class ELEMENT>
EultingSphereProblem<ELEMENT>::EultingSphereProblem()
{ 
 // Setup mesh  -don't forget to include the timestepping in the mesh build
 //------------------------------------------------------------------------

 // pi definition
 const double pi = MathematicalConstants::Pi;
     
 // # of elements in r-direction
 unsigned n_r=10;

 // # of elements in theta-direction
 unsigned n_theta=10;

 // Radius of inner sphere
 double R_inner = 1.0;

 // Radius of outer sphere
 double R_outer=30.0;
 
  // Build and assign mesh
 Problem::mesh_pt() = 
  new RectangularQuadMesh<ELEMENT>(n_r,n_theta,R_inner,R_outer,
                                   0.0,pi);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 
 // Pin the concentration on boundaries 1 and 3
 for(unsigned ibound=1;ibound<num_bound;ibound = ibound + 2)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  } // end loop over boundaries 1 and 3
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

   //Set the Reynolds number
   el_pt->pe_pt() = &Global_Physical_Variables::Pe;
   
   //Set the wind function
   el_pt->wind_fct_pt() = Wind_Function::get_wind;
  } // end loop over elements


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
} // end_of_constructor


//=========start of set_boundary_conditions===============================
/// Set the boundary conditions so that the inner sphere has
/// a constant angular rotation of angular velocity one.
//========================================================================
template<class ELEMENT>
void EultingSphereProblem<ELEMENT>::set_boundary_conditions()
{
 //Set velocity for boundary 1 - outer wall (zero)
 unsigned ibound=1;
 
 // Loop over the nodes on boundary
 unsigned num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   nod_pt->set_value(0,0.0);
  }
 
 //Setting for boundary 3 (inner sphere that elutes the chemical)
 ibound=3;
 
 // Loop over the nodes on boundary
 num_nod=mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
   nod_pt->set_value(0,1.0);
  }
  
  
} // end of actions_before_implicit_timestep


//==start_of_destructor===================================================
/// Destructor for EultingSphere problem 
//========================================================================
template<class ELEMENT>
EultingSphereProblem<ELEMENT>::~EultingSphereProblem()
{ 
 //Delete the error estimator
 //delete mesh_pt()->spatial_error_estimator_pt();

 //Clean up the memory allocated for the mesh
 delete Problem::mesh_pt();

} // end_of_destructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void EultingSphereProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 
 
 
 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for Pe =" << Global_Physical_Variables::Pe 
      << std::endl;
 cout << "=================================================" << std::endl;


 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%g.dat",doc_info.directory().c_str(),
         Global_Physical_Variables::Pe);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 
 some_file.close();
} // end_of_doc_solution


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for EultingSphere test problem using a 
//  Crouzeix-Raviart interpolation.
//=====================================================================
int main()
{
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");

 // ---------------
 // end of Set up doc info

  EultingSphereProblem<
  QSphericalAdvectionDiffusionElement<3> > problem;

  //Set the boundary conditions
  problem.set_boundary_conditions();
  
  // Solve the problem
  problem.steady_newton_solve();
  
  problem.doc_solution(doc_info);

} // end_of_main
