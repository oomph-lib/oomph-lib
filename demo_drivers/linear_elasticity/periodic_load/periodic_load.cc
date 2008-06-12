//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
// Driver for a periodically loaded elastic body

// The oomphlib headers
#include "generic.h"
#include "linear_elasticity.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===start_of_namespace=================================================
/// Namespace for global parameters
//======================================================================
namespace Global_Parameters
{
 const double pi = 4.0*atan(1.0);

 double Amplitude = 1.0;

 double Lambda = 0.7;

 double Mu = 1.0;

 double Lx = 1.0;

 double Ly = 5.0;

 /// The elasticity tensor
 IsotropicElasticityTensor E(Lambda,Mu);

 /// The exact solution for infinite depth case
 void exact_solution(const Vector<double> &x,
                     Vector<double> &u)
 {
  u[0] = -cos(2.0*pi*x[0])*exp(2.0*pi*(x[1]-Ly))/(4.0*Mu*pi);
  u[1] = -sin(2.0*pi*x[0])*exp(2.0*pi*(x[1]-Ly))/(4.0*Mu*pi);
 }

 /// The traction function
 void periodic_traction(const Vector<double> &x,
                        const Vector<double> &n,
                        Vector<double> &result)
 {
  //The sinusoidal but
  double t22 = -Amplitude*sin(2.0*pi*x[0]);
  double t12 = - Amplitude*cos(2.0*pi*x[0]);

  //It's a vertical traction
  result[0] = t12*n[1];
  result[1] = t22*n[1];
 }

} // end of namespace


//===start_of_problem_class=============================================
/// Rayleigh-type problem: 2D channel whose upper
/// wall oscillates periodically.
//======================================================================
template<class ELEMENT>
class PeriodicLoadProblem : public Problem
{
public:

 /// Constructor: Pass number of elements in x and y directions and 
 /// lengths
 PeriodicLoadProblem(const unsigned &nx, const unsigned &ny, 
                     const double &lx, const double &ly);

 //Update before solve is empty
 void actions_before_newton_solve() {}

 /// \short Update after solve is empty
 void actions_after_newton_solve() {}

 /// \short Allocate traction elements on the top surface
 void assign_traction_elements(const unsigned &bound);
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Pointer to the bulk mesh
 Mesh* Bulk_mesh_pt;

 /// Pointer to the mesh of traction elements
 Mesh* Surface_mesh_pt;

}; // end of problem_class


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT>
PeriodicLoadProblem<ELEMENT>::PeriodicLoadProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& ly)
{
 //Initialise the surface mesh to zero
 Surface_mesh_pt=0;
 //Now create the mesh with periodic boundary conditions in x direction
 bool periodic_in_x=true;
 Bulk_mesh_pt = 
  new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,periodic_in_x);

 //Now make the surface traction elements
 assign_traction_elements(2);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 unsigned num_bound=Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pinned in y at the bottom
     if(ibound==0)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);

       //Pin only at corners
       /*Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
         if((nod_pt->is_on_boundary(1)) || (nod_pt->is_on_boundary(3)))
         {
         nod_pt->pin(1);
         }*/
       
      }
     // No horizontal motion on left or right
     /*else if ((ibound==1)||(ibound==3))
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       }*/
    }
  } // end loop over boundaries

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->elasticity_tensor_pt() = &Global_Parameters::E;
  }

 //Loop over the traction elements
 unsigned n_traction =  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_traction;e++)
  {
   //Cast to a fluid element
   LinearElasticityTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<LinearElasticityTractionElement<ELEMENT>* >
    (Surface_mesh_pt->element_pt(e));
   
   //Set the Reynolds number, etc
   el_pt->traction_fct_pt() = &Global_Parameters::periodic_traction;
  }



 //Add the submeshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 //Now build the global mesh
 build_global_mesh();

 //Assgn equation numbers
 cout << assign_eqn_numbers() << std::endl; 
} // end of constructor


//===================================================================
/// Make traction elements along the top boundary of the bulk mesh
//===================================================================
template<class ELEMENT>
void PeriodicLoadProblem<ELEMENT>::assign_traction_elements(
 const unsigned &bound)
{
 //Clear the existing surface mesh and flush it
 if(Surface_mesh_pt!=0) 
  {
   Surface_mesh_pt->flush_element_and_node_storage();
   delete Surface_mesh_pt;
  }

 //Now make a new mesh
 Surface_mesh_pt = new Mesh;

 //Now set up the neighbouring element schemes
 unsigned n_neigh = Bulk_mesh_pt->nboundary_element(bound);
 

 //Now loop over neighbours and create the free surface elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   //Create the free surface element
   FiniteElement *traction_element_pt 
    = new LinearElasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
   //Push it back onto the stack
   Surface_mesh_pt->add_element_pt(traction_element_pt);
  }
 
 cout << Surface_mesh_pt->nelement() 
      << " free surface elements assigned" << std::endl;
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PeriodicLoadProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Doc error
 //----------
 /*double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,
                          Global_Parameters::exact_solution,
                          error,norm); 
                          some_file.close();*/

} // end_of_doc_solution   


//===start_of_main======================================================
/// Driver code for PeriodicLoad channel problem
//======================================================================
int main(int argc, char* argv[]) 
{
 // Number of elements in x-direction
 unsigned nx=10;

 // Number of elements in y-direction
 unsigned ny=10;

 // Solve with Crouzeix-Raviart elements
 {
  // Set up doc info
  DocInfo doc_info;
  doc_info.number()=0;
  doc_info.set_directory("RESLT_fixed");
  
  //Set up problem
  PeriodicLoadProblem<QLinearElasticityElement<2,3> > 
   problem(nx,ny,Global_Parameters::Lx, Global_Parameters::Ly);
  
  // Run the unsteady simulation
  problem.newton_solve();
  doc_info.number()++;
  //Print the answer
  problem.doc_solution(doc_info);
 }

} // end of main
