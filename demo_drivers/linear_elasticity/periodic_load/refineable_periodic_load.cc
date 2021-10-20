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
// Driver for a periodically loaded elastic body -- spatially adaptive version

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
 /// Amplitude of traction applied
 double Amplitude = 1.0;

 /// Specify problem to be solved (boundary conditons for finite or
 /// infinite domain).
 bool Finite=false;

 /// Define Poisson coefficient Nu
 double Nu = 0.3;

 /// Length of domain in x direction
 double Lx = 1.0;

 /// Length of domain in y direction
 double Ly = 2.0;

 /// The elasticity tensor
 IsotropicElasticityTensor E(Nu);

 /// The exact solution for infinite depth case
 void exact_solution(const Vector<double> &x,
                     Vector<double> &u)
 {
  u[0] = -Amplitude*cos(2.0*MathematicalConstants::Pi*x[0]/Lx)*
	exp(2.0*MathematicalConstants::Pi*(x[1]-Ly))/
	(2.0/(1.0+Nu)*MathematicalConstants::Pi);
  u[1] = -Amplitude*sin(2.0*MathematicalConstants::Pi*x[0]/Lx)*
	exp(2.0*MathematicalConstants::Pi*(x[1]-Ly))/
	(2.0/(1.0+Nu)*MathematicalConstants::Pi);
 }

 /// The traction function
void periodic_traction(const double &time,
                       const Vector<double> &x,
                       const Vector<double> &n,
                       Vector<double> &result)
 {
  result[0] = -Amplitude*cos(2.0*MathematicalConstants::Pi*x[0]/Lx);
  result[1] = -Amplitude*sin(2.0*MathematicalConstants::Pi*x[0]/Lx);
 }
} // end_of_namespace


//===start_of_problem_class=============================================
/// Periodic loading problem
//======================================================================
template<class ELEMENT>
class RefineablePeriodicLoadProblem : public Problem
{
public:

 ///  Constructor: Pass number of elements in x and y directions 
 /// and lengths.
 RefineablePeriodicLoadProblem(const unsigned &nx, const unsigned &ny, 
                               const double &lx, const double &ly);

 /// Update before solve is empty
 void actions_before_newton_solve() {}

 /// Update after solve is empty
 void actions_after_newton_solve() {}

 /// Actions before adapt: Wipe the mesh of traction elements
 void actions_before_adapt()
  {
   // Kill the traction elements and wipe surface mesh
   delete_traction_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }

 /// Actions after adapt: Rebuild the mesh of traction elements
 void actions_after_adapt()
  {
   // Create traction elements
   assign_traction_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Allocate traction elements on the top surface
 void assign_traction_elements();
 
 /// Kill traction elements on the top surface
 void delete_traction_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_mesh_pt->nelement();
   
   // Loop over the traction elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_mesh_pt->flush_element_and_node_storage();
  }
 
 /// Pointer to the (refineable!) bulk mesh
 TreeBasedRefineableMeshBase* Bulk_mesh_pt;

 /// Pointer to the mesh of traction elements
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class


//===start_of_constructor=============================================
/// Problem constructor: Pass number of elements in the coordinate
/// directions and the domain sizes.
//====================================================================
template<class ELEMENT>
RefineablePeriodicLoadProblem<ELEMENT>::RefineablePeriodicLoadProblem
(const unsigned &nx, const unsigned &ny,
 const double &lx, const double& ly)
{
 // Create the mesh 
 Bulk_mesh_pt = new RefineableRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly);
 
 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Make the mesh periodic in the x-direction by setting the nodes on
 // right boundary (boundary 1) to be the periodic counterparts of
 // those on the left one (boundary 3).
 unsigned n_node = Bulk_mesh_pt->nboundary_node(1);
 for(unsigned n=0;n<n_node;n++)
  {
   Bulk_mesh_pt->boundary_node_pt(1,n)
    ->make_periodic(Bulk_mesh_pt->boundary_node_pt(3,n));
  }
 

 // Now establish the new neighbours (created by "wrapping around"
 // the domain) in the TreeForst representation of the mesh

 // Get pointers to tree roots associated with elements on the 
 // left and right boundaries
  Vector<TreeRoot*> left_root_pt(ny);
  Vector<TreeRoot*> right_root_pt(ny);
  for(unsigned i=0;i<ny;i++) 
   {
    left_root_pt[i] = 
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i*nx))->
     tree_pt()->root_pt();
    right_root_pt[i] = 
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(nx-1+i*nx))->
     tree_pt()->root_pt();
   }

  // Switch on QuadTreeNames for enumeration of directions
   using namespace QuadTreeNames;

  //Set the neighbour and periodicity
  for(unsigned i=0;i<ny;i++) 
   {
    // The western neighbours of the elements on the left
    // boundary are those on the right
    left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
    left_root_pt[i]->set_neighbour_periodic(W); 
    
    // The eastern neighbours of the elements on the right
    // boundary are those on the left
    right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
    right_root_pt[i]->set_neighbour_periodic(E);     
   } // done


 //Create the surface mesh of traction elements
 Surface_mesh_pt=new Mesh;
 assign_traction_elements();
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin & set the ones that have Dirichlet 
 // conditions here
 unsigned ibound=0;
 unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Get pointer to node
   Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

   // Pinned in x & y at the bottom and set value
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Check which boundary conditions to set and set them
   if (Global_Parameters::Finite)
     {
      // Set both displacements to zero
      nod_pt->set_value(0,0);
      nod_pt->set_value(1,0);
     }
   else
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);

      // Compute the value of the exact solution at the nodal point
      Vector<double> u(2);
      Global_Parameters::exact_solution(x,u);

      // Assign these values to the nodal values at this node
      nod_pt->set_value(0,u[0]);
      nod_pt->set_value(1,u[1]);
     };
  } // end_loop_over_boundary_nodes

 // Complete the problem setup to make the elements fully functional

 // Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   // Cast to a bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the elasticity tensor
   el_pt->elasticity_tensor_pt() = &Global_Parameters::E;
  }// end loop over elements


 // Do selective refinement of one element so that we can test
 // whether periodic hanging nodes work: Choose a single element 
 // (the zero-th one) as the to-be-refined element.
 // This creates a hanging node on the periodic boundary
 Vector<unsigned> refine_pattern(1,0);
 Bulk_mesh_pt->refine_selected_elements(refine_pattern);
  
 // Add the submeshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Now build the global mesh
 build_global_mesh();

 // Assign equation numbers
 cout << assign_eqn_numbers() << " equations assigned" << std::endl; 
} // end of constructor


//===start_of_traction===============================================
/// Make traction elements along the top boundary of the bulk mesh
//===================================================================
template<class ELEMENT>
void RefineablePeriodicLoadProblem<ELEMENT>::assign_traction_elements()
{

 // How many bulk elements are next to boundary 2 (the top boundary)?
 unsigned bound=2;
 unsigned n_neigh = Bulk_mesh_pt->nboundary_element(bound); 
 
 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   FiniteElement *traction_element_pt 
    = new LinearElasticityTractionElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
 
   // Add to mesh
   Surface_mesh_pt->add_element_pt(traction_element_pt);
  }
 
 // Now set function pointer to applied traction
 unsigned n_traction =  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_traction;e++)
  {
   // Cast to a surface element
   LinearElasticityTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<LinearElasticityTractionElement<ELEMENT>* >
    (Surface_mesh_pt->element_pt(e));
   
   // Set the applied traction
   el_pt->traction_fct_pt() = &Global_Parameters::periodic_traction;
  }// end loop over traction elements
 
 
} // end of assign_traction_elements

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineablePeriodicLoadProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 sprintf(filename,"%s/exact_soln.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,
                          Global_Parameters::exact_solution); 
 some_file.close();

 // Doc error
 double error=0.0;
 double norm=0.0;
 sprintf(filename,"%s/error.dat",doc_info.directory().c_str());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             Global_Parameters::exact_solution, 
                             error,norm);
 some_file.close();

// Doc error norm:
 cout << "\nNorm of error    " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;


} // end_of_doc_solution   


//===start_of_main======================================================
/// Driver code for PeriodicLoad linearly elastic problem
//======================================================================
int main(int argc, char* argv[]) 
{
 // Number of elements in x-direction
 unsigned nx=2; 
 
 // Number of elements in y-direction (for (approximately) square elements)
 unsigned ny=unsigned(double(nx)*Global_Parameters::Ly/Global_Parameters::Lx);
 
 // Set up doc info
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Set up problem
 RefineablePeriodicLoadProblem<RefineableQLinearElasticityElement<2,3> >
  problem(nx,ny,Global_Parameters::Lx, Global_Parameters::Ly);
 
 // Run the simulation, allowing for up to 4 spatial adaptations
 unsigned max_adapt=4; 
 problem.newton_solve(max_adapt);
 
 // Output the solution
 problem.doc_solution(doc_info);
  
} // end_of_main
