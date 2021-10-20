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
//Driver for 2D moving block

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "poisson.h"

// Header
#include "meshes/quad_from_triangle_mesh.h"

using namespace std;
using namespace oomph;
 

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=50.0;

} // end_of_namespace


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_problem_class============================================
/// Moving block problem
//====================================================================
template<class ELEMENT,class MESH>
class MovingBlockProblem : public Problem
{

public:

 /// Constructor
 MovingBlockProblem();

 /// Destructor (empty)
 ~MovingBlockProblem(){}

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e,const unsigned &pdof,const double &pvalue)
  {
   // Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->fix_pressure(pdof,pvalue);
  } // end of fix_pressure

 ///  Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  {   
   // No flow along the boundaries
   unsigned num_bound=mesh_pt()->nboundary();
   for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod=mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
     for (unsigned i=0;i<2;i++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
     }
    }
   }
   
   // Setup vertical flow along boundary 1:
   unsigned ibound=1; 
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Vertical flow
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,1.0);
   }
   
  } // end_of_actions_before_newton_solve
 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}
 
 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to the "bulk" mesh
 MESH* Bulk_mesh_pt;
 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for MovingBlock problem 
//========================================================================
template<class ELEMENT,class MESH>
MovingBlockProblem<ELEMENT,MESH>::MovingBlockProblem()
{ 
 // Set the maximum residuals value
 Problem::Max_residuals=1000.0;

 // Convert arguments to strings that specify the input file names
 string node_file_name("box_hole.1.node");
 string element_file_name("box_hole.1.ele");
 string poly_file_name("box_hole.1.poly");

 // Create the bulk mesh
 Bulk_mesh_pt=new MESH(node_file_name,element_file_name,poly_file_name);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=4e-07;
 Bulk_mesh_pt->max_permitted_error()=1e-06;
 
 // Create the main mesh
 add_sub_mesh(Bulk_mesh_pt);
 
 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Output the mesh (just to make sure everything looks good!)
 unsigned npts=2;
 std::ofstream outfile;
 outfile.open("RESLT/mesh.dat");
 mesh_pt()->output(outfile,npts);
 outfile.close();
 mesh_pt()->output_boundaries("RESLT/mesh_boundaries.dat");
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Loop over values (u and v velocities)
   for (unsigned i=0;i<2;i++)
   {
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
   }
  }
 } // end loop over boundaries


 // Complete the build of all elements so they are fully functional:

 // Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to the present element
  ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

  // Set the Reynolds number
  el_pt->re_pt() = &Global_Physical_Variables::Re;
 } // end loop over elements

 // Now set the first pressure value in element 0 to 0.0
 fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;
 
} // end_of_constructor


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT,class MESH>
void MovingBlockProblem<ELEMENT,MESH>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=2; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
	 doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for MovingBlock test problem -- test drive
/// with two different types of element.
//=====================================================================
int main()
{
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;
 
 // Typedef the ELEMENT and MESH type
 typedef RefineableQCrouzeixRaviartElement<2> ELEMENT;
 typedef RefineableQuadFromTriangleMesh<ELEMENT> MESH;
 
 // Build the problem with QCrouzeixRaviartElements
 MovingBlockProblem<ELEMENT,MESH> problem;
 
 std::cout << "Doing QCrouzeixRaviartElement<2>" << std::endl;

 // Maximum number of adaptations
 unsigned max_adapt=2;
 
 // Solve the problem
 problem.newton_solve(max_adapt);
 
 // Output the solution
 problem.doc_solution(doc_info);
 
} // end_of_main

