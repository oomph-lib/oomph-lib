//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
//Driver for flow past rectangular box -- meshed with triangle

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

// The mesh
#include "meshes/triangle_mesh.h"

using namespace std;

using namespace oomph; 

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=10.0; 

} // end_of_namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Flow past a box in a channel
//====================================================================
template<class ELEMENT>
class FlowPastBoxProblem : public Problem
{

public:


 /// Constructor: Pass filenames for mesh
 FlowPastBoxProblem(const string& node_file_name,
                    const string& element_file_name,
                    const string& poly_file_name);

 /// Destructor (empty)
 ~FlowPastBoxProblem(){}

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve. 
 /// Re-set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
 {
  // Find max. and min y-coordinate at inflow
  double y_min=1.0e20;
  double y_max=-1.0e20;
  unsigned ibound=0;
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
    if (y>y_max)
     {
      y_max=y;
     }
    if (y<y_min)
     {
      y_min=y;
     }
   }

  // Loop over all boundaries
  for (unsigned ibound=0;ibound<mesh_pt()->nboundary();ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      // Flow towards the north east on outer boundary (boundary 1)
      if (ibound==1)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,-1.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,-1.0);
       }
      // Zero flow elsewhere 
      else
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
       }
     }
   }
 } // end_of_actions_before_newton_solve


 /// Access function for the specific mesh
 TriangleMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<TriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for FlowPastBox problem. Pass filenames for mesh
//========================================================================
template<class ELEMENT>
FlowPastBoxProblem<ELEMENT>::FlowPastBoxProblem(
 const string& node_file_name,
 const string& element_file_name,
 const string& poly_file_name)
{ 

 Problem::Max_residuals=1000.0;

 //Create mesh
 Problem::mesh_pt() = new TriangleMesh<ELEMENT>(node_file_name,
                                                element_file_name,
                                                poly_file_name);
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin velocity everywhere
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
    }
  } // end loop over boundaries


 // Pin the zero-th pressure dof in the zero-th element 
 // and set it to zero
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->
                          fix_pressure(0,0.0);


 // Complete the build of all elements so they are fully functional

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
  } // end loop over elements


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FlowPastBoxProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
} // end_of_doc_solution





/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////







//==start_of_main======================================================
/// Driver for FlowPastBox test problem
//=====================================================================
int main(int argc, char* argv[])
{

// Store command line arguments
 CommandLineArgs::setup(argc,argv);


 // Check number of command line arguments: Need exactly two.
 if (argc!=4)
  {
   std::string error_message =
    "Wrong number of command line arguments.\n";
   error_message +=
    "Must specify the following file names  \n";
   error_message += 
    "filename.node then filename.ele then filename.poly\n";
   
   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }


 // Convert arguments to strings that specify the input file names
 string node_file_name(argv[1]);
 string element_file_name(argv[2]);
 string poly_file_name(argv[3]);


 // Set up doc info
 // ---------------

 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;


// Doing QTaylorHoodElements
 {
  // Build the problem with TTaylorHoodElements
  FlowPastBoxProblem<TTaylorHoodElement<2> > problem(node_file_name,
                                                     element_file_name,
                                                     poly_file_name);
  // Output boundaries 
  problem.mesh_pt()->output_boundaries("RESLT/boundaries.dat");

  // Outpt the solution
  problem.doc_solution(doc_info);
  
  // Step number
  doc_info.number()++;
  
  // Solve the problem
  problem.newton_solve();
  
  // Outpt the solution
  problem.doc_solution(doc_info);
  
  // Step number
  doc_info.number()++;
  
 } // end of QTaylorHoodElements
 

} // end_of_main










