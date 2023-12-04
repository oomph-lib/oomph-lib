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
/// Driver for a 3D navier stokes flow with tetgen

//Generic routines
#include "generic.h"
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

using namespace std;

using namespace oomph;

//=start_of_namespace================================================
/// Namespace for physical parameters
//===================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=10.0;
} // end_of_namespace




//=start_of_problem_class=============================================
/// Entry flow problem in quarter tube domain
//====================================================================
template<class ELEMENT>
class NavierStokesProblem : public Problem
{

public:

 /// Constructor: Pass DocInfo object and file names
 NavierStokesProblem(DocInfo& doc_info, 
                  const string& node_file_name,
                  const string& element_file_name,
                  const string& face_file_name);

 /// Destructor (empty)
 ~NavierStokesProblem() {}

 /// Doc the solution after solve
 void actions_after_newton_solve() 
  {
   // Doc solution after solving
   doc_solution();

   // Increment label for output files
   Doc_info.number()++;
  }

 /// Update the problem specs before solve 
 void actions_before_newton_solve();

 /// Doc the solution
 void doc_solution();

 //Access function for the specific mesh
TetgenMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<TetgenMesh<ELEMENT>*>(Problem::mesh_pt());
  }


private:

 
 /// Doc info object
 DocInfo Doc_info;

}; // end_of_problem_class




//=start_of_constructor===================================================
/// Constructor: Pass DocInfo object mesh files
//========================================================================
template<class ELEMENT>
NavierStokesProblem<ELEMENT>::NavierStokesProblem(DocInfo& doc_info,
                                            const string& node_file_name,
                                            const string& element_file_name,
                                            const string& face_file_name)
 : Doc_info(doc_info)
{ 
 //Create mesh
 Problem::mesh_pt() = new TetgenMesh<ELEMENT>(node_file_name,
                                              element_file_name,
                                              face_file_name);

 //Doc the boundaries
 ofstream some_file;
 char filename[100];
 sprintf(filename,"boundaries.dat");
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();
 
 
 // Set the boundary conditions for this problem: All nodal values are
 // free by default -- just pin the ones that have Dirichlet conditions

 // Pin transverse velocities on outer walls -- the outer boundary
 // behaves like a channel (open along the x-axis) with slippery walls 
 // on the transverse boundaries
 {
  unsigned ibound=0;
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(2);
    }
   }
 }


 // Pin all velocity components on internal block -- behaves like 
 // a rigid body
 {
  unsigned ibound=1;
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(2);
    }
   }
 }
   

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  }


 //Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=start_of_actions_before_newton_solve==========================================
/// Set the inflow boundary conditions
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::actions_before_newton_solve()
{

 // Apply conditions on obstacle
 unsigned ibound=1;
 unsigned num_nod= mesh_pt()->nboundary_node(ibound); 
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Block moves in z-direction
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(2,1.0);
   
  }

} // end_of_actions_before_newton_solve


//=start_of_doc_solution==================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::doc_solution()
{ 
 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution

 


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=start_of_main=======================================================
/// 3D Navier Stokes on an unstructured mesh
//=====================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Check number of command line arguments: Need exactly three.
 if (argc!=4)
  {
   std::string error_message =
    "Wrong number of command line arguments.\n";
   error_message +=
    "Must specify the following file names  \n";
   error_message += 
    "filename.node then filename.ele then filename.face\n";

   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Convert arguments to strings that specify the input file names
 string node_file_name(argv[1]);
 string element_file_name(argv[2]);
 string face_file_name(argv[3]);

 // Set up doc info
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Build problem
 NavierStokesProblem<TTaylorHoodElement<3> > 
  problem(doc_info,node_file_name,element_file_name,face_file_name);
 
 // Solve the problem 
 problem.newton_solve();


} // end_of_main


