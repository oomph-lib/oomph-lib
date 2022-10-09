//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Driver for solution of 2D Poisson equation in fish-shaped domain 
 
// Generic oomph-lib headers
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The fish mesh 
#include "meshes/fish_mesh.h"

using namespace std;

using namespace oomph;


//============ start_of_namespace=====================================
/// Namespace for const source term in Poisson equation
//====================================================================
namespace ConstSourceForPoisson
{ 
 
 /// Strength of source function: default value -1.0
 double Strength=-1.0;

/// Const source function
 void source_function(const Vector<double>& x, double& source)
 {
  source = Strength;
 }

} // end of namespace




//======start_of_problem_class========================================
///  Poisson problem in fish-shaped domain.
/// Template parameter identifies the element type.
//====================================================================
template<class ELEMENT>
class FishPoissonProblem : public Problem
{

public:

 /// Constructor
 FishPoissonProblem();

 /// Destructor: Empty
 virtual ~FishPoissonProblem(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 FishMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<FishMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Doc the solution. Output directory and labels are specified 
 /// by DocInfo object
 void doc_solution(DocInfo& doc_info);

}; // end of problem class





//===========start_of_constructor=========================================
/// Constructor for  Poisson problem in fish-shaped
/// domain.
//========================================================================
template<class ELEMENT>
FishPoissonProblem<ELEMENT>::FishPoissonProblem()
{ 
    
 // Build fish mesh -- this is a coarse base mesh consisting 
 // of four elements.
 Problem::mesh_pt()=new FishMesh<ELEMENT>;
   
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. Since the boundary values are never changed, we set
 // them here rather than in actions_before_newton_solve(). 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     // Pin the single scalar value at this node
     mesh_pt()->boundary_node_pt(i,n)->pin(0); 

     // Assign the homogenous boundary condition for the one and only
     // nodal value
     mesh_pt()->boundary_node_pt(i,n)->set_value(0,0.0); 
    }
  } 

 // Loop over elements and set pointers to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::source_function;
  }

 // Setup the equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//=======start_of_doc=====================================================
/// Doc the solution in tecplot format.
//========================================================================
template<class ELEMENT>
void FishPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points in each coordinate direction.
 unsigned npts;
 npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output solution 
 sprintf(filename,"%s/soln_nodes%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,4);
 some_file.close();

 // Output solution 
 sprintf(filename,"%s/soln_fine%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,20*npts);
 some_file.close();


 // Output boundaries
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();
 
} // end of doc

 






//=================start_of_main==========================================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain.
//========================================================================
int main()
{


  //Set up the problem with nine-node  Poisson elements
  FishPoissonProblem<QPoissonElement<2,3> > problem;
  
  // Setup labels for output
  //------------------------
  DocInfo doc_info;
  
  // Set output directory
  doc_info.set_directory("RESLT"); 
  
  // Step number
  doc_info.number()=0;


  
  // Solve/doc the problem
  //----------------------

  // Solve the problem
  problem.newton_solve();
  
  //Output solution
  problem.doc_solution(doc_info);   

} // end of main


