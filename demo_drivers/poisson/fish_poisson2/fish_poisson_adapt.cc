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
// Driver for solution of 2D Poisson equation in fish-shaped domain with
// adaptivity
 
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
 void get_source(const Vector<double>& x, double& source)
 {
  source = Strength;
 }

} // end of namespace




//======start_of_problem_class========================================
/// Refineable Poisson problem in fish-shaped domain.
/// Template parameter identifies the element type.
//====================================================================
template<class ELEMENT>
class RefineableFishPoissonProblem : public Problem
{

public:

 /// Constructor
 RefineableFishPoissonProblem();

 /// Destructor: Empty 
 virtual ~RefineableFishPoissonProblem(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 ///  Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableFishMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableFishMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 ///  Doc the solution. Output directory and labels are specified 
 /// by DocInfo object
 void doc_solution(DocInfo& doc_info);

}; // end of problem class





//===========start_of_constructor=========================================
/// Constructor for adaptive Poisson problem in fish-shaped
/// domain.
//========================================================================
template<class ELEMENT>
RefineableFishPoissonProblem<ELEMENT>::RefineableFishPoissonProblem()
{ 
    
 // Build fish mesh -- this is a coarse base mesh consisting 
 // of four elements. We'll refine/adapt the mesh later.
 Problem::mesh_pt()=new RefineableFishMesh<ELEMENT>;
 
 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. Since the boundary values are never changed, we set
 // them here rather than in actions_before_newton_solve(). 
 unsigned num_bound = mesh_pt()->nboundary();

 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound); 
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin the single scalar value at this node
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 

     // Assign the homogenous boundary condition to the one
     // and only nodal value
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0); 
    }
  }

 // Loop over elements and set pointers to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::get_source;
  }

 // Setup the equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//=======start_of_doc=====================================================
/// Doc the solution in tecplot format.
//========================================================================
template<class ELEMENT>
void RefineableFishPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 
} // end of doc

 




//=====================start_of_incremental===============================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain with mesh adaptation. First we solve on the original
/// coarse mesh. Next we do a few uniform refinement steps and resolve.
/// Finally, we enter into an automatic adapation loop.
//========================================================================
void solve_with_incremental_adaptation()
{
 
 //Set up the problem with 4 node refineable Poisson elements
 RefineableFishPoissonProblem<RefineableQPoissonElement<2,2> > problem;
 
 // Setup labels for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT"); 
 
 // Step number
 doc_info.number()=0;
  
  
 // Doc (default) refinement targets
 //----------------------------------
 problem.mesh_pt()->doc_adaptivity_targets(cout);
 
 // Solve/doc the problem on the initial, very coarse mesh
 //-------------------------------------------------------
 
 // Solve the problem
 problem.newton_solve();
 
 //Output solution
 problem.doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;


 // Do three rounds of uniform mesh refinement and re-solve
 //--------------------------------------------------------
 unsigned n_uniform=3;
 for (unsigned isolve=0;isolve<n_uniform;isolve++)
  {
   
   // Refine the problem uniformly
   problem.refine_uniformly();
   
   // Solve the problem 
   problem.newton_solve();
   
   //Output solution
   problem.doc_solution(doc_info);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }

 
 // Now do (up to) four rounds of fully automatic adapation in response to 
 //-----------------------------------------------------------------------
 // error estimate
 //---------------
 unsigned max_solve=4;
 for (unsigned isolve=0;isolve<max_solve;isolve++)
  {
   // Adapt problem/mesh
   problem.adapt(); 
         
   // Re-solve the problem if the adaptation has changed anything
   if ((problem.mesh_pt()->nrefined()  !=0)||
       (problem.mesh_pt()->nunrefined()!=0))
    {
     problem.newton_solve();
    }
   else
    {
     cout << "Mesh wasn't adapted --> we'll stop here" << std::endl;
     break;
    }
   
   //Output solution
   problem.doc_solution(doc_info);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }
 
 
} // end of incremental


//=================start_of_main==========================================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain with mesh adaptation.
//========================================================================
int main()
{
 // Solve with adaptation, docing the intermediate steps
 solve_with_incremental_adaptation();
 
} // end of main

