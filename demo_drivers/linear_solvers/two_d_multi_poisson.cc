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
//Driver for a simple 2D multi_poisson problem

#include<fenv.h>

//Generic routines
#include "generic.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

// The elements
#include "multi_poisson_elements.h"

// My own simple preconditioners
#include "multi_poisson_block_preconditioners.h"

using namespace std;

using namespace oomph;

// Number of fields
#define N_FIELD 5

//===== start_of_namespace=============================================
/// Namespace for exact solution for MultiPoisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForMultiPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=5.0;

 /// Parameter for angle Phi of "step" (45 degrees)
 double TanPhi=1.0;

 /// Interaction parameter
 double Beta=1.0;

 /// Scalar solution
 double scalar_soln(const Vector<double>& x)
 {
  return tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Scalar source
 double scalar_source(const Vector<double>& x)
 {
  return 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }


 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  double u_scalar=scalar_soln(x);
  for (unsigned i=0;i<N_FIELD;i++)
   {
    u[i]=double(i+1)*u_scalar;
   }
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x, Vector<double>& source)
 {
  double scalar_f=scalar_source(x);
  double scalar_u=scalar_soln(x);
  for (unsigned i=0;i<N_FIELD;i++)
   {
    source[i] = double(i+1)*scalar_f;
    for (unsigned k=0;k<N_FIELD;k++)
     {
      source[i]+=Beta*double(k+1)*scalar_u;
     }
    source[i]*=double(i+1);
   }
 }
 
} // end of namespace




//====== start_of_problem_class=======================================
/// 2D MultiPoisson problem on rectangular domain, discretised with
/// 2D QMultiPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class MultiPoissonProblem : public Problem
{

public:

 /// Constructor
 MultiPoissonProblem();

 /// Destructor (empty)
 ~MultiPoissonProblem(){}

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);


}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for MultiPoisson problem
//========================================================================
template<class ELEMENT>
MultiPoissonProblem<ELEMENT>::MultiPoissonProblem()
{ 
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     for (unsigned ii=0;ii<N_FIELD;ii++)
      {
       mesh_pt()->boundary_node_pt(i,n)->pin(ii);
      } 
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = &TanhSolnForMultiPoisson::source_function;

   //Set the pointer to the interaction parameter
   el_pt->beta_pt() = &TanhSolnForMultiPoisson::Beta;
  }

 // Create oomph-lib iterative linear solver
 GMRES<CRDoubleMatrix>* solver_pt=new GMRES<CRDoubleMatrix>;
 
 solver_pt->tolerance()=1.0e-12;

 // Set linear solver
 linear_solver_pt() = solver_pt;
 
 // Doc convergence
 solver_pt->open_convergence_history_file_stream
  ("RESLT/iterative_solver_convergence.dat");

 // Our own simple block diagonal block preconditioner
 if (CommandLineArgs::command_line_flag_has_been_set("--simple"))
  {
   oomph_info << "Using simple block triangular preconditioner\n";

   // Create preconditioner
   Diagonal<CRDoubleMatrix>* block_prec_pt=
    new Diagonal<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
 else if (CommandLineArgs::command_line_flag_has_been_set("--upper_triangular"))
  {
   oomph_info << "Using upper triangular preconditioner\n";

   // Create preconditioner
   UpperTriangular<CRDoubleMatrix>* block_prec_pt=
    new UpperTriangular<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());
   
   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  } 
else if (CommandLineArgs::
         command_line_flag_has_been_set("--two_plus_three"))
  {
   oomph_info << "Using 2+3 block diagonal preconditioner\n";

   // Create preconditioner
   TwoPlusThree<CRDoubleMatrix>* 
    block_prec_pt = new 
    TwoPlusThree<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());
   
   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set("--two_plus_three_upper_triangular"))
  {
   oomph_info << "Using 2+3 upper triangular preconditioner\n";

   // Create preconditioner
   TwoPlusThreeUpperTriangular<CRDoubleMatrix>* 
    block_prec_pt = new 
    TwoPlusThreeUpperTriangular<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set(
          "--two_plus_three_upper_triangular_with_sub"))
  {
   oomph_info << "Using 2+3 upper triangular with one level subsidiary"
              << " preconditioner\n";

   // Create preconditioner
   TwoPlusThreeUpperTriangularWithOneLevelSubsidiary<CRDoubleMatrix>* 
    block_prec_pt = new 
    TwoPlusThreeUpperTriangularWithOneLevelSubsidiary<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set(
          "--two_plus_three_upper_triangular_with_two_sub"))
  {
   oomph_info << "Using 2+3 block diagonal preconditioner with 2 levels of"
              << " subsidiaries.\n";

   // Create preconditioner
   TwoPlusThreeUpperTriangularWithTwoLevelSubsidiary<CRDoubleMatrix>* 
    block_prec_pt = new 
    TwoPlusThreeUpperTriangularWithTwoLevelSubsidiary<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set(
          "--two_plus_three_upper_triangular_with_replace"))
  {
   oomph_info << "Using 2+3 block diagonal preconditioner with replacement\n";

   // Create preconditioner
   TwoPlusThreeUpperTriangularWithReplace<CRDoubleMatrix>* 
    block_prec_pt = new 
    TwoPlusThreeUpperTriangularWithReplace<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set(
          "--coarse_two_plus_two_plus_one"))
  {
   oomph_info << "Using 2+2+1 block diagonal preconditioner with coarsening.\n";

   // Create preconditioner
   CoarseTwoPlusTwoPlusOne<CRDoubleMatrix>* 
    block_prec_pt = new 
    CoarseTwoPlusTwoPlusOne<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else if (CommandLineArgs::
         command_line_flag_has_been_set(
          "--one_plus_four_with_two_coarse"))
  {
   oomph_info << "Using 1+4 block diagonal preconditioner with \n"
              << "two level coarsening.\n";

   // Create preconditioner
   OnePlusFourWithTwoCoarse<CRDoubleMatrix>* 
    block_prec_pt = new 
    OnePlusFourWithTwoCoarse<CRDoubleMatrix>;
   
   // Set mesh
   block_prec_pt->set_multi_poisson_mesh(mesh_pt());

   // Set preconditioner
   solver_pt->preconditioner_pt()=block_prec_pt;
  }
else
  {
   // Create preconditioner
   MatrixBasedDiagPreconditioner* prec_pt=
    new MatrixBasedDiagPreconditioner;
      
   // Set preconditioner
   solver_pt->preconditioner_pt()=prec_pt;

   oomph_info << "Using diagonal preconditioner on entire system.\n";
  }
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

 // Block diagonal with one plus 4 repackaging of dofs into blocks
 // and block diagonal preconditioner for 4x4 sub-block
} // end of constructor




//=================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void MultiPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = mesh_pt()->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(i,n);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(N_FIELD);
     TanhSolnForMultiPoisson::get_exact_u(x,u);

     // Assign the value to nodal values at this node
     for (unsigned ii=0;ii<N_FIELD;ii++)
      {
       nod_pt->set_value(ii,u[ii]);
      }
    }
  } 
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void MultiPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForMultiPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForMultiPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc

 




//===== start_of_main=====================================================
/// Driver code for 2D MultiPoisson problem
//========================================================================
int main(int argc, char* argv[])
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 CommandLineArgs::specify_command_line_flag("--simple");

 CommandLineArgs::specify_command_line_flag
  ("--upper_triangular");

 CommandLineArgs::specify_command_line_flag
  ("--two_plus_three");

 CommandLineArgs::specify_command_line_flag
  ("--two_plus_three_upper_triangular");

 CommandLineArgs::specify_command_line_flag
  ("--two_plus_three_upper_triangular_with_sub");

 CommandLineArgs::specify_command_line_flag
  ("--two_plus_three_upper_triangular_with_two_sub");

 CommandLineArgs::specify_command_line_flag
  ("--two_plus_three_upper_triangular_with_replace");

 CommandLineArgs::specify_command_line_flag
  ("--coarse_two_plus_two_plus_one");

 CommandLineArgs::specify_command_line_flag
  ("--one_plus_four_with_two_coarse");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create the problem with 2D nine-node elements from the
 // QMultiPoissonElement family.
 MultiPoissonProblem<QMultiPoissonElement<2,3,N_FIELD> > problem;

 // Create label for output
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;
 
 oomph_info << "\n\nSolving for TanhSolnForMultiPoisson::Alpha="
            << TanhSolnForMultiPoisson::Alpha << std::endl << std::endl;
 
 // Solve the problem
 problem.newton_solve();
 
 //Output the solution
 problem.doc_solution(doc_info);

} //end of main









