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
//Driver using simple 2D poisson problem

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;
using namespace oomph;


//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step"
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=10.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }
 
} // end of namespace





//====== start_of_problem_class=======================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~PoissonProblem()
  {
   delete mesh_pt()->spatial_error_estimator_pt();
   delete mesh_pt();
  }

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 
 /// Bump up counter for number of Newton iterations
 void actions_before_newton_convergence_check()
  {
   Newton_count++;
  }

 /// Report number of Newton iterations taken
 unsigned newton_count(){return Newton_count-1;}

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

 /// Number of Newton iterations
 unsigned Newton_count;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
 :  Source_fct_pt(source_fct_pt)
{
 // The problem is linear
 //Problem_is_nonlinear = false;

 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x=10;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt() =
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here.
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     mesh_pt()->boundary_node_pt(i,n)->pin(0);
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
   el_pt->source_fct_pt() = Source_fct_pt;
  }


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl;

} // end of constructor




//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: Set nodal values to zero and
/// (Re-)set boundary conditions to the values from the exact solution
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::actions_before_newton_solve()
{

 // Reset counter
 Newton_count=0;

 // Loop over nodes and set values to zero
 unsigned n_node=mesh_pt()->nnode();
 for (unsigned j=0; j<n_node; j++)
  {
   Node* node_pt=mesh_pt()->node_pt(j);
   node_pt->set_value(0,0.0);
  }
  
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
     Vector<double> u(1);
     TanhSolnForPoisson::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  } 
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln_trilinos_%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << endl;
 cout << "Norm of solution: " << sqrt(norm) << endl << endl;

} // end of doc






//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


  // Note: there is a memory leak in trilinos itself as 
  // demonstrated by this bit of code!
  if (false)
   {
    // set up a parameter list
    Teuchos::ParameterList ml_parameters;
    ML_Epetra::SetDefaults("SA", ml_parameters);
    
    std::cout << "hello" << std::endl;
    exit(0);
   }

 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;

 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0;
  
 // Create the problem with 2D nine-node elements from the
 // QPoissonElement family. Pass pointer to source function.
 PoissonProblem<RefineableQPoissonElement<2,2> >
  problem(&TanhSolnForPoisson::source_function);   


 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // File to report the number of Newton iterations
 ofstream conv_file;
 char filename[100];
 sprintf(filename,"%s/conv.dat",doc_info.directory().c_str());
 conv_file.open(filename);
 

 // Trilinos solvers with the ML preconditioner
 //--------------------------------------------
 {

  cout << "==================================="
       << "=================================" << endl;
  cout << "Testing Trilinos solver preconditioned with "
       << "TrilinosMLPreconditioner" << endl;
  cout << "======================================="
       << "=============================" << endl;

  // Create a Trilinos Solver 
  TrilinosAztecOOSolver* linear_solver_pt = new TrilinosAztecOOSolver;
  
  // Create the Trilinos ML preconditioner
  TrilinosMLPreconditioner* preconditioner_pt = new TrilinosMLPreconditioner;

  // Set the preconditioner pointer
  linear_solver_pt->preconditioner_pt() = preconditioner_pt;

  // Set linear solver
  problem.linear_solver_pt() = linear_solver_pt;
  
  linear_solver_pt->enable_doc_time();
  
  // Loop over specific built-in Trilinos solvers
  for (unsigned s=0; s<=2; s++)
   {
    // Select the linear solver
    switch(s)
     {
     case 0: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::CG;
      cout << "Using CG as solver" << endl;
      cout << "==================" << endl;
      break;
     case 1: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
      cout << "Using GMRES as solver" << endl;
      cout << "=====================" << endl;
      break;
     case 2: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::BiCGStab;
      cout << "Using BiCGStab as solver" << endl;
      cout << "========================" << endl;
      break;
     default: cout << "Error selecting solver" << endl;
     } 
    
    // Solve the problem
     problem.newton_solve(); 
    
     // Doc
     problem.doc_solution(doc_info);
     doc_info.number()++;
     oomph_info << "Number of linear solver iterations: " 
                << dynamic_cast<TrilinosAztecOOSolver*>(
                 linear_solver_pt)->iterations() << std::endl << std::endl;
     conv_file << "Number of linear solver iterations: " 
               << dynamic_cast<TrilinosAztecOOSolver*>(
                linear_solver_pt)->iterations() << std::endl;

   } // end of loop over preconditioner
 
  // delete the Trilinos solver
  delete linear_solver_pt;
  delete preconditioner_pt;  
 }
 
 

 // Trilinos solvers with the IFPACK preconditioner
 //------------------------------------------------
 {

  cout << "================================="
       << "=======================================" << endl;
  cout << "Testing Trilinos solver preconditioned" 
       << " with TrilinosIFPACKPreconditioner" << endl;
  cout << "==========================================="
       << "=============================" << endl;

    // Create a Trilinos Solver 
  TrilinosAztecOOSolver* linear_solver_pt = new TrilinosAztecOOSolver;
  
  // create the Trilinos IFPACK preconditioner
  TrilinosIFPACKPreconditioner* preconditioner_pt = 
   new TrilinosIFPACKPreconditioner;

  // set the preconditioner pointer
  linear_solver_pt->preconditioner_pt() = preconditioner_pt;
  
  // Set linear solver
  problem.linear_solver_pt() = linear_solver_pt;
  
  linear_solver_pt->enable_doc_time();
  
  // Loop over specific built-in Trilinos solvers
  for (unsigned s=0; s<=2; s++)
   {
    // Select the linear solver
    switch(s)
     {
     case 0: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::CG;
      cout << "Using CG as solver" << endl;
      cout << "==================" << endl;
      break;
     case 1: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
      cout << "Using GMRES as solver" << endl;
      cout << "=====================" << endl;
      break;
     case 2: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::BiCGStab;
      cout << "Using BiCGStab as solver" << endl;
      cout << "========================" << endl;
      break;
     default: cout << "Error selecting solver" << endl;
     } 
    
    // Solve the problem
     problem.newton_solve();

     // Doc
     problem.doc_solution(doc_info);
     doc_info.number()++;
     oomph_info << "Number of linear solver iterations: " 
                << dynamic_cast<TrilinosAztecOOSolver*>(
                 linear_solver_pt)->iterations() << std::endl << std::endl;
     conv_file << "Number of linear solver iterations: " 
               << dynamic_cast<TrilinosAztecOOSolver*>(
                linear_solver_pt)->iterations() << std::endl;
    
   } // end of loop over preconditioner
 
  // delete the Trilinos solver
  delete linear_solver_pt;
  delete preconditioner_pt;
 }
 
 



 // Test the oomph-lib iterative linear solvers with Trilinos preconditioners
 //--------------------------------------------------------------------------
 
 cout
  << "======================================" 
  << "==================================" 
  << endl;
 cout 
  << "Testing oomph-lib iterative linear solvers " 
  << "with Trilinos preconditioners"
  << endl;
 cout
  << "=======================================" 
  << "=================================" 
  << endl  << endl  << endl;

 

 // Pointer to linear solver
 LinearSolver* linear_solver_pt=0;

 // Loop over two Trilonos preconditioners
 for (unsigned i_prec=0;i_prec<2;i_prec++)
  {

   // Pointer to preconditioner
   Preconditioner* preconditioner_pt=0;

   if (i_prec==0)
    {
     // Create a IFPACK preconditioner
     preconditioner_pt = new TrilinosIFPACKPreconditioner;
     cout << "Using  TrilinosIFPACKPreconditioner as preconditioner" << endl;
     cout << "=====================================================" << endl;
    }
   else
    {
     // Create a ML preconditioner
     preconditioner_pt = new TrilinosMLPreconditioner;
     cout << "Using  TrilinosMLPreconditioner as preconditioner" << endl;
     cout << "=================================================" << endl;
    }
 

   // Test CG
   //--------
   cout << "CG iterative solver" << endl;
   cout << "-------------------" << endl;
   linear_solver_pt = new CG<CRDoubleMatrix>;
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->tolerance()=1.0e-10;
   problem.linear_solver_pt() = linear_solver_pt;
   
   // set the preconditioner in the iterative solver
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->
    preconditioner_pt()= preconditioner_pt;
   
   // Solve the problem
   problem.newton_solve(); 
  

   // Doc
   problem.doc_solution(doc_info);
   doc_info.number()++;
   oomph_info << "Number of linear solver iterations: " 
              << dynamic_cast<CG<CRDoubleMatrix>*>(
               linear_solver_pt)->iterations() << std::endl << std::endl;
   conv_file <<  "Number of linear solver iterations: " 
             << dynamic_cast<CG<CRDoubleMatrix>*>(
              linear_solver_pt)->iterations() << std::endl;

   delete linear_solver_pt;
   
   
   // Test BiCGStab
   //--------------
   cout << "BiCGStab iterative solver" << endl;
   cout << "-------------------------" << endl;
   linear_solver_pt = new BiCGStab<CRDoubleMatrix>;
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->tolerance()=1.0e-10;
   problem.linear_solver_pt() = linear_solver_pt;
   
   // set the preconditioner in the iterative solver
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->preconditioner_pt()
    = preconditioner_pt;
   
   // Solve the problem
   problem.newton_solve(); 
  
   // Doc
   problem.doc_solution(doc_info);
   doc_info.number()++;
   oomph_info << "Number of linear solver iterations: " 
              << dynamic_cast<BiCGStab<CRDoubleMatrix>*>(
               linear_solver_pt)->iterations() << std::endl << std::endl;
   conv_file << "Number of linear solver iterations: " 
             << dynamic_cast<BiCGStab<CRDoubleMatrix>*>(
              linear_solver_pt)->iterations() << std::endl;

   delete linear_solver_pt;

   
   // Test GMRES
   //-----------
   cout << "GMRES iterative solver" << endl;
   cout << "----------------------" << endl;
   linear_solver_pt = new GMRES<CRDoubleMatrix>;
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->tolerance()=1.0e-10;
   problem.linear_solver_pt() = linear_solver_pt;
     
   // set the preconditioner in the iterative solver
   static_cast<IterativeLinearSolver*>(linear_solver_pt)->preconditioner_pt()
    = preconditioner_pt;
  
   // Solve the problem
   problem.newton_solve();

   // Doc
   problem.doc_solution(doc_info);
   doc_info.number()++;
   oomph_info << "Number of linear solver iterations: " 
              << dynamic_cast<GMRES<CRDoubleMatrix>*>(
               linear_solver_pt)->iterations() << std::endl << std::endl;
   conv_file << "Number of linear solver iterations: " 
             << dynamic_cast<GMRES<CRDoubleMatrix>*>(
              linear_solver_pt)->iterations() << std::endl;

   delete linear_solver_pt;

   // Kill preconditioner
   delete preconditioner_pt;
   preconditioner_pt=0;

  } // end of loop over preconditioners


 // Trilinos solvers with the OOMPH-lib
 //------------------------------------------------
 {

  cout << "==================================" 
       << "=======================================" << endl;
  cout << "Testing Trilinos solver preconditioned with" 
       << " MatrixBasedDiagPreconditioner" << endl;
  cout << "=======================================" 
       << "==================================" << endl;

  // Create a Trilinos Solver 
  TrilinosAztecOOSolver* linear_solver_pt = new TrilinosAztecOOSolver;
  
  // create the oomphlib diagonal preconditioner
  Preconditioner* preconditioner_pt = 
   new MatrixBasedDiagPreconditioner;

  // set the preconditioner pointer
  linear_solver_pt->preconditioner_pt() = preconditioner_pt;

  // Set linear solver
  problem.linear_solver_pt() = linear_solver_pt;
  
  linear_solver_pt->enable_doc_time();
  
  // Loop over specific built-in Trilinos solvers
  for (unsigned s=0; s<=2; s++)
   {
        // Select the linear solver
    switch(s)
     {
     case 0: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::CG;
      cout << "Using CG as solver" << endl;
      cout << "==================" << endl;
      break;
     case 1: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
      cout << "Using GMRES as solver" << endl;
      cout << "=====================" << endl;
      break;
     case 2: 
      linear_solver_pt->solver_type() = TrilinosAztecOOSolver::BiCGStab;
      cout << "Using BiCGStab as solver" << endl;
      cout << "========================" << endl;
      break;
     default: cout << "Error selecting solver" << endl;
     } 
    
    // Solve the problem
    problem.newton_solve();
    
    // Doc
    problem.doc_solution(doc_info);
    doc_info.number()++;
    oomph_info << "Number of linear solver iterations: " 
               << dynamic_cast<TrilinosAztecOOSolver*>(
                linear_solver_pt)->iterations() << std::endl << std::endl;
    conv_file << "Number of linear solver iterations: " 
              << dynamic_cast<TrilinosAztecOOSolver*>(
               linear_solver_pt)->iterations() << std::endl;
    
   } // end of loop over preconditioner
 
  // delete the Trilinos solver
  delete linear_solver_pt;

  // Delete oomph-lib preconditioner
  delete preconditioner_pt;  
 }

 conv_file.close();

 // This code is not used, just compiled -- don't delete as it's
 // extracted by doxygen for a tutorial
 bool never_get_here=true;
 if (!never_get_here)
  {
   // Create oomph-lib linear solver
   IterativeLinearSolver* linear_solver_pt=new GMRES<CRDoubleMatrix>;

   // Create Trilinos IFPACK preconditioner as oomph-lib Preconditioner
   Preconditioner* preconditioner_pt=new TrilinosIFPACKPreconditioner;

   // Pass pointer to preconditioner to oomph-lib IterativeLinearSolver
   linear_solver_pt->preconditioner_pt()=preconditioner_pt;
  }

#ifdef OOMPH_HAS_MPI
MPI_Helpers::finalize();
#endif


} //end of main









