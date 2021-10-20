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
//Driver for a simple 1D poisson problem


// Generic oomph-lib routines
#include "generic.h"

// Include Poisson elements/equations
#include "poisson.h"

// Include the mesh
#include "meshes/one_d_mesh.h"


using namespace std;

using namespace oomph;


//==start_of_namespace================================================
/// Namespace for fish-shaped solution of 1D Poisson equation
//====================================================================
namespace FishSolnOneDPoisson
{

 ///  Sign of the source function 
 /// (- gives the upper half of the fish, + the lower half)
 int Sign=-1;


 /// Exact, fish-shaped solution as a 1D vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = double(Sign)*((sin(sqrt(30.0))-1.0)*x[0]-sin(sqrt(30.0)*x[0]));
 }


 /// Source function required to make the fish shape an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = double(Sign)*30.0*sin(sqrt(30.0)*x[0]);
 }

} // end of namespace







//==start_of_problem_class============================================
/// 1D Poisson problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 OneDPoissonProblem(const unsigned& n_element, 
                    PoissonEquations<1>::PoissonSourceFctPt source_fct_pt);

 /// Destructor - delete the mesh
 ~OneDPoissonProblem()
  {
   delete Mesh_pt;
  };

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// zero the dofs
 void zero_dofs()
  {
   unsigned ndof = Problem::ndof();
   for (unsigned i = 0; i < ndof; i++)
    {
     *Dof_pt[i] = 0.0;
    }
  }

 ///  Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(DocInfo& doc);

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;

 /// the mesh for this problem
 Mesh* Mesh_pt;
}; // end of problem class





//=====start_of_constructor===============================================
///  Constructor for 1D Poisson problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDPoissonProblem<ELEMENT>::OneDPoissonProblem(const unsigned& n_element,
 PoissonEquations<1>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 // Set domain length 
 double L=1.0;

 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = Mesh_pt = new OneDMesh<ELEMENT>(n_element,L);

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 mesh_pt()->boundary_node_pt(1,0)->pin(0);

 // Complete the setup of the 1D Poisson problem:

 // Loop over elements and set pointers to source function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer
   elem_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering scheme
 assign_eqn_numbers();
} // end of constructor




//===start_of_actions_before_newton_solve======================================
///  Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//=============================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 
 // Assign boundary values for this problem by reading them out
 // from the exact solution.

 // Left boundary is node 0 in the mesh:
 Node* left_node_pt=mesh_pt()->node_pt(0);

 // Determine the position of the boundary node (the exact solution
 // requires the coordinate in a 1D vector!)
 Vector<double> x(1);
 x[0]=left_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 Vector<double> u(1);
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 left_node_pt->set_value(0,u[0]);

 // Right boundary is last node in the mesh:
 unsigned last_node=mesh_pt()->nnode()-1;
 Node* right_node_pt=mesh_pt()->node_pt(last_node);

 // Determine the position of the boundary node
 x[0]=right_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 right_node_pt->set_value(0,u[0]);

 
} // end of actions before solve



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution with specified number of plot points per element
 sprintf(filename,"%s/soln_direct_solver_%i.dat",doc.directory().c_str(),doc.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 sprintf(filename,"%s/exact_soln%i.dat",doc.directory().c_str(),doc.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,20*npts,FishSolnOneDPoisson::get_exact_u); 
 some_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc.directory().c_str(),doc.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,FishSolnOneDPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end of doc

 

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


//======start_of_main==========================================================
/// Driver for 1D Poisson problem
//=============================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 // Setup mpi but don't make a copy of mpi_comm_world because
 // mumps wants to work with the real thing.
 bool make_copy_of_mpi_comm_world=false;
 MPI_Helpers::init(argc,argv,make_copy_of_mpi_comm_world);
#endif

 //Number of elements
 unsigned n_element=40;

 // global distribution
 // LinearAlgebraDistribution 
 //  distributed_dist(problem.communicator_pt(),ndof,true);

 // the doc info
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

#ifdef OOMPH_HAS_MPI
 //////////////////////////////////////////////////////////////////////////////
 // SuperLU_dist Test
 //////////////////////////////////////////////////////////////////////////////

 // MATRIX BASED SOLVE / CRDoubleMatrix (global)
 {
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl;
  oomph_info 
   << "TESTING: SuperLU_dist matrix based solve w/ global CRDoubleMatrix"
   << std::endl;
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl << std::endl;
  OneDPoissonProblem<QPoissonElement<1,4> >
    problem(n_element,FishSolnOneDPoisson::source_function);
  unsigned ndof = problem.ndof();
  LinearAlgebraDistribution global_dist(problem.communicator_pt(),ndof,false);
  CRDoubleMatrix A(&global_dist);
  DoubleVector b(&global_dist,0.0);
  DoubleVector x(&global_dist,0.0);
  SuperLUSolver solver;
  solver.set_solver_type(SuperLUSolver::Distributed);
  problem.get_jacobian(b,A);
  solver.enable_doc_time();
  solver.solve(&A,b,x);
  x.output("RESLT/SuperLU_dist_CRDoubleMatrix_global.dat");
 }

 // MATRIX BASED SOLVE / CRDoubleMatrix (distributed)
 {
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl;
  oomph_info 
   << "TESTING: SuperLU_dist matrix based solve w/ dist CRDoubleMatrix"
   << std::endl;
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl << std::endl;
  OneDPoissonProblem<QPoissonElement<1,4> >
    problem(n_element,FishSolnOneDPoisson::source_function);
  unsigned ndof = problem.ndof();
  LinearAlgebraDistribution distributed_dist(problem.communicator_pt(),
					     ndof,true);
  CRDoubleMatrix A(&distributed_dist);
  DoubleVector b(&distributed_dist,0.0);
  DoubleVector x(&distributed_dist,0.0);
  SuperLUSolver solver;
  solver.set_solver_type(SuperLUSolver::Distributed);
  problem.get_jacobian(b,A);
  solver.enable_doc_time();
  solver.solve(&A,b,x);
  x.output("RESLT/SuperLU_dist_CRDoubleMatrix_distributed.dat");
 }

 // PROBLEM BASED SOLVE (global)
 {
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl;
  oomph_info << "TESTING: SuperLU_dist global problem based solve"
             << std::endl;
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl << std::endl;
  OneDPoissonProblem<QPoissonElement<1,4> >
    problem(n_element,FishSolnOneDPoisson::source_function);
  DoubleVector x;
  SuperLUSolver solver;
  solver.set_solver_type(SuperLUSolver::Distributed);
  solver.use_global_solve_in_superlu_dist();
  problem.linear_solver_pt() = &solver;
  problem.newton_solve();
  doc_info.number()++;
  problem.doc_solution(doc_info);
  problem.zero_dofs();
 }

 // PROBLEM BASED SOLVE (distributed)
 {
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl;
  oomph_info << "TESTING: SuperLU_dist distributed solve"
             << std::endl;
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl << std::endl;
  OneDPoissonProblem<QPoissonElement<1,4> >
    problem(n_element,FishSolnOneDPoisson::source_function);
  DoubleVector x;
  SuperLUSolver solver;
  solver.set_solver_type(SuperLUSolver::Distributed);
  solver.use_distributed_solve_in_superlu_dist();
  problem.linear_solver_pt() = &solver;
  problem.newton_solve();
  doc_info.number()++;
  problem.doc_solution(doc_info);
  problem.zero_dofs();
 }
#endif


#ifdef OOMPH_HAS_MUMPS

 // PROBLEM BASED SOLVE (distributed, using mumps)
 {
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl;
  oomph_info << "TESTING: MUMPS global problem based solve"
             << std::endl;
  oomph_info 
   << "///////////////////////////////////////////////////////////////////////"
   << std::endl << std::endl;
  OneDPoissonProblem<QPoissonElement<1,4> >
    problem(n_element,FishSolnOneDPoisson::source_function);
  DoubleVector x;
  MumpsSolver solver;
  problem.linear_solver_pt() = &solver;
  problem.newton_solve();
  doc_info.number()++;
  problem.doc_solution(doc_info);
  problem.zero_dofs();
 }

#else

 ofstream some_file;
 char filename[100];

 // Output solution with specified number of plot points per element
 sprintf(filename,"RESLT/dummy_mumps.dat");
 some_file.open(filename);
 some_file << "dummy data for missing mumps\n";
 some_file.close();

#endif


#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
} // end of main









