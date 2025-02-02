//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Driver for a simple 2D poisson problem -- compare various linear solvers

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=1.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
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

 /// Constructor: Pass pointer to source function and flag to indicate
 /// if order of elements should be messed up. nel_1d is the number of 
 /// elements along the 1d mesh edges -- total number of elements is 
 /// nel_1d x nel_1d. 
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                const unsigned& nel_1d,
                const bool& mess_up_order);

 /// Destructor (empty)
 ~PoissonProblem(){}

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function
/// and a flag that specifies if the order of the elements should
/// be messed up. nel_1d is the number of elements along
/// the 1d mesh edges -- total number of elements is nel_1d x nel_1d.
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
       const unsigned& nel_1d, const bool& mess_up_order)
       :  Source_fct_pt(source_fct_pt)
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=nel_1d;

 // # of elements in y-direction
 unsigned n_y=nel_1d;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);



 // Mess up order of elements in the Problem's mesh:
 //-------------------------------------------------
 if (mess_up_order)
  {
   unsigned n_element=mesh_pt()->nelement();
   
   // Make copy
   Vector<ELEMENT*> tmp_element_pt(n_element);
   for (unsigned e=0;e<n_element;e++)
    {
     tmp_element_pt[e]=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
    }
   
   // Reorder
   unsigned e_half=unsigned(0.5*double(n_element));
   unsigned e_lo=e_half-3;
   unsigned e_up=n_element-e_lo;
   unsigned count=0;
   for (unsigned e=0;e<e_lo;e++)
    {
     mesh_pt()->element_pt(count)=tmp_element_pt[e];
     count++;
     mesh_pt()->element_pt(count)=tmp_element_pt[n_element-e-1];
     count++;
    }
   for (unsigned e=e_lo;e<e_up;e++)
    {
     mesh_pt()->element_pt(count)=tmp_element_pt[e];
     count++;
    }
  }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
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
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned num_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // How many nodes are there on this boundary?
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);

   // Loop over the nodes on boundary
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

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
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc

 




//===== start_of_run======================================================
/// Build and run problem with specified linear solver. Also pass
/// flag to specify if order of elements should be messed up.
/// nel_1d is the number of elements along the 1D mesh edge. 
/// Total number of elements is nel_1d x nel_1d.
//========================================================================
void run(const string& dir_name,
         LinearSolver* linear_solver_pt,
         const unsigned nel_1d,
         bool mess_up_order)
{

 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node elements from the
 // QPoissonElement family. Pass pointer to source function
 // and flag to specify if order of elements should be messed up.
 PoissonProblem<QPoissonElement<2,3> > 
  problem(&TanhSolnForPoisson::get_source,nel_1d,mess_up_order);
 

 /// Set linear solver: 
 //--------------------
 problem.linear_solver_pt() = linear_solver_pt;


 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory(dir_name);

 // Step number
 doc_info.number()=0;

 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test:";
 if (problem.self_test()==0) 
  {
   cout << "passed: Problem can be solved." << std::endl;
  }
 else 
  {
   throw OomphLibError("Self test failed",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0; 

 // Solve the problem
 problem.newton_solve();
 
 //Output the solution
 problem.doc_solution(doc_info);
 

} //end of run





//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem -- compare different solvers
//========================================================================
int main()
{

 // Pointer to linear solver
 LinearSolver* linear_solver_pt;

 // Result directory
 string dir_name;

 // Cpu start/end times
 clock_t cpu_start,cpu_end;

 // Storage for cpu times with and without messed up order
 Vector<double> superlu_cr_cpu(2);
 Vector<double> superlu_cc_cpu(2);
 Vector<double> hsl_ma42_cpu(2);
 Vector<double> hsl_ma42_reordered_cpu(2);
 Vector<double> dense_lu_cpu(2);
 Vector<double> fd_lu_cpu(2);

 // Flag to indicate if order should be messed up:
 bool mess_up_order;
 
 // Number of elements along 1D edge of mesh; total number of elements
 // is nel_1d x nel_1d
 unsigned nel_1d;

 // Do run with and without messed up elements
 for (unsigned mess=0;mess<2;mess++)
  {

   // Mess up order?
   if (mess==0)
    {
     mess_up_order=true;
    }
   else
    {
     mess_up_order=false;
    }

   /// Use SuperLU with compressed row storage
   //-----------------------------------------
   
   cout << std::endl;
   cout << " Use SuperLU with compressed row storage: " << std::endl;
   cout << "========================================= " << std::endl;
   cout << std::endl;

   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new SuperLUSolver;
   
   /// Use compressed row storage
   static_cast<SuperLUSolver*>(linear_solver_pt)
    ->use_compressed_row_for_superlu_serial();
   
   /// Switch on full doc
   static_cast<SuperLUSolver*>(linear_solver_pt)->enable_doc_stats();
   
   // Choose result directory
   dir_name="RESLT_cr";

   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=20;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
   
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   superlu_cr_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   
   
   
   
   /// Use SuperLU with compressed column storage
   //--------------------------------------------


   cout << std::endl;
   cout << " Use SuperLU with compressed column storage: " << std::endl;
   cout << "============================================ " << std::endl;
   cout << std::endl;

   
   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new SuperLUSolver;
   
   /// Use compressed row storage
   static_cast<SuperLUSolver*>(linear_solver_pt)
    ->use_compressed_column_for_superlu_serial();
   
   /// Switch on full doc
   static_cast<SuperLUSolver*>(linear_solver_pt)->enable_doc_stats();

   // Choose result directory
   dir_name="RESLT_cc";

   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=20;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
   
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   superlu_cc_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   
   
#ifdef HAVE_HSL_SOURCES
   
   /// Use HSL frontal solver MA42 without element re-ordering
   //---------------------------------------------------------
   
   cout << std::endl;
   cout << " Use HSL frontal solver MA42 without element re-ordering: " << std::endl;
   cout << "========================================================= " << std::endl;
   cout << std::endl;

   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new HSL_MA42;
   
   /// Switch on full doc
   static_cast<HSL_MA42*>(linear_solver_pt)->enable_doc_stats();
   
   // Choose result directory
   dir_name="RESLT_frontal";


   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=20;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
   
   
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   hsl_ma42_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   
   
   
   /// Use HSL frontal solver MA42 with element reordering by Sloan's algorithm
   //--------------------------------------------------------------------------
   
   cout << std::endl;
   cout << " Use HSL frontal solver MA42 with element re-ordering: " << std::endl;
   cout << "====================================================== " << std::endl;
   cout << std::endl;

   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new HSL_MA42;
   
   /// Switch on full doc
   static_cast<HSL_MA42*>(linear_solver_pt)->enable_doc_stats();
   
   
   /// Switch on re-ordering
   static_cast<HSL_MA42*>(linear_solver_pt)->enable_reordering();
   
   // Choose result directory
   dir_name="RESLT_frontal_reordered";


   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=20;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
   
   
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   hsl_ma42_reordered_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   
 
#endif

   /// Use dense matrix LU decomposition
   //-----------------------------------
   
   cout << std::endl;
   cout << " Use dense matrix LU decomposition: " << std::endl;
   cout << "=================================== " << std::endl;
   cout << std::endl;

   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new DenseLU;
   
   // Choose result directory
   dir_name="RESLT_dense_LU";

   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=4;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
   
   
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   dense_lu_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   



   /// Use dense matrix LU decomposition
   //-----------------------------------
   
   cout << std::endl;
   cout << " Use dense FD-ed matrix LU decomposition: " << std::endl;
   cout << "========================================= " << std::endl;
   cout << std::endl;

   // Start cpu clock
   cpu_start=clock();
   
   // Build linear solver
   linear_solver_pt = new FD_LU;
   
   // Choose result directory
   dir_name="RESLT_FD_LU";

   // Number of elements along 1D edge of mesh; total number of elements
   // is nel_1d x nel_1d
   nel_1d=4;

   // Run it
   run(dir_name,linear_solver_pt,nel_1d,mess_up_order);
   
   // Note: solver does not have to be deleted -- it's killed automatically
   // in the problem destructor.
      
   // End cpu clock
   cpu_end=clock();
   
   // Timing
   fd_lu_cpu[mess]=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   
  }


 // Doc timings with and without messed up elements
 for (unsigned mess=0;mess<2;mess++)
  {
   cout << std::endl   << std::endl   << std::endl ;
   if (mess==0)
    {
     cout << "TIMINGS WITH MESSED UP ORDERING OF ELEMENTS: " << std::endl;
     cout << "============================================ " << std::endl;

    }
   else
    {
     cout << "TIMINGS WITHOUT MESSED UP ORDERING OF ELEMENTS: " << std::endl;
     cout << "=============================================== " << std::endl;
    }

   cout << "CPU time with SuperLU compressed row                 : " 
        <<  superlu_cr_cpu[mess] << std::endl;
   cout << "CPU time with SuperLU compressed col                 : " 
        <<  superlu_cc_cpu[mess] << std::endl;
#ifdef HAVE_HSL_SOURCES
   cout << "CPU time with MA42 frontal solver                    : " 
        <<  hsl_ma42_cpu[mess] << std::endl;
   cout << "CPU time with MA42 frontal solver (incl. reordering) : " 
        <<  hsl_ma42_reordered_cpu[mess] << std::endl;
#endif
   cout << "CPU time with dense LU solver (fewer elements!)      : " 
        <<  dense_lu_cpu[mess] << std::endl;  
   cout << "CPU time with dense LU solver & FD (fewer elements!) : " 
        <<  fd_lu_cpu[mess] << std::endl;  
   cout << std::endl   << std::endl   << std::endl ;
  }

 


} //end of main





