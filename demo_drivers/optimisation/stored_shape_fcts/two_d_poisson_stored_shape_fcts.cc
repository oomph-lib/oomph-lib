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
// Driver code to document the cost and benefits of using
// stored shape functions.

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

 ///  Constructor: Pass pointer to source function
 /// Mesh has 2^|n_power| x 2^|n_power| elements.
 PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                const unsigned& n_power);

 /// Destructor (empty)
 ~PoissonProblem(){}

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Treat the problem as being nonlinear
 void set_problem_is_nonlinear() {Problem::Problem_is_nonlinear = true;}

 /// Treat the problem as being linear
 void set_problem_is_linear() {Problem::Problem_is_nonlinear = false;}

 ///  Return the flag to determine whether the problem is being
 /// treated as linear or nonlinear
 bool is_problem_nonlinear() const {return Problem::Problem_is_nonlinear;}

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);


 ///  Run the probl
 void run_it(DocInfo& doc_info);

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
/// Mesh has 2^|n_power| x 2^|n_power| elements.
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
      PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                     const unsigned& n_power)
       :  Source_fct_pt(source_fct_pt)
{ 
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=unsigned(pow(2.0,int(n_power)));
 
 // # of elements in y-direction
 unsigned n_y=unsigned(pow(2.0,int(n_power)));

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;


 cout << " Building a " << n_x << " x " << n_y << " mesh." 
      << std::endl;

 //  Initialise timers
 clock_t t_start = clock();

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Finish/doc timing
 clock_t t_end = clock();
 double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
 cout << std::endl;
 cout << "======================================================= " << std::endl;
 cout << "Total time for Mesh setup [sec]: " << total_time << std::endl;
 cout << "======================================================= " << std::endl;
 cout << std::endl;

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

 // Do output
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 FILE* file_pt = fopen(filename,"w");
 mesh_pt()->output(file_pt,npts);
 fclose(file_pt);

 // Doc error and return of the square of the L2 error
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

 



//======start_of_run_it===============================================
/// Run the problem
//====================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::run_it(DocInfo& doc_info)
{ 

 // Repeated assembly of Jacobian -- stored shape functions only
 // pay of on subsequent solves.
 DoubleVector residuals;
 CRDoubleMatrix jacobian;
 
 unsigned n_assemble=10;
 for (unsigned i=0;i<n_assemble;i++)
  {
   
   //  Initialise timers
   clock_t t_start = clock();
   
   // Assemble Jacobian and residual vector
   get_jacobian(residuals,jacobian);
   
   // Finish/doc timing
   clock_t t_end = clock();
   double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
   cout << "    ======================================================= " 
        << std::endl;
   cout << "    Time for " << i << "th jac/residual computation [sec]: " 
        << total_time << std::endl;
   cout << "    ======================================================= " 
        << std::endl;
  }



 // Solve as linear problem only
 this->set_problem_is_linear();
 
 //  Initialise timers
 clock_t t_start = clock();
 
 // Solve the problem
 newton_solve();
 
 // Finish/doc timing
 clock_t t_end = clock();
 double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
 cout << "======================================================= " << std::endl;
 cout << "Total time for Newton solve [sec]: " << total_time << std::endl;
 cout << "======================================================= " << std::endl;
 
 //Output the solution
 doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Give user a chance to check the memory usage
 cout << std::endl << std::endl;
 cout << "Execution is paused while problem is in core" << std::endl;
 pause("Have a look at the memory usage now"); 
 

}


//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main(int argc, char* argv[])
{
 // Number of uniform refinements relative to a 2x2 base mesh
 unsigned n_refine;
 
 // Get number of refinement levels from command line or use default
 if (argc==1)
  {
   n_refine=6; //5
  }
 else if (argc==2)
  {
   n_refine=atoi(argv[1]);
  }
 else
  {
   std::string error_message =
    "Wrong number of input arguments. The options are: \n";
   error_message +=
    "No args: Default number of refinements\n";
   error_message +=
    "One arg: Required number of refinements\n";

   throw OomphLibError(error_message,
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }
 

  cout << sizeof(Shape) << std::endl;
  cout << sizeof(DShape) << std::endl;
 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0; 




 //Set up the problem with normal Poisson elements
 //-----------------------------------------------
 {
  //  Initialise timers
  clock_t t_start = clock();
  
  // Create the problem with 2D nine-node elements from the
  // QPoissonElement family. Pass pointer to source function
  // and number of refinements required
  PoissonProblem<QPoissonElement<2,3> > 
   problem(&TanhSolnForPoisson::source_function,n_refine);
  
  // Finish/doc timing
  clock_t t_end = clock();
  double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << std::endl;
  cout << "======================================================= " << std::endl;
  cout << "Total time for Problem setup [sec]: " << total_time << std::endl;
  cout << "======================================================= " << std::endl;
  cout << std::endl;
  
  // Run problem
  problem.run_it(doc_info);
 }



 //Set up the problem with Poisson elements with storable shape fcts
 //-----------------------------------------------------------------
 {
  //  Initialise timers
  clock_t t_start = clock();
  
  typedef StorableShapeElement<QPoissonElement<2,3> >  ELEMENT;
  // Create the problem with 2D nine-node elements from the
  // QPoissonElement family. Pass pointer to source function
  // and number of refinements required
  PoissonProblem<ELEMENT> 
   problem(&TanhSolnForPoisson::source_function,n_refine);
  
  // Finish/doc timing
  clock_t t_end = clock();
  double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << "Total time for Problem setup with storable shape fcts [sec]: " 
       << total_time << std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << std::endl;
  
  // Run problem with nothing stored (the default)
  //----------------------------------------------
  cout << "================================================================= "
       << std::endl;
  cout << "No shape functions actually stored                                "
       << std::endl;
  cout << "================================================================= "
       << std::endl << std::endl;
  problem.run_it(doc_info);



  // Store the derivatives w.r.t. to the Eulerian coordinates
  //-------------------------------------------------------------
  //  Initialise timers
  t_start = clock();
  unsigned nelem=problem.mesh_pt()->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    dynamic_cast<ELEMENT*>(problem.mesh_pt()->element_pt(e))
     ->pre_compute_dshape_eulerian_at_knots();
   }

  // Finish/doc timing
  t_end = clock();
  total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << "Total time for pre-calculation of shape functions and derivatives "
       << std::endl; 
  cout << total_time << " sec" <<  std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << std::endl;
  //Now run the problem
  problem.run_it(doc_info);
  
  // Set them all to refer to the first elements
  //-------------------------------------------------------------
  //  Initialise timers
  t_start = clock();
  for (unsigned e=1;e<nelem;e++)
   {
    dynamic_cast<ELEMENT*>(problem.mesh_pt()->element_pt(e))
     ->set_dshape_eulerian_stored_from_element(
      dynamic_cast<ELEMENT*>(problem.mesh_pt()->element_pt(0)));
   }
  // Finish/doc timing
  t_end = clock();
  total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << "Total time for assignment of shape functions and derivatives from "
       << std::endl << "the first element: " << total_time << " sec" <<  std::endl;
  cout << "================================================================= " 
       << std::endl;
  cout << std::endl;

  problem.run_it(doc_info);
 }


} //end of main









