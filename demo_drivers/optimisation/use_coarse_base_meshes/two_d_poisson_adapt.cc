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
// Driver code to document the relative cost of 
// (a) building a coarse base mesh and refining it afterwards to 
// (b) building a fine mesh directly.
// The latter is much more expensive so the take-home message is:
// always build a your initial mesh as coarsely as possible.

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

// We need to include the refineable_quad_mesh's 
// templated source here to allow the build of
// all templates. 
//#include "generic/refineable_quad_mesh.template.cc"


using namespace std;

using namespace oomph;

//==============================start_of_mesh======================
/// Refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
//=================================================================
template<class ELEMENT>
class SimpleRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 /// \short  Pass number of elements in the horizontal 
 /// and vertical directions, and the corresponding dimensions.
 /// Timestepper defaults to Static.
 SimpleRefineableRectangularQuadMesh(const unsigned &Nx,
                                     const unsigned &Ny, 
                                     const double &Lx, const double &Ly,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly,time_stepper_pt)
  {


   //  Initialise timers
   clock_t t_start = clock();
   
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();


   // Finish/doc timing
   clock_t t_end = clock();
   double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
   cout << std::endl;
   cout << "======================================================= " << std::endl;
   cout << "Total time for black-box setup of quadtree " << std::endl;
   cout << "forest                              [sec]: " << total_time << std::endl;
   cout << "======================================================= " << std::endl;
   cout << std::endl;

  } // end of constructor
 

 /// Destructor: Empty
 virtual ~SimpleRefineableRectangularQuadMesh(){}

}; // end of mesh



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
/// refineable 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineablePoissonProblem : public Problem
{

public:

 /// \short Constructor: Pass pointer to source function and number of 
 /// refinements for the base mesh. Mesh has 2^|n_power| x 2^|n_power| 
 /// elements. If n_power>0 this is achieved by building a base mesh
 /// with that number of elements; if  n_power<0, we apply the
 /// required number of uniform refinements to a 2x2 base mesh.
 RefineablePoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
                        const int& n_power);

 /// Destructor (empty)
 ~RefineablePoissonProblem(){};

 /// \short Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Treat the problem as being nonlinear
 void set_problem_is_nonlinear() {Problem::Problem_is_nonlinear = true;}

 /// Treat the problem as being linear
 void set_problem_is_linear() {Problem::Problem_is_nonlinear = false;}

 /// \short Return the flag to determine whether the problem is being
 /// treated as linear or nonlinear
 bool is_problem_nonlinear() const {return Problem::Problem_is_nonlinear;}

 /// \short Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// \short Overloaded version of the Problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 SimpleRefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<SimpleRefineableRectangularQuadMesh<ELEMENT>*>
    (Problem::mesh_pt());
  }

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function
/// and power for number of elements in base mesh.
//========================================================================
template<class ELEMENT>
RefineablePoissonProblem<ELEMENT>::
      RefineablePoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt, const int& n_power)
       :  Source_fct_pt(source_fct_pt)
{ 

 // Setup mesh
 unsigned n_x,n_y;

 if (n_power>0)
  {
   // # of elements in x-direction
   n_x=unsigned(pow(2.0,int(n_power)));
   
   // # of elements in y-direction
   n_y=unsigned(pow(2.0,int(n_power)));
  }
 else
  {
   // # of elements in x-direction
   n_x=2;
   
   // # of elements in y-direction
   n_y=2;
  }

 cout << " Building a " << n_x << " x " << n_y << " refineable base mesh." 
      << std::endl;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;


 //  Initialise timers
 clock_t t_start = clock();

 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);


 // Do uniform refinement if required
 if (n_power<0)
  {
   for (int r=0;r<std::abs(n_power);r++)
    {
     cout << "Doing one round of uniform mesh refinement" << std::endl;
     mesh_pt()->refine_uniformly();
    }
  }

 // Finish/doc timing
 clock_t t_end = clock();
 double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
 cout << std::endl;
 cout << "======================================================= " << std::endl;
 cout << "Total time for Mesh setup [sec]: " << total_time << std::endl;
 cout << "======================================================= " << std::endl;
 cout << std::endl;


 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
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
void RefineablePoissonProblem<ELEMENT>::actions_before_newton_solve()
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
void RefineablePoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution (C++ style) 
 //----------------------------
 {
  // Start timer
  clock_t t_start = clock();
  
  // Do output
  sprintf(filename,"%s/soln%i.cpp_style.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();

  // Finish/doc timing
  clock_t t_end = clock();
  double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << "------------------------------------------------------- " << std::endl;
  cout << "Total time for C++ style output [sec]: " << total_time << std::endl;
  cout << "------------------------------------------------------- " << std::endl;
 }


 // Output solution (C style) 
 //--------------------------
 {
  // Start timer
  clock_t t_start = clock();
  
  // Do output
  sprintf(filename,"%s/soln%i.c_style.dat",doc_info.directory().c_str(),
          doc_info.number());
  FILE* file_pt = fopen(filename,"w");
  mesh_pt()->output(file_pt,npts);
  fclose(file_pt);

  // Finish/doc timing
  clock_t t_end = clock();
  double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << "------------------------------------------------------- " << std::endl;
  cout << "Total time for C style output [sec]:   " << total_time << std::endl;
  cout << "------------------------------------------------------- " << std::endl;
 }



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
   n_refine=5;
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
 
 
 // Loop over two versions of refinement: Build a fine base mesh
 // or refine coarse base mesh uniformly. The former version is
 // much more expensive because the black-box quadtree forest generator
 // creates the quadtree forest in the base mesh.
 for (unsigned i=0;i<2;i++)
  {

   // Exponent for number of elements in base mesh
   int n_power;
   if (i==0)
    {
     cout << std::endl;
     cout << "/////////////////////////////////////////////////////////" 
          << std::endl;
     cout << "Building problem with a coarse base mesh" << std::endl;
     cout << "followed by uniform refinements." << std::endl;
     cout << "/////////////////////////////////////////////////////////" 
          << std::endl;
     cout << std::endl;
     n_power=-n_refine;
    }
   else
    {
     cout << std::endl;
     cout << "/////////////////////////////////////////////////////////" 
          << std::endl;
     cout << "Building problem with a fine base mesh" << std::endl;
     cout << "/////////////////////////////////////////////////////////" 
          << std::endl;
     cout << std::endl;
     n_power=n_refine+1;
    }
   
   //Set up the problem
   //------------------
   
   //  Initialise timers
   clock_t t_start = clock();
   
   // Create the problem with 2D nine-node refineable elements from the
   // RefineableQuadPoissonElement family. Pass pointer to source function. 
   RefineablePoissonProblem<RefineableQPoissonElement<2,3> > 
    problem(&TanhSolnForPoisson::get_source,n_power);
   
   // Finish/doc timing
   clock_t t_end = clock();
   double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
   cout << std::endl;
   cout << "======================================================= " << std::endl;
   cout << "Total time for Problem setup [sec]: " << total_time << std::endl;
   cout << "======================================================= " << std::endl;
   cout << std::endl;
   
   
   // Create label for output
   //------------------------
   DocInfo doc_info;
   
   // Set output directory
   doc_info.set_directory("RESLT");
   
   // Step number
   doc_info.number()=0;
      
   // Set the orientation of the "step" to 45 degrees
   TanhSolnForPoisson::TanPhi=1.0;
   
   // Choose a large value for the steepness of the "step"
   TanhSolnForPoisson::Alpha=1.0; 
   
   
   // Solve as linear and nonlinear problem
   //--------------------------------------
   problem.set_problem_is_nonlinear();
   
   unsigned nstep=2;
   for (unsigned istep=0;istep<nstep;istep++)
    {
     
     //  Initialise timers
     clock_t t_start = clock();
     
     if (problem.is_problem_nonlinear())
      {
       cout << std::endl << std::endl;
       cout << "============================ " << std::endl;
       cout << "Solving as nonlinear problem " << std::endl;
       cout << "============================ " << std::endl;
       cout << std::endl << std::endl;
      }
     else
      {
       cout << std::endl << std::endl;
       cout << "============================ " << std::endl;
       cout << "Solving as linear problem " << std::endl;
       cout << "============================ " << std::endl;
       cout << std::endl << std::endl;
      }
     
     // Solve the problem
     problem.newton_solve();
     
     // Finish/doc timing
     clock_t t_end = clock();
     double total_time=double(t_end-t_start)/CLOCKS_PER_SEC;
     cout << "======================================================= " 
          << std::endl;
     cout << "Total time for Newton solve [sec]: " << total_time << std::endl;
     cout << "======================================================= " 
          << std::endl;
     
     //Output the solution
     problem.doc_solution(doc_info);
     
     //Increment counter for solutions 
     doc_info.number()++;
     
     // Next time around solve as linear problem
     problem.set_problem_is_linear();
     
    }

   // Give user a chance to check the memory usage
   cout << std::endl << std::endl;
   cout << "Execution is paused while problem is in core" << std::endl;
   pause("Have a look at the memory usage now"); 
   
  }
 
} //end of main









