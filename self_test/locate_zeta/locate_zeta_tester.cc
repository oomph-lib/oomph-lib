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
//Driver for a simple 2D poisson problem with adaptive mesh refinement

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"



using namespace std;

using namespace oomph;

#include "locate_zeta_tester.h"


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

 ///   Pass number of elements in the horizontal 
 /// and vertical directions, and the corresponding dimensions.
 /// Timestepper defaults to Static.
 SimpleRefineableRectangularQuadMesh(const unsigned &Nx,
                                     const unsigned &Ny, 
                                     const double &Lx, const double &Ly,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly,time_stepper_pt)
  {
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();

  } // end of constructor
 

 /// Destructor: Empty
 virtual ~SimpleRefineableRectangularQuadMesh() {}

}; // end of mesh



//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace GlobalParameters
{

 /// Function to distort mesh
 double distorted_y(const double& y_normalised)
 {
  return 0.5*y_normalised*(3.0-y_normalised);
 }

 /// Parameter for steepness of "step"
 double Alpha=5.0;

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

 /// Constructor: Pass pointer to source function
 RefineablePoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                          source_fct_pt);

 /// Destructor (empty)
 ~RefineablePoissonProblem()
  {
   delete mesh_pt()->spatial_error_estimator_pt();
   delete Problem::mesh_pt();

  }

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 ///  Overloaded version of the Problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 SimpleRefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<SimpleRefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
RefineablePoissonProblem<ELEMENT>::
      RefineablePoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                               source_fct_pt)
       :  Source_fct_pt(source_fct_pt)
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
 Problem::mesh_pt() = 
  new SimpleRefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Distort the mesh
 if (CommandLineArgs::command_line_flag_has_been_set("--distort_mesh"))
  {
   unsigned nnod=mesh_pt()->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     double y_orig=mesh_pt()->node_pt(j)->x(1);
     
     mesh_pt()->node_pt(j)->x(1)=
      l_y*GlobalParameters::distorted_y(y_orig/l_y);
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




//=================================start_of_actions_before_newton_solve===
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
     GlobalParameters::get_exact_u(x,u);

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
 mesh_pt()->output_fct(some_file,npts,GlobalParameters::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,GlobalParameters::get_exact_u,
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

#ifdef OOMPH_HAS_CGAL
 oomph_info << "have it" << std::endl;
#else
 oomph_info << "don't have it" << std::endl;
#endif


 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Validation of projection?
 CommandLineArgs::specify_command_line_flag("--distort_mesh");

 // Number of biased refinements for locate zeta test
 unsigned nrefine_for_locate_zeta_test=2; // hierher
 CommandLineArgs::specify_command_line_flag("--nrefine_for_locate_zeta_test",
                                            &nrefine_for_locate_zeta_test);

 // Number of genuine refinements for locate zeta test
 unsigned nrefine_genuine=2; // hierher
 CommandLineArgs::specify_command_line_flag("--nrefine_genuine",
                                            &nrefine_genuine);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node refineable elements from the
 // RefineableQuadPoissonElement family. Pass pointer to source function. 
 RefineablePoissonProblem<RefineableQPoissonElement<2,3> > 
  problem(&GlobalParameters::get_source);
 
 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;
 
 // Set the orientation of the "step" to 45 degrees
 GlobalParameters::TanPhi=1.0;
 
 // Choose a large value for the steepness of the "step"
 GlobalParameters::Alpha=50.0; 

 // Solve the problem, performing adaptive refinements
 problem.newton_solve(nrefine_genuine);
 
 //Output the solution
 problem.doc_solution(doc_info);
 doc_info.number()++;

 // Setup fake error estimator: Boundaries of refinement region
 Vector<double> lower_left(2);
 Vector<double> upper_right(2);

 double box_size=0.7;
 lower_left[0]=0.2;
 lower_left[1]=0.9;
 upper_right[0]=lower_left[0]+box_size;
 upper_right[1]=lower_left[1]+box_size;
 unsigned central_node_number=8;

 // Stop unrefinement and allow for infinite refinement
 problem.mesh_pt()->min_permitted_error()=-1.0;
 problem.mesh_pt()->max_refinement_level()=UINT_MAX;

 // Adapt lots of times
 for (unsigned i=0;i<nrefine_for_locate_zeta_test;i++)
  {
   // Clean up
   delete problem.mesh_pt()->spatial_error_estimator_pt();

   // Make a new error estimator
   problem.mesh_pt()->spatial_error_estimator_pt()=
    new DummyErrorEstimator(problem.mesh_pt(),lower_left,upper_right,
                            central_node_number);  
   
   // Adapt
   problem.adapt();

   // Shrink box
   box_size*=0.5;
   upper_right[0]=lower_left[0]+box_size;
   upper_right[1]=lower_left[1]+box_size;
  }

 //Output the solution
 problem.doc_solution(doc_info);

 // Check it out!
 check_locate_zeta(problem.mesh_pt());
   
 
} //end of main









