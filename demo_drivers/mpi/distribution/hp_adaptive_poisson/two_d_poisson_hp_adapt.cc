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
//Driver for a simple 2D poisson problem with adaptive mesh refinement

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
 double Alpha=5.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Exact gradient as a Vector
 void get_exact_gradient(const Vector<double>& x, Vector<double>& dudx)
 {
  //dudx[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
  double CoshTerm = cosh(1.0-Alpha*(TanPhi*x[0]-x[1]))
                   *cosh(1.0-Alpha*(TanPhi*x[0]-x[1]));
  dudx[0] = -Alpha*TanPhi/CoshTerm;
  dudx[1] = Alpha/CoshTerm;
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



//==============================start_of_mesh======================
/// p-refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
//=================================================================
template<class ELEMENT>
class SimpleRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 ///  Pass number of elements in the horizontal 
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



//====== start_of_problem_class=======================================
/// hp-refineable 2D Poisson problem on rectangular domain, discretised
/// with hp-refineable 2D QPoisson elements. The specific type of
/// element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableTwoDPoissonProblem : public Problem
{
public:
    Vector<Integral*> integration_scheme;

public:

 /// Constructor: Pass pointer to source function
 RefineableTwoDPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                          source_fct_pt);

 /// Destructor (empty)
 ~RefineableTwoDPoissonProblem()
  {
   delete mesh_pt()->spatial_error_estimator_pt();
   delete Problem::mesh_pt();
  }

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// Overloaded version of the Problem's access function to 
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
RefineableTwoDPoissonProblem<ELEMENT>::
      RefineableTwoDPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt 
                               source_fct_pt)
       :  Source_fct_pt(source_fct_pt)
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=8;

 // # of elements in y-direction
 unsigned n_y=8;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new SimpleRefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

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




//=================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void RefineableTwoDPoissonProblem<ELEMENT>::actions_before_newton_solve()
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
     
     // Pin the node
     nod_pt->pin(0);

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
void RefineableTwoDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Get current process rank
 int my_rank=this->communicator_pt()->my_rank();

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc



//===== start_of_main=====================================================
/// Driver code for hp-adaptive solution to the 2D Poisson problem
//========================================================================
int main(int argc, char** argv)
{

#ifdef OOMPH_HAS_MPI

 // Initialise MPI
 MPI_Helpers::init(argc,argv);

#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 //Set up the problem
 //------------------
 
 // Create the problem with 2D hp-refineable elements from the
 // PRefineableQPoissonElement family. Pass pointer to source function. 
 RefineableTwoDPoissonProblem<PRefineableQPoissonElement<2> > 
  problem(&TanhSolnForPoisson::get_source);
 
 //Are there command-line arguments?
 if(CommandLineArgs::Argc==1)
  {
 
#ifdef OOMPH_HAS_MPI
 
   // Provide storage for each element's partition number
   const unsigned n_element=problem.mesh_pt()->nelement();
   Vector<unsigned> out_element_partition(n_element);
 
   // Distribute the problem
   bool report_stats=true;
   out_element_partition=problem.distribute(report_stats);
 
   // Write partition to disk
   std::ofstream output_file;
   char filename[100];
   sprintf(filename,"out_hp_adaptive_poisson_partition.dat");
   output_file.open(filename);
   for (unsigned e=0;e<n_element;e++)
    {
     output_file << out_element_partition[e] << std::endl;
    }
 
   // Check halo schemes (optional)
   problem.check_halo_schemes(doc_info);

#endif
   
  }
 else // Validation run - read in partition from file
  {

#ifdef OOMPH_HAS_MPI

   DocInfo mesh_doc_info;
   mesh_doc_info.set_directory("RESLT_MESH");
   mesh_doc_info.number()=0;
   std::ifstream input_file;
   char filename[100];

   // Get the partition to be used from file
   const unsigned n_element=problem.mesh_pt()->nelement();
   Vector<unsigned> element_partition(n_element);
   sprintf(filename,"hp_adaptive_poisson_partition.dat");
   input_file.open(filename);
   std::string input_string;
   for (unsigned e=0;e<n_element;e++)
    {
     getline(input_file,input_string,'\n');
     element_partition[e]=atoi(input_string.c_str());
    }

   // Now perform the distribution
   bool report_stats=true;
   problem.distribute(element_partition,mesh_doc_info,report_stats);

#endif

  } // end of mesh distribution


 // Check if we're ready to go:
 //----------------------------
 cout << "\n\n\nProblem self-test ";
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
 
 // Choose a large value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=20.0;
 
 // Initially uniformly refine the mesh
 for (unsigned p=0; p<0; p++)
  {
   cout << "p-refining:" << endl;
   problem.p_refine_uniformly();
  }
 for (unsigned h=0; h<0; h++)
  {
   cout << "h-refining:" << endl;
   problem.refine_uniformly();
  }
 
 problem.newton_solve();

 //Output the solution
 doc_info.number()=1;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=2;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=3;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=4;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=5;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=6;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=7;
 problem.doc_solution(doc_info);
 
 problem.p_adapt();
 problem.newton_solve();
 doc_info.number()=8;
 problem.doc_solution(doc_info);
 
 problem.adapt();
 problem.newton_solve();
 doc_info.number()=9;
 problem.doc_solution(doc_info);

 
 // Count hanging nodes
 unsigned hang_no=0;
 for (unsigned n=0; n<problem.mesh_pt()->nnode(); n++)
  {
   if (problem.mesh_pt()->node_pt(n)->is_hanging()) hang_no++;
  }
 cout << hang_no << " hanging nodes in mesh." << endl;
 
 //Output the solution
 doc_info.number()=0;
 problem.doc_solution(doc_info);


// Finalise MPI
#ifdef OOMPH_HAS_MPI

 MPI_Helpers::finalize();

#endif
 
} //end of main

