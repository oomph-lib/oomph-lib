// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

#ifdef OOMPH_HAS_MPI
// MPI header
#include "mpi.h"
#endif

// Error-catching
//#include "fenv.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;
} // end_of_namespace



//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class RefineableDrivenCavityProblem : public Problem
{ // perhaps this should have an ifdef OOMPH_HAS_MPI so it works in serial too?

public:

 /// Constructor
 RefineableDrivenCavityProblem();

 /// Destructor: Empty -- all memory gets cleaned up in base destructor
 ~RefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along all other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
  } // end_of_actions_before_newton_solve


 /// After adaptation: Unpin pressure and pin redundant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   // The "first element" contains the point (0,0)
   // need it here for distributed problem
   unsigned nnod=mesh_pt()->nnode();
   for (unsigned j=0; j<nnod; j++)
    {
     if (mesh_pt()->node_pt(j)->x(0)==0 && 
         mesh_pt()->node_pt(j)->x(1)==0) // 2d problem only
      {
         fix_pressure(0,0,0.0);
      }
    }

  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 
private:

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
RefineableDrivenCavityProblem<ELEMENT>::RefineableDrivenCavityProblem()
{ 
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
  // Now set the first pressure dof in the first element to 0.0
  // The "first element" contains the point (0,0)
  // no real need to do it here because problem not yet distributed
  unsigned nnod=mesh_pt()->nnode();
  for (unsigned j=0; j<nnod; j++)
   {
    if (mesh_pt()->node_pt(j)->x(0)==0 && 
        mesh_pt()->node_pt(j)->x(1)==0) // 2d problem only
        {
         fix_pressure(0,0,0.0);
        }
   }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),MPI_Helpers::My_rank);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main(int argc, char **argv)
{
// feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // initialise MPI and setup MPI_Helpers
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // output
// char filename[100];
 std::ofstream eqn_file;
 std::ofstream internal_file;
 std::ofstream jacobian_file;

 // Set max. number of black-box adaptation
 unsigned max_adapt=3;

 // Solve problem with Taylor Hood elements
 //---------------------------------------
 {
  oomph_info << " " << std::endl
             << "------------------------------------------" << std::endl
             << "-- Solving MPI problem with TH elements --" << std::endl
             << "------------------------------------------" << std::endl;

  //Build problem
  RefineableDrivenCavityProblem<RefineableQTaylorHoodElement<2> > problem;

#ifdef OOMPH_HAS_MPI
  // Setup solver
  problem.linear_solver_pt() = new SuperLU_dist;
  static_cast<SuperLU_dist*>(problem.linear_solver_pt())->
   enable_distributed_solve();
#endif

  // Initial solve
  problem.newton_solve();

  // Initial uniform refinement
  problem.refine_uniformly();
  
  // Set output directory
  DocInfo doc_info;
  doc_info.set_directory("RESLT_TH");
  DocInfo mesh_doc_info;
  mesh_doc_info.set_directory("RESLT_TH_MESH");

  doc_info.number()=0;
  problem.doc_solution(doc_info);

  // Distribute
#ifdef OOMPH_HAS_MPI
  mesh_doc_info.number()=0;

  std::ifstream input_file;
  std::ofstream output_file;
  char filename[100];

  // Get the partition to be used from file
  unsigned n_partition=problem.mesh_pt()->nelement();
  Vector<unsigned> element_partition(n_partition);
  sprintf(filename,"adaptive_cavity_1_partition.dat");
  input_file.open(filename);
  std::string input_string;
  for (unsigned e=0;e<n_partition;e++)
   {
    getline(input_file,input_string,'\n');
    element_partition[e]=atoi(input_string.c_str());
   }

//  Vector<unsigned> out_element_partition;
  bool report_stats=false;
  problem.distribute(mesh_doc_info,report_stats,element_partition);
//                     out_element_partition);

//   sprintf(filename,"out_adaptive_cavity_1_partition.dat");
//   output_file.open(filename);
//   for (unsigned e=0;e<n_partition;e++)
//    {
//     output_file << out_element_partition[e] << std::endl;
//    }

  // Check halos
  problem.check_halo_schemes(mesh_doc_info);
#endif

  // Re-solve with adaptation
  problem.newton_solve(max_adapt);

  // change doc_info number
  doc_info.number()=1;

  //Output solution
  problem.doc_solution(doc_info);

  mesh_doc_info.number()=1;
  problem.mesh_pt()->doc_mesh_distribution(mesh_doc_info);

 } // end of Taylor Hood elements
 

 // Solve problem with Crouzeix Raviart elements
 //--------------------------------------------
 {
  oomph_info << " " << std::endl
             << "------------------------------------------" << std::endl
             << "-- Solving MPI problem with CR elements --" << std::endl
             << "------------------------------------------" << std::endl;

  // Build problem
  RefineableDrivenCavityProblem<RefineableQCrouzeixRaviartElement<2> > problem;

#ifdef OOMPH_HAS_MPI
  // Setup solver
  problem.linear_solver_pt() = new SuperLU_dist;
  static_cast<SuperLU_dist*>(problem.linear_solver_pt())->
   enable_distributed_solve();
#endif

  // Initial solve
  problem.newton_solve();

  // Initial refine
  problem.refine_uniformly();

  // Set output directory
  DocInfo doc_info;
  doc_info.set_directory("RESLT_CR");
  DocInfo mesh_doc_info;
  mesh_doc_info.set_directory("RESLT_CR_MESH");

  doc_info.number()=0;
  problem.doc_solution(doc_info);

  // Distribute
#ifdef OOMPH_HAS_MPI
  mesh_doc_info.number()=0;

  std::ifstream input_file;
  std::ofstream output_file;
  char filename[100];

  // Get the partition to be used from file
  unsigned n_partition=problem.mesh_pt()->nelement();
  Vector<unsigned> element_partition(n_partition);
  sprintf(filename,"adaptive_cavity_2_partition.dat");
  input_file.open(filename);
  std::string input_string;
  for (unsigned e=0;e<n_partition;e++)
   {
    getline(input_file,input_string,'\n');
    element_partition[e]=atoi(input_string.c_str());
   }

//  Vector<unsigned> out_element_partition;
  bool report_stats=false;
  problem.distribute(mesh_doc_info,report_stats,element_partition);
//                     out_element_partition);

//   sprintf(filename,"out_adaptive_cavity_2_partition.dat");
//   output_file.open(filename);
//   for (unsigned e=0;e<n_partition;e++)
//    {
//     output_file << out_element_partition[e] << std::endl;
//    }

  // Check halos
  problem.check_halo_schemes(mesh_doc_info);
#endif

  // Re-solve with adaptation
  problem.newton_solve(max_adapt);

  // change doc_info number
  doc_info.number()=1;

  //Output solution
  problem.doc_solution(doc_info);

  mesh_doc_info.number()=1;
  problem.mesh_pt()->doc_mesh_distribution(mesh_doc_info);

 } // end of Crouzeix Raviart elements

#ifdef OOMPH_HAS_MPI
// Finalise MPI
 MPI_Helpers::finalize();
#endif

} // end_of_main

