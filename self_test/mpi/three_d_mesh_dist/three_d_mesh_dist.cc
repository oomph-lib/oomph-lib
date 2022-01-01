//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
//Driver for a 3D poisson problem

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

using namespace std;

using namespace oomph;


//=============start_of_namespace=====================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of step
 double Alpha=1;

 /// Orientation (non-normalised x-component of unit vector in direction
 /// of step plane)
 double N_x=-1.0;

 /// Orientation (non-normalised y-component of unit vector in direction
 /// of step plane)
 double N_y=-1.0;

 /// Orientation (non-normalised z-component of unit vector in direction
 /// of step plane)
 double N_z=1.0;


 /// Orientation (x-coordinate of step plane) 
 double X_0=0.0;

 /// Orientation (y-coordinate of step plane) 
 double Y_0=0.0;

 /// Orientation (z-coordinate of step plane) 
 double Z_0=0.0;


 // Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
                     N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
                     N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  u = tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
                     N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
                     N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }


 /// Source function to make it an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {

  double s1,s2,s3,s4;

  s1 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_x*N_x/(N_x*N_x+N_y*N_y+N_z*N_z);
      s3 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_y*N_y/(N_x*N_x+N_y*N_y+N_z*N_z);
      s4 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_z*N_z/(N_x*N_x+N_y*N_y+N_z*N_z);
      s2 = s3+s4;
      source = s1+s2;
 }


} // end of namespace




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_class_definition====================================
/// Poisson problem
//====================================================================
template<class ELEMENT>
class ThreeDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 ThreeDPoissonProblem(
  PoissonEquations<3>::PoissonSourceFctPt source_fct_pt);

 /// Destructor: Empty
 ~ThreeDPoissonProblem(){}

 /// Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 RefineableSimpleCubicMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableSimpleCubicMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()  {}

 /// Update the problem specs before solve: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_newton_solve()
 {
  //Loop over the boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
   {
    // Loop over the nodes on boundary
    unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     double u;
     Vector<double> x(3);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     x[2]=nod_pt->x(2);
     TanhSolnForPoisson::get_exact_u(x,u);
     nod_pt->set_value(0,u);
    }
   }
 }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to source function
 PoissonEquations<3>::PoissonSourceFctPt Source_fct_pt;


}; // end of class definition





//====================start_of_constructor================================
/// Constructor for Poisson problem 
//========================================================================
template<class ELEMENT>
ThreeDPoissonProblem<ELEMENT>::ThreeDPoissonProblem(
   PoissonEquations<3>::PoissonSourceFctPt source_fct_pt) : 
         Source_fct_pt(source_fct_pt)

{ 

 // Setup parameters for exact tanh solution
 // Steepness of step
 TanhSolnForPoisson::Alpha=5.0;

 /// 2x2x2 elements in a 5x5x5 box
 Problem::mesh_pt() = new RefineableSimpleCubicMesh<ELEMENT>(
  2,2,2,5.0,5.0,5.0);

 //Doc the mesh boundaries
 ofstream some_file;
 some_file.open("boundaries.dat");
 mesh_pt()->output_boundaries(some_file);
 some_file.close();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here (all the nodes on the boundary)
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  } // end of pinning


 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering 
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ThreeDPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 
 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();
 
 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 
//  // Output exact solution 
//  //----------------------
//  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
//  some_file.close();

 doc_info.number()++;

} // end of doc


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////




//===== start_of_main=====================================================
/// Test mesh distribution in parallel, starting with n_uniform_first 
/// uniform refinments before starting the distribution. Then
/// refine and unrefine uniformly as specified. Re-distribute
/// after every uniform refinement if requested by the boolean flag.
//========================================================================
void parallel_test(const unsigned& n_refine_first, 
                   const unsigned& n_refine_once_distributed, 
                   const bool& redistribute,
                   const bool& solve,
                   const bool& doc)
{

 oomph_info 
  << std::endl  << std::endl
  << "Running with " << n_refine_first << " initial serial refinements and \n"
  << "with " << n_refine_once_distributed 
  << " subsequent distributed refinements"
  << std::endl << std::endl;

 double av_efficiency; 
 double max_efficiency;
 double min_efficiency;
 double av_number_halo_nodes;
 unsigned max_number_halo_nodes;
 unsigned min_number_halo_nodes;
 
 char filename[200];
 std::ofstream some_file;

 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 DocInfo quiet_doc_info;
 quiet_doc_info.set_directory("RESLT");

 // Create the problem
 //-------------------
 ThreeDPoissonProblem<RefineableQPoissonElement<3,3> >* problem_pt=
  new ThreeDPoissonProblem<RefineableQPoissonElement<3,3> >
  (&TanhSolnForPoisson::get_source);

 // Storage for number of processors and current processor
 int n_proc=problem_pt->communicator_pt()->nproc();
 int my_rank=problem_pt->communicator_pt()->my_rank(); 

 // Initial number of dofs in base mesh
 unsigned n_dof_orig=problem_pt->ndof();
 unsigned n_dof=problem_pt->ndof();
 

 // Start with a few rounds of uniform refinement
 //----------------------------------------------
 clock_t t_start = clock();
 for (unsigned i=0;i<n_refine_first;i++)
  {
   // Refine uniformly
   problem_pt->refine_uniformly();
  } 
 clock_t t_end = clock();
 

 // Doc stats
 if (my_rank==0)
  {
   double t_tot=double(t_end-t_start)/CLOCKS_PER_SEC;
   oomph_info << "Time for initial adaptation [ndof=" << problem_pt->ndof() 
              << "]: " << t_tot << std::endl;
 
   n_dof=problem_pt->ndof();
   sprintf(
    filename,
    "%s/three_d_initial_refinement_nrefinement_%i_ndofbefore_%i_ndofafter_%i.dat",
    doc_info.directory().c_str(),
    n_refine_first,
    n_dof_orig,
    n_dof);
   some_file.open(filename);
   some_file << n_refine_first << " " 
             << n_dof_orig << " " 
             << n_dof << " " 
             << t_tot << " "
             << std::endl;
   some_file.close();    
 }


 
 // Now distribute the problem
 //---------------------------
 t_start = clock();
 if (doc)
  {
   std::ifstream input_file;
   std::ofstream output_file;

   // Get partition from file
   unsigned n_partition=problem_pt->mesh_pt()->nelement();
   Vector<unsigned> element_partition(n_partition,0);
   sprintf(filename,"three_d_mesh_dist_partition.dat");
   input_file.open(filename);
   std::string input_string;
   for (unsigned e=0;e<n_partition;e++)
    {
     getline(input_file,input_string,'\n');
     element_partition[e]=atoi(input_string.c_str());
    }

//   Vector<unsigned> out_element_partition;
   bool report_stats=false;
   problem_pt->distribute(element_partition,doc_info,report_stats);
//                           out_element_partition);

//    sprintf(filename,"out_three_d_mesh_dist_partition.dat");
//    output_file.open(filename);
//    for (unsigned e=0;e<n_partition;e++)
//     {
//      output_file << out_element_partition[e] << std::endl;
//     }

//   problem_pt->distribute (doc_info,report_stats);
  }
 else
  {
   problem_pt->distribute();
  }
 t_end = clock();
 

 // Assess the quality of the distribution
 //---------------------------------------
 problem_pt->mesh_pt()->get_halo_node_stats(av_number_halo_nodes,
                                            max_number_halo_nodes,
                                            min_number_halo_nodes);
 problem_pt->mesh_pt()->get_efficiency_of_mesh_distribution
  (av_efficiency,max_efficiency,min_efficiency);
 

 // Doc stats
 if (my_rank==0)
  {
   double t_tot=double(t_end-t_start)/CLOCKS_PER_SEC;
   oomph_info << "Time for problem distribution [ndof=" << problem_pt->ndof() 
              << "]: " << t_tot << std::endl;

   sprintf(
    filename,
    "%s/three_d_initial_distr_np_%i_nrefinement_%i_ndof_%i.dat",
    doc_info.directory().c_str(),
    n_proc,
    n_refine_first,
    n_dof);
   some_file.open(filename);
   some_file << n_proc << " " 
             << n_dof << " " 
             << t_tot << " "
             << av_number_halo_nodes << " " 
             << max_number_halo_nodes << " " 
             << min_number_halo_nodes  << " " 
             << av_efficiency << " " 
             << max_efficiency << " " 
             << min_efficiency << " " 
             << std::endl;
   some_file.close();    
  }


 // Check things
 //-------------
 problem_pt->check_halo_schemes(quiet_doc_info);


 // Solve?
 //-------
 if (solve)
  {

#ifdef OOMPH_HAS_TRILINOS

  // Create a Trilinos solver
  TrilinosAztecOOSolver* linear_solver_pt = new TrilinosAztecOOSolver;

  // Create the Trilinos ML preconditioner
  TrilinosMLPreconditioner* preconditioner_pt = new TrilinosMLPreconditioner;

  // Set the preconditioner pointer
  linear_solver_pt->preconditioner_pt() = preconditioner_pt;

  // Choose solver type (GMRES)
  linear_solver_pt->solver_type()=TrilinosAztecOOSolver::GMRES;

  // Set linear solver
  problem_pt->linear_solver_pt() = linear_solver_pt;
  linear_solver_pt->disable_doc_time();
    
#endif

   // Solve
   problem_pt->newton_solve();
   
   if (doc)
    {
     // Doc the solution
     problem_pt->doc_solution(doc_info);
    }
  }


   // distributed... so now check halo schemes, 
   // as there seems to be a problem after the uniform refinement...
   // andy trying to find bug
//   DocInfo tmp_doc_info;

   // Check files into DIST directory
//   tmp_doc_info.set_directory("DIST");
//   problem_pt->check_halo_schemes(tmp_doc_info); // andy checking for bug 
//   problem_pt->mesh_pt()->doc_mesh_distribution(tmp_doc_info);
 

 // Now refine the distributed problem uniformly
 //---------------------------------------------
 for (unsigned i=0;i<n_refine_once_distributed;i++)
  {
   // Number of dofs before uniform refinement
   n_dof_orig=problem_pt->ndof();

   // Refine uniformly
   //-----------------
   t_start = clock();   
   problem_pt->refine_uniformly();
   t_end = clock();
 

   // Assess efficiency
   //------------------
   problem_pt->mesh_pt()->get_halo_node_stats(av_number_halo_nodes,
                                              max_number_halo_nodes,
                                              min_number_halo_nodes);
   problem_pt->mesh_pt()->get_efficiency_of_mesh_distribution
    (av_efficiency,max_efficiency,min_efficiency);
 

   if (my_rank==0)
    {
     double t_tot=double(t_end-t_start)/CLOCKS_PER_SEC;
     oomph_info << "Time for distributed adaptation [ndof=" 
                << problem_pt->ndof() 
                << "]: " << t_tot << std::endl;

     // Number of dofs achieved after uniform refinement
     n_dof=problem_pt->ndof();

     // Snippet indicating if  redistribution of haloes is used or not
     string snippet="";
     if (!redistribute) snippet="out";

     sprintf(
      filename,
      "%s/three_d_distributed_refinement_np_%i_ninitialrefinement_%i_ntotalrefinement_%i_ndof_%i_with%sredistribution.dat",
      doc_info.directory().c_str(),
      n_proc,
      n_refine_first,
      n_refine_first+i+1,
      n_dof,
      snippet.c_str());
     some_file.open(filename);
     some_file << n_proc << " " 
               << n_dof_orig << " " 
               << n_dof << " " 
               << t_tot << " "
               << av_number_halo_nodes << " " 
               << max_number_halo_nodes << " " 
               << min_number_halo_nodes  << " " 
               << av_efficiency << " " 
               << max_efficiency << " " 
               << min_efficiency << " " 
               << std::endl;
     some_file.close();   

    }


 
   // Redistribute (if required)
   //---------------------------
   if (redistribute)
    {
      t_start = clock();   
      bool report_stats=false;
      doc_info.disable_doc();
      problem_pt->prune_halo_elements_and_nodes(doc_info, 
                                                report_stats);
      t_end = clock();
      double t_tot=double(t_end-t_start)/CLOCKS_PER_SEC;
     

      if (my_rank==0)
       {

        oomph_info << "Av., min. max. efficiency before redistribution " 
                   << av_efficiency << " " 
                   << min_efficiency << " " 
                   << max_efficiency << " " 
                   << std::endl;
        sprintf(
         filename,
         "%s/three_d_redistribution_np_%i_ninitialrefinement_%i_ntotalrefinement_%i_ndof_%i.dat",
         doc_info.directory().c_str(),
         n_proc,
         n_refine_first,
         n_refine_first+i+1,
         n_dof);
        some_file.open(filename);
        some_file << n_proc << " " 
                  << n_dof << " " 
                  << t_tot << " "
                  << av_number_halo_nodes << " " 
                  << max_number_halo_nodes << " " 
                  << min_number_halo_nodes  << " " 
                  << av_efficiency << " " 
                  << max_efficiency << " " 
                  << min_efficiency << " " 
                  << std::endl;
       }


      // Re-assess efficiency
      //----------------------
      problem_pt->mesh_pt()->get_halo_node_stats(av_number_halo_nodes,
                                                 max_number_halo_nodes,
                                                 min_number_halo_nodes);
      problem_pt->mesh_pt()->
       get_efficiency_of_mesh_distribution(av_efficiency,
                                           max_efficiency,
                                           min_efficiency);
      
      
      if (my_rank==0)
       {
        oomph_info << "Time for redistribution  [ndof=" 
                   << problem_pt->ndof() 
                   << "]: " << t_tot << std::endl;


        some_file << n_proc << " " 
                  << n_dof << " " 
                  << t_tot << " "
                  << av_number_halo_nodes << " " 
                  << max_number_halo_nodes << " " 
                  << min_number_halo_nodes  << " " 
                  << av_efficiency << " " 
                  << max_efficiency << " " 
                  << min_efficiency << " " 
                  << std::endl;
        some_file.close();   
       }
            
      problem_pt->mesh_pt()->
       get_efficiency_of_mesh_distribution(av_efficiency,
                                           max_efficiency,
                                           min_efficiency);
      oomph_info << "Av., min., max. efficiency after redistribution " 
                 << av_efficiency << " " 
                 << min_efficiency << " " 
                 << max_efficiency << " " 
                 << std::endl;
    }
   
     
   // Check things
   //-------------
   problem_pt->check_halo_schemes(quiet_doc_info);


   // Solve again
   //------------
   if (solve)
    {
     // Solve again
     problem_pt->newton_solve();
     
     if (doc)
      {
       // Doc the solution
       problem_pt->doc_solution(doc_info);
      }
    }

 }

 
 // Unrefine the distributed problem uniformly (only if not redistributed)
 //-----------------------------------------------------------------------
 if (!redistribute)
  {
   // Number of dofs before unrefinement
   n_dof_orig=problem_pt->ndof();
   
   
   t_start = clock();
   problem_pt->unrefine_uniformly();
   t_end = clock();
   
   if (my_rank==0)
    {
     double t_tot=double(t_end-t_start)/CLOCKS_PER_SEC;
     oomph_info << "Time for distributed unrefinement [ndof=" 
                << problem_pt->ndof() 
                << "]: " << t_tot << std::endl;
     
     // Number of dofs after unrefinement
     n_dof=problem_pt->ndof();
     
     sprintf(
      filename,
      "%s/three_d_distributed_unrefinement_np_%i_ninitialrefinement_%i_ntotalrefinement_%i_ndof_%i.dat",
      doc_info.directory().c_str(),
      n_proc,
      n_refine_first,
      n_refine_first+n_refine_once_distributed,
      n_dof);
     some_file.open(filename);
     some_file << n_proc << " " 
               << n_dof_orig << " " 
               << n_dof << " " 
               << t_tot << " "
               << av_number_halo_nodes << " " 
               << max_number_halo_nodes << " " 
               << min_number_halo_nodes  << " " 
               << av_efficiency << " " 
               << max_efficiency << " " 
               << min_efficiency << " " 
               << std::endl;
     some_file.close();   
     
    }
   
      
   // Check things
   problem_pt->check_halo_schemes(quiet_doc_info);
   
   // Solve again
   if (solve)
    {
     // Solve again
     problem_pt->newton_solve();
     
     
     if (doc)
      {
       // Doc the solution
       problem_pt->doc_solution(doc_info);
      }
    }
  }

}






//===== start_of_main=====================================================
/// Driver code for 3D Poisson problem
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 // Setup MPI helpers
 MPI_Helpers::init(argc,argv);
#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 

 // Default min. number of refinements before we distribute
 unsigned n_refine_first_min=1;

 // Default max. number of refinements before we distribute
 unsigned n_refine_first_max=1;

 // Default min. number of uniform refinements after distribution
 unsigned n_uniform_min=1; 

 // Default max. number of uniform refinements after distribution
 unsigned n_uniform_max=2; 

 // Redistribute after uniform refinement?
 bool redistribute=true;
                 

        
 if (CommandLineArgs::Argc==1)
  {
   std::cout << "Using default settings for flags" << std::endl;
  }
 else if (CommandLineArgs::Argc==6)
  {

    // Min. number of refinements before we distribute
    n_refine_first_min=atoi(argv[1]);

    // Max. number of refinements before we distribute
    n_refine_first_max=atoi(argv[2]);

    // Min. number of uniform refinements after distribution
    n_uniform_min=atoi(argv[3]); 

    // Max. number of uniform refinements after distribution
    n_uniform_max=atoi(argv[4]); 

    // Redistribute after uniform refinement?
    redistribute=atoi(argv[5]); 
  }
 else
  {
   oomph_info 
    << "Wrong number of input args. Specify none or five " 
    << "(n_refine_first_min/max, n_uniform_min/max, redistributon flag)" 
    << std::endl;
   exit(1);
  }



 // Tell us what we're going to do
 oomph_info << "Run with  n_refine_first_min/max, n_uniform_min/max, redistribute: " 
            << n_refine_first_min << " "
            << n_refine_first_max << " "
            << n_uniform_min   << " " 
            << n_uniform_max   << " " 
            << redistribute  << " " 
            << std::endl;
 

 // switch off output for all but one processor
// if (my_rank>0)
//  {
//   oomph_info.stream_pt() = &oomph_nullstream;
//  }


 // Solve and self test?
 bool solve=true; // (andy) only testing mesh distribution here
 bool doc=true;

 // Write headers for analysis
 ofstream some_file;
 some_file.open("RESLT/three_d_initial_refinement_header.dat");
 some_file << "VARIABLES= \"n<SUB>refine first</SUB>\"," 
           << "\"n<SUB>dof orig</SUB>\","
           << "\"n<SUB>dof</SUB>\","
           << "\"t<SUB>serial adapt</SUB>\""
           << std::endl;
 some_file.close();    


 some_file.open("RESLT/three_d_initial_distr_header.dat");
 some_file << "VARIABLES= \"n_<SUB>p</SUB>\","
           << "\"n<SUB>dof</SUB>\","
           << "\"t<SUB>distr</SUB>\","
           << "\"average n<SUB>halo nodes</SUB>\","
           << "\"max. n<SUB>halo nodes</SUB>\","
           << "\"min. n<SUB>halo nodes</SUB>\","
           << "\"average efficiency\","
           << "\"max. efficiency\","
           << "\"min. efficiency\""
           << std::endl;
 some_file.close();    


 some_file.open("RESLT/three_d_distributed_refinement_header.dat");
 some_file << "VARIABLES= \"n_<SUB>p</SUB>\","
           << "\"n<SUB>dof orig</SUB>\","
           << "\"n<SUB>dof</SUB>\","
           << "\"t<SUB>distributed refinement</SUB>\","
           << "\"average n<SUB>halo nodes</SUB>\","
           << "\"max. n<SUB>halo nodes</SUB>\","
           << "\"min. n<SUB>halo nodes</SUB>\","
           << "\"average efficiency\","
           << "\"max. efficiency\","
           << "\"min. efficiency\""
           << std::endl;
 some_file.close();   


 some_file.open("RESLT/three_d_redistribution_header.dat");
 some_file << "VARIABLES= \"n_<SUB>p</SUB>\","
           << "\"n<SUB>dof</SUB>\","
           << "\"t<SUB>redistribution</SUB>\","
           << "\"average n<SUB>halo nodes</SUB>\","
           << "\"max. n<SUB>halo nodes</SUB>\","
           << "\"min. n<SUB>halo nodes</SUB>\","
           << "\"average efficiency\","
           << "\"max. efficiency\","
           << "\"min. efficiency\""
           << std::endl;
 some_file.close();   



 some_file.open("RESLT/three_d_distributed_unrefinement_header.dat");
 some_file << "VARIABLES= \"n_<SUB>p</SUB>\","
           << "\"n<SUB>dof orig</SUB>\","
           << "\"n<SUB>dof</SUB>\","
           << "\"t<SUB>distributed unrefinement</SUB>\","
           << "\"average n<SUB>halo nodes</SUB>\","
           << "\"max. n<SUB>halo nodes</SUB>\","
           << "\"min. n<SUB>halo nodes</SUB>\","
           << "\"average efficiency\","
           << "\"max. efficiency\","
           << "\"min. efficiency\""
           << std::endl;
 some_file.close();   


 // Test the procedures in parallel mode
 for (unsigned i=n_refine_first_min;i<=n_refine_first_max;i++)

  {
   for (unsigned j=n_uniform_min;j<=n_uniform_max;j++)
    {
     parallel_test(i,j,redistribute,solve,doc);
    }
  }

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

}
