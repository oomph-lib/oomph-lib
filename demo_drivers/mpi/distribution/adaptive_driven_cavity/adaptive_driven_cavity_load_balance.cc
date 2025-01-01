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
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

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
{

public:

 /// Constructor
 RefineableDrivenCavityProblem();

 /// Destructor: Empty
 ~RefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve. 
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


 /// After adaptation: Unpin pressures and pin redundant pressure dofs
 /// and pressure at origin
 void actions_after_adapt()
  {
   pin_only_pressure_at_origin();
  } 

 /// Build the mesh
 void build_mesh();

 /// Doc the solution
 void doc_solution(const std::string& header);

private:

 /// Doc info object
 DocInfo Doc_info;

 /// Pin redundant pressures and pressure at origin
 void pin_only_pressure_at_origin()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   
   // Loop over all elements
   const unsigned n_element=mesh_pt()->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     // If the lower left node of this element is (0,0), then fix the 
     // pressure dof in this element to zero
     if (mesh_pt()->finite_element_pt(e)->node_pt(0)->x(0)==0.0 && 
         mesh_pt()->finite_element_pt(e)->node_pt(0)->x(1)==0.0) // 2d problem
      {
       // Fix the pressure in element e at pdof=0 to 0.0
       unsigned pdof=0;
       double pvalue=0.0;

       //Cast to proper element and fix pressure
       dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
        fix_pressure(pdof,pvalue);
      }
    }
  }

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
RefineableDrivenCavityProblem<ELEMENT>::RefineableDrivenCavityProblem()
{ 

 // Set output directory
 Doc_info.set_directory("RESLT_LOAD_BALANCE");

 // Build the mesh
 build_mesh();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor


//==start_of_build_mesh===================================================
/// Build the mesh for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::build_mesh()
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
 
 // Fine tune error targets to get "interesting" refinement pattern
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  max_permitted_error()=1.0e-5;
 
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  min_permitted_error()=1.0e-6;
 
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
 const unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin pressure at origin no matter which processor contains that node
 pin_only_pressure_at_origin();

}// end_of_build_mesh


//==start_of_doc_solution================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::doc_solution(
 const std::string& header)
{ 

 ofstream some_file;
 ofstream some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Get current process rank
 int my_rank=this->communicator_pt()->my_rank();
 
 // Output solution 
 sprintf(filename,"%s/soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file << "TEXT X=25,Y=93,F=HELV,HU=POINT,C=BLUE,H=26,T=\"" 
           << header << "\"" << std::endl; 
 some_file.close();


 // Output solution 
 sprintf(filename,"%s/soln_with_haloes%i_on_proc%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->enable_output_of_halo_elements();
 mesh_pt()->output(some_file,npts);
 some_file << "TEXT X=25,Y=93,F=HELV,HU=POINT,C=BLUE,H=26,T=\"" 
           << header << "\"" << std::endl; 
 mesh_pt()->disable_output_of_halo_elements();
 some_file.close();
 
// Output mesh distribution
 mesh_pt()->doc_mesh_distribution(Doc_info);

 // Output boundaries
 sprintf(filename,"%s/boundaries%i_on_proc%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number(),my_rank);
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file << "TEXT X=25,Y=93,F=HELV,HU=POINT,C=BLUE,H=26,T=\"" 
           << header << "\"" << std::endl; 
 some_file.close();


 {
  unsigned count_non_halo=0;
  unsigned nel=mesh_pt()->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    if (!(mesh_pt()->element_pt(e)->is_halo())) count_non_halo++;
   }
  unsigned total_count_non_halo=0;
  MPI_Allreduce(&count_non_halo,&total_count_non_halo,1,MPI_UNSIGNED,MPI_SUM,
                communicator_pt()->mpi_comm());
  
  oomph_info << "Number of non halo elements: " 
             << total_count_non_halo << " " 
             << count_non_halo << "\n"; 
 }
 
 {
  
  sprintf(filename,"%s/soln_non_halo_nodes%i_on_proc%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number(),my_rank);
  some_file.open(filename);
  
  sprintf(filename,"%s/ndof%i_on_proc%i.dat",
          Doc_info.directory().c_str(),
          Doc_info.number(),my_rank);
  some_file2.open(filename);
  unsigned count_non_halo=0;
  unsigned nnod=mesh_pt()->nnode();
  for (unsigned j=0;j<nnod;j++)
   {
    Node* nod_pt=mesh_pt()->node_pt(j);
    if (!(nod_pt->is_halo()))
     {
      count_non_halo++;
      some_file << nod_pt->x(0) << " " 
                << nod_pt->x(1) << "\n"; 
      unsigned nval=nod_pt->nvalue();
      unsigned cnt=0;
      for (unsigned i=0;i<nval;i++)
       {
        if (!(nod_pt->is_pinned(i))) cnt++;
       }
      some_file2 << nod_pt->x(0) << " " 
                 << nod_pt->x(1) << " " 
                 << cnt << "\n"; 
     }
   }
  some_file2.close();
  some_file.close();
  unsigned total_count_non_halo=0;
  MPI_Allreduce(&count_non_halo,&total_count_non_halo,1,MPI_UNSIGNED,MPI_SUM,
                communicator_pt()->mpi_comm());
  
  oomph_info << "Number of non halo nodes: " 
             << total_count_non_halo << " " 
             << count_non_halo << "\n"; 
 }
 
 
 Doc_info.number()++;
 
} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main(int argc, char **argv)
{

#ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
#endif

  // Switch off output modifier
 oomph_info.output_modifier_pt() = &default_output_modifier;

 // Define processor-labeled output file for all on-screen stuff
 std::ofstream output_stream;
 char filename[100];
 sprintf(filename,"OUTPUT.%i",MPI_Helpers::communicator_pt()->my_rank());
 output_stream.open(filename);
 oomph_info.stream_pt() = &output_stream;
 OomphLibWarning::set_stream_pt(&output_stream);
 OomphLibError::set_stream_pt(&output_stream);   

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Use manual distribution of elements rather than metis
 CommandLineArgs::specify_command_line_flag("--validate");
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Build problem
 RefineableDrivenCavityProblem<RefineableQTaylorHoodElement<2> > problem;

 // Tell us what you're doing
 bool report_stats=true;

 // Metis or manual distribution?
 if (!CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   // Distribute with metis
   problem.distribute(report_stats);
  }
 else
  {    
   // Make up some pre-determined distribution
   unsigned n_element=problem.mesh_pt()->nelement();
   Vector<unsigned> element_partition(n_element);
   std::string input_string;
   for (unsigned e=0;e<n_element;e++)
    {
     unsigned target_part=0;
     if ( (e<unsigned(0.25*double(n_element))) ||
          ( (e>unsigned(0.5*double(n_element))) &&
            (e<unsigned(0.75*double(n_element))) ) )
      {
       target_part=1;
      }
     element_partition[e]=target_part;
    }
   problem.distribute(element_partition,report_stats);
   
   // This is a validation run: use the default partition when load balancing
   problem.set_default_partition_in_load_balance();    
  }
 
 //Output initial mesh
 problem.doc_solution("Initial mesh");
 
 // Solve the problem with automatic adaptation
 problem.newton_solve(1);
 
 //Output solution
 problem.doc_solution("First adapted solution");
 
 // Prune
 problem.prune_halo_elements_and_nodes(report_stats);
 
 //Output solution
 problem.doc_solution("After first prune");
 
 // Solve the problem again with automatic adaptation 
 problem.newton_solve(1);
 
 //Output solution
 problem.doc_solution("Second adapted solution after pruning");
 
 // Do load balancing
 problem.load_balance(report_stats);
  
 //Output solution
 problem.doc_solution("After load balancing");

 // Prune
 problem.prune_halo_elements_and_nodes(report_stats);
 
 //Output solution
 problem.doc_solution("After second pruning");
 
 // Solve the problem again with automatic adaptation
 problem.newton_solve(1);
 
 //Output solution
 problem.doc_solution("Third adapted solution after pruning");
 
 // Do a uniform refinement
 problem.refine_uniformly();
 
 //Output solution
 problem.doc_solution("After uniform refinement");
 
 // Prune
 problem.prune_halo_elements_and_nodes(report_stats);
 
 //Output solution
 problem.doc_solution("After third pruning");
 
 // Solve the problem again with automatic adaptation
 problem.newton_solve(1);
 
 //Output solution
 problem.doc_solution("Fourth adapted solution after pruning");
 
 // Do load balancing
 problem.load_balance(report_stats);
 
 //Output solution
 problem.doc_solution("After load balancing");
  
 // Solve the problem again with automatic adaptation
 problem.newton_solve(1);
 
 //Output solution
 problem.doc_solution("Fifth adapted solution after load balancing");
 
 oomph_info << "done\n";
 std::cout << "done\n";
 
  // Finalise MPI
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif
  
  
} // end_of_main

