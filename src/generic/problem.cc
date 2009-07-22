//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include<list>
#include<algorithm>

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include "oomph_utilities.h"
#include "problem.h"
#include "timesteppers.h"
#include "explicit_timesteppers.h"
#include "refineable_mesh.h"
#include "linear_solver.h"
#include "eigen_solver.h"
#include "assembly_handler.h"
#include "dg_elements.h"
#include "partitioning.h"

namespace oomph
{


//////////////////////////////////////////////////////////////////
//Non-inline functions for the problem class
//////////////////////////////////////////////////////////////////

//================================================================
/// Constructor: Allocate space for one time stepper
/// and set all pointers to NULL and set defaults for all
/// parameters.
//===============================================================
 Problem::Problem() : 
  Mesh_pt(0), Time_pt(0), Explicit_time_stepper_pt(0), Saved_dof_pt(0), 
#ifdef OOMPH_HAS_MPI
  Must_recompute_load_balance_for_assembly(true),
#endif
  Newton_solver_tolerance(1.0e-8), 
  Max_newton_iterations(10), Max_residuals(10.0),
  Jacobian_reuse_is_enabled(false), Jacobian_has_been_computed(false),
  Problem_is_nonlinear(true),
  Pause_at_end_of_sparse_assembly(false),
  Sparse_assembly_method(Perform_assembly_using_vectors_of_pairs), 
  Sparse_assemble_with_arrays_initial_allocation(400),
  Sparse_assemble_with_arrays_allocation_increment(150),
  Numerical_zero_for_sparse_assembly(0.0),
  Mass_matrix_reuse_is_enabled(false), Mass_matrix_has_been_computed(false),
  Discontinuous_element_formulation(false),
  Minimum_dt(1.0e-12), Maximum_dt(1.0e12), DTSF_max_increase(4.0),
  Scale_arc_length(true), Desired_proportion_of_arc_length(0.5),
  Theta_squared(1.0), Sign_of_jacobian(0), Continuation_direction(1.0), 
  Parameter_derivative(1.0), Parameter_current(0.0),
  Ds_current(0.0), Desired_newton_iterations_ds(5), 
  Minimum_ds(1.0e-10), Bifurcation_detection(false), 
  First_jacobian_sign_change(false),Arc_length_step_taken(false),
#ifdef OOMPH_HAS_MPI
  Dist_problem_matrix_distribution(Uniform_matrix_distribution),
  Parallel_sparse_assemble_previous_allocation(0),
  Problem_has_been_distributed(false),
  Max_permitted_error_for_halo_check(1.0e-14),
#endif
  Shut_up_in_newton_solve(false)
 {
  // By default no submeshes:
  Sub_mesh_pt.resize(0);
  // No timesteppers
  Time_stepper_pt.resize(0);
  
  //Set the linear solvers, eigensolver and assembly handler
  Linear_solver_pt = Default_linear_solver_pt = new SuperLUSolver;
  
  Eigen_solver_pt = Default_eigen_solver_pt = new ARPACK;
  
  Assembly_handler_pt = Default_assembly_handler_pt = new AssemblyHandler;

  // setup the communicator
#ifdef OOMPH_HAS_MPI
  if (MPI_Helpers::MPI_has_been_initialised)
   {
    Communicator_pt = new OomphCommunicator(MPI_Helpers::Communicator_pt);
   }
  else
   {
    Communicator_pt = new OomphCommunicator();
   }
#else
  Communicator_pt = new OomphCommunicator();
#endif

  // just create an empty linear algebra distribution for the 
  // DOFs
  // this is setup when assign_eqn_numbers(...) is called.
  Dof_distribution_pt = new LinearAlgebraDistribution;
 }

//================================================================
/// Destructor to clean up memory
//================================================================
 Problem::~Problem()
 {

  // Delete the memory assigned for the global time
  // (it's created on the fly in Problem::add_time_stepper_pt()
  // so we are entitled to delete it.
  if (Time_pt!=0)
   {
    delete Time_pt; 
    Time_pt = 0;
   }

  // We're not using the default linear solver,
  // somebody else must have built it, so that person 
  // must be in charge of killing it. 
  // We can safely delete the defaults, however
  delete Default_linear_solver_pt;

  delete Default_eigen_solver_pt;
  delete Default_assembly_handler_pt;
  delete Communicator_pt;
  delete Dof_distribution_pt;

 // Delete any copies of the problem that have been created for
 // use in adaptive bifurcation tracking.
 // ALH: This will eventually go
 unsigned n_copies = Copy_of_problem_pt.size();
 for(unsigned c=0;c<n_copies;c++)
  {
   delete Copy_of_problem_pt[c];
  }

 // if this problem has sub meshes then we must delete the Mesh_pt
 if (Sub_mesh_pt.size() != 0)
  {
   Mesh_pt->flush_element_and_node_storage();
   delete Mesh_pt;
  }
 }


#ifdef OOMPH_HAS_MPI

 //==================================================================
 /// Distribute the problem without doc; report stats if required
 //==================================================================
 Vector<unsigned>& Problem::distribute(const bool& report_stats)
  {
   // Set dummy doc paramemters
   DocInfo doc_info;
   doc_info.doc_flag()=false;
   // Set the sizes of the input and output vectors
   unsigned n_element=mesh_pt()->nelement();
   Vector<unsigned> element_partition(n_element,0);
   Element_partition.resize(n_element);
   // Distribute
   Element_partition=distribute(element_partition,doc_info,report_stats);
   // Return partition that was used
   return Element_partition;
  }

 //==================================================================
 /// Distribute the problem according to specified partition
 //==================================================================
 Vector<unsigned>& Problem::distribute
 (const Vector<unsigned>& element_partition, const bool& report_stats)
  {
   // Set dummy doc paramemters
   DocInfo doc_info;
   doc_info.doc_flag()=false;
   // Set the size of the output vector
   unsigned n_element=element_partition.size();
   Element_partition.resize(n_element);
   // Distribute
   Element_partition=distribute(element_partition,doc_info,report_stats);
   // Return the partition that was used
   return Element_partition;
  }

 //==================================================================
 /// Distribute the problem and doc to specified DocInfo
 //==================================================================
 Vector<unsigned>& Problem::distribute(DocInfo& doc_info,
                                       const bool& report_stats)
  {
   // Set the sizes of the input and output vectors
   unsigned n_element=mesh_pt()->nelement();
   // Dummy input vector
   Vector<unsigned> element_partition(n_element,0);
   // Set size of output vector
   Element_partition.resize(n_element);
   // Distribute
   Element_partition=distribute(element_partition,doc_info,report_stats);
   // Return partition that was used
   return Element_partition;
  }

 //==================================================================
 /// Distribute the problem and doc
 //==================================================================
 Vector<unsigned>& Problem::distribute
 (const Vector<unsigned>& element_partition,
  DocInfo& doc_info, const bool& report_stats)
 {
  // Storage for number of processors and number of elements in global mesh
  int n_proc=this->communicator_pt()->nproc();
  int n_element=mesh_pt()->nelement();

  // Vector to be returned
  Element_partition.resize(n_element);

  // Buffer extreme cases
  if (n_proc==1) // single-process job - don't do anything
   {
    if (report_stats)
     {
      std::ostringstream warn_message;
      warn_message << "WARNING: You've tried to distribute a problem over\n"
                   << "only one processor: this would make METIS crash.\n"
                   << "Ignoring your request for distribution.\n";
      OomphLibWarning(warn_message.str(),
                      "Problem::distribute()",
                      OOMPH_EXCEPTION_LOCATION);
     }
   }
  else if (n_proc>n_element) // more processors than elements
   {
    // Throw an error
    std::ostringstream error_stream;
    error_stream << "You have tried to distribute a problem\n"
                 << "but there are less elements than processors.\n"  
                 << "Please re-run with more elements!\n"
                 << "Please also ensure that actions_before_distribute().\n"
                 << "and actions_after_distribute() are correctly set up.\n"
                 << std::endl;
    throw OomphLibError(error_stream.str(),
                        "Problem::distribute()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  else
   {
    // We only distribute uniformly-refined meshes; buffer the case where
    // either mesh is not uniformly refined
    bool a_mesh_is_not_uniformly_refined=false;
    unsigned n_mesh=nsub_mesh();
    if (n_mesh==0)
     {
      // Check refinement levels 
      if (RefineableMeshBase* mmesh_pt = 
          dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
       {
        unsigned min_ref_level=0;
        unsigned max_ref_level=0;
        mmesh_pt->get_refinement_levels(min_ref_level,max_ref_level);
        // If they are not the same
        if (max_ref_level!=min_ref_level)
         {
          a_mesh_is_not_uniformly_refined=true;
         }
       }
     }
    else
     {
      for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
       {
        // Check refinement levels for each mesh individually
        // (one mesh is allowed to be "more uniformly refined" than another)
        if (RefineableMeshBase* mmesh_pt = 
            dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
         {
          unsigned min_ref_level=0;
          unsigned max_ref_level=0;
          mmesh_pt->get_refinement_levels(min_ref_level,max_ref_level);
          // If they are not the same
          if (max_ref_level!=min_ref_level)
           {
            a_mesh_is_not_uniformly_refined=true;
           }
         }
       }
     }

    // If any mesh is not uniformly refined
    if (a_mesh_is_not_uniformly_refined)
     {
      // Again it may make more sense to throw an error here as the user
      // will probably not be running a problem that is small enough to
      // fit the whole of on each processor
      std::ostringstream error_stream;
      error_stream << "You have tried to distribute a problem\n"
                   << "but at least one of your meshes is no longer\n"  
                   << "uniformly refined.  In order to preserve the Tree\n"
                   << "and TreeForest structure, Problem::distribute() can\n"
                   << "only be called while meshes are uniformly refined.\n"
                   << std::endl;
      throw OomphLibError(error_stream.str(),
                          "Problem::distribute()",
                          OOMPH_EXCEPTION_LOCATION);
     }
    else
     {
      // Is there any global data?  If so, distributing the problem won't work
      if (nglobal_data() > 0)
       {
        std::ostringstream error_stream;
        error_stream << "You have tried to distribute a problem\n"
                     << "and there is some global data.\n"  
                     << "This is not likely to work...\n"
                     << std::endl;
        throw OomphLibError(error_stream.str(),
                            "Problem::distribute()",
                            OOMPH_EXCEPTION_LOCATION);
       }

      // Need to partition the global mesh before distributing
      Mesh* global_mesh_pt = mesh_pt();

      // Vector listing the affiliation of each element
      unsigned nelem=global_mesh_pt->nelement();
      Vector<unsigned> element_domain(nelem);

      // Partition the mesh, unless the partition has already been passed in
      // If it hasn't then the sum of all the entries of the vector should be 0
      unsigned sum_element_partition=0;
      for (unsigned e=0;e<nelem;e++)
       {
        sum_element_partition+=element_partition[e];
       }
      if (sum_element_partition==0)
       {
        partition_global_mesh(global_mesh_pt,doc_info,element_domain);
       }
      else
       {
        if (report_stats)
         {
          oomph_info << "INFO: using pre-set partition of elements"
                     << std::endl;
         }
        element_domain=element_partition;
       }

      // Set the GLOBAL Mesh_has_been_distributed flag
      global_mesh_pt->mesh_has_been_distributed()=true;

      // Set the returned vector
      Element_partition=element_domain;

      // Prepare vector of vectors for submesh element domains
      Vector<Vector<unsigned> > submesh_element_domain(n_mesh);
  
      // The submeshes need to know their own element domains
      if (n_mesh!=0)
       {
        unsigned count=0;
        for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
         {
          unsigned nsub_elem=mesh_pt(i_mesh)->nelement();
          submesh_element_domain[i_mesh].resize(nsub_elem);
          for (unsigned e=0; e<nsub_elem; e++)
           {
            submesh_element_domain[i_mesh][e]=element_domain[count];
            count++;
           }
         }
       }

      // Partitioning complete; call actions before distribute
      actions_before_distribute();

      // Distribute the (sub)meshes (i.e. sort out their halo lookup schemes)
      n_mesh=nsub_mesh();
      if (n_mesh==0)
       {
        mesh_pt()->distribute(this->communicator_pt(),
                              element_domain,doc_info,report_stats);
       }
      else // There are submeshes, "distribute" each one separately
       {
        for (unsigned i_mesh=0; i_mesh<n_mesh; i_mesh++)
         {
          if (report_stats)
           {
            oomph_info << "Distributing submesh " << i_mesh << std::endl
                       << "--------------------" << std::endl;
           }
          // Set the doc_info number to reflect the submesh
          doc_info.number()=i_mesh;
          mesh_pt(i_mesh)->distribute(this->communicator_pt(),
                                      submesh_element_domain[i_mesh],
                                      doc_info,report_stats);
         }
        // Rebuild the global mesh
        rebuild_global_mesh();
       }
      // Now the problem has been distributed
      Problem_has_been_distributed=true;

      // Call actions after distribute
      actions_after_distribute();

      // Re-assign the equation numbers (incl synchronisation if reqd)
      oomph_info << "Number of equations: " << assign_eqn_numbers()
                 << std::endl;

     } // end if to check for uniformly refined mesh(es)

   } // end if to check number of processors vs. number of elements etc.

  // Return the partition vector used in the distribution
  return Element_partition;

 }

 //==================================================================
 /// Partition the global mesh, return vector specifying the processor
 /// number for each element. Virtual so that it can be overloaded by
 /// any user; the default is to use METIS to perform the partitioning
 /// (with a bit of cleaning up afterwards to sort out "special cases").
 //==================================================================
 void Problem::partition_global_mesh(Mesh* &global_mesh_pt, DocInfo& doc_info,
                                     Vector<unsigned>& element_domain,
                                     const bool& report_stats)
 {
  // Storage for number of processors and current processor
  int n_proc=this->communicator_pt()->nproc();
  int rank=this->communicator_pt()->my_rank();

  char filename[100];
  std::ofstream some_file;

  // Doc the original mesh on proc 0
  //--------------------------------
  if (doc_info.doc_flag())
   {
    if (rank==0)
     {
      sprintf(filename,"%s/complete_mesh%i.dat",doc_info.directory().c_str(),
              doc_info.number());
      global_mesh_pt->output(filename,5);
     }
   }

//  oomph_info << "Partioning mesh....";
//  clock_t t_start=clock();

  // Partition the mesh
  //-------------------
  // METIS Objective (0: minimise edge cut; 1: minimise total comm volume)
  unsigned objective=0;

  // Do the partitioning
  METIS::partition_mesh(this->communicator_pt(),global_mesh_pt,
                        n_proc,objective,element_domain);

  // On very coarse meshes with larger numbers of processors, METIS 
  // occasionally returns an element_domain Vector for which a particular 
  // processor has no elements affiliated to it; the following fixes this

  // Convert element_domain to integer storage
  unsigned nelem=element_domain.size();
  Vector<int> int_element_domain(nelem);
  for (unsigned e=0;e<nelem;e++)
   {
    int_element_domain[e]=element_domain[e];
   }

  // Global storage for number of elements on each process
  int my_number_of_elements=0;
  Vector<int> number_of_elements(n_proc,0);

  for (unsigned e=0;e<nelem;e++)
   {
    if (int_element_domain[e]==rank)
     {
      my_number_of_elements++;
     }
   }

  // Communicate the correct value for each single process into
  // the global storage vector
  MPI_Allgather(&my_number_of_elements,1,MPI_INT,
                &number_of_elements[0],1,MPI_INT,
                this->communicator_pt()->mpi_comm());

  // If a process has no elements then switch an element with the
  // process with the largest number of elements, assuming
  // that it still has enough elements left to share
  int max_number_of_elements=0;
  int process_with_max_elements=0;
  for (int d=0;d<n_proc;d++)
   {
    if (number_of_elements[d]==0)
     {
      // Find the process with maximum number of elements
      if (max_number_of_elements<=1)
       {
        for (int dd=0;dd<n_proc;dd++)
         {
          if (number_of_elements[dd]>max_number_of_elements)
           {
            max_number_of_elements=number_of_elements[dd];
            process_with_max_elements=dd;
           }
         }
       }

      // Check that this number of elements is okay for sharing
      if (max_number_of_elements<=1)
       {
        // Throw an error if elements can't be shared
        std::ostringstream error_stream;
        error_stream << "No process has more than 1 element, and\n"
                     << "at least one process has no elements!\n"
                     << "Suggest rerunning with more refinement.\n"
                     << std::endl;
        throw OomphLibError(error_stream.str(),
                            "Problem::partition_global_mesh()",
                            OOMPH_EXCEPTION_LOCATION);

       }

      // Loop over the element domain vector and switch
      // one value for process "process_with_max_elements" with d
      for (unsigned e=0;e<nelem;e++)
       {
        if (int_element_domain[e]==process_with_max_elements)
         {
          int_element_domain[e]=d;
          // Change the numbers associated with these processes
          number_of_elements[d]++;
          number_of_elements[process_with_max_elements]--;
          // Reduce the number of elements available on "max" process
          max_number_of_elements--;
          // Inform the user that a switch has taken place
          if (report_stats)
           {
            oomph_info << "INFO: Switched element domain at position " << e 
                       << std::endl 
                       << "from process " << process_with_max_elements 
                       << " to process " << d 
                       << std::endl
                       << "which was given no elements by METIS partition" 
                       << std::endl;
           }
          // Only need to do this once for this element loop, otherwise
          // this will take all the elements from "max" process and put them
          // in process d, thus leaving essentially the same problem!
          break;
         }
       }
     }

   }

  // Reassign new values to the element_domain vector
  for (unsigned e=0;e<nelem;e++)
   {
    element_domain[e]=int_element_domain[e];
   }

  unsigned count_elements=0;
  for (unsigned e=0; e<nelem; e++)
   {
    if(int(element_domain[e])==rank)
     {
      count_elements++;
     }
   }

  if (report_stats)
   { 
    oomph_info << "I have " << count_elements
               << " elements from this partition" << std::endl << std::endl;
   }
 }

 //==================================================================
 /// (Irreversibly) prune halo(ed) elements and nodes, usually
 /// after another round of refinement, to get rid of
 /// excessively wide halo layers. Note that the current
 /// mesh will be now regarded as the base mesh and no unrefinement
 /// relative to it will be possible once this function 
 /// has been called.
 //==================================================================
 void Problem::prune_halo_elements_and_nodes(DocInfo& doc_info, 
                                             const bool& report_stats)
 {  
  // Storage for number of processors and current processor
  int n_proc=this->communicator_pt()->nproc();

  // Has the problem been distributed yet?
  if (!Problem_has_been_distributed)
   {
    oomph_info 
     << "WARNING: Problem::prune_halo_elements_and_nodes() was called on a "
     << "non-distributed Problem!" << std::endl;
    oomph_info << "Ignoring your request..." << std::endl;
   }
  else
   {
    // There are no halo layers to prune if it's a single-process job   
    if (n_proc==1)
     {
      oomph_info << "WARNING: You've tried to prune halo layers on a problem\n"
                 << "with only one processor: this is unnecessary.\n" 
                 << "Ignoring your request."
                 << std::endl << std::endl;
     }
    else
     {
      // Call actions before distribute
      actions_before_distribute();

      // Prune the halo elements and nodes of the mesh(es)
      unsigned n_mesh=nsub_mesh();
      if (n_mesh==0)
       {
        // Prune halo elements and nodes for the (single) global mesh
        mesh_pt()->prune_halo_elements_and_nodes
         (this->communicator_pt(),doc_info,report_stats);
       }
      else
       {
        // Loop over individual submeshes and prune separately
        for (unsigned i_mesh=0; i_mesh<n_mesh; i_mesh++)
         {
          mesh_pt(i_mesh)->prune_halo_elements_and_nodes
           (this->communicator_pt(),doc_info,report_stats);
         }

        // Rebuild the global mesh
        rebuild_global_mesh();
       }

      // Call actions after distribute
      actions_after_distribute();
  
      // Re-assign the equation numbers (incl synchronisation if reqd)
      oomph_info << "No. of equations: " << assign_eqn_numbers() << std::endl;
     }
   }

 }



 
#endif


//===================================================================
/// Build a single (global) mesh from a number
/// of submeshes which are passed as a vector of pointers to the
/// submeshes. The ordering is not necessarily optimal.
//==============================================================
 void Problem::build_global_mesh()
 {
#ifdef PARANOID
  //Has a global mesh already been built
  if(Mesh_pt!=0)
   {
    std::string error_message = 
     "Problem::build_global_mesh() called,\n";
    error_message += " but a global mesh has already been built:\n";
    error_message += "Problem::Mesh_pt is not zero!\n";

    throw OomphLibError(error_message,
                        "Problem::build_global_mesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  //Check that there are submeshes
  if(Sub_mesh_pt.size()==0)
   {
    std::string error_message = 
     "Problem::build_global_mesh() called,\n";
    error_message += " but there are no submeshes:\n";
    error_message += "Problem::Sub_mesh_pt has no entries\n";
   
    throw OomphLibError(error_message,
                        "Problem::build_global_mesh()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //Create an empty mesh
  Mesh_pt = new Mesh();

  //Call the rebuild function to construct the mesh
  rebuild_global_mesh();
 }

//====================================================================
/// If one of the submeshes has changed (e.g. by
/// mesh adaptation) we need to update the global mesh.
/// \b Note: The nodes boundary information refers to the
/// boundary numbers within the submesh!
/// N.B. This is essentially the same function as the Mesh constructor
/// that assembles a single global mesh from submeshes
//=====================================================================
 void Problem::rebuild_global_mesh()
 {

  //Number of submeshes
  unsigned nsub_mesh=Sub_mesh_pt.size();

  // Initialise element, node and boundary counters for global mesh
  unsigned long n_element=0;
  unsigned long n_node=0;
  unsigned n_bound=0;

  // Loop over submeshes and get total number of elements, nodes and
  // boundaries
  for(unsigned imesh=0;imesh<nsub_mesh;imesh++)
   {
    n_element += Sub_mesh_pt[imesh]->nelement();
    n_node += Sub_mesh_pt[imesh]->nnode();
    n_bound += Sub_mesh_pt[imesh]->nboundary();
   }  
 
  // Reserve storage for element and node pointers 
  Mesh_pt->Element_pt.clear();
  Mesh_pt->Element_pt.reserve(n_element);
  Mesh_pt->Node_pt.clear();
  Mesh_pt->Node_pt.reserve(n_node);
  //Resize vector of vectors of nodes
  Mesh_pt->Boundary_node_pt.clear();
  Mesh_pt->Boundary_node_pt.resize(n_bound);
 

  // Sets of pointers to elements and nodes (to exlude duplicates -- they
  // shouldn't occur anyway but if they do, they must only be added
  // once in the global mesh to avoid trouble in the timestepping)
  std::set<GeneralisedElement*> element_set_pt;
  std::set<Node*> node_set_pt;

  //Counter for total number of boundaries in all the submeshes
  unsigned ibound_global=0;   
  //Loop over the number of submeshes 
  for(unsigned imesh=0;imesh<nsub_mesh;imesh++)
   {
    //Loop over the elements of the submesh and add to vector
    //duplicates are ignored
    unsigned nel_before=0;
    unsigned long n_element=Sub_mesh_pt[imesh]->nelement();
    for (unsigned long e=0;e<n_element;e++)
     {
      GeneralisedElement* el_pt=Sub_mesh_pt[imesh]->element_pt(e);
      element_set_pt.insert(el_pt);
      // Was it a duplicate?
      unsigned nel_now=element_set_pt.size();
      if (nel_now==nel_before)
       {
        std::ostringstream warning_stream;
        warning_stream  <<"WARNING: " << std::endl
                        <<"Element " << e << " in submesh " << imesh 
                        <<" is a duplicate \n and was ignored when assembling " 
                        <<"global mesh." << std::endl;
        OomphLibWarning(warning_stream.str(),
                        "Problem::rebuild_global_mesh()",
                        OOMPH_EXCEPTION_LOCATION);
       }
      else
       {
        Mesh_pt->Element_pt.push_back(el_pt);
       }
      nel_before=nel_now;
     }

    //Loop over the nodes of the submesh and add to vector
    //duplicates are ignored
    unsigned nnod_before=0;
    unsigned long n_node=Sub_mesh_pt[imesh]->nnode();
    for (unsigned long n=0;n<n_node;n++)
     {
      Node* nod_pt=Sub_mesh_pt[imesh]->node_pt(n);
      node_set_pt.insert(nod_pt);
      // Was it a duplicate?
      unsigned nnod_now=node_set_pt.size();
      if (nnod_now==nnod_before)
       {
        std::ostringstream warning_stream;
        warning_stream << "WARNING: " << std::endl
                       << "Node " << n << " in submesh " << imesh 
                       << " is a duplicate \n and was ignored when assembling " 
                       << "global mesh." << std::endl;
        OomphLibWarning(warning_stream.str(),
                        "Problem::rebuild_global_mesh()",
                        OOMPH_EXCEPTION_LOCATION);
       }
      else
       {
        Mesh_pt->Node_pt.push_back(nod_pt);
       }
      nnod_before=nnod_now;
     }

    //Loop over the boundaries of the submesh
    unsigned n_bound=Sub_mesh_pt[imesh]->nboundary();
    for (unsigned ibound=0;ibound<n_bound;ibound++)
     {
      //Loop over the number of nodes on the boundary and add to the 
      //global vector
      unsigned long n_bound_node=Sub_mesh_pt[imesh]->nboundary_node(ibound);
      for (unsigned long n=0;n<n_bound_node;n++)
       {
        Mesh_pt->Boundary_node_pt[ibound_global].push_back(
         Sub_mesh_pt[imesh]->boundary_node_pt(ibound,n));
       }
      //Increase the number of the global boundary counter
      ibound_global++;
     }
   } //End of loop over submeshes

 }





//================================================================
///  Add a timestepper to the problem. The function will automatically
/// create or resize the Time object so that it contains the appropriate
/// number of levels of storage.
//================================================================
 void Problem::add_time_stepper_pt(TimeStepper* const &time_stepper_pt) 
 {
  //Add the timestepper to the vector
  Time_stepper_pt.push_back(time_stepper_pt);
  //Find the number of timesteps required by the timestepper
  unsigned ndt = time_stepper_pt->ndt();
  //If time has not been allocated, create time object with the 
  //required number of time steps
  if(Time_pt==0) 
   {
    Time_pt = new Time(ndt);
    oomph_info << "Created Time with " << ndt << " timesteps" << std::endl;
   }
  else
   {
    //If the required number of time steps is greater than currently stored
    //resize the time storage
    if(ndt > Time_pt->ndt()) 
     {
      Time_pt->resize(ndt);
      oomph_info << "Resized Time to include " << ndt << " timesteps" 
                 << std::endl;
     }
    //Otherwise report that we are OK
    else
     {
      oomph_info << "Time object already has storage for " << ndt 
                 << " timesteps" << std::endl;
     }
   }
 
  //Pass the pointer to time to the timestepper
  time_stepper_pt->time_pt() = Time_pt;
 }

//================================================================
/// Set the explicit time stepper for the problem and also
/// ensure that a time object has been created.
//================================================================
 void Problem::set_explicit_time_stepper_pt(ExplicitTimeStepper* const 
                                            &explicit_time_stepper_pt) 
 {
  //Set the explicit time stepper
  Explicit_time_stepper_pt = explicit_time_stepper_pt;

  //If time has not been allocated, create time object with the 
  //required number of time steps
  if(Time_pt==0) 
   {
    Time_pt = new Time(0);
    oomph_info << "Created Time with storage for no previous timestep" 
               << std::endl;
   }
  else
   {
    oomph_info << "Time object already exists " << std::endl;
   }

 }


#ifdef OOMPH_HAS_MPI

//================================================================
/// Set default first and last elements for parallel assembly
/// of non-distributed problem.
//================================================================
 void Problem::set_default_first_and_last_element_for_assembly()
 {
  // Resize and make default assignments
  int n_proc=this->communicator_pt()->nproc();
  unsigned n_elements=Mesh_pt->nelement();    
  First_el_for_assembly.resize(n_proc,0);
  Last_el_for_assembly.resize(n_proc,n_elements-1);
  
  // In the absence of any better knowledge distribute work evenly 
  // over elements
  unsigned long range = 
   static_cast<unsigned long>(double(n_elements)/double(n_proc));
  for (int p=0;p<n_proc;p++)
   {
    First_el_for_assembly[p] = p*range;
    Last_el_for_assembly[p] = (p+1)*range-1;
   }
  
  // Last one needs to incorporate any dangling elements
  Last_el_for_assembly[n_proc-1] = n_elements-1;
 
  // Doc
  if (n_proc>1)
   {
    oomph_info << "\nProblem is not distributed. Parallel assembly of "
               << "Jacobian uses default partitioning: "<< std::endl;
    for (int p=0;p<n_proc;p++)
     {
      oomph_info << "Proc " << p << " assembles from element " 
                 <<  First_el_for_assembly[p] << " to " 
                 <<  Last_el_for_assembly[p] << " \n"; 
     }
   }
 }
 

//=======================================================================
/// Helper function to re-assign the first and last elements to be 
/// assembled by each processor during parallel assembly for 
/// non-distributed problem. On each processor the vector 
/// elemental_assembly_time must contain the timings for the assembly
/// of the elements. Each processor ONLY fills in the timings for
/// elements it's in charge of when using the default distribution
/// which is re-assigned every time assign_eqn_numbers() is called.
/// First and last elements are then re-assigned to load-balance
/// any subsequent assemblies.
//=======================================================================
 void Problem::recompute_load_balanced_assembly(
  Vector<double>& elemental_assembly_time)
 {
  // Wait until all processes have completed/timed their assembly
  MPI_Barrier(this->communicator_pt()->mpi_comm());

  // Storage for number of processors and current processor
  int n_proc=this->communicator_pt()->nproc();
  int rank=this->communicator_pt()->my_rank();
  
  // Setup vectors storing the number of element timings to be sent
  // and the offset in the final vector
  Vector<int> receive_count(n_proc);
  Vector<int> displacement(n_proc);
  int offset=0;
  for (int p=0;p<n_proc;p++)
   {
    // Default distribution of labour
    unsigned el_lo = First_el_for_assembly[p];
    unsigned el_hi = Last_el_for_assembly[p];;
     
    // Number of timings to be sent and offset from start in
    // final vector
    receive_count[p]=el_hi-el_lo+1;
    displacement[p]=offset;
    offset+=el_hi-el_lo+1;
   }
      
  // Gather timings on root processor:
  unsigned n_e = Last_el_for_assembly[rank]-First_el_for_assembly[rank]+1;
  double* el_ass_time = new double[n_e];
  for (unsigned i = 0; i < n_e; i++)
   {
    el_ass_time[i] = elemental_assembly_time[i];
   }
  MPI_Gatherv(
   el_ass_time,n_e,MPI_DOUBLE,
   &elemental_assembly_time[0],&receive_count[0],&displacement[0],
   MPI_DOUBLE,0,this->communicator_pt()->mpi_comm());
  delete[] el_ass_time;

  // We have determined load balancing for current setup.
  // This can remain the same until assign_eqn_numbers() is called
  // again -- the flag is re-set to true there.
  Must_recompute_load_balance_for_assembly=false;
   
  // Vector of first and last elements for each processor
  Vector<int> first_and_last_element(2);
   
  // Re-distribute work
  if (rank==0)
   {
    oomph_info 
     << std::endl
     << "Re-assigning distribution of element assembly over processors:" 
     << std::endl;
     
    // Get total assembly time
    double total=0.0;
    unsigned n_elements=Mesh_pt->nelement();    
    for (unsigned e=0;e<n_elements;e++)
     {
      total+=elemental_assembly_time[e];
     }
     
    // Target load per processor
    double target_load=total/double(n_proc);
     
    // We're on the root processor: Always start with the first element
    int proc=0;
    First_el_for_assembly[0]=0;
     
    // Initialise total work allocated
    total=0.0;
    for (unsigned e=0;e<n_elements;e++)
     {
      total+=elemental_assembly_time[e];
       
      if (total>target_load)
       {
         
        // Keep it local on root
        if (proc==0)
         {
          // Last element for current processor
          Last_el_for_assembly[0]=e;
           
          // First element for next one
          first_and_last_element[0]=e+1;
           
          // Doc
          oomph_info 
           << "Processor " << 0 << " assembles Jacobians" 
           <<  " from elements " << First_el_for_assembly[0] << " to " 
           <<  Last_el_for_assembly[0] << " " 
           << std::endl;
         }
        else if (proc<(n_proc-1))
         {
          // Last element for current processor
          first_and_last_element[1]=e;

          // Send two ints to processor p:
          MPI_Send(&first_and_last_element[0],2,MPI_INT,proc,0,
                   this->communicator_pt()->mpi_comm());

          // Doc
          oomph_info 
           << "Processor " << proc << " assembles Jacobians" 
           <<  " from elements " << first_and_last_element[0] << " to " 
           <<  first_and_last_element[1] << " " 
           << std::endl;


          // Set first element for next one
          first_and_last_element[0]=e+1;
         }
         
        // Move on to the next processor
        proc++;
         
        // Re-initialise
        total=0.0;

       } // end of test for "total exceeds target"
     }
     
    // Last one
    first_and_last_element[1]=n_elements-1;

    // Send two ints to processor p:
    MPI_Send(&first_and_last_element[0],2,MPI_INT,n_proc-1,
             0,this->communicator_pt()->mpi_comm());
     
    // Doc
    oomph_info 
     << "Processor " << n_proc-1 << " assembles Jacobians" 
     <<  " from elements " << first_and_last_element[0] << " to " 
     <<  first_and_last_element[1] << " " 
     << std::endl;
    
   }
  // Receive first and last element from root on non-master processors
  else
   {
    Vector<int> aux(2);
    MPI_Status status;
    MPI_Recv(&aux[0],2,MPI_INT,0,0,
             this->communicator_pt()->mpi_comm(),&status);
    First_el_for_assembly[rank]=aux[0];
    Last_el_for_assembly[rank]=aux[1];
   }

  // Wipe all others
  for (int p=0;p<n_proc;p++)
   {
    if (p!=rank)
     {
      First_el_for_assembly[p]=0;
      Last_el_for_assembly[p]=0;
     }
   }

 }
 
#endif
//================================================================
/// Assign all equation numbers for problem: Deals with global
/// data (= data that isn't attached to any elements) and then
/// does the equation numbering for the elements.  Bool argument
/// can be set to false to ignore assigning local equation numbers 
/// (necessary in the parallel implementation of locate_zeta 
/// between multiple meshes).
//================================================================
unsigned long Problem::assign_eqn_numbers(const bool& assign_local_eqn_numbers)
{

#ifdef OOMPH_HAS_MPI

 // Storage for number of processors
 int n_proc=this->communicator_pt()->nproc();

 // If the problem has been distributed we first have to 
 // classify any potentially newly created nodes as
 // halo or haloed (or neither)
 if (Problem_has_been_distributed)
  {
   // Classify any so-far unclassified nodes
   // Perform at submesh level
   unsigned nmesh=nsub_mesh();

#ifdef PARANOID
   bool report_stats=true;
#endif
   if (nmesh==0)
    {
     // Cast to a refineable mesh and call synchronisation routine
     if(RefineableMeshBase* mmesh_pt = 
        dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
      {
       unsigned ncont_interpolated_values=dynamic_cast<RefineableElement*>
           (mmesh_pt->element_pt(0))->ncont_interpolated_values();

       mmesh_pt->synchronise_hanging_nodes(this->communicator_pt(),
                                           ncont_interpolated_values);
      }

     // Now classify all the halo and haloed nodes
#ifdef PARANOID 
     mesh_pt()->
      classify_halo_and_haloed_nodes(this->communicator_pt(),report_stats);
#else
     mesh_pt()->classify_halo_and_haloed_nodes(this->communicator_pt());
#endif

     // The degrees of freedom must be synchronised 
     // before calling check_halo_schemes
     synchronise_dofs(mesh_pt());

#ifdef PARANOID
     // Check that the halo schemes are okay if all eqn numbers 
     // are being assigned
     if (assign_local_eqn_numbers)
      {
       oomph_info << "Calling check_halo_schemes() from assign_eqn_numbers()"
                  << std::endl << std::endl;
       check_halo_schemes();
      }
#endif

    }
   else // nmesh!=0
    {
     // Check the synchronicity of hanging nodes for a refineable mesh
     // - In a multi-physics case this should be called for every mesh
     // - It is also possible that a single mesh contains different elements 
     //   (with different values of ncont_interpolated_values);
     //   in this instance, this routine needs a rethink.
     for (unsigned imesh=0;imesh<nmesh;imesh++)
      {
       // Cast to a refineable mesh and call synchronisation routine
       if(RefineableMeshBase* mmesh_pt = 
          dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
         {
          unsigned n_int_values=0;
          // Ensure the mesh has elements on this processor
          if (mesh_pt(imesh)->nelement()>0)
           {
            n_int_values=dynamic_cast<RefineableElement*>
             (mmesh_pt->element_pt(0))->ncont_interpolated_values();
           }

          unsigned ncont_interpolated_values=0;
          // Need to use the largest value of n_int_values 
          // when calling the routine
          MPI_Allreduce(&n_int_values,&ncont_interpolated_values,1,
                        MPI_INT,MPI_MAX,this->communicator_pt()->mpi_comm());

          // All processes call the routine to ensure correct communications
          mmesh_pt->synchronise_hanging_nodes(this->communicator_pt(),
                                              ncont_interpolated_values);
         }
      }

     // Now classify all halo/haloed nodes on each mesh in turn
     for (unsigned imesh=0;imesh<nmesh;imesh++)
      {
#ifdef PARANOID
       mesh_pt(imesh)->
        classify_halo_and_haloed_nodes(this->communicator_pt(),report_stats);
#else
       mesh_pt(imesh)->classify_halo_and_haloed_nodes(this->communicator_pt());
#endif

       // Degrees of freedom must be synchronised before
       // calling check_halo_schemes
       synchronise_dofs(mesh_pt(imesh));
      }

#ifdef PARANOID
     // Check that the halo schemes are okay if all eqn numbers 
     // are being assigned
     if (assign_local_eqn_numbers)
      {
       oomph_info << "Calling check_halo_schemes() from assign_eqn_numbers()"
                  << std::endl << std::endl;
       check_halo_schemes();
      }
#endif

    }
  }
  // Re-distribution of elements over processors during assembly
  // must be recomputed
  else
   {
    if (n_proc>1)
     {
      // Force re-analysis of time spent on assembly each
      // elemental Jacobian
      Must_recompute_load_balance_for_assembly=true;
     }
    else
     {
      Must_recompute_load_balance_for_assembly=false;
     }

    // Set default first and last elements for parallel assembly
    // of non-distributed problem.
    set_default_first_and_last_element_for_assembly();
     
   }

#endif

  //(Re)-set the dof pointer to zero length because entries are 
  //pushed back onto it -- if it's not reset here then we get into
  //trouble during mesh refinement when we reassign all dofs
  Dof_pt.resize(0);
  unsigned long n_dof=0;

  //Reset the equation number
  unsigned long equation_number=0;

  // Loop over all elements in the mesh and set up any additional 
  // dependencies that they may have (e.g. storing the geometric
  // Data, i.e. Data that affects an element's shape in elements
  // with algebraic node-update functions
  unsigned nel=Mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    Mesh_pt->element_pt(e)->complete_setup_of_dependencies();
   }

#ifdef OOMPH_HAS_MPI
 // Complete setup of dependencies for external halo elements too
 unsigned n_mesh=nsub_mesh();
 for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
  {
   for (int iproc=0;iproc<n_proc;iproc++)
    {
     unsigned n_ext_halo_el=mesh_pt(i_mesh)->nexternal_halo_element(iproc);
     for (unsigned e=0;e<n_ext_halo_el;e++)
      {
       mesh_pt(i_mesh)->external_halo_element_pt(iproc,e)
        ->complete_setup_of_dependencies();
      }
    }
  }
#endif

 //Now set equation numbers for the global Data
 unsigned Nglobal_data = nglobal_data();
 for(unsigned i=0;i<Nglobal_data;i++)
  {Global_data_pt[i]->assign_eqn_numbers(equation_number,Dof_pt);}


 //Check that the Mesh_pt has been assigned
 if(Mesh_pt==0)
  {
   std::string error_message =
    "(Global) Mesh_pt must be assigned before calling ";
   error_message += " assign_eqn_numbers()\n";

   throw OomphLibError(error_message,
                       "Problem::assign_eqn_numbers()",
                       OOMPH_EXCEPTION_LOCATION);
  }
   
 // Loop over the submeshes: Note we need to call the submeshes' own
 // assign_*_eqn_number() otherwise we miss additional functionality
 // that is implemented (e.g.) in SolidMeshes!
 unsigned nsub_mesh=Sub_mesh_pt.size();
 if (nsub_mesh==0)
  {
   n_dof=Mesh_pt->assign_global_eqn_numbers(Dof_pt);
   if (assign_local_eqn_numbers)
    {
     Mesh_pt->assign_local_eqn_numbers();
    }
   //Clear out the temporary global storage for hijacked equation numbers
   //HijackedElementBase::reset_hijacked_data_pt();
  }
 else
  {
   //Assign global equation numbers first
   for (unsigned i=0;i<nsub_mesh;i++)
    {
     Sub_mesh_pt[i]->assign_global_eqn_numbers(Dof_pt);
     n_dof=Dof_pt.size();
    }

   //Assign local equation numbers if required
   if (assign_local_eqn_numbers)
    {
     for (unsigned i=0;i<nsub_mesh;i++)
      {
       Sub_mesh_pt[i]->assign_local_eqn_numbers();
      }
    }
   //Clear out the temporaray global storage for hijacked equation numbers
   //HijackedElementBase::reset_hijacked_data_pt();
   n_dof=Dof_pt.size();
  }


#ifdef OOMPH_HAS_MPI
  // reset previous allocation
  Parallel_sparse_assemble_previous_allocation = 0;
 
  // Only synchronise if the problem has actually been
  // distributed.
  if (Problem_has_been_distributed)
   {
   
    // Synchronise the equation numbers and return the total
    // number of degrees of freedom in the overall problem
    n_dof=synchronise_eqn_numbers(assign_local_eqn_numbers);
   }
  // just setup the Dof_distribution_pt
  // NOTE - this is setup by synchronise_eqn_numbers(...)
  // if Problem_has_been_distributed
  else
#endif
   {
    Dof_distribution_pt->rebuild(Communicator_pt,n_dof,false);
   }

  // and return the total number of DOFs
  return n_dof;
 }


//================================================================
/// Get the vector of dofs, i.e. a vector containing the current
/// values of all unknowns.
//================================================================
 void Problem::get_dofs(DoubleVector& dofs)
 {
  //Find number of dofs
  const unsigned long n_dof = ndof();

  //Resize the vector
  dofs.build(Dof_distribution_pt,0.0);

  //Copy dofs into vector
  for(unsigned long l=0;l<n_dof;l++)
   {
    dofs[l] = *Dof_pt[l];
   }
 }

//=======================================================================
/// Function that sets the values of the dofs in the object
//======================================================================
 void Problem::set_dofs(const DoubleVector &dofs)
 {
  const unsigned long n_dof = this->ndof();
#ifdef PARANOID
  if(n_dof != dofs.nrow())
   {
    std::ostringstream error_stream;
    error_stream << "Number of degrees of freedom in vector argument "
                 << dofs.nrow() << "\n"
                 << "does not equal number of degrees of freedom in problem "
                 << n_dof;
    throw OomphLibError(error_stream.str(),
                        "Problem::set_dofs()",
                        OOMPH_EXCEPTION_LOCATION);  
   } 
#endif 
  for(unsigned long l=0;l<n_dof;l++)
   {
    *Dof_pt[l] = dofs[l];
   }
 }

//===================================================================
///Function that adds the values to the dofs
//==================================================================
 void Problem::add_to_dofs(const double &lambda,
                           const DoubleVector &increment_dofs)
 {
  const unsigned long n_dof = this->ndof();
  for(unsigned long l=0;l<n_dof;l++)
   {
    *Dof_pt[l] += lambda*increment_dofs[l];
   }
 }



//=========================================================================
///Return the residual vector multiplied by the inverse mass matrix
///Virtual so that it can be overloaded for mpi problems
//=========================================================================
 void Problem::get_inverse_mass_matrix_times_residuals(DoubleVector &Mres)
 {
  //Find the number of degrees of freedom in the problem
  const unsigned n_dof = this->ndof();

  //Resize the vector
  LinearAlgebraDistribution dist(this->communicator_pt(),n_dof,false);
  Mres.build(&dist,0.0);
 
  //If we have discontinuous formulation
  if(Discontinuous_element_formulation)
   {
    //Loop over the elements and get their residuals
    const unsigned n_element = Problem::mesh_pt()->nelement();
    Vector<double> element_Mres;
    for(unsigned e=0;e<n_element;e++)
     {
      //Cache the element
      DGElement* const elem_pt =
       dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));
     
      const unsigned n_el_dofs = elem_pt->ndof();
      elem_pt->get_inverse_mass_matrix_times_residuals(element_Mres);
     
      for(unsigned i=0;i<n_el_dofs;i++)
       {
        Mres[elem_pt->eqn_number(i)] = element_Mres[i];
       }
     }
   }
  //Otherwise it's continous
  else
   {
    //Now do the linear solve -- recycling Mass matrix if requested
    //If we already have the factorised mass matrix, then resolve
    if(Mass_matrix_reuse_is_enabled && Mass_matrix_has_been_computed)
     {     
      if(!Shut_up_in_newton_solve) 
       {
        oomph_info << "Not recomputing Mass Matrix " << std::endl;
       }
     
      //Get the residuals
      DoubleVector residuals(&dist,0.0);
      this->get_residuals(residuals);
     
      // Resolve the linear system
      this->linear_solver_pt()->resolve(residuals,Mres);
     }
    //Otherwise solve for the first time
    else
     {
      //If we wish to reuse the mass matrix, then enable resolve
      if(Mass_matrix_reuse_is_enabled)
       {
        if(!Shut_up_in_newton_solve) 
         {
          oomph_info << "Enabling resolve in explicit timestep" << std::endl;
         }
        this->linear_solver_pt()->enable_resolve();
       }
     
      //Store the old assembly handler
      AssemblyHandler* old_assembly_handler_pt = this->assembly_handler_pt();
      //Set the assembly handler to the explicit timestep handler
      this->assembly_handler_pt() = new ExplicitTimeStepHandler;
     
      //Solve the linear system
      this->linear_solver_pt()->solve(this,Mres);
      //The mass matrix has now been computed
      Mass_matrix_has_been_computed=true;
     
      //Delete the Explicit Timestep handler
      delete this->assembly_handler_pt();
      //Reset the assembly handler to the original handler
      this->assembly_handler_pt() = old_assembly_handler_pt;
     }
   }
 }


//================================================================
/// Get the total residuals Vector for the problem
//================================================================
 void Problem::get_residuals(DoubleVector &residuals)
 {
  // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true 
  // this means MPI_Helpers::init() has been called.  This could happen on a
  // code compiled with MPI but run serially; in this instance the
  // get_residuals function still works on one processor.
  //
  // Secondly, if a code has been compiled with MPI, but MPI_Helpers::init()
  // has not been called, then MPI_Helpers::MPI_has_been_initialised=false 
  // and the code calls...
  //
  // Thirdly, the serial version (compiled by all, but only run when compiled
  // with MPI if MPI_Helpers::MPI_has_been_initialised=false
 
  // start by setting the distribution of the residuals vector if it is not 
  // setup
  if (!residuals.distribution_setup())
   {
    int n_dof=ndof();
    LinearAlgebraDistribution dist(Communicator_pt,n_dof,false);
    residuals.build(&dist,0.0);
   }
  // otherwise just zero the residuals
  else
   {
#ifdef PARANOID
    // PARANOID check - if the residuals are distributed then this method
    // cannot be used, a distributed residuals can only be assembled by
    // get_jacobian(...) for CRDoubleMatrices
    if (residuals.distributed())
     {
      throw OomphLibError
       ("this method can only assemble a non-distributed residuals vector",
        "Problem::get_residuals()",OOMPH_EXCEPTION_LOCATION);
     }
#endif

    // and zero
    residuals.initialise(0.0);
   }


#ifdef OOMPH_HAS_MPI
  if (MPI_Helpers::MPI_has_been_initialised)
   {
    // Storage for the number of processors
    int n_proc=this->communicator_pt()->nproc();
 
    // cache the number of local rows
    unsigned nrow = residuals.nrow();

    // get a pointer to the underlying values
    double* residuals_pt;
    if (n_proc>1)
     {
      residuals_pt = new double[nrow];
      for (unsigned i = 0; i < residuals.nrow(); i++)
       {
        residuals_pt[i] = 0.0;  
       }
     }
    else
     {
      residuals_pt = residuals.values_pt();
     }

    // number of elements
    int n_el=mesh_pt()->nelement();
 
    // Default assignments for distributed problem
    unsigned j_lo=0;
    unsigned j_hi=n_el;
 
    // Otherwise just loop over fractions of the elements
    if (!Problem_has_been_distributed)
     {
      // Distribute work evenly
      unsigned range=unsigned(double(n_el)/double(Communicator_pt->nproc()));
      j_lo=Communicator_pt->my_rank()*range;
      j_hi=(Communicator_pt->my_rank()+1)*range;
    
      // Last one needs to incorporate any dangling elements
      if (Communicator_pt->my_rank() == Communicator_pt->nproc()-1)
       {
        j_hi=n_el;
       }
     }

    // Assemble the partial residual vector: Note that this
    // is a full-length vector but only contributions from 
    // a sub-set of elements are filled in, other values
    // are set to zero
 
    //Loop over all the elements for this processor
    for(unsigned long e=j_lo;e<j_hi;e++)
     {
      // Get element
      GeneralisedElement* el_pt=mesh_pt()->element_pt(e);
   
      // Is it a halo?
      if (!el_pt->is_halo())
       {
        //Find number of dofs in the element
        const unsigned n_el_dofs = el_pt->ndof();
     
        //Set up a Vector
        Vector<double> element_residuals(n_el_dofs);
     
        //Fill the array
        el_pt->get_residuals(element_residuals);

        //Now loop over the dofs and assign values to global Vector
        for(unsigned l=0;l<n_el_dofs;l++)
         {
          residuals_pt[el_pt->eqn_number(l)]+=element_residuals[l];
         }
       }
     }
 
    // Receive from the other processors and assemble if required
    if (n_proc>1)
     {
      // clear and resize residuals
      MPI_Allreduce(residuals_pt,residuals.values_pt(),nrow,
                    MPI_DOUBLE,MPI_SUM,Communicator_pt->mpi_comm());
      delete[] residuals_pt;
     }
   
   
   }
  else // !MPI_Helpers::MPI_has_been_initialised
#endif // OOMPH_HAS_MPI
   {

    //Locally cache pointer to assembly handler
    AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

    //Loop over all the elements
    unsigned long Element_pt_range = Mesh_pt->nelement();
    for(unsigned long e=0;e<Element_pt_range;e++)
     {
      //Get the pointer to the element
      GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);
      //Find number of dofs in the element
      unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt); 
      //Set up an array
      Vector<double> element_residuals(n_element_dofs);
      //Fill the array
      assembly_handler_pt->get_residuals(elem_pt,element_residuals);
      //Now loop over the dofs and assign values to global Vector
      for(unsigned l=0;l<n_element_dofs;l++)
       {
        residuals[assembly_handler_pt->eqn_number(elem_pt,l)] 
         += element_residuals[l];
       }
     }
   }
 }

//=============================================================================
/// Get the fully assembled residual vector and Jacobian matrix
/// in dense storage. The DoubleVector residuals returned will be 
/// non-distributed. If on calling this method the DoubleVector residuals is
/// setup then it must be non-distributed and of the correct length. \n
/// The matrix type DenseDoubleMatrix is not distributable and therefore 
/// the residual vector is also assumed to be non distributable.
//=============================================================================
 void Problem::get_jacobian(DoubleVector &residuals, 
                            DenseDoubleMatrix& jacobian)
 {
  // get the number of degrees of freedom
  unsigned n_dof=ndof(); 

#ifdef PARANOID
  // PARANOID checks : if the distribution of residuals is setup then it must
  // must not be distributed, have the right number of rows, and the same 
  // communicator as the problem
  if (residuals.distribution_pt()->setup())
   {
    if (residuals.distribution_pt()->distributed())
     {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must not "
                   << "be distributed.";
      throw OomphLibError(error_stream.str(),
                          "Problem::get_jacobian(...)",
                          OOMPH_EXCEPTION_LOCATION);
     }
    if (residuals.distribution_pt()->nrow() != n_dof)
     {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must have"
                   << " the correct number of rows";
      throw OomphLibError(error_stream.str(),
                          "Problem::get_jacobian(...)",
                          OOMPH_EXCEPTION_LOCATION);
     }
    if (!(*Communicator_pt == *residuals.distribution_pt()->communicator_pt()))
     {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must have"
                   << " the same communicator as the problem.";
      throw OomphLibError(error_stream.str(),
                          "Problem::get_jacobian(...)",
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif

  // set the residuals distribution if it is not setup
  if (!residuals.distribution_pt()->setup())
   {
    LinearAlgebraDistribution dist(Communicator_pt,n_dof,false);
    residuals.build(&dist,0.0);
   }
  // else just zero the residuals
  else
   {
    residuals.initialise(0.0);
   }

  // Resize the matrices -- this cannot always be done externally
  // because get_jacobian exists in many different versions for
  // different storage formats -- resizing a CC or CR matrix doesn't
  // make sense.

  // resize the jacobian
  jacobian.resize(n_dof,n_dof);
  jacobian.initialise(0.0);
 
  //Locally cache pointer to assembly handler
  AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

  //Loop over all the elements
  unsigned long n_element = Mesh_pt->nelement();
  for(unsigned long e=0;e<n_element;e++)
   {
    //Get the pointer to the element
    GeneralisedElement* elem_pt = Mesh_pt->element_pt(e);  
    //Find number of dofs in the element
    unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt); 
    //Set up an array
    Vector<double> element_residuals(n_element_dofs);
    //Set up a matrix
    DenseMatrix<double> element_jacobian(n_element_dofs);
    //Fill the array
    assembly_handler_pt->get_jacobian(elem_pt,
                                      element_residuals,element_jacobian);
    //Now loop over the dofs and assign values to global Vector
    for(unsigned l=0;l<n_element_dofs;l++)
     {
      unsigned long eqn_number = assembly_handler_pt->eqn_number(elem_pt,l);
      residuals[eqn_number] += element_residuals[l];
      for(unsigned l2=0;l2<n_element_dofs;l2++)
       {
        jacobian(eqn_number ,
                 assembly_handler_pt->eqn_number(elem_pt,l2)) += 
         element_jacobian(l,l2);
       }
     }
   }
 }

//=============================================================================
/// Return the fully-assembled Jacobian and residuals for the problem,
/// in the case where the Jacobian matrix is in a distributable 
/// row compressed storage format. \n
/// 1. If the distribution of the jacobian and residuals is setup then, they 
/// will be returned with that distribution.
/// Note. the jacobian and residuals must have the same distribution. \n
/// 2. If the distribution of the jacobian and residuals are not setup then
/// their distribution will computed based on:
/// Distributed_problem_matrix_distribution.
//=============================================================================
 void Problem::get_jacobian(DoubleVector &residuals, CRDoubleMatrix &jacobian)
 {

  // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true 
  // this means MPI_Helpers::setup() has been called.  This could happen on a
  // code compiled with MPI but run serially; in this instance the
  // get_residuals function still works on one processor.
  //
  // Secondly, if a code has been compiled with MPI, but MPI_Helpers::setup()
  // has not been called, then MPI_Helpers::MPI_has_been_initialised=false 
  // and the code calls...
  //
  // Thirdly, the serial version (compiled by all, but only run when compiled
  // with MPI if MPI_Helpers::MPI_has_been_initialised=false
  //
  // The only case where an MPI code cannot run serially at present
  // is one where the distribute function is used (i.e. METIS is called)

  //Allocate storage for the matrix entries
  //The generalised Vector<Vector<>> structure is required
  //for the most general interface to sparse_assemble() which allows
  //the assembly of multiple matrices at once.
  Vector<int* > column_index(1);
  Vector<int* > row_start(1);
  Vector<double* > value(1); 
  Vector<unsigned> nnz(1);

#ifdef PARANOID
  // PARANOID checks that the distribution of the jacobian matches that of the
  // residuals (if they are setup) and that they have the right number of rows
  if (residuals.distribution_pt()->setup() && 
      jacobian.distribution_pt()->setup())
   {
    if (!(*residuals.distribution_pt() == *jacobian.distribution_pt()))
     {                                    
      std::ostringstream error_stream;
      error_stream << "If the distribution of the residuals must "
                   << "be the same as the distribution of the jacobian."; 
      throw OomphLibError(error_stream.str(),              
                          "Problem::get_jacobian(...)", 
                          OOMPH_EXCEPTION_LOCATION); 
     }                                                 
    if (jacobian.distribution_pt()->nrow() != this->ndof())
     {
      std::ostringstream error_stream;       
      error_stream << "The distribution of the jacobian and residuals does not"
                   << "have the correct number of global rows.";
      throw OomphLibError(error_stream.str(),                      
                          "Problem::get_jacobian(...)",       
                          OOMPH_EXCEPTION_LOCATION);              
     }                 
   }
  else if (residuals.distribution_pt()->setup() != 
           jacobian.distribution_pt()->setup())
   {
    std::ostringstream error_stream; 
    error_stream << "The distribution of the jacobian and residuals must "
                 << "both be setup or both not setup";
    throw OomphLibError(error_stream.str(),                           
                        "Problem::get_jacobian(...)",              
                        OOMPH_EXCEPTION_LOCATION); 
   } 
#endif 


  //Allocate generalised storage format for passing to sparse_assemble()
  Vector<double* > res(1);
 
  // number of rows
  unsigned nrow = this->ndof();

  // determine the distribution for the jacobian.
  // IF the jacobian has distribution setup then use that
  // ELSE determine the distribution based on the 
  // distributed_matrix_distribution enum
  LinearAlgebraDistribution* dist_pt=0;
  if (jacobian.distribution_pt()->setup())
   {
    dist_pt = new LinearAlgebraDistribution(jacobian.distribution_pt());
   }
  else
   {
#ifdef OOMPH_HAS_MPI
    // if problem is only one one processor
    if (Communicator_pt->nproc() == 1)
     {
      dist_pt = new LinearAlgebraDistribution(Communicator_pt,nrow,false);
     }
    // if the problem is not distributed then assemble the jacobian with
    // a uniform distributed distribution
    else if (!Problem_has_been_distributed)
     {
      dist_pt = new LinearAlgebraDistribution(Communicator_pt,nrow,true);
     }
    // otherwise the problem is a distributed problem
    else
     {
      switch (Dist_problem_matrix_distribution)
       {
       case Uniform_matrix_distribution:
        dist_pt = new LinearAlgebraDistribution(Communicator_pt,nrow,true);
        break;
       case Problem_matrix_distribution:
        dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
        break;
       case Default_matrix_distribution:
        LinearAlgebraDistribution* uniform_dist_pt = 
         new LinearAlgebraDistribution(Communicator_pt,nrow,true);
        bool use_problem_dist = true;
        unsigned nproc = Communicator_pt->nproc();
        for (unsigned p = 0; p < nproc; p++)
         {
          if ((double)Dof_distribution_pt->nrow_local(p) > 
              ((double)uniform_dist_pt->nrow_local(p))*1.1)
           {
            use_problem_dist = false;
           }
         }
        if (use_problem_dist)
         {
          dist_pt = new LinearAlgebraDistribution(Dof_distribution_pt);
         }
        else
         {
          dist_pt = new LinearAlgebraDistribution(uniform_dist_pt);
         }
        delete uniform_dist_pt;
        break;
       }
     }
#else
    dist_pt = new LinearAlgebraDistribution(Communicator_pt,nrow,false);
#endif
   }


  //The matrix is in compressed row format
  bool compressed_row_flag=true;
  
#ifdef OOMPH_HAS_MPI
  //
  if (Communicator_pt->nproc() == 1)
   {
#endif
    sparse_assemble_row_or_column_compressed(column_index,
                                             row_start,
                                             value,
                                             nnz,
                                             res,
                                             compressed_row_flag);
    jacobian.build(dist_pt);
    jacobian.build_matrix_without_copy(dist_pt->nrow(),nnz[0],
                                       value[0],column_index[0],row_start[0]);
    residuals.build(dist_pt,0.0);
    residuals.set_external_values(res[0],true);
#ifdef OOMPH_HAS_MPI
   }
  else
   {
    if (dist_pt->distributed())
     {
      parallel_sparse_assemble(dist_pt,
                               column_index,
                               row_start,
                               value,
                               nnz,
                               res);
      jacobian.build(dist_pt);
      jacobian.build_matrix_without_copy(dist_pt->nrow(),nnz[0],
                                         value[0],column_index[0],
                                         row_start[0]); 
      residuals.build(dist_pt,0.0);
      residuals.set_external_values(res[0],true);   
     }
    else
     {
      LinearAlgebraDistribution* temp_dist_pt = 
       new LinearAlgebraDistribution(Communicator_pt,dist_pt->nrow(),true);
      parallel_sparse_assemble(temp_dist_pt,
                               column_index,
                               row_start,
                               value,
                               nnz,
                               res);
      jacobian.build(temp_dist_pt);
      jacobian.build_matrix_without_copy(dist_pt->nrow(),nnz[0],
                                         value[0],column_index[0],
                                         row_start[0]); 
      jacobian.redistribute(dist_pt);
      residuals.build(temp_dist_pt,0.0);
      residuals.set_external_values(res[0],true);   
      residuals.redistribute(dist_pt);
      delete temp_dist_pt;
     }
   }
#endif

  // clean up dist_pt and residuals_vector pt
  delete dist_pt;
 }

//=============================================================================
/// Return the fully-assembled Jacobian and residuals for the problem,
/// in the case when the jacobian matrix is in column-compressed storage
/// format.
//=============================================================================
 void Problem::get_jacobian(DoubleVector &residuals, CCDoubleMatrix &jacobian)
 {
  // Three different cases; if MPI_Helpers::MPI_has_been_initialised=true 
  // this means MPI_Helpers::setup() has been called.  This could happen on a
  // code compiled with MPI but run serially; in this instance the
  // get_residuals function still works on one processor.
  //
  // Secondly, if a code has been compiled with MPI, but MPI_Helpers::setup()
  // has not been called, then MPI_Helpers::MPI_has_been_initialised=false 
  // and the code calls...
  //
  // Thirdly, the serial version (compiled by all, but only run when compiled
  // with MPI if MPI_Helpers::MPI_has_been_5Binitialised=false
  //
  // The only case where an MPI code cannot run serially at present
  // is one where the distribute function is used (i.e. METIS is called)

  // get the number of degrees of freedom  
  unsigned n_dof=ndof();   

#ifdef PARANOID   
  // PARANOID checks : if the distribution of residuals is setup then it must
  // must not be distributed, have the right number of rows, and the same 
  // communicator as the problem   
  if (residuals.distribution_pt()->setup())
   {                                                            
    if (residuals.distribution_pt()->distributed())  
     {                  
      std::ostringstream error_stream;                
      error_stream 
       << "If the DoubleVector residuals is setup then it must not "
       << "be distributed.";         
      throw OomphLibError(error_stream.str(),      
                          "Problem::get_jacobian(...)",      
                          OOMPH_EXCEPTION_LOCATION);       
     }                                        
    if (residuals.distribution_pt()->nrow() != n_dof)     
     {                               
      std::ostringstream error_stream;                 
      error_stream 
       << "If the DoubleVector residuals is setup then it must have"
       << " the correct number of rows";     
      throw OomphLibError(error_stream.str(),                          
                          "Problem::get_jacobian(...)",     
                          OOMPH_EXCEPTION_LOCATION);        
     }                                 
    if (!(*Communicator_pt == *residuals.distribution_pt()->communicator_pt()))
     {                  
      std::ostringstream error_stream;      
      error_stream 
       << "If the DoubleVector residuals is setup then it must have"
       << " the same communicator as the problem.";  
      throw OomphLibError(error_stream.str(),                  
                          "Problem::get_jacobian(...)",     
                          OOMPH_EXCEPTION_LOCATION);       
     }                                                
   }                                           
#endif  

  //Allocate storage for the matrix entries
  //The generalised Vector<Vector<>> structure is required
  //for the most general interface to sparse_assemble() which allows
  //the assembly of multiple matrices at once.
  Vector<int* > row_index(1);
  Vector<int* > column_start(1);
  Vector<double* > value(1); 

  //Allocate generalised storage format for passing to sparse_assemble()
  Vector<double* > res(1);

  // allocate storage for the number of non-zeros in each matrix
  Vector<unsigned> nnz(1);
 
  //The matrix is in compressed column format
  bool compressed_row_flag=false;
  
  // get the distribution for the residuals
  LinearAlgebraDistribution* dist_pt;
  if (!residuals.distribution_setup())
   {
    dist_pt 
     = new LinearAlgebraDistribution(Communicator_pt,this->ndof(),false);
   }
  else
   {
    dist_pt = new LinearAlgebraDistribution(residuals.distribution_pt());
   }

#ifdef OOMPH_HAS_MPI
  if (communicator_pt()->nproc() == 1)
   {
#endif
    sparse_assemble_row_or_column_compressed(row_index,
                                             column_start,
                                             value,
                                             nnz,
                                             res,
                                             compressed_row_flag);    
    jacobian.build_without_copy(value[0],row_index[0],column_start[0],nnz[0],
                                n_dof,n_dof); 
    residuals.build(dist_pt,0.0);
    residuals.set_external_values(res[0],true);
#ifdef OOMPH_HAS_MPI
   }
  else
   {
   std::ostringstream error_stream;
   error_stream
    << "Cannot assemble a CCDoubleMatrix Jacobian on more "
    << "than one processor."; 
   throw OomphLibError(error_stream.str(),
                       "Problem::get_jacobian(...)",
                       OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // clean up
  delete dist_pt;
 }


//===================================================================
/// \short Set all pinned values to zero.
/// Used to set boundary conditions to be homogeneous in the copy
/// of the problem  used in adaptive bifurcation tracking 
/// (ALH: TEMPORARY HACK, WILL BE FIXED)
//==================================================================
void Problem::set_pinned_values_to_zero()
{
 //NOTE THIS DOES NOT ZERO ANY SPINE DATA, but otherwise everything else 
 //should be zeroed

 //Zero any pinned global Data
 const unsigned n_global_data = nglobal_data();
 for(unsigned i=0;i<n_global_data;i++)
  {
   Data* const local_data_pt = Global_data_pt[i];
   const unsigned n_value = local_data_pt->nvalue();
   for(unsigned j=0;j<n_value;j++)
    {
     //If the data value is pinned set the value to zero
     if(local_data_pt->is_pinned(j)) {local_data_pt->set_value(j,0.0);}
    }
  }

 // Loop over the submeshes:
 const unsigned n_sub_mesh=Sub_mesh_pt.size();
 if(n_sub_mesh==0)
  {
   //Loop over the nodes in the element
   const unsigned n_node = Mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* const local_node_pt = Mesh_pt->node_pt(n);
     const unsigned n_value = local_node_pt->nvalue();
     for(unsigned j=0;j<n_value;j++)
      {
       //If the data value is pinned set the value to zero
       if(local_node_pt->is_pinned(j)) {local_node_pt->set_value(j,0.0);}
      }
     
     //Try to cast to a solid node
     SolidNode* const local_solid_node_pt = 
      dynamic_cast<SolidNode*>(local_node_pt);
     //If we are successful
     if(local_solid_node_pt)
      {
       //Find the dimension of the node
       const unsigned n_dim = local_solid_node_pt->ndim();
       //Find number of positions
       const unsigned n_position_type = local_solid_node_pt->nposition_type();
       
       for(unsigned k=0;k<n_position_type;k++)
        {
         for(unsigned i=0;i<n_dim;i++)
          {
           //If the generalised position is pinned,
           //set the value to zero
           if(local_solid_node_pt->position_is_pinned(k,i))
            {
             local_solid_node_pt->x_gen(k,i) = 0.0;
            }
          }
        }
      }
    }

   //Now loop over the element's and zero the internal data
   const unsigned n_element = Mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     GeneralisedElement* const local_element_pt = Mesh_pt->element_pt(e);
     const unsigned n_internal = local_element_pt->ninternal_data();
     for(unsigned i=0;i<n_internal;i++)
      {
       Data* const local_data_pt = local_element_pt->internal_data_pt(i);
       const unsigned n_value = local_data_pt->nvalue();
       for(unsigned j=0;j<n_value;j++)
        {
         //If the data value is pinned set the value to zero
         if(local_data_pt->is_pinned(j)) {local_data_pt->set_value(j,0.0);}
        }
      }
    } //End of loop over elements
  }
 else
  {
   //Alternatively loop over all sub meshes
   for (unsigned m=0;m<n_sub_mesh;m++)
    {
     //Loop over the nodes in the element
     const unsigned n_node = Sub_mesh_pt[m]->nnode();
     for(unsigned n=0;n<n_node;n++)
      {
       Node* const local_node_pt = Sub_mesh_pt[m]->node_pt(n);
       const unsigned n_value = local_node_pt->nvalue();
       for(unsigned j=0;j<n_value;j++)
        {
         //If the data value is pinned set the value to zero
         if(local_node_pt->is_pinned(j)) {local_node_pt->set_value(j,0.0);}
        }
       
       //Try to cast to a solid node
       SolidNode* const local_solid_node_pt = 
        dynamic_cast<SolidNode*>(local_node_pt);
       //If we are successful
       if(local_solid_node_pt)
        {
         //Find the dimension of the node
         const unsigned n_dim = local_solid_node_pt->ndim();
         //Find number of positions
         const unsigned n_position_type = 
          local_solid_node_pt->nposition_type();
         
         for(unsigned k=0;k<n_position_type;k++)
          {
           for(unsigned i=0;i<n_dim;i++)
            {
             //If the generalised position is pinned,
             //set the value to zero
             if(local_solid_node_pt->position_is_pinned(k,i))
              {
               local_solid_node_pt->x_gen(k,i) = 0.0;
              }
            }
          }
        }
      }

     //Now loop over the element's and zero the internal data
     const unsigned n_element = Sub_mesh_pt[m]->nelement();
     for(unsigned e=0;e<n_element;e++)
      {
       GeneralisedElement* const local_element_pt = 
        Sub_mesh_pt[m]->element_pt(e);
       const unsigned n_internal = local_element_pt->ninternal_data();
       for(unsigned i=0;i<n_internal;i++)
        {
         Data* const local_data_pt = local_element_pt->internal_data_pt(i);
         const unsigned n_value = local_data_pt->nvalue();
         for(unsigned j=0;j<n_value;j++)
          {
           //If the data value is pinned set the value to zero
           if(local_data_pt->is_pinned(j)) {local_data_pt->set_value(j,0.0);}
          }
        }
      } //End of loop over elements
    }
  }
}




//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors.
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by changing the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
/// We provide four different assembly methods, each with different
/// memory requirements/execution speeds. The method is set by
/// the public flag Problem::Sparse_assembly_method.
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value,
 Vector<unsigned > &nnz, 
 Vector<double* > &residuals,
 bool compressed_row_flag)
{

 // Choose the actual method
 switch(Sparse_assembly_method)
  {

  case Perform_assembly_using_vectors_of_pairs:

   sparse_assemble_row_or_column_compressed_with_vectors_of_pairs(
    column_or_row_index, 
    row_or_column_start, 
    value, 
    nnz,
    residuals,
    compressed_row_flag);

   break;

  case Perform_assembly_using_two_vectors:

   sparse_assemble_row_or_column_compressed_with_two_vectors(
    column_or_row_index, 
    row_or_column_start, 
    value,
    nnz,
    residuals,
    compressed_row_flag);

   break;

  case Perform_assembly_using_maps:

   sparse_assemble_row_or_column_compressed_with_maps(
    column_or_row_index, 
    row_or_column_start, 
    value,
    nnz,
    residuals,
    compressed_row_flag);

   break;

  case Perform_assembly_using_lists:

   sparse_assemble_row_or_column_compressed_with_lists(
    column_or_row_index, 
    row_or_column_start, 
    value,
    nnz,
    residuals,
    compressed_row_flag);

   break;

  case Perform_assembly_using_two_arrays:

  sparse_assemble_row_or_column_compressed_with_two_arrays(
    column_or_row_index, 
    row_or_column_start, 
    value,
    nnz,
    residuals,
    compressed_row_flag);

   break;

  default:

   std::ostringstream error_stream;
   error_stream
    << "Error: Incorrect value for Problem::Sparse_assembly_method" 
    << Sparse_assembly_method << std::endl
    << "It should be one of the enumeration Problem::Assembly_method" 
    << std::endl;
   throw OomphLibError(error_stream.str(),
                       "Problem::sparse_assemble_row_or_column_compressed",
                       OOMPH_EXCEPTION_LOCATION);
  }
}





//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors, using maps
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by chaging the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed_with_maps(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value,
 Vector<unsigned> &nnz,
 Vector<double* > &residuals,
 bool compressed_row_flag)
{
 //Total number of elements
 const unsigned long n_elements = mesh_pt()->nelement();

 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;

#ifdef OOMPH_HAS_MPI
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[Communicator_pt->my_rank()];
   el_hi=Last_el_for_assembly[Communicator_pt->my_rank()];
  }
#endif 

 // number of dofs
 unsigned ndof = this->ndof();

 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();

 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_or_row_index.size();

 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

//Error check dimensions
#ifdef PARANOID
 if(row_or_column_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_or_column_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_or_row_index.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_maps",
    OOMPH_EXCEPTION_LOCATION);
  }

 if(value.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error in Problem::sparse_assemble_row_or_column_compressed " 
    << std::endl
    << "value.size() " << value.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_or_row_index.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_maps",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif


//The idea behind this sparse assembly routine is to use a vector of
//maps for the entries in each row or column of the complete matrix.
//The key for each map is the global row or column number and
//the default comparison operator for integers means that each map 
//is ordered by the global row or column number. Thus, we need not
//sort the maps, that happens at each insertion of a new entry. The 
//price we pay  is that for large maps, inseration is not a 
//cheap operation. Hash maps can be used to increase the speed, but then
//the ordering is lost and we would have to sort anyway. The solution if
//speed is required is to use lists, see below.
 
 
//Set up a vector of vectors of maps of entries of each  matrix,
//indexed by either the column or row. The entries of the vector for
//each matrix correspond to all the rows or columns of that matrix. 
//The use of the map storage
//scheme, with its implicit ordering on the first index, gives
//a sparse ordered list of the entries in the given row or column. 
 Vector<Vector<std::map<unsigned,double> > > matrix_data_map(n_matrix);
 //Loop over the number of matrices being assembled and resize
 //each vector of maps to the number of rows or columns of the matrix
 for(unsigned m=0;m<n_matrix;m++) {matrix_data_map[m].resize(ndof);}
 
 //Resize the residuals vectors
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals[v] = new double[ndof];
   for (unsigned i = 0; i < ndof; i++)
    {
     residuals[v][i] = 0;
    }
  }


#ifdef OOMPH_HAS_MPI

 
 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

#endif
   
 //----------------Assemble and populate the maps-------------------------
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element.
  //This means that the storage is only allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);

  //Loop over the elements for this processor
  for(unsigned long e=el_lo;e<=el_hi;e++)
   {

#ifdef OOMPH_HAS_MPI
    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }
#endif

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
    //Ignore halo elements
    if (!elem_pt->is_halo())
     {
#endif

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

      //Resize the storage for elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}
    
      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);

      //---------------Insert the values into the maps--------------

      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {
        //Get the local equation number
        unsigned eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);
 
        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in each residuals vector
          residuals[v][eqn_number] += el_residuals[v][i];
         }

        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);

          //Loop over the matrices
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only bother to add to the map if it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              //If it's compressed row storage, then our vector of maps
              //is indexed by row (equation number)
              if(compressed_row_flag)
               {
                //Add the data into the map using the unknown as the map key
                matrix_data_map[m][eqn_number][unknown] += value;
               }
              //Otherwise it's compressed column storage and our vector is
              //indexed by column (the unknown)
              else
               {
                //Add the data into the map using the eqn_numbe as the map key
                matrix_data_map[m][unknown][eqn_number] += value;
               }
             }
           } //End of loop over matrices
         }
       }

#ifdef OOMPH_HAS_MPI
     } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
  // Time it?
  if ((!Problem_has_been_distributed)&&
      Must_recompute_load_balance_for_assembly)
   {
    elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
   }
#endif

   } //End of loop over the elements

 } //End of map assembly
 

#ifdef OOMPH_HAS_MPI
 
 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

#endif


 //-----------Finally we need to convert the beautiful map storage scheme
 //------------------------to the containers required by SuperLU
 
 //Loop over the number of matrices
 for(unsigned m=0;m<n_matrix;m++)
  {
   //Set the number of rows or columns
   row_or_column_start[m] = new int[ndof+1];
   //Counter for the total number of entries in the storage scheme
   unsigned long entry_count=0;
   row_or_column_start[m][0] = entry_count;
   
   // first we compute the number of non-zeros
   nnz[m] = 0;
   for(unsigned long i_global=0;i_global<ndof;i_global++)
    {
     nnz[m] += matrix_data_map[m][i_global].size();
    }

   // and then resize the storage
   column_or_row_index[m] = new int[nnz[m]];
   value[m] = new double[nnz[m]];

   //Now we merely loop over the number of rows or columns 
   for(unsigned long i_global=0;i_global<ndof;i_global++)
    {
     //Start index for the present row
     row_or_column_start[m][i_global] = entry_count;
     //If there are no entries in the map then skip the rest of the loop
     if(matrix_data_map[m][i_global].empty()) {continue;}

     //Loop over all the entries in the map corresponding to the given
     //row or column. It will be ordered
     
     for(std::map<unsigned,double>::iterator 
          it = matrix_data_map[m][i_global].begin();
         it!=matrix_data_map[m][i_global].end();++it)
      {
       //The first value is the column or row index
       column_or_row_index[m][entry_count] = it->first;
       //The second value is the actual data value
       value[m][entry_count] = it->second;
       //Increase the value of the counter
       entry_count++;
      }
    }
     
   //Final entry in the row/column start vector
   row_or_column_start[m][ndof] = entry_count;
  } //End of the loop over the matrices

 if (Pause_at_end_of_sparse_assembly)
  {
   oomph_info << "Pausing at end of sparse assembly." << std::endl;
   pause("Check memory usage now.");
  }
} 






//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors using lists
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by chaging the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed_with_lists(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value,
 Vector<unsigned> &nnz,
 Vector<double* > &residuals,
 bool compressed_row_flag)
{
 //Total number of elements
 const unsigned long n_elements = mesh_pt()->nelement();

 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;

#ifdef OOMPH_HAS_MPI
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[Communicator_pt->my_rank()];
   el_hi=Last_el_for_assembly[Communicator_pt->my_rank()];
  }
#endif 

 // number of dofs
 unsigned ndof = this->ndof();

 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();

 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_or_row_index.size();

 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

//Error check dimensions
#ifdef PARANOID
 if(row_or_column_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_or_column_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_or_row_index.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_lists",
    OOMPH_EXCEPTION_LOCATION);
  }

 if(value.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error in Problem::sparse_assemble_row_or_column_compressed " 
    << std::endl
    << "value.size() " << value.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_or_row_index.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_lists",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

//The idea behind this sparse assembly routine is to use a vector of
//lists for the entries in each row or column of the complete matrix.
//The lists contain pairs of entries (global row/column number, value).
//All non-zero contributions from each element are added to the lists.
//We then sort each list by global row/column number and then combine
//the entries corresponding to each row/column before adding to the
//vectors column_or_row_index and value.
 
//Note the trade off for "fast assembly" is that we will require
//more memory during the assembly phase. Then again, if we can
//only just assemble the sparse matrix, we're in real trouble.

//Set up a vector of lists of paired entries of 
//(row/column index, jacobian matrix entry). 
//The entries of the vector correspond to all the rows or columns. 
//The use of the list storage scheme, should give fast insertion 
//and fast sorts later.
 Vector<Vector<std::list<std::pair<unsigned,double> > > > 
  matrix_data_list(n_matrix);
 //Loop over the number of matrices and resize
 for(unsigned m=0;m<n_matrix;m++) {matrix_data_list[m].resize(ndof);}

 //Resize the residuals vectors
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals[v] = new double[ndof];
   for (unsigned i = 0; i < ndof; i++)
    {
     residuals[v][i] = 0;
    }
  }

#ifdef OOMPH_HAS_MPI


 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

#endif

 //------------Assemble and populate the lists-----------------------
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element.
  //This means that the stored is only allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);


  //Pointer to a single list to be used during the assembly
  std::list<std::pair<unsigned,double> > *list_pt;
    
  //Loop over the all elements
  for(unsigned long e=el_lo;e<=el_hi;e++)
   {

#ifdef OOMPH_HAS_MPI
    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }
#endif

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
    //Ignore halo elements
    if (!elem_pt->is_halo())
     {
#endif

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);
    
      //Resize the storage for the elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}
    
      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);
    
      //---------------- Insert the values into the lists -----------
      
      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {
      //Get the local equation number
        unsigned eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);

        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in the residuals vector
          residuals[v][eqn_number] += el_residuals[v][i];
         }

        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);
          
          //Loop over the matrices
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only add to theif it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              //If it's compressed row storage, then our vector is indexed
              //by row (the equation number)
              if(compressed_row_flag)
               {
                //Find the list that corresponds to the desired row
                list_pt = &matrix_data_list[m][eqn_number];
                //Insert the data into the list, the first entry
                //in the pair is the unknown (column index),
                //the second is the value itself.
                list_pt->
                 insert(list_pt->end(),std::make_pair(unknown,value));
               }
              //Otherwise it's compressed column storage, and our
              //vector is indexed by column (the unknown)
              else
               {
                //Find the list that corresponds to the desired column
                list_pt = &matrix_data_list[m][unknown];
                //Insert the data into the list, the first entry
                //in the pair is the equation number (row index),
                //the second is the value itself.
                list_pt->
                 insert(list_pt->end(),std::make_pair(eqn_number,value));
               }
             }
           }
         }
       }

#ifdef OOMPH_HAS_MPI
     } // endif halo element
#endif      


#ifdef OOMPH_HAS_MPI
  // Time it?
  if ((!Problem_has_been_distributed)&&
      Must_recompute_load_balance_for_assembly)
   {
    elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
   }
#endif

   } //End of loop over the elements
    
 } //list_pt goes out of scope


#ifdef OOMPH_HAS_MPI
 
 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

#endif


 //----Finally we need to convert the beautiful list storage scheme---
 //----------to the containers required by SuperLU--------------------

 //Loop over the number of matrices
 for(unsigned m=0;m<n_matrix;m++)
  {
   //Set the number of rows or columns
   row_or_column_start[m] = new int[ndof+1];
   //Counter for the total number of entries in the storage scheme
   unsigned long entry_count=0;
   //The first entry is 0
   row_or_column_start[m][0] = entry_count;

   // first we compute the number of non-zeros
   nnz[m] = 0;
   for(unsigned long i_global=0;i_global<ndof;i_global++)
    {
     nnz[m] += matrix_data_list[m][i_global].size();
    }

   // and then resize the storage
   column_or_row_index[m] = new int[nnz[m]];
   value[m] = new double[nnz[m]];
     
   //Now we merely loop over the number of rows or columns 
   for(unsigned long i_global=0;i_global<ndof;i_global++)
    {
     //Start index for the present row is the number of entries so far
     row_or_column_start[m][i_global] = entry_count;
     //If there are no entries in the list then skip the loop
     if(matrix_data_list[m][i_global].empty()) {continue;}

     //Sort the list corresponding to this row or column by the
     //column or row index (first entry in the pair). 
     //This might be inefficient, but we only have to do the sort ONCE
     //for each list. This is faster than using a map storage scheme, where
     //we are sorting for every insertion (although the map structure 
     //is cleaner and more memory efficient)
     matrix_data_list[m][i_global].sort();
       
     //Set up an iterator for start of the list
     std::list<std::pair<unsigned,double> >::iterator it
      = matrix_data_list[m][i_global].begin();
       
     //Get the first row or column index in the list...
     unsigned current_index = it->first;
     //...and the corresponding value
     double current_value = it->second;
       
     //Loop over all the entries in the sorted list 
     //Increase the iterator so that we start at the second entry
     for(++it;it!=matrix_data_list[m][i_global].end();++it)
      {
       //If the index has not changed, then we must add the contribution
       //of the present entry to the value. 
       //Additionally check that the entry is non-zero
       if((it->first == current_index) && 
          (std::abs(it->second) > Numerical_zero_for_sparse_assembly))
        {
         current_value += it->second;
        }
       //Otherwise, we have added all the contributions to the index
       //to current_value, so add it to the SuperLU data structure
       else
        {
         //Add the row or column index to the vector
         column_or_row_index[m][entry_count] = current_index;
         //Add the actual value to the vector
         value[m][entry_count] = current_value;
         //Increase the counter for the number of entries in each vector
         entry_count++;
           
         //Set the index and value to be those of the current entry in the 
         //list
         current_index = it->first;
         current_value = it->second;
        }
      } //End of loop over all list entries for this global row or column
       
     //There are TWO special cases to consider. 
     //If there is only one equation number in the list, then it
     //will NOT have been added. We test this case by comparing the
     //number of entries with those stored in row_or_column_start[i_global]
     //Otherwise
     //If the final entry in the list has the same index as the penultimate
     //entry, then it will NOT have been added to the SuperLU storage scheme
     //Check this by comparing the current_index with the final index
     //stored in the SuperLU scheme. If they are not the same, then
     //add the current_index and value.
       
     //If single equation number in list
     if((static_cast<int>(entry_count) == row_or_column_start[m][i_global])
        //If we have a single equation number, this will not be evaluated.
        //If we don't then we do the test to check that the final
        //entry is added
        ||(static_cast<int>(current_index) != 
           column_or_row_index[m][entry_count-1]))
      {
       //Add the row or column index to the vector
       column_or_row_index[m][entry_count] = current_index;
       //Add the actual value to the vector
       value[m][entry_count] = current_value;
       //Increase the counter for the number of entries in each vector
       entry_count++;
      }
       
    } //End of loop over the rows or columns of the entire matrix
     
   //Final entry in the row/column start vector
   row_or_column_start[m][ndof] = entry_count;
  } //End of loop over matrices

 if (Pause_at_end_of_sparse_assembly)
  {
   oomph_info << "Pausing at end of sparse assembly." << std::endl;
   pause("Check memory usage now.");
  }

}



//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors using vectors of pairs
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by chaging the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed_with_vectors_of_pairs(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value, 
 Vector<unsigned> &nnz,
 Vector<double*> &residuals,
 bool compressed_row_flag)
{
 //Total number of elements
 const unsigned long n_elements = mesh_pt()->nelement();

 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;
 
#ifdef OOMPH_HAS_MPI
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[Communicator_pt->my_rank()];
   el_hi=Last_el_for_assembly[Communicator_pt->my_rank()];
  }
#endif 

 // number of local eqns
 unsigned ndof = this->ndof();

 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();
 
 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_or_row_index.size();
 
 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;
 
//Error check dimensions
#ifdef PARANOID
 if(row_or_column_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_or_column_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_or_row_index.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_vectors_of_pairs",
    OOMPH_EXCEPTION_LOCATION);
  }

 if(value.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " 
    << std::endl
    << "value.size() " << value.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_or_row_index.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_vectors_of_pairs",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif


// The idea behind this sparse assembly routine is to use a Vector of
// Vectors of pairs for each complete matrix.
// Each inner Vector stores pairs and holds the row (or column) index
// and the value of the matrix entry.
 
// Set up Vector of Vectors to store the entries of each matrix,
// indexed by either the column or row.
 Vector<Vector< Vector<std::pair<unsigned,double> > > > matrix_data(n_matrix);
 
 //Loop over the number of matrices being assembled and resize
 //each Vector of Vectors to the number of rows or columns of the matrix
 for(unsigned m=0;m<n_matrix;m++) {matrix_data[m].resize(ndof);} 
 
 //Resize the residuals vectors
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals[v] = new double[ndof];
   for (unsigned i = 0; i < ndof; i++)
    {
     residuals[v][i] = 0;
    }
  }

#ifdef OOMPH_HAS_MPI

 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

#endif
 
 //----------------Assemble and populate the vector storage scheme--------
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element
  //This means that the storage is only allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);

  //Loop over the elements
  for(unsigned long e=el_lo;e<=el_hi;e++) 
   {

#ifdef OOMPH_HAS_MPI
    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }
#endif

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);
    
#ifdef OOMPH_HAS_MPI
    //Ignore halo elements
    if (!elem_pt->is_halo())
     {
#endif

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);

      //Resize the storage for elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}

      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);
    
      //---------------Insert the values into the vectors--------------
    
      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {

        //Get the local equation number
        unsigned eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);

        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in each residuals vector
          residuals[v][eqn_number] += el_residuals[v][i];
         }
      
        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);
        
          //Loop over the matrices
          //If it's compressed row storage, then our vector of maps
          //is indexed by row (equation number)
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only bother to add to the vector if it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              //If it's compressed row storage, then our vector of maps
              //is indexed by row (equation number)
              if(compressed_row_flag)
               {
                //Find the correct position and add the data into the vectors
                const unsigned size = matrix_data[m][eqn_number].size();
                for(unsigned k=0; k<=size; k++)
                 {
                  if(k==size)
                   {
                    matrix_data[m][eqn_number].push_back(
                     std::make_pair(unknown,value));
                    break;
                   }
                  else if(matrix_data[m][eqn_number][k].first == unknown)
                   {
                    matrix_data[m][eqn_number][k].second += value;
                    break;
                   }
                 }
               }
              //Otherwise it's compressed column storage and our vector is
              //indexed by column (the unknown)
              else
               {
                //Add the data into the vectors in the correct position
                const unsigned size = matrix_data[m][unknown].size();
                for(unsigned k=0; k<=size; k++)
                 {
                  if(k==size)
                   {
                    matrix_data[m][unknown].push_back(
                     std::make_pair(eqn_number,value));
                    break;
                   }
                  else if(matrix_data[m][unknown][k].first == eqn_number)
                   {
                    matrix_data[m][unknown][k].second += value;
                    break;
                   }
                 }
               }
             }
           } //End of loop over matrices
         }
       }

#ifdef OOMPH_HAS_MPI
     } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
  // Time it?
  if ((!Problem_has_been_distributed)&&
      Must_recompute_load_balance_for_assembly)
   {
    elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
   }
#endif
  
   } //End of loop over the elements
  
  
 } //End of vector assembly



#ifdef OOMPH_HAS_MPI
 
 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

#endif


 //-----------Finally we need to convert this vector storage scheme
 //------------------------to the containers required by SuperLU
 
 //Loop over the number of matrices
 for(unsigned m=0;m<n_matrix;m++)
  {
   //Set the number of rows or columns
   row_or_column_start[m]  = new int[ndof+1];
   
   // fill row_or_column_start and find the number of entries
   row_or_column_start[m][0] = 0;
   for(unsigned long i=0;i<ndof;i++)
    {
     row_or_column_start[m][i+1] = row_or_column_start[m][i]
      + matrix_data[m][i].size();
    }
   const unsigned entries = row_or_column_start[m][ndof];
   
   // resize vectors
   column_or_row_index[m] = new int[entries];
   value[m] = new double[entries];
   nnz[m] = entries;
   
   //Now we merely loop over the number of rows or columns
   for(unsigned long i_global=0;i_global<ndof;i_global++)
    {
     //If there are no entries in the vector then skip the rest of the loop
     if(matrix_data[m][i_global].empty()) {continue;}
     
     //Loop over all the entries in the vectors corresponding to the given
     //row or column. It will NOT be ordered
     unsigned p = 0;
     for(int j=row_or_column_start[m][i_global];
         j<row_or_column_start[m][i_global+1]; j++)
      {
       column_or_row_index[m][j] = matrix_data[m][i_global][p].first;
       value[m][j] = matrix_data[m][i_global][p].second;
       ++p;
      }
    }
  } //End of the loop over the matrices

 if (Pause_at_end_of_sparse_assembly)
  {
   oomph_info << "Pausing at end of sparse assembly." << std::endl;
   pause("Check memory usage now.");
  }
} 








//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors using two vectors.
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by chaging the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed_with_two_vectors(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value, 
 Vector<unsigned> &nnz,
 Vector<double* > &residuals,
 bool compressed_row_flag)
{
 //Total number of elements
 const unsigned long  n_elements = mesh_pt()->nelement();

 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;
 

#ifdef OOMPH_HAS_MPI
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[Communicator_pt->my_rank()];
   el_hi=Last_el_for_assembly[Communicator_pt->my_rank()];
  }
#endif 

 // number of local eqns
 unsigned ndof = this->ndof();

 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();

 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_or_row_index.size();

 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

//Error check dimensions
#ifdef PARANOID
 if(row_or_column_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_or_column_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_or_row_index.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_lists",
    OOMPH_EXCEPTION_LOCATION);
  }

 if(value.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " 
    << std::endl
    << "value.size() " << value.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_or_row_index.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_two_vectors",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

// The idea behind this sparse assembly routine is to use Vectors of
// Vectors for the entries in each complete matrix. And a second
// Vector of Vectors stores the global row (or column) indeces. This
// will not have the memory overheads associated with the methods using
// lists or maps, but insertion will be more costly.
 
// Set up two vector of vectors to store the entries of each  matrix,
// indexed by either the column or row. The entries of the vector for
// each matrix correspond to all the rows or columns of that matrix.
 Vector<Vector<Vector<unsigned> > > matrix_row_or_col_indices(n_matrix);
 Vector<Vector<Vector<double> > > matrix_values(n_matrix);
 
 //Loop over the number of matrices being assembled and resize
 //each vector of vectors to the number of rows or columns of the matrix
 for(unsigned m=0;m<n_matrix;m++)
  {
   matrix_row_or_col_indices[m].resize(ndof);
   matrix_values[m].resize(ndof);
  }
 
 //Resize the residuals vectors
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals[v] = new double[ndof];
   for (unsigned i = 0; i < ndof; i++)
    {
     residuals[v][i] = 0;
    }
  }

#ifdef OOMPH_HAS_MPI

 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

#endif

 
 //----------------Assemble and populate the vector storage scheme-------
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element
  //This means that the storage will only be allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);
  
  //Loop over the elements
  for(unsigned long e=el_lo;e<=el_hi;e++)
   {

#ifdef OOMPH_HAS_MPI
    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }
#endif

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
    //Ignore halo elements
    if (!elem_pt->is_halo())
     {
#endif

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);
    
      //Resize the storage for elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}

      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);
    
      //---------------Insert the values into the vectors--------------
    
      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {
        //Get the local equation number
        unsigned eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);
       
        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in each residuals vector
          residuals[v][eqn_number] += el_residuals[v][i];
         }
      
        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);
        
          //Loop over the matrices
          //If it's compressed row storage, then our vector of maps
          //is indexed by row (equation number)
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only bother to add to the vector if it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              //If it's compressed row storage, then our vector of maps
              //is indexed by row (equation number)
              if(compressed_row_flag)
               {
                //Find the correct position and add the data into the vectors
                const unsigned size = 
                 matrix_row_or_col_indices[m][eqn_number].size();
              
                for(unsigned k=0; k<=size; k++)
                 {
                  if(k==size)
                   {
                    matrix_row_or_col_indices[m][eqn_number].
                     push_back(unknown);
                    matrix_values[m][eqn_number].push_back(value);
                    break;
                   }
                  else if(matrix_row_or_col_indices[m][eqn_number][k] == 
                          unknown)
                   {
                    matrix_values[m][eqn_number][k] += value;
                    break;
                   }
                 }
               }
              //Otherwise it's compressed column storage and our vector is
              //indexed by column (the unknown)
              else
               {
                //Add the data into the vectors in the correct position
                const unsigned size = 
                 matrix_row_or_col_indices[m][unknown].size();
                for (unsigned k=0; k<=size; k++)
                 {
                  if (k==size)
                   {
                    matrix_row_or_col_indices[m][unknown].
                     push_back(eqn_number);
                    matrix_values[m][unknown].push_back(value);
                    break;
                   }
                  else if (matrix_row_or_col_indices[m][unknown][k] == 
                           eqn_number)
                   {
                    matrix_values[m][unknown][k] += value;
                    break;
                   }
                 }
               }
             }
           } //End of loop over matrices
         }
       }

#ifdef OOMPH_HAS_MPI
     } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
  // Time it?
  if ((!Problem_has_been_distributed)&&
      Must_recompute_load_balance_for_assembly)
   {
    elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
   }
#endif

   } //End of loop over the elements
  
 } //End of vector assembly



#ifdef OOMPH_HAS_MPI
 
 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

#endif
 
   //-----------Finally we need to convert this lousy vector storage scheme
   //------------------------to the containers required by SuperLU
 
   //Loop over the number of matrices
 for(unsigned m=0;m<n_matrix;m++)
  {
   //Set the number of rows or columns
   row_or_column_start[m] = new int[ndof+1];
   
   // fill row_or_column_start and find the number of entries
   row_or_column_start[m][0] = 0;
   for (unsigned long i=0;i<ndof;i++)
    {
     row_or_column_start[m][i+1] = row_or_column_start[m][i]
      + matrix_values[m][i].size();
    }
   const unsigned entries = row_or_column_start[m][ndof];
   
   // resize vectors
   column_or_row_index[m] = new int[entries];
   value[m] = new double[entries];
   nnz[m] = entries;
   
   //Now we merely loop over the number of rows or columns
   for(unsigned long  i_global=0;i_global<ndof;i_global++)
    {
     //If there are no entries in the vector then skip the rest of the loop
     if(matrix_values[m][i_global].empty()) {continue;}
     
     //Loop over all the entries in the vectors corresponding to the given
     //row or column. It will NOT be ordered
     unsigned p = 0;
     for(int j=row_or_column_start[m][i_global];
         j<row_or_column_start[m][i_global+1]; j++)
      {
       column_or_row_index[m][j] = matrix_row_or_col_indices[m][i_global][p];
       value[m][j] = matrix_values[m][i_global][p];
       ++p;
      }
    }
  } //End of the loop over the matrices

 if (Pause_at_end_of_sparse_assembly)
  {
   oomph_info << "Pausing at end of sparse assembly." << std::endl;
   pause("Check memory usage now.");
  }

}


//=====================================================================
/// This is a (private) helper function that is used to assemble system
/// matrices in compressed row or column format
/// and compute residual vectors using two vectors.
/// The default action is to assemble the jacobian matrix and 
/// residuals for the Newton method. The action can be
/// overloaded at an elemental level by chaging the default
/// behaviour of the function Element::get_all_vectors_and_matrices().
/// column_or_row_index: Column [or row] index of given entry
/// row_or_column_start: Index of first entry for given row [or column]
/// value              : Vector of nonzero entries
/// residuals          : Residual vector
/// compressed_row_flag: Bool flag to indicate if storage format is 
///                      compressed row [if false interpretation of
///                      arguments is as stated in square brackets].
//=====================================================================
void Problem::sparse_assemble_row_or_column_compressed_with_two_arrays(
 Vector<int* > &column_or_row_index, 
 Vector<int* > &row_or_column_start, 
 Vector<double* > &value, 
 Vector<unsigned> &nnz,
 Vector<double* > &residuals,
 bool compressed_row_flag)
{
 
 //Total number of elements
 const unsigned long  n_elements = mesh_pt()->nelement();

 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;
 

#ifdef OOMPH_HAS_MPI
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[Communicator_pt->my_rank()];
   el_hi=Last_el_for_assembly[Communicator_pt->my_rank()];
  }
#endif 

 // number of local eqns
 unsigned ndof = this->ndof();

 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();

 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_or_row_index.size();

 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

//Error check dimensions
#ifdef PARANOID
 if(row_or_column_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_or_column_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_or_row_index.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_lists",
    OOMPH_EXCEPTION_LOCATION);
  }

 if(value.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " 
    << std::endl
    << "value.size() " << value.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_or_row_index.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::sparse_assemble_row_or_column_compressed_with_two_vectors",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

// The idea behind this sparse assembly routine is to use Vectors of
// Vectors for the entries in each complete matrix. And a second
// Vector of Vectors stores the global row (or column) indeces. This
// will not have the memory overheads associated with the methods using
// lists or maps, but insertion will be more costly.
 
// Set up two vector of vectors to store the entries of each  matrix,
// indexed by either the column or row. The entries of the vector for
// each matrix correspond to all the rows or columns of that matrix.
 Vector<unsigned**> matrix_row_or_col_indices(n_matrix);
 Vector<double**> matrix_values(n_matrix);
 
 //Loop over the number of matrices being assembled and resize
 //each vector of vectors to the number of rows or columns of the matrix
 for(unsigned m=0;m<n_matrix;m++)
  {
   matrix_row_or_col_indices[m] = new unsigned*[ndof];
   matrix_values[m] = new double*[ndof];
  }
 
 //Resize the residuals vectors
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals[v] = new double[ndof];
   for (unsigned i = 0; i < ndof; i++)
    {
     residuals[v][i] = 0;
    }
  }

#ifdef OOMPH_HAS_MPI

 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

#endif

 // number of coefficients in each row
 Vector<Vector<unsigned> > ncoef(n_matrix);
 for (unsigned m = 0; m < n_matrix; m++)
  {
   ncoef[m].resize(ndof,0);
  }

 if (Sparse_assemble_with_arrays_previous_allocation.size() == 0)
  {
   Sparse_assemble_with_arrays_previous_allocation.resize(n_matrix);
   for (unsigned m = 0; m < n_matrix; m++)
    {
     Sparse_assemble_with_arrays_previous_allocation[m].resize(ndof,0);
    }
  }

 //----------------Assemble and populate the vector storage scheme-------
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element
  //This means that the storage will only be allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);
  
  //Loop over the elements
  for(unsigned long e=el_lo;e<=el_hi;e++)
   {

#ifdef OOMPH_HAS_MPI
    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }
#endif

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

#ifdef OOMPH_HAS_MPI
    //Ignore halo elements
    if (!elem_pt->is_halo())
     {
#endif

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);
    
      //Resize the storage for elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}

      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);
    
      //---------------Insert the values into the vectors--------------
    
      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {
      //Get the local equation number
        unsigned eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);
     
        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in each residuals vector
          residuals[v][eqn_number] += el_residuals[v][i];
         }
      
        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);
        
          //Loop over the matrices
          //If it's compressed row storage, then our vector of maps
          //is indexed by row (equation number)
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only bother to add to the vector if it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              // number of entrys in this row
              const unsigned size = ncoef[m][eqn_number];
              
              // if no data has been allocated for this row then allocate
              if (size == 0)
               {
                // do we have previous allocation data
                if (Sparse_assemble_with_arrays_previous_allocation
                    [m][eqn_number] != 0)
                 {
                  matrix_row_or_col_indices[m][eqn_number] = 
                   new unsigned
                   [Sparse_assemble_with_arrays_previous_allocation[m]
                    [eqn_number]];
                  matrix_values[m][eqn_number] = 
                   new double
                   [Sparse_assemble_with_arrays_previous_allocation[m]
                    [eqn_number]];
                 }
                else
                 {
                  matrix_row_or_col_indices[m][eqn_number] = 
                   new unsigned
                   [Sparse_assemble_with_arrays_initial_allocation];
                  matrix_values[m][eqn_number] = 
                   new double
                   [Sparse_assemble_with_arrays_initial_allocation];
                  Sparse_assemble_with_arrays_previous_allocation[m]
                   [eqn_number]=
                   Sparse_assemble_with_arrays_initial_allocation;
                 }
               }
              
              //If it's compressed row storage, then our vector of maps
              //is indexed by row (equation number)
              if(compressed_row_flag)
               {
                // next add the data
                for(unsigned k=0; k<=size; k++)
                 {
                  if(k==size)
                   {
                    // do we need to allocate more storage
                    if (Sparse_assemble_with_arrays_previous_allocation
                        [m][eqn_number] == ncoef[m][eqn_number])
                     {
                      unsigned new_allocation = ncoef[m][eqn_number]+
                       Sparse_assemble_with_arrays_allocation_increment;
                      double* new_values = new double[new_allocation];
                      unsigned* new_indices = new unsigned[new_allocation];
                      for (unsigned c = 0; c < ncoef[m][eqn_number]; c++)
                       {
                        new_values[c] = matrix_values[m][eqn_number][c];
                        new_indices[c] = 
                         matrix_row_or_col_indices[m][eqn_number][c];
                       }
                      delete[] matrix_values[m][eqn_number];
                      delete[] matrix_row_or_col_indices[m][eqn_number];
                      matrix_values[m][eqn_number]=new_values;
                      matrix_row_or_col_indices[m][eqn_number]=new_indices;
                      Sparse_assemble_with_arrays_previous_allocation
                       [m][eqn_number] = new_allocation;
                     }
                    // and now add the data
                    unsigned entry = ncoef[m][eqn_number];
                    ncoef[m][eqn_number]++;
                    matrix_row_or_col_indices[m][eqn_number][entry] = unknown;
                    matrix_values[m][eqn_number][entry] = value;
                    break;
                   }
                  else if(matrix_row_or_col_indices[m][eqn_number][k] == 
                          unknown)
                   {
                    matrix_values[m][eqn_number][k] += value;
                    break;
                   }
                 }
               }
              //Otherwise it's compressed column storage and our vector is
              //indexed by column (the unknown)
              else
               {
                //Add the data into the vectors in the correct position
                for(unsigned k=0; k<=size; k++)
                 {
                  if(k==size)
                   {
                    // do we need to allocate more storage
                    if (Sparse_assemble_with_arrays_previous_allocation
                        [m][unknown] == ncoef[m][unknown])
                     {
                      unsigned new_allocation = ncoef[m][unknown]+
                       Sparse_assemble_with_arrays_allocation_increment;
                      double* new_values = new double[new_allocation];
                      unsigned* new_indices = new unsigned[new_allocation];
                      for (unsigned c = 0; c < ncoef[m][unknown]; c++)
                       {
                        new_values[c] = matrix_values[m][unknown][c];
                        new_indices[c] = 
                         matrix_row_or_col_indices[m][unknown][c];
                       }
                      delete[] matrix_values[m][unknown];
                      delete[] matrix_row_or_col_indices[m][unknown];
                      Sparse_assemble_with_arrays_previous_allocation
                       [m][unknown] = new_allocation;
                     }
                    // and now add the data
                    unsigned entry = ncoef[m][unknown];
                    ncoef[m][unknown]++;
                    matrix_row_or_col_indices[m][unknown][entry] = eqn_number;
                    matrix_values[m][unknown][entry] = value;
                    break;
                   }
                  else if(matrix_row_or_col_indices[m][unknown][k] == 
                          eqn_number)
                   {
                    matrix_values[m][unknown][k] += value;
                    break;
                   }
                 }
               }
             }
           } //End of loop over matrices
         }
       }

#ifdef OOMPH_HAS_MPI
     } // endif halo element
#endif


#ifdef OOMPH_HAS_MPI
  // Time it?
  if ((!Problem_has_been_distributed)&&
      Must_recompute_load_balance_for_assembly)
   {
    elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
   }
#endif

   } //End of loop over the elements
  
 } //End of vector assembly



#ifdef OOMPH_HAS_MPI
 
 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

#endif
 
   //-----------Finally we need to convert this lousy vector storage scheme
   //------------------------to the containers required by SuperLU
 
   //Loop over the number of matrices
 for(unsigned m=0;m<n_matrix;m++)
  {
   //Set the number of rows or columns
   row_or_column_start[m] = new int[ndof+1];
   
   // fill row_or_column_start and find the number of entries
   row_or_column_start[m][0] = 0;
   for (unsigned long i=0;i<ndof;i++)
    {
     row_or_column_start[m][i+1] = row_or_column_start[m][i]
      + ncoef[m][i];
     Sparse_assemble_with_arrays_previous_allocation[m][i] = ncoef[m][i];
    }
   const unsigned entries = row_or_column_start[m][ndof];
   
   // resize vectors
   column_or_row_index[m] = new int[entries];
   value[m] = new double[entries];
   nnz[m] = entries;
   
   //Now we merely loop over the number of rows or columns
   for(unsigned long  i_global=0;i_global<ndof;i_global++)
    {
     //If there are no entries in the vector then skip the rest of the loop
     if(ncoef[m][i_global]==0) {continue;}
     
     //Loop over all the entries in the vectors corresponding to the given
     //row or column. It will NOT be ordered
     unsigned p = 0;
     for(int j=row_or_column_start[m][i_global];
         j<row_or_column_start[m][i_global+1]; j++)
      {
       column_or_row_index[m][j] = matrix_row_or_col_indices[m][i_global][p];
       value[m][j] = matrix_values[m][i_global][p];
       ++p;
      }

     // and delete
     delete[] matrix_row_or_col_indices[m][i_global];
     delete[] matrix_values[m][i_global];
    }

   //
   delete[] matrix_row_or_col_indices[m];
   delete[] matrix_values[m];
  } //End of the loop over the matrices

 if (Pause_at_end_of_sparse_assembly)
  {
   oomph_info << "Pausing at end of sparse assembly." << std::endl;
   pause("Check memory usage now.");
  }
}


#ifdef OOMPH_HAS_MPI
//=============================================================================
/// \short Helper method to assemble CRDoubleMatrices from distributed
/// on multiple processors.
//=============================================================================
void Problem::parallel_sparse_assemble
(const LinearAlgebraDistribution* const& target_dist_pt,
 Vector<int* > &column_indices, 
 Vector<int* > &row_start, 
 Vector<double* > &values,
 Vector<unsigned > &nnz, 
 Vector<double* > &residuals)
{

 // my rank and nproc
 unsigned my_rank = Communicator_pt->my_rank();
 unsigned nproc = Communicator_pt->nproc();
 
 // start by assembling the sorted set of equations which this processor
 // contributes to
 //======================================================================
 Vector<unsigned> my_eqns;
 unsigned dof_first_row = Dof_distribution_pt->first_row();
 unsigned dof_nrow_local = Dof_distribution_pt->nrow_local();
 unsigned dof_last_row = dof_first_row + dof_nrow_local;

 // first the local eqn numbers
 for (unsigned i = dof_first_row; i < dof_last_row; i++)
  {
   my_eqns.push_back(i);
  }

 // Each processor assembles equations from each mesh... I have no idea
 // how this worked for any distributed FSI problem... assuming it did... ?
 // This current method won't pick up any halo equation numbers at all as the
 // global mesh doesn't store any halo objects at all; they are all stored
 // by the relevant submesh

 unsigned nmesh=nsub_mesh();

 // If there are no submeshes then the current method works fine
 if (nmesh==0)
  {
   // hierher: I suspect that for a single-mesh problem there's no need
   //          to look for any "external" eqn numbers (whether halo or not)

//    // hang on, there are "local" external eqn numbers too
//    unsigned n_ext_node=mesh_pt()->nexternal_node();
//    for (unsigned j=0;j<n_ext_node;j++)
//     {
//      Node* ext_nod_pt=mesh_pt()->external_node_pt(j);
//      unsigned nval=ext_nod_pt->nvalue();
//      for (unsigned ival=0;ival<nval;ival++)
//       {
//        int eqn_num = ext_nod_pt->eqn_number(ival);
//        if (eqn_num>=0)
//         {
//          my_eqns.push_back(eqn_num);
//         }
//       }
//     }

   // then the halo equation numbers
   for (unsigned p = 0; p < nproc; p++)
    {
     if (p != my_rank)
      {       
       unsigned n_halo_node=mesh_pt()->nhalo_node(p);
       for (unsigned j=0;j<n_halo_node;j++)
        {
         Node* halo_node_pt=mesh_pt()->halo_node_pt(p,j);
         unsigned nval=halo_node_pt->nvalue();
         for (unsigned ival=0;ival<nval;ival++)
          {
           int eqn_num = halo_node_pt->eqn_number(ival);
           if (eqn_num>=0)
            {
             my_eqns.push_back(eqn_num);
            }
          }
         SolidNode* solid_node_pt=dynamic_cast<SolidNode*>(halo_node_pt);
         if (solid_node_pt!=0)
          {
           unsigned nval=solid_node_pt->variable_position_pt()->nvalue();
           for (unsigned ival=0; ival<nval; ival++)
            {
             int eqn_num = 
              solid_node_pt->variable_position_pt()->eqn_number(ival);
             if (eqn_num>=0)
              {
               my_eqns.push_back(eqn_num);
              }
            }
          }
        }
      }
    }

//    // and finally the external halo equation numbers
//    for (unsigned p = 0; p < nproc; p++)
//     {
//      if (p != my_rank)
//       {       
//        unsigned n_ext_halo_node=mesh_pt()->nexternal_halo_node(p);
//        for (unsigned j=0;j<n_ext_halo_node;j++)
//         {
//          Node* ext_halo_node_pt = mesh_pt()->external_halo_node_pt(p,j);
//          unsigned nval = ext_halo_node_pt->nvalue();
//          for (unsigned ival=0;ival<nval;ival++)
//           {
//            int eqn_num = ext_halo_node_pt->eqn_number(ival);
//            if (eqn_num>=0)
//             {
//              my_eqns.push_back(eqn_num);
//             }
//           }
//          SolidNode* solid_node_pt=dynamic_cast<SolidNode*>(ext_halo_node_pt);
//          if (solid_node_pt!=0)
//           {
//            unsigned nval=solid_node_pt->variable_position_pt()->nvalue();
//            for (unsigned ival=0; ival<nval; ival++)
//             {
//              int eqn_num = 
//               solid_node_pt->variable_position_pt()->eqn_number(ival);
//              if (eqn_num>=0)
//               {
//                my_eqns.push_back(eqn_num);
//               }
//             }
//           }
//         }
//       }
//     }  

//    // finally loop over the internal data in external halo elements
//    for (unsigned p = 0; p < nproc; p++)
//     {
//      if (p!=my_rank)
//       {
//        unsigned n_ext_halo_element = mesh_pt()->nexternal_halo_element(p);
//        for (unsigned j = 0; j < n_ext_halo_element; j++)
//         {
//          FiniteElement* ext_halo_element_pt = 
//           mesh_pt()->external_halo_element_pt(p,j);
//          unsigned ndata = ext_halo_element_pt->ninternal_data();
//          for (unsigned i = 0; i < ndata; i++)
//           {
//            unsigned nvalue = 
//             ext_halo_element_pt->internal_data_pt(i)->nvalue();
//            for (unsigned k = 0; k < nvalue; k++)
//             {
//              int eqn_num 
//               = ext_halo_element_pt->internal_data_pt(i)->eqn_number(k);
//              if (eqn_num>=0)
//               {
//                my_eqns.push_back(eqn_num);
//               }
//             }
//           }
//         }
//       }
//     }

//    // finally loop over the internal data in halo elements
//    for (unsigned p = 0; p < nproc; p++)
//     {
//      if (p!=my_rank)
//       {
//        Vector<FiniteElement*> halo_element_pt = mesh_pt()->halo_element_pt(p);
//        unsigned n_halo_element = halo_element_pt.size();
//        for (unsigned j = 0; j < n_halo_element; j++)
//         {
//          unsigned ndata = halo_element_pt[j]->ninternal_data();
//          for (unsigned i = 0; i < ndata; i++)
//           {
//            unsigned nvalue = halo_element_pt[j]->internal_data_pt(i)->nvalue();
//            for (unsigned k = 0; k < nvalue; k++)
//             {
//              int eqn_num 
//               = halo_element_pt[j]->internal_data_pt(i)->eqn_number(k);
//              if (eqn_num>=0)
//               {
//                my_eqns.push_back(eqn_num);
//               }
//             }
//           }
//         }
//       }
//     }
  }
 else // there are submeshes
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<nmesh;imesh++)
    {
     // There are "local" external eqn numbers
     unsigned n_ext_node=mesh_pt(imesh)->nexternal_node();
     for (unsigned j=0;j<n_ext_node;j++)
      {
       Node* ext_nod_pt=mesh_pt(imesh)->external_node_pt(j);
       unsigned nval=ext_nod_pt->nvalue();
       for (unsigned ival=0;ival<nval;ival++)
        {
         int eqn_num = ext_nod_pt->eqn_number(ival);
         if (eqn_num>=0)
          {
           my_eqns.push_back(eqn_num);
          }
        }
      }

     // then the halo equation numbers
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p != my_rank)
        {       
         unsigned n_halo_node=mesh_pt(imesh)->nhalo_node(p);
         for (unsigned j=0;j<n_halo_node;j++)
          {
           Node* halo_node_pt=mesh_pt(imesh)->halo_node_pt(p,j);
           unsigned nval=halo_node_pt->nvalue();
           for (unsigned ival=0;ival<nval;ival++)
            {
             int eqn_num = halo_node_pt->eqn_number(ival);
             if (eqn_num>=0)
              {
               my_eqns.push_back(eqn_num);
              }
            }
           SolidNode* solid_node_pt=dynamic_cast<SolidNode*>(halo_node_pt);
           if (solid_node_pt!=0)
            {
             unsigned nval=solid_node_pt->variable_position_pt()->nvalue();
             for (unsigned ival=0; ival<nval; ival++)
              {
               int eqn_num = 
                solid_node_pt->variable_position_pt()->eqn_number(ival);
               if (eqn_num>=0)
                {
                 my_eqns.push_back(eqn_num);
                }
              }
            }
          }
        }
      }

     // and finally the external halo equation numbers
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p != my_rank)
        {       
         unsigned n_ext_halo_node=mesh_pt(imesh)->nexternal_halo_node(p);
         for (unsigned j=0;j<n_ext_halo_node;j++)
          {
           Node* ext_halo_node_pt = mesh_pt(imesh)->external_halo_node_pt(p,j);
           unsigned nval = ext_halo_node_pt->nvalue();
           for (unsigned ival=0;ival<nval;ival++)
            {
             int eqn_num = ext_halo_node_pt->eqn_number(ival);
             if (eqn_num>=0)
              {
               my_eqns.push_back(eqn_num);
              }
            }
           SolidNode* solid_node_pt=dynamic_cast<SolidNode*>(ext_halo_node_pt);
           if (solid_node_pt!=0)
            {
             unsigned nval=solid_node_pt->variable_position_pt()->nvalue();
             for (unsigned ival=0; ival<nval; ival++)
              {
               int eqn_num = 
                solid_node_pt->variable_position_pt()->eqn_number(ival);
               if (eqn_num>=0)
                {
                 my_eqns.push_back(eqn_num);
                }
              }
            }
          }
        }
      }  

     // finally loop over the internal data in external halo elements
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p!=my_rank)
        {
         unsigned n_ext_halo_element = 
          mesh_pt(imesh)->nexternal_halo_element(p);
         for (unsigned j = 0; j < n_ext_halo_element; j++)
          {
           FiniteElement* ext_halo_element_pt = 
            mesh_pt(imesh)->external_halo_element_pt(p,j);
           unsigned ndata = ext_halo_element_pt->ninternal_data();
           for (unsigned i = 0; i < ndata; i++)
            {
             unsigned nvalue = 
              ext_halo_element_pt->internal_data_pt(i)->nvalue();
             for (unsigned k = 0; k < nvalue; k++)
              {
               int eqn_num 
                = ext_halo_element_pt->internal_data_pt(i)->eqn_number(k);
               if (eqn_num>=0)
                {
                 my_eqns.push_back(eqn_num);
                }
              }
            }
          }
        }
      }

     // finally loop over the internal data in halo elements
     for (unsigned p = 0; p < nproc; p++)
      {
       if (p!=my_rank)
        {
         Vector<FiniteElement*> halo_element_pt = 
          mesh_pt(imesh)->halo_element_pt(p);
         unsigned n_halo_element = halo_element_pt.size();
         for (unsigned j = 0; j < n_halo_element; j++)
          {
           unsigned ndata = halo_element_pt[j]->ninternal_data();
           for (unsigned i = 0; i < ndata; i++)
            {
             unsigned nvalue = 
              halo_element_pt[j]->internal_data_pt(i)->nvalue();
             for (unsigned k = 0; k < nvalue; k++)
              {
               int eqn_num 
                = halo_element_pt[j]->internal_data_pt(i)->eqn_number(k);
               if (eqn_num>=0)
                {
                 my_eqns.push_back(eqn_num);
                }
              }
            }
          }
        }
      }
    }
  }

 // and sort and remove duplicate entries
 std::sort(my_eqns.begin(),my_eqns.end());
 Vector<unsigned>::iterator it = std::unique(my_eqns.begin(),my_eqns.end());
 my_eqns.resize(it-my_eqns.begin());

 // number of equations
 unsigned my_n_eqn = my_eqns.size();

 // next we assemble the data into an array of arrays
 // =================================================

 //Total number of elements
 const unsigned long  n_elements = mesh_pt()->nelement();
   
 // Default range of elements for distributed problems
 unsigned long el_lo=0;
 unsigned long el_hi=n_elements-1;
   
 // Otherwise just loop over a fraction of the elements
 // (This will either have been initialised in
 // Problem::set_default_first_and_last_element_for_assembly() or
 // will have been re-assigned during a previous assembly loop
 // Note that following the re-assignment only the entries
 // for the current processor are relevant.
 if (!Problem_has_been_distributed)
  {
   el_lo=First_el_for_assembly[my_rank];
   el_hi=Last_el_for_assembly[my_rank];
  }
      
 //Find the number of vectors to be assembled
 const unsigned n_vector = residuals.size();
   
 //Find the number of matrices to be assembled
 const unsigned n_matrix = column_indices.size();
   
 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;
   
//Error check dimensions
#ifdef PARANOID
 if(row_start.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " << std::endl
    << "row_or_column_start.size() " << row_start.size() 
    << " does not equal "
    << "column_or_row_index.size() " 
    <<  column_indices.size() << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::parallel_sparse_assemble()",
    OOMPH_EXCEPTION_LOCATION);
  }
   
 if(values.size() != n_matrix)
  {
   std::ostringstream error_stream;
   error_stream
    << "Error: " 
    << std::endl
    << "value.size() " << values.size() << " does not equal "
    << "column_or_row_index.size() " 
    << column_indices.size() << std::endl<< std::endl
    << std::endl;
   throw OomphLibError(
    error_stream.str(),
    "Problem::parallel_sparse_assemble()",
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

// The idea behind this sparse assembly routine is to use an array of
// arrays for the entries in each complete matrix. And a second
// array of arrays stores the global row (or column) indeces.
 
// Set up two vector of vectors to store the entries of each  matrix,
// indexed by either the column or row. The entries of the vector for
// each matrix correspond to all the rows or columns of that matrix.
 Vector<unsigned**> matrix_col_indices(n_matrix);
 Vector<double**> matrix_values(n_matrix);
 
 //Loop over the number of matrices being assembled and resize
 //each vector of vectors to the number of rows or columns of the matrix
 for(unsigned m=0;m<n_matrix;m++)
  {
   matrix_col_indices[m] = new unsigned*[my_n_eqn];
   matrix_values[m] = new double*[my_n_eqn];
   for (unsigned i = 0; i < my_n_eqn; i++)
    {
     matrix_col_indices[m][i]=0;
     matrix_values[m][i]=0;
    }
  }
 
 //Resize the residuals vectors
 Vector<double*> residuals_data(n_vector);
 for(unsigned v=0;v<n_vector;v++) 
  { 
   residuals_data[v] = new double[my_n_eqn];
   for (unsigned i = 0; i < my_n_eqn; i++)
    {
     residuals_data[v][i] = 0;
    }
  }

 // Storage for assembly time for elements
 double t_assemble_start=0.0;
 
 // Storage for assembly times
 Vector<double> elemental_assembly_time;
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   elemental_assembly_time.resize(n_elements);
  }

 // number of coefficients in each row
 Vector<Vector<unsigned> > ncoef(n_matrix);
 for (unsigned m = 0; m < n_matrix; m++)
  {
   ncoef[m].resize(my_n_eqn,0);
  }

 // Sparse_assemble_with_arrays_previous_allocation stores the number of
 // coefs in each row.
 // if a matrix of this size has not been assembled before then resize this
 // storage 
 if (Sparse_assemble_with_arrays_previous_allocation.size() == 0)
  {
   Sparse_assemble_with_arrays_previous_allocation.resize(n_matrix);
   for (unsigned m = 0; m < n_matrix; m++)
    {
     Sparse_assemble_with_arrays_previous_allocation[m].resize(my_n_eqn,0);
    }
  }

 // assemble and populate an array based storage scheme
 {
  //Allocate local storage for the element's contribution to the
  //residuals vectors and system matrices of the size of the maximum
  //number of dofs in any element
  //This means that the storage will only be allocated (and deleted) once
  Vector<Vector<double> > el_residuals(n_vector);
  Vector<DenseMatrix<double> > el_jacobian(n_matrix);
  
  //Loop over the elements
  for(unsigned long e=el_lo;e<=el_hi;e++)
   {

    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      t_assemble_start=TimingHelpers::timer();
     }

    //Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

    //Ignore halo elements
    if (!elem_pt->is_halo())
     {

      //Find number of degrees of freedom in the element
      const unsigned nvar = assembly_handler_pt->ndof(elem_pt);
    
      //Resize the storage for elemental jacobian and residuals
      for(unsigned v=0;v<n_vector;v++) {el_residuals[v].resize(nvar);}
      for(unsigned m=0;m<n_matrix;m++) {el_jacobian[m].resize(nvar);}

      //Now get the residuals and jacobian for the element
      assembly_handler_pt->
       get_all_vectors_and_matrices(elem_pt,el_residuals, el_jacobian);
    
      //---------------Insert the values into the vectors--------------
    
      //Loop over the first index of local variables
      for(unsigned i=0;i<nvar;i++)
       {
        //Get the local equation number
        unsigned global_eqn_number 
         = assembly_handler_pt->eqn_number(elem_pt,i);
        
        // determine the element number in my set of eqns using the 
        // bisection method
        int left = 0;
        int right = my_n_eqn-1;
        int eqn_number = right/2;
        while (my_eqns[eqn_number] != global_eqn_number)
         {
          if (left == right)
           {
            std::ostringstream error_stream;
            error_stream
             << "Internal Error: " 
             << std::endl
             << "Could not find global equation number "
             << global_eqn_number << "in my_eqns vector of equation numbers.";
            throw OomphLibError(
             error_stream.str(),
             "Problem::parallel_sparse_assemble()",
             OOMPH_EXCEPTION_LOCATION);
           }
          if (my_eqns[eqn_number] > global_eqn_number)
           {
            right = std::max(eqn_number-1,left);
           }
          else
           {
            left = std::min(eqn_number+1,right);
           }
          eqn_number = (right+left)/2;
         }

        //Add the contribution to the residuals
        for(unsigned v=0;v<n_vector;v++)
         {
          //Fill in each residuals vector
          residuals_data[v][eqn_number] += el_residuals[v][i];
         }
      
        //Now loop over the other index
        for(unsigned j=0;j<nvar;j++)
         {
          //Get the number of the unknown
          unsigned unknown = assembly_handler_pt->eqn_number(elem_pt,j);
        
          //Loop over the matrices
          //If it's compressed row storage, then our vector of maps
          //is indexed by row (equation number)
          for(unsigned m=0;m<n_matrix;m++)
           {
            //Get the value of the matrix at this point
            double value = el_jacobian[m](i,j);
            //Only bother to add to the vector if it's non-zero
            if(std::abs(value) > Numerical_zero_for_sparse_assembly)
             {
              // number of entrys in this row
              const unsigned size = ncoef[m][eqn_number];
              
              // if no data has been allocated for this row then allocate
              if (size == 0)
               {
                // do we have previous allocation data
                if (Sparse_assemble_with_arrays_previous_allocation
                    [m][eqn_number] != 0)
                 {
                  matrix_col_indices[m][eqn_number] = 
                   new unsigned
                   [Sparse_assemble_with_arrays_previous_allocation[m]
                    [eqn_number]];
                  matrix_values[m][eqn_number] = 
                   new double
                   [Sparse_assemble_with_arrays_previous_allocation[m]
                    [eqn_number]];
                 }
                else
                 {
                  matrix_col_indices[m][eqn_number] = 
                   new unsigned
                   [Sparse_assemble_with_arrays_initial_allocation];
                  matrix_values[m][eqn_number] = 
                   new double
                   [Sparse_assemble_with_arrays_initial_allocation];
                  Sparse_assemble_with_arrays_previous_allocation[m]
                   [eqn_number]=
                   Sparse_assemble_with_arrays_initial_allocation;
                 }
               }

              // next add the data
              for(unsigned k=0; k<=size; k++)
               {
                if(k==size)
                 {
                  // do we need to allocate more storage
                  if (Sparse_assemble_with_arrays_previous_allocation
                      [m][eqn_number] == ncoef[m][eqn_number])
                   {
                    unsigned new_allocation = ncoef[m][eqn_number]+
                     Sparse_assemble_with_arrays_allocation_increment;
                    double* new_values = new double[new_allocation];
                    unsigned* new_indices = new unsigned[new_allocation];
                    for (unsigned c = 0; c < ncoef[m][eqn_number]; c++)
                     {
                      new_values[c] = matrix_values[m][eqn_number][c];
                      new_indices[c] = matrix_col_indices[m][eqn_number][c];
                     }
                    delete[] matrix_values[m][eqn_number];
                    delete[] matrix_col_indices[m][eqn_number];
                    matrix_values[m][eqn_number] = new_values;
                    matrix_col_indices[m][eqn_number] = new_indices;
                    Sparse_assemble_with_arrays_previous_allocation
                     [m][eqn_number] = new_allocation;
                   }
                  // and now add the data
                  unsigned entry = ncoef[m][eqn_number];
                  ncoef[m][eqn_number]++;
                  matrix_col_indices[m][eqn_number][entry] = unknown;
                  matrix_values[m][eqn_number][entry] = value;
                  break;
                 }
                else if(matrix_col_indices[m][eqn_number][k] == unknown)
                 {
                  matrix_values[m][eqn_number][k] += value;
                  break;
                 }
               }
             }//numerical zero check
           } //End of loop over matrices
         }
       }
     } // endif halo element

    // Time it?
    if ((!Problem_has_been_distributed)&&
        Must_recompute_load_balance_for_assembly)
     {
      elemental_assembly_time[e]=TimingHelpers::timer()-t_assemble_start;
     }
   } //End of loop over the elements
 } //End of vector assembly

 // Postprocess timing information and re-allocate distribution of
 // elements during subsequent assemblies.
 if ((!Problem_has_been_distributed)&&
     Must_recompute_load_balance_for_assembly)
  {
   recompute_load_balanced_assembly(elemental_assembly_time);
  }

 // next we compute the number of equations and number of non-zeros to be
 // sent to each processor, and send/recv that information
 // =====================================================================
 
 // determine the number of eqns to be sent to each processor
 Vector<unsigned> n_eqn_for_proc(nproc,0);
 Vector<unsigned> first_eqn_element_for_proc(nproc);
 unsigned current_p = target_dist_pt->rank_of_global_row(my_eqns[0]);
 first_eqn_element_for_proc[current_p] = 0;
 n_eqn_for_proc[current_p] = 1;
 for (unsigned i = 1; i < my_n_eqn; i++)
  {
   unsigned next_p = target_dist_pt->rank_of_global_row(my_eqns[i]);
   if (next_p != current_p)
    {
     current_p = next_p;
     first_eqn_element_for_proc[current_p] = i;
    }
   n_eqn_for_proc[current_p]++;
  }

 // determine the number of non-zeros to be sent to each processor for each 
 // matrix
 DenseMatrix<unsigned> nnz_for_proc(nproc,n_matrix,0);
 for (unsigned p = 0; p < nproc; p++)
  {
   int first_eqn_element = first_eqn_element_for_proc[p];
   int last_eqn_element = (int)(first_eqn_element + n_eqn_for_proc[p]) - 1;
   for (unsigned m = 0; m < n_matrix; m++)
    {
     for (int i = first_eqn_element; i <= last_eqn_element; i++)
      {
       nnz_for_proc(p,m) += ncoef[m][i];
      }
    }
  }

 // next post the sends and recvs to the corresponding processors
 Vector<unsigned*> temp_send_storage(nproc);
 Vector<unsigned*> temp_recv_storage(nproc);
 Vector<MPI_Request> send_nnz_reqs;
 Vector<MPI_Request> recv_nnz_reqs;
 for (unsigned p = 0; p < nproc; p++)
  {
   if (p != my_rank)
    {
   temp_send_storage[p] = new unsigned[n_matrix+1];
   temp_send_storage[p][0] = n_eqn_for_proc[p];
   for (unsigned m = 0; m < n_matrix; m++)
    {
     temp_send_storage[p][m+1] = nnz_for_proc(p,m);
    }
   MPI_Request sreq;
   MPI_Isend(temp_send_storage[p],n_matrix+1,MPI_UNSIGNED,p,0,
             Communicator_pt->mpi_comm(),&sreq);
   send_nnz_reqs.push_back(sreq);
   temp_recv_storage[p] = new unsigned[n_matrix+1];
   MPI_Request rreq;
   MPI_Irecv(temp_recv_storage[p],n_matrix+1,MPI_UNSIGNED,p,0,
             Communicator_pt->mpi_comm(),&rreq);
   recv_nnz_reqs.push_back(rreq);
  }
  }

 // assemble the data to be sent to each processor
 // ==============================================

 // storage
 Vector<unsigned*> eqns_for_proc(nproc);
 DenseMatrix<double*> residuals_for_proc(nproc,n_vector);
 DenseMatrix<unsigned*> row_start_for_proc(nproc,n_matrix);
 DenseMatrix<unsigned*> column_indices_for_proc(nproc,n_matrix);
 DenseMatrix<double*> values_for_proc(nproc,n_matrix);
   
 // equation numbers
 for (unsigned p = 0; p < nproc; p++)
  {
   unsigned first_eqn_element = first_eqn_element_for_proc[p];
   unsigned first_row = target_dist_pt->first_row(p);
   unsigned n_eqns_p = n_eqn_for_proc[p];
   if (n_eqns_p > 0)
    {
     eqns_for_proc[p] = new unsigned[n_eqns_p];
     for (unsigned i = 0; i < n_eqns_p; i++)
      {
       eqns_for_proc[p][i] = my_eqns[i+first_eqn_element]-first_row;
      }
    }
  }

 // residuals for p
 for (unsigned v = 0; v < n_vector; v++)
  {
   for (unsigned p = 0; p < nproc; p++)
    {
     unsigned n_eqns_p = n_eqn_for_proc[p];
     if (n_eqns_p > 0)
      {
//       oomph_info << "allocating: " << v << " " << p << std::endl;
       unsigned first_eqn_element = first_eqn_element_for_proc[p];
       residuals_for_proc(p,v) = new double[n_eqns_p];
       for (unsigned i = 0; i < n_eqns_p; i++)
        {
         residuals_for_proc(p,v)[i] = residuals_data[v][first_eqn_element+i];
        }
      }
    }
   delete[] residuals_data[v];
  }

 // matrices for p
 for (unsigned m = 0; m < n_matrix; m++)
  {
   for (unsigned p = 0; p < nproc; p++)
    {
     unsigned n_eqns_p = n_eqn_for_proc[p];
     if (n_eqns_p > 0)
      {
       unsigned first_eqn_element = first_eqn_element_for_proc[p];     
       row_start_for_proc(p,m) = new unsigned[n_eqns_p+1];
       column_indices_for_proc(p,m) = new unsigned[nnz_for_proc(p,m)];
       values_for_proc(p,m) = new double[nnz_for_proc(p,m)];
       unsigned entry = 0;
       for (unsigned i = 0; i < n_eqns_p; i++)
        {
         row_start_for_proc(p,m)[i] = entry;
         unsigned n_coef_in_row = ncoef[m][first_eqn_element+i];
         for (unsigned j = 0; j < n_coef_in_row; j++)
          {
           column_indices_for_proc(p,m)[entry] 
            = matrix_col_indices[m][i+first_eqn_element][j];
           values_for_proc(p,m)[entry] 
            = matrix_values[m][i+first_eqn_element][j];
           entry++;
          }
        }
       row_start_for_proc(p,m)[n_eqns_p]=entry;
      }
    }
   for (unsigned i = 0; i < my_n_eqn; i++)
    {
     delete[] matrix_col_indices[m][i];
     delete[] matrix_values[m][i];
    }
   delete[] matrix_col_indices[m];
   delete[] matrix_values[m];
  }
  
 // need to wait for the recv nnzs to complete
 // before we can allocate storage for the matrix recvs
 // ===================================================

 // recv and copy the datafrom the recv storage to 
 // + nnz_from_proc
 // + n_eqn_from_proc
 Vector<MPI_Status> recv_nnz_stat(nproc-1);
 MPI_Waitall(nproc-1,&recv_nnz_reqs[0],&recv_nnz_stat[0]);
 Vector<unsigned> n_eqn_from_proc(nproc-1);
 DenseMatrix<unsigned> nnz_from_proc(nproc,n_matrix);
 for (unsigned p = 0; p < nproc; p++)
  {
   if (p != my_rank)
    {
   n_eqn_from_proc[p] = temp_recv_storage[p][0];
   for (unsigned m = 0; m < n_matrix; m++)
    {
     nnz_from_proc(p,m) = temp_recv_storage[p][m+1];
    }
   delete[] temp_recv_storage[p];
  }
   else
    {
     n_eqn_from_proc[p] = n_eqn_for_proc[p];
     for (unsigned m = 0; m < n_matrix; m++)
      {
       nnz_from_proc(p,m) = nnz_for_proc(p,m);
      }
    }
  }
 recv_nnz_stat.clear();
 recv_nnz_reqs.clear();

 // allocate the storage for the data to be recv and post the sends recvs
 // =====================================================================

 // storage
 Vector<unsigned*> eqns_from_proc(nproc);
 DenseMatrix<double*> residuals_from_proc(nproc,n_vector);
 DenseMatrix<unsigned*> row_start_from_proc(nproc,n_matrix);
 DenseMatrix<unsigned*> column_indices_from_proc(nproc,n_matrix);
 DenseMatrix<double*> values_from_proc(nproc,n_matrix);

 // allocate and post sends and recvs
 double base;
 MPI_Aint communication_base;
 MPI_Address(&base,&communication_base);
 unsigned n_comm_types = 1 + 1*n_vector + 3*n_matrix;
 Vector<MPI_Request> recv_reqs;
 Vector<MPI_Request> send_reqs;
 for (unsigned p = 0; p < nproc; p++)
  {
   if (p != my_rank)
    {

     // allocate
     if (n_eqn_from_proc[p] > 0)
      {
       eqns_from_proc[p] = new unsigned[n_eqn_from_proc[p]];
       for (unsigned v = 0; v < n_vector; v++)
        {
         residuals_from_proc(p,v) = new double[n_eqn_from_proc[p]];
        }
       for (unsigned m = 0; m < n_matrix; m++)
        {
         row_start_from_proc(p,m) = new unsigned[n_eqn_from_proc[p]+1];
         column_indices_from_proc(p,m) = new unsigned[nnz_from_proc(p,m)];
         values_from_proc(p,m) = new double[nnz_from_proc(p,m)];
        }
      }
     
     // recv
     if (n_eqn_from_proc[p] > 0)
      {
       MPI_Datatype types[n_comm_types];
       MPI_Aint offsets[n_comm_types];
       int count[n_comm_types];
       int pt = 0;
       
       // equations
       count[pt] = 1;
       MPI_Address(eqns_from_proc[p],&offsets[pt]);
       offsets[pt] -= communication_base;
       MPI_Type_contiguous(n_eqn_from_proc[p],MPI_UNSIGNED,&types[pt]);
       MPI_Type_commit(&types[pt]);
       pt++;
       
       // vectors
       for (unsigned v = 0; v < n_vector; v++)
        {
         count[pt] = 1;
         MPI_Address(residuals_from_proc(p,v),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(n_eqn_from_proc[p],MPI_DOUBLE,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;
        }

       // matrices
       for (unsigned m = 0; m < n_matrix; m++)
        {
         // row start
         count[pt] = 1;
         MPI_Address(row_start_from_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(n_eqn_from_proc[p]+1,MPI_UNSIGNED,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;         
         

         // column indices
         count[pt] = 1;
         MPI_Address(column_indices_from_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(nnz_from_proc(p,m),MPI_UNSIGNED,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;   

         // values
         count[pt] = 1;
         MPI_Address(values_from_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(nnz_from_proc(p,m),MPI_DOUBLE,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;
        }

       // build the combined type
       MPI_Datatype recv_type;
       MPI_Type_struct(n_comm_types,count,offsets,types,&recv_type);
       MPI_Type_commit(&recv_type);
       for (unsigned t = 0; t < n_comm_types; t++)
        {
         MPI_Type_free(&types[t]);
        }
       MPI_Request req;
       MPI_Irecv(&base,1,recv_type,p,1,
                 Communicator_pt->mpi_comm(),&req);
       MPI_Type_free(&recv_type);
       recv_reqs.push_back(req);
      }
     
     // send
     if (n_eqn_for_proc[p] > 0)
      {
       MPI_Datatype types[n_comm_types];
       MPI_Aint offsets[n_comm_types];
       int count[n_comm_types];
       int pt = 0;
       
       // equations
       count[pt] = 1;
       MPI_Address(eqns_for_proc[p],&offsets[pt]);
       offsets[pt] -= communication_base;
       MPI_Type_contiguous(n_eqn_for_proc[p],MPI_UNSIGNED,&types[pt]);
       MPI_Type_commit(&types[pt]);
       pt++;
       
       // vectors
       for (unsigned v = 0; v < n_vector; v++)
        {
         count[pt] = 1;
         MPI_Address(residuals_for_proc(p,v),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(n_eqn_for_proc[p],MPI_DOUBLE,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;
        }

       // matrices
       for (unsigned m = 0; m < n_matrix; m++)
        {
         // row start
         count[pt] = 1;
         MPI_Address(row_start_for_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(n_eqn_for_proc[p]+1,MPI_UNSIGNED,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;         
         

         // column indices
         count[pt] = 1;
         MPI_Address(column_indices_for_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(nnz_for_proc(p,m),MPI_UNSIGNED,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;   

         // values
         count[pt] = 1;
         MPI_Address(values_for_proc(p,m),&offsets[pt]);
         offsets[pt] -= communication_base;
         MPI_Type_contiguous(nnz_for_proc(p,m),MPI_DOUBLE,&types[pt]);
         MPI_Type_commit(&types[pt]);
         pt++;
        }

       // build the combined type
       MPI_Datatype send_type;
       MPI_Type_struct(n_comm_types,count,offsets,types,&send_type);
       MPI_Type_commit(&send_type);
       for (unsigned t = 0; t < n_comm_types; t++)
        {
         MPI_Type_free(&types[t]);
        }
       MPI_Request req;
       MPI_Isend(&base,1,send_type,p,1,
                 Communicator_pt->mpi_comm(),&req);
       MPI_Type_free(&send_type);
       send_reqs.push_back(req);
      }
    }
   // otherwise send to self
   else
    {
     eqns_from_proc[p] = eqns_for_proc[p];
     for (unsigned v = 0; v < n_vector; v++)
      {
       residuals_from_proc(p,v) = residuals_for_proc(p,v);
      }
     for (unsigned m = 0; m < n_matrix; m++)
      {
       row_start_from_proc(p,m) = row_start_for_proc(p,m);
       column_indices_from_proc(p,m) = column_indices_for_proc(p,m);
       values_from_proc(p,m) = values_for_proc(p,m);
      }
    }
  }

 // wait for the recvs to complete
 unsigned n_recv_req = recv_reqs.size();
 if (n_recv_req > 0)
  {
   Vector<MPI_Status> recv_stat(n_recv_req);
   MPI_Waitall(n_recv_req,&recv_reqs[0],&recv_stat[0]);
  }

 // ==============================================
 unsigned target_nrow_local = target_dist_pt->nrow_local();

 // loop over the matrices
 for (unsigned m = 0; m < n_matrix; m++)
  {

   // allocate row_start
   row_start[m] = new int[target_nrow_local+1];
   row_start[m][0] = 0;

   // initially allocate storage based on the maximum number of non-zeros
   // from any one processor
   unsigned nnz_allocation = Parallel_sparse_assemble_previous_allocation;
   for (unsigned p = 0; p < nproc; p++)
    {
     nnz_allocation = std::max(nnz_allocation,nnz_from_proc(p,m));
    }
   Vector<double*> values_chunk(1);
   values_chunk[0] = new double[nnz_allocation];
   Vector<int*> column_indices_chunk(1);
   column_indices_chunk[0] = new int[nnz_allocation];
   Vector<unsigned> ncoef_in_chunk(1,0);
   Vector<unsigned> size_of_chunk(1,0);
   size_of_chunk[0] = nnz_allocation;
   unsigned current_chunk = 0;

   // for each row on this processor
   for (unsigned i = 0; i < target_nrow_local; i++)
    {
     row_start[m][i]=0;

     // determine the processors that this row is on
     Vector<int> row_on_proc(nproc,-1);
     for (unsigned p = 0; p < nproc; p++)
      {
       if (n_eqn_from_proc[p]==0)
        {
         row_on_proc[p]=-1;
        }
       else
        {
       int left = 0;
       int right = n_eqn_from_proc[p]-1;
       int midpoint = right/2;
       bool complete = false;
       while (!complete)
        {
         midpoint = (right+left)/2;
         if (midpoint > right) { midpoint = right; }
         if (midpoint < left) { midpoint = left; }
         if (left==right)
          {
           if (eqns_from_proc[p][midpoint] == i)
            {
             midpoint=left;
            }
           else
            {
             midpoint = -1;
            }
           complete=true;
          }
         else if (eqns_from_proc[p][midpoint] == i)
          {
           complete = true;
          }
         else if (eqns_from_proc[p][midpoint] > i)
          {
           right = std::max(midpoint-1,left);
          }
         else
          {
           left = std::min(midpoint+1,right);
          }
        }
       row_on_proc[p] = midpoint;
        }
      }

     // for each processor build this row of the matrix
     unsigned check_first = ncoef_in_chunk[current_chunk];
     unsigned check_last = check_first;
     for (unsigned p = 0; p < nproc; p++)
      {
       if (row_on_proc[p] != -1)
        {
         int row = row_on_proc[p];
         unsigned first = row_start_from_proc(p,m)[row];
         unsigned last = row_start_from_proc(p,m)[row+1];
         for (unsigned l = first; l < last; l++)
          {
           bool done = false;
           for (unsigned j = check_first; j <=check_last && !done; j++)
            {
             if (j==check_last)
              {
               // is this temp array full, do we need to allocate
               // a new temp array
               if (ncoef_in_chunk[current_chunk] ==
                   size_of_chunk[current_chunk])
                {

                 // number of chunks allocated
                 unsigned n_chunk = values_chunk.size();

                 // determine the number of non-zeros added so far
                 // (excluding the current row)
                 unsigned nnz_so_far = 0;
                 for (unsigned c = 0; c < n_chunk; c++)
                  {
                   nnz_so_far += ncoef_in_chunk[c];
                  }
                 nnz_so_far -= row_start[m][i];

                 // average number of non-zeros per row
                 unsigned avg_nnz = nnz_so_far/(i+1);

                 // number of rows left +1
                 unsigned nrows_left = target_nrow_local-i;

                 // allocation for next chunk
                 unsigned next_chunk_size = avg_nnz*nrows_left+row_start[m][i];

                 // allocate storage in next chunk
                 current_chunk++;
                 n_chunk++;
                 values_chunk.resize(n_chunk);
                 values_chunk[current_chunk] = new double[next_chunk_size];
                 column_indices_chunk.resize(n_chunk);
                 column_indices_chunk[current_chunk] 
                  = new int[next_chunk_size];
                 size_of_chunk.resize(n_chunk);
                 size_of_chunk[current_chunk] = next_chunk_size;
                 ncoef_in_chunk.resize(n_chunk);

                 // copy current row from previous chunk to new chunk
                 for (unsigned k = check_first; k < check_last; k++)
                  {
                   values_chunk[current_chunk][k-check_first] =
                    values_chunk[current_chunk-1][k];
                   column_indices_chunk[current_chunk][k-check_first] = 
                    column_indices_chunk[current_chunk-1][k];
                  }
                 ncoef_in_chunk[current_chunk-1]-=row_start[m][i];
                 ncoef_in_chunk[current_chunk]=row_start[m][i];
          
                 // update first_check and last_check
                 check_first=0;
                 check_last=row_start[m][i];
                 j=check_last;
                }

               // add the coefficient             
               values_chunk[current_chunk][j]=values_from_proc(p,m)[l];
               column_indices_chunk[current_chunk][j]
                =column_indices_from_proc(p,m)[l];
               ncoef_in_chunk[current_chunk]++;
               row_start[m][i]++;
               check_last++;
               done=true;
              }
             else if (column_indices_chunk[current_chunk][j] == 
                      (int)column_indices_from_proc(p,m)[l])
              {
               values_chunk[current_chunk][j] += values_from_proc(p,m)[l];
               done=true;
              }
            }
          }
        }
      }
    }

   // delete recv data for this matrix
   for (unsigned p = 0; p < nproc; p++)
    {
     if (n_eqn_from_proc[p] > 0)
      {
       delete[] row_start_from_proc(p,m);
       delete[] column_indices_from_proc(p,m);
       delete[] values_from_proc(p,m);
      }
    }

   // next we take the chunk base storage of the column indices and values
   // and copy into a single contiguous block of memory
   // ====================================================================
   unsigned n_chunk = values_chunk.size();
   nnz[m]=0;
   for (unsigned c = 0; c < n_chunk; c++)
    {
     nnz[m]+=ncoef_in_chunk[c];
    }
   Parallel_sparse_assemble_previous_allocation=nnz[m];

   // allocate
   values[m] = new double[nnz[m]];
   column_indices[m] = new int[nnz[m]];
                               
   // copy
   unsigned pt = 0;
   for (unsigned c = 0; c < n_chunk; c++)
    {
     unsigned nc = ncoef_in_chunk[c];
     for (unsigned i = 0; i < nc; i++)
      {
       values[m][pt+i]=values_chunk[c][i];
       column_indices[m][pt+i]=column_indices_chunk[c][i];
      }
     pt += nc;
     delete[] values_chunk[c];
     delete[] column_indices_chunk[c];
    }

   // the row_start vector currently contains the number of coefs in each
   // row. Update
   // ===================================================================
   unsigned g = row_start[m][0];
   row_start[m][0] = 0;
   for (unsigned i = 1; i < target_nrow_local; i++)
    {
     unsigned h = g+row_start[m][i];
     row_start[m][i] = g;
     g=h;
    }
   row_start[m][target_nrow_local]=g;
  }

 // next accumulate the residuals
 for (unsigned v = 0; v < n_vector; v++)
  {
   residuals[v] = new double[target_nrow_local];
   for (unsigned i = 0; i < target_nrow_local; i++)
    {
     residuals[v][i] = 0;
    }
   for (unsigned p = 0; p < nproc; p++)
    {
     if (n_eqn_from_proc[p] > 0)
      {
       unsigned n_eqn_p = n_eqn_from_proc[p];
       for (unsigned i = 0; i < n_eqn_p; i++)
        {
         residuals[v][eqns_from_proc[p][i]]+=residuals_from_proc(p,v)[i];
        }
       delete[] residuals_from_proc(p,v);
      }
    }
  }

 // delete list of eqns from proc
 for (unsigned p = 0; p < nproc; p++)
  {
   if (n_eqn_from_proc[p] > 0)
    {
     delete[] eqns_from_proc[p];
    }
  }
 
 // and wait for sends to complete
 Vector<MPI_Status> send_nnz_stat(nproc-1);
 MPI_Waitall(nproc-1,&send_nnz_reqs[0],&send_nnz_stat[0]);
 for (unsigned p = 0; p < nproc; p++)
  {
   if (p != my_rank)
    {
   delete[] temp_send_storage[p];
  }
  }
 send_nnz_stat.clear();
 send_nnz_reqs.clear();

 // wait for the matrix data sends to complete and delete the data
 unsigned n_send_reqs = send_reqs.size();
 if (n_send_reqs > 0)
  {
   Vector<MPI_Status> send_stat(n_send_reqs);
   MPI_Waitall(n_send_reqs,&send_reqs[0],&send_stat[0]);
   for (unsigned p = 0; p < nproc; p++)
    {
     if (p != my_rank)
      {
       if (n_eqn_for_proc[p])
        {
         delete[] eqns_for_proc[p];
         for (unsigned m = 0; m < n_matrix; m++)
          {
           delete[] row_start_for_proc(p,m);
           delete[] column_indices_for_proc(p,m);
           delete[] values_for_proc(p,m);
          }
         for (unsigned v = 0; v < n_vector; v++)
          {
           delete[] residuals_for_proc(p,v);
          }
        }
      }
    }
  }
}
#endif


//================================================================
/// \short Get the full Jacobian by finite differencing
//================================================================
void Problem::get_fd_jacobian(DoubleVector &residuals, 
                              DenseMatrix<double> &jacobian)
{

#ifdef OOMPH_HAS_MPI

 if (Problem_has_been_distributed)
  {
   OomphLibWarning("This is unlikely to work with a distributed problem",
                   " Problem::get_fd_jacobian()",
                   OOMPH_EXCEPTION_LOCATION);
  }
#endif


 //Find number of dofs
 const unsigned long n_dof = ndof();
 
// Advanced residuals
 DoubleVector residuals_pls;

 // Get reference residuals
 get_residuals(residuals);

 const double FD_step=1.0e-8;
 //Loop over all dofs
 for(unsigned long jdof=0;jdof<n_dof;jdof++) 
  {
   double backup=*Dof_pt[jdof];
   *Dof_pt[jdof]+=FD_step;

   // We're checking if the new values for Dof_pt[] actually
   // solve the entire problem --> update as if problem had
   // been solved 
   actions_before_newton_solve();
   actions_before_newton_convergence_check();
   actions_after_newton_solve();

   // Get advanced residuals
   get_residuals(residuals_pls);

   for (unsigned long ieqn=0;ieqn<n_dof;ieqn++) 
    {
     jacobian(ieqn,jdof)=(residuals_pls[ieqn]-residuals[ieqn])/FD_step;
    }
   *Dof_pt[jdof]=backup;
  }

 // Reset problem to state it was in
 actions_before_newton_solve();
 actions_before_newton_convergence_check();
 actions_after_newton_solve();

}

//======================================================================
/// \short Get derivative of the residuals vector wrt a global parameter
/// This is required in continuation problems
//=======================================================================
void Problem::get_derivative_wrt_global_parameter(double* const &parameter_pt,
                                                  DoubleVector &result)
{

#ifdef OOMPH_HAS_MPI

 if (Problem_has_been_distributed)
  {
   OomphLibWarning("This is unlikely to work with a distributed problem",
                   "Problem::get_derivative_wrt_global_parameter()",
                   OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Find the number of degrees of freedom in the problem
 const unsigned long n_dof = this->ndof();

 // create the linear algebra distribution for this solver
 // currently only global (non-distributed) distributions are allowed
 LinearAlgebraDistribution* dist_pt = new 
  LinearAlgebraDistribution(Communicator_pt,n_dof,false);

 // rebuild the result
 result.build(dist_pt,0.0);

 //Initialise the result vector to zero 
 result.initialise(0.0);

 //Get the residuals and store in the result vector
 get_residuals(result);

 //Storage for the new residuals
 DoubleVector newres(dist_pt,0.0);
  
 //Increase the global parameter
 const double FD_step = 1.0e-8;
 
 //Store the current value of the parameter
 double param_value = *parameter_pt;

 //Increase the parameter
 *parameter_pt += FD_step;

 //Do any possible updates
 actions_after_change_in_global_parameter();

 //Get the new residuals
 get_residuals(newres);

 //Do the finite differencing
 for(unsigned long n=0;n<n_dof;++n)
  {
   result[n] = (newres[n] - result[n])/FD_step;
  }

 //Reset the value of the parameter
 *parameter_pt = param_value;

 //Do any possible updates
 actions_after_change_in_global_parameter();
}

//================================================================
/// \short Get derivative of an element in the problem wrt a global
/// parameter, to be used in continuation problems
//================================================================
void Problem::get_derivative_wrt_global_parameter(
 double* const &parameter_pt,
 GeneralisedElement* const &elem_pt,
 Vector<double> &result)
{

#ifdef OOMPH_HAS_MPI

 if (Problem_has_been_distributed)
  {
   OomphLibWarning("This is unlikely to work with a distributed problem",
                   "Problem::get_derivative_wrt_global_parameter()",
                   OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Locally cache pointer to assembly handler
 AssemblyHandler* const assembly_handler_pt = Assembly_handler_pt;

 //Should definitely give this a more global scope
 double FD_Jstep = 1.0e-8;
 
 //Find the number of variables in the element, e
 unsigned nvar = assembly_handler_pt->ndof(elem_pt);
 //Create storage for residuals
 Vector<double> residuals(nvar), newres(nvar);

 //Get the "original" residuals
 assembly_handler_pt->get_residuals(elem_pt,residuals);

 //Save the old value of the global parameter
 double old_var = *parameter_pt;

 //Increment the value
 *parameter_pt += FD_Jstep;
 
 //Now do any possible updates
 actions_after_change_in_global_parameter();

 //Get the "new" residuals
 assembly_handler_pt->get_residuals(elem_pt,newres);

 //Do the finite differences
 for(unsigned m=0;m<nvar;m++)
  {
   result[m] = (newres[m] - residuals[m])/FD_Jstep;
  }
 
 //Reset value of the global parameter
 *parameter_pt = old_var;

 //Now do any possible updates
 actions_after_change_in_global_parameter();
}

//==================================================================
/// Solve the eigenproblem
//==================================================================
void Problem::solve_eigenproblem(const unsigned &n_eval,
                                 Vector<std::complex<double> > &eigenvalue,
                                 Vector<DoubleVector> &eigenvector,
                                 const bool &steady)
{
 //If the boolean flag is steady, then make all the timesteppers steady
 //before solving the eigenproblem. This will "switch off" the time-derivative
 //terms in the jacobian matrix
 if(steady)
  {
   //Find out how many timesteppers there are
   const unsigned n_time_steppers = ntime_stepper();
   
   // Vector of bools to store the is_steady status of the various
   // timesteppers when we came in here
   std::vector<bool> was_steady(n_time_steppers);
   
   //Loop over them all and make them (temporarily) static
   for(unsigned i=0;i<n_time_steppers;i++)
    {
     was_steady[i]=time_stepper_pt(i)->is_steady();
     time_stepper_pt(i)->make_steady();
    }
   
   //Call the Eigenproblem for the eigensolver
   Eigen_solver_pt->solve_eigenproblem(this,n_eval,eigenvalue,eigenvector);
   
   // Reset the is_steady status of all timesteppers that
   // weren't already steady when we came in here and reset their
   // weights
   for(unsigned i=0;i<n_time_steppers;i++)
    {
     if (!was_steady[i])
      {
       time_stepper_pt(i)->undo_make_steady();
      }
    }
  }
 //Otherwise if we don't want to make the problem steady, just
 //assemble and solve the eigensystem
 else
  {
   //Call the Eigenproblem for the eigensolver
   Eigen_solver_pt->solve_eigenproblem(this,n_eval,eigenvalue,eigenvector);
  }
}




//===================================================================
/// Get the matrices required to solve an eigenproblem
/// WARNING: temporarily this method only works with non-distributed 
/// matrices
//===================================================================
void Problem::get_eigenproblem_matrices(CRDoubleMatrix &mass_matrix,
                                        CRDoubleMatrix &main_matrix,
                                        const double &shift)
{
#ifdef PARANOID
 if (mass_matrix.distribution_setup())
  {
   if (mass_matrix.nrow() != this->ndof())
    {   
     std::ostringstream error_stream;
     error_stream
      << "mass_matrix has a distribution, but the number of rows is not "
      << "equal to the number of degrees of freedom in the problem.";
     throw OomphLibError(error_stream.str(),
                         "Problem::get_eigen_problem_matrices(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (mass_matrix.distribution_pt()->distributed())
    {
     std::ostringstream error_stream;
     error_stream
      << "mass_matrix cannot be distributed";
     throw OomphLibError(error_stream.str(),
                         "Problem::get_eigen_problem_matrices(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
 if (main_matrix.distribution_setup())
  {
   if (main_matrix.nrow() != this->ndof())
    {
     std::ostringstream error_stream;
     error_stream
      << "main_matrix has a distribution, but the number of rows is not "
      << "equal to the number of degrees of freedom in the problem.";
     throw OomphLibError(error_stream.str(),
                         "Problem::get_eigen_problem_matrices(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (main_matrix.distribution_pt()->distributed())
    {
     std::ostringstream error_stream;
     error_stream
      << "main_matrix cannot be distributed";
     throw OomphLibError(error_stream.str(),
                         "Problem::get_eigen_problem_matrices(...)",
                         OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif 

 // if the matrices are not setup then build them
 if (!main_matrix.distribution_setup() || !main_matrix.distribution_setup())
  {
   LinearAlgebraDistribution dist(this->communicator_pt(),this->ndof(),false);
   if (!main_matrix.distribution_setup())
    {
     main_matrix.build(&dist);
    }
   if (!mass_matrix.distribution_setup())
    {
     mass_matrix.build(&dist);
    }
  }

 //Store the old assembly handler
 AssemblyHandler* old_assembly_handler_pt = Assembly_handler_pt;
 //Now setup the eigenproblem handler, pass in the value of the shift
 Assembly_handler_pt = new EigenProblemHandler(shift);

 //Prepare the storage formats.
 Vector<int* > column_or_row_index(2);
 Vector<int* > row_or_column_start(2);
 Vector<double* > value(2); 
 Vector<unsigned> nnz(2);
 //Allocate pointer to residuals, although not used in these problems
 Vector<double* > residuals_vectors(0);

 bool compressed_row_flag=true;

 sparse_assemble_row_or_column_compressed(column_or_row_index,
                                          row_or_column_start,
                                          value,
                                          nnz,
                                          residuals_vectors,
                                          compressed_row_flag);

 //Get the number of dofs (size of matrices)
 unsigned long n_dof = ndof();

 //The main matrix is the first entry
 main_matrix.build_matrix_without_copy(n_dof,nnz[0],value[0],
                                       column_or_row_index[0],
                                       row_or_column_start[0]);
 //The mass matrix is the second entry
 mass_matrix.build_matrix_without_copy(n_dof,nnz[1],value[1],
                                       column_or_row_index[1],
                                       row_or_column_start[1]);   

 //Delete the eigenproblem handler
 delete Assembly_handler_pt;
 //Reset the assembly handler to the original handler
 Assembly_handler_pt = old_assembly_handler_pt;
}


//=======================================================================
/// Stored the current values of the dofs
//=======================================================================
void Problem::store_current_dof_values()
{
 //Find the number of dofs
 unsigned long n_dof = ndof();
 //If memory has not been allocated, then allocated memory for the saved
 //dofs
 if(Saved_dof_pt==0) {Saved_dof_pt = new Vector<double>;}
 //Resize the vector
 Saved_dof_pt->resize(n_dof);

 //Transfer the values over
 for(unsigned long n=0;n<n_dof;n++) {(*Saved_dof_pt)[n] = dof(n);}
}

//====================================================================
/// Restore the saved dofs
void Problem::restore_dof_values()
{
 //Check that we can do this
 if(Saved_dof_pt==0) 
  {
   throw OomphLibError(
    "There are no stored values, use store_current_dof_values()\n",
    "Problem::restore_dof_values()",
    OOMPH_EXCEPTION_LOCATION);
  }

 //Find the number of dofs
 unsigned long n_dof = ndof();

 if(Saved_dof_pt->size() != n_dof)
  {
   throw OomphLibError(
    "The number of stored values is not equal to the current number of dofs\n",
    "Problem::restore_dof_values()",
    OOMPH_EXCEPTION_LOCATION);
  }

 //Transfer the values over
 for(unsigned long n=0;n<n_dof;n++) {dof(n) = (*Saved_dof_pt)[n];}

 //Delete the memory
 delete Saved_dof_pt;
 Saved_dof_pt = 0;
}

//======================================================================
/// Assign the eigenvector passed to the function to the dofs
//======================================================================
void Problem::assign_eigenvector_to_dofs(DoubleVector &eigenvector) 
{
 unsigned long n_dof = ndof();
 //Check that the eigenvector has the correct size
 if(eigenvector.nrow() != n_dof)
  {
   std::ostringstream error_message;
   error_message << "Eigenvector has size " << eigenvector.nrow() 
                 << ", not equal to the number of dofs in the problem," 
                 << n_dof << std::endl;

   throw OomphLibError(error_message.str(),
                       "Problem::assign_eigenvector_to_dofs()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Loop over the dofs and assign the eigenvector
 for(unsigned long n=0;n<n_dof;n++)
  {
   dof(n) = eigenvector[n];
  }
}



//======================================================================
/// Add the eigenvector passed to the function to the dofs with 
/// magnitude epsilon
//======================================================================
void Problem::add_eigenvector_to_dofs(const double &epsilon,
                                      DoubleVector &eigenvector) 
{
 unsigned long n_dof = ndof();
 //Check that the eigenvector has the correct size
 if(eigenvector.nrow() != n_dof)
  {
   std::ostringstream error_message;
   error_message << "Eigenvector has size " << eigenvector.nrow() 
                 << ", not equal to the number of dofs in the problem," 
                 << n_dof << std::endl;

   throw OomphLibError(error_message.str(),
                       "Problem::assign_eigenvector_to_dofs()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Loop over the dofs and add the eigenvector
 for(unsigned long n=0;n<n_dof;n++)
  {
   dof(n) += epsilon*eigenvector[n];
  }
}


//================================================================
/// General Newton solver. Requires only a convergence tolerance. 
/// The linear solver takes a pointer to the problem (which defines
/// the Jacobian \b J and the residual Vector \b r) and returns 
/// the solution \b x of the system 
/// \f[ {\bf J} {\bf x} = - \bf{r} \f]. 
//================================================================
void Problem::newton_solve()
{

 // Initialise timers
 double total_linear_solver_time=0.0;
 double t_start = TimingHelpers::timer();

 //Find total number of dofs
 unsigned long n_dofs = ndof();

 //Set up the Vector to hold the solution
 DoubleVector dx;

 //Set the counter
 unsigned count=0;
 //Set the loop flag
 unsigned LOOP_FLAG=1;

 //Update anything that needs updating
 actions_before_newton_solve();

 //Now do the Newton loop
 do
  {
   count++;
   
   //Do any updates that are required 
   actions_before_newton_step();


   // No degrees of freedom? What are you solving for?
   if (n_dofs==0)
    {
     oomph_info << std::endl << std::endl << std::endl;
     oomph_info << "This is a bit bizarre: The problem has no dofs." 
                << std::endl;
     oomph_info 
      << "I'll just return from the Newton solver without doing anything." 
      << std::endl;

     // Do any updates that would have been performed
     actions_before_newton_convergence_check();
     actions_after_newton_step();
     actions_before_newton_convergence_check();
     actions_after_newton_solve();

     oomph_info << "I hope this is what you intended me to do..." << std::endl;
     oomph_info 
      << std::endl << "Note: All actions_...() functions were called" 
      << std::endl;
     oomph_info << std::endl << "      before returning." << std::endl;
     oomph_info << std::endl << std::endl << std::endl;
     return;
    }

   //Calculate initial residuals
   if(count==1)
    {
     // Is the problem nonlinear? If not ignore the pre-iteration 
     // convergence check.
     if (Problem_is_nonlinear)
      {
#ifdef OOMPH_HAS_MPI
       // Synchronise the solution on different processors (on each submesh)
       unsigned nmesh=nsub_mesh();
       if (nmesh==0)
        {
         synchronise_dofs(mesh_pt());
        }
       else
        {
         // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
         for (unsigned imesh=0; imesh<nmesh; imesh++)
          {
           synchronise_dofs(mesh_pt(imesh));
          }
         for (unsigned imesh=0; imesh<nmesh; imesh++)
          {
           synchronise_external_dofs(mesh_pt(imesh));
          }
        }
#endif

       actions_before_newton_convergence_check();
       dx.clear();
       get_residuals(dx);

       //Get maximum residuals
       double maxres = dx.max();

       if (!Shut_up_in_newton_solve) 
        {
         oomph_info << "Initial Maximum residuals " << maxres << std::endl;
        }
       if(maxres < Newton_solver_tolerance) {LOOP_FLAG=0; continue;}
      }
     else
      {
       if (!Shut_up_in_newton_solve) 
        {
         oomph_info 
          << "Linear problem -- convergence in one iteration assumed." 
          << std::endl;
        }
      }
    }

   // Initialise timer for linear solver
   double t_solver_start = TimingHelpers::timer();
   
   //Now do the linear solve -- recycling Jacobian if requested
   if (Jacobian_reuse_is_enabled&&Jacobian_has_been_computed)
    {     
     if (!Shut_up_in_newton_solve) 
      {
       oomph_info << "Not recomputing Jacobian! " << std::endl;
      }
     
     // If we're doing the first iteration and the problem is nonlinear, 
     // the residuals have already been computed above during the
     // initial convergence check. Otherwise compute them here.
     if ((count!=1)||(!Problem_is_nonlinear)) get_residuals(dx);

     // Backup residuals
     DoubleVector resid(dx);
     
     // Resolve
     Linear_solver_pt->resolve(resid,dx);
    }
   else
    {
     if (Jacobian_reuse_is_enabled)
      {
       if (!Shut_up_in_newton_solve) 
        {
         oomph_info << "Enabling resolve" << std::endl;
        }
       Linear_solver_pt->enable_resolve();
      }  
     Linear_solver_pt->solve(this,dx);
     Jacobian_has_been_computed=true;
    }

   // End of linear solver
   double t_solver_end = TimingHelpers::timer();
   total_linear_solver_time+=t_solver_end-t_solver_start;

   if (!Shut_up_in_newton_solve) 
    {
     oomph_info << std::endl;
     oomph_info << "Time for linear solver (ndof="
                << n_dofs << ") [sec]: " 
                << t_solver_end-t_solver_start 
                << std::endl << std::endl;
    }

   //Subtract the new values from the true dofs
   dx.redistribute(Dof_distribution_pt);
   double* dx_pt = dx.values_pt();
   unsigned ndof_local = Dof_distribution_pt->nrow_local();
   for (unsigned l = 0; l < ndof_local; l++)
    {
     *Dof_pt[l] -= dx_pt[l];
    }

#ifdef OOMPH_HAS_MPI
   // Synchronise the solution on different processors (on each submesh)
   unsigned nmesh=nsub_mesh();
   if (nmesh==0)
    {
     synchronise_dofs(mesh_pt());
    }
   else
    {
     // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
     for (unsigned imesh=0; imesh<nmesh; imesh++)
      {
       synchronise_dofs(mesh_pt(imesh));
      }
     for (unsigned imesh=0; imesh<nmesh; imesh++)
      {
       synchronise_external_dofs(mesh_pt(imesh));
      }
    }
#endif

   // Do any updates that are required 
   actions_after_newton_step();
   actions_before_newton_convergence_check();

   // Maximum residuals
   double maxres=0.0;
   // If the user has declared that the Problem is linear
   // we ignore the convergence check
   if (Problem_is_nonlinear)
    {
     //Get the maximum residuals
     //maxres = std::abs(*std::max_element(dx.begin(),dx.end(),
     //                                    AbsCmp<double>()));
     //std::cout << "Maxres correction " << maxres << "\n";

     //Calculate the new residuals
     dx.clear();
     get_residuals(dx);

     //Get the maximum residuals
     maxres = dx.max();

     if (!Shut_up_in_newton_solve) 
      {
       oomph_info << "Newton Step " << count << ": Maximum residuals "
                  << maxres << std::endl;
      }
    }

   //If we have converged jump straight to the test at the end of the loop
   if(maxres < Newton_solver_tolerance) {LOOP_FLAG=0; continue;} 

   //This section will not be reached if we have converged already
   //If the maximum number of residuals is too high or the maximum number
   //of iterations has been reached
   if((maxres > Max_residuals) || (count == Max_newton_iterations))
    {
     // Print a warning -- regardless of what the throw does
     if (maxres > Max_residuals) 
      {
       oomph_info << "Max. residual (" << Max_residuals 
                  << ") has been exceeded in Newton solver." << std::endl;
      }
     if (count == Max_newton_iterations)
      {
       oomph_info << "Reached max. number of iterations (" 
                  << Max_newton_iterations
                  << ") in Newton solver." << std::endl;
      }
     // Now throw...
     throw NewtonSolverError(count,maxres);
    }

  }
 while(LOOP_FLAG);

 //Now update anything that needs updating
 actions_after_newton_solve();

 // Finalise/doc timings
 if (!Shut_up_in_newton_solve) 
  {
   oomph_info << std::endl;
   oomph_info << "Total time for linear solver (ndof="<< n_dofs << ") [sec]: " 
              << total_linear_solver_time << std::endl;
  }

 double t_end = TimingHelpers::timer();
 double total_time=t_end-t_start;

 if (!Shut_up_in_newton_solve) 
  {
   oomph_info << "Total time for Newton solver (ndof="<< n_dofs << ") [sec]: " 
              << total_time << std::endl;
  }
 if (total_time>0.0)
  {
   if (!Shut_up_in_newton_solve) 
    {
     oomph_info << "Time outside linear solver        : "
                << (total_time-total_linear_solver_time)/total_time*100.0
                << " %"
                << std::endl;
    }
  }
 else
  {
   if (!Shut_up_in_newton_solve) 
    {
     oomph_info << "Time outside linear solver        : "
                << (total_time-total_linear_solver_time)/total_time*100.0
                << "[too fast]"
                << std::endl;
    }
  }
 if (!Shut_up_in_newton_solve) oomph_info << std::endl;
}  


//========================================================================
/// Solve a steady problem, in the context of an overall unsteady problem.
/// This is achieved by setting the weights in the timesteppers to be zero
/// which has the effect of rendering them steady timesteppers
/// The optional argument max_adapt specifies the max. number of 
/// adaptations of all refineable submeshes are performed to 
/// achieve the the error targets specified in the refineable submeshes.
//========================================================================
void Problem::steady_newton_solve(unsigned const &max_adapt)
{
 //Find out how many timesteppers there are
 unsigned n_time_steppers = ntime_stepper();

 // Vector of bools to store the is_steady status of the various
 // timesteppers when we came in here
 std::vector<bool> was_steady(n_time_steppers);

 //Loop over them all and make them (temporarily) static
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   was_steady[i]=time_stepper_pt(i)->is_steady();
   time_stepper_pt(i)->make_steady();
  }

 try
  {
   //Solve the non-linear problem with Newton's method
   if (max_adapt==0)
    {
     newton_solve();
    }
   else
    {
     newton_solve(max_adapt);
    }
  }
 //Catch any exceptions thrown in the Newton solver
 catch(NewtonSolverError &error)
  {
   oomph_info << std::endl << "USER-DEFINED ERROR IN NEWTON SOLVER " 
              << std::endl;
   //Check whether it's the linear solver
   if(error.linear_solver_error)
    {
     oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
    }
   //Check to see whether we have reached Max_iterations
   else if(error.iterations==Max_newton_iterations)
    {
     oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" << error.iterations << 
      ") REACHED WITHOUT CONVERGENCE " << std::endl;
    }
   //If not, it must be that we have exceeded the maximum residuals
   else
    {
     oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                << " EXCEEDS PREDEFINED MAXIMUM " << Max_residuals
                << std::endl;
    }

   //Die horribly!!
   std::ostringstream error_stream;
   error_stream << "Error occured in Newton solver. "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       "Problem::steady_newton_solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }


 // Reset the is_steady status of all timesteppers that
 // weren't already steady when we came in here and reset their
 // weights
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   if (!was_steady[i])
    {
     time_stepper_pt(i)->undo_make_steady();
    }
  }

 // Since we performed a steady solve, the history values 
 // now have to be set as if we had performed an impulsive start from
 // the current solution. This ensures that the time-derivatives
 // evaluate to zero even now that the timesteppers have been
 // reactivated. 
 assign_initial_values_impulsive();

}

//===========================================================================
///Perform a basic continuation step using Newton's method. The governing
///parameter of the problem is passed as a pointer to the routine. The
///number of Newton steps taken is returned
//==========================================================================
unsigned Problem::
newton_solve_continuation(double* const &parameter_pt)
{
 //Set up memory for z 
 unsigned long n_dofs = ndof();
 LinearAlgebraDistribution dist(Communicator_pt,n_dofs,false);
 DoubleVector z(&dist,0.0);
 //Call the solver
 return newton_solve_continuation(parameter_pt,z);
}


//===================================================================
/// This function performs a basic continuation step using the Newton method. 
/// The number of Newton steps taken is returned, to be used in any 
/// external step-size control routines. 
/// The governing parameter of the problem is passed as a pointer to the
/// routine, as is the sign of the Jacobian and a Vector in which 
/// to store the derivatives wrt the parameter, if required.
//==================================================================
unsigned Problem::
newton_solve_continuation(double* const &parameter_pt,
                          DoubleVector &z)
{

 //Find the total number of dofs
 unsigned long n_dofs = ndof();

 // create the distribution (not distributed)
 LinearAlgebraDistribution dist(this->communicator_pt(),n_dofs,false);

 //Assign memory for solutions of the equations
 DoubleVector y(&dist,0.0);

 //Assign memory for the dot products of the uderivatives and y and z
 double uderiv_dot_y = 0.0, uderiv_dot_z = 0.0;
 //Set and initialise the counter
 unsigned count=0;
 //Set the loop flag
 unsigned LOOP_FLAG=1;
     
 //Update anything that needs updating
 actions_before_newton_solve();

 //Check the arc-length constraint
 double arc_length_constraint_residual=0.0;

 //Are we storing the matrix in the linear solve
 bool enable_resolve = Linear_solver_pt->resolve_is_enabled();
 
 //For this problem, we must store the residuals
 Linear_solver_pt->enable_resolve();
 
 //Now do the Newton loop
 do
  {
   count++;
       
   //Do any updates that are required 
   actions_before_newton_step();
   
   //Calculate initial residuals
   if(count==1)
    {
#ifdef OOMPH_HAS_MPI
     // Synchronise the solution on different processors (on each submesh)
     unsigned nmesh=nsub_mesh();
     if (nmesh==0)
      {
       synchronise_dofs(mesh_pt());
      }
     else
      {
       // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
       for (unsigned imesh=0; imesh<nmesh; imesh++)
        {
         synchronise_dofs(mesh_pt(imesh));
        }
       for (unsigned imesh=0; imesh<nmesh; imesh++)
        {
         synchronise_external_dofs(mesh_pt(imesh));
        }
      }
#endif

     actions_before_newton_convergence_check();

     get_residuals(y);
     //Get maximum residuals, using our own abscmp function
     double maxres = y.max();

     //Assemble the residuals for the arc-length step
     arc_length_constraint_residual = 0.0; 
     //Add the variables
     for(unsigned long l=0;l<n_dofs;l++)
      {
       arc_length_constraint_residual +=
        Dof_derivatives[l]*(*Dof_pt[l] - Dofs_current[l]);
      }
     arc_length_constraint_residual *= Theta_squared;
     arc_length_constraint_residual += 
      Parameter_derivative*(*parameter_pt - Parameter_current) - Ds_current;

     //Is it the max
     if(std::abs(arc_length_constraint_residual) > maxres)
      {
       maxres = std::abs(arc_length_constraint_residual);
      }

     //Find the max

     //If we are below the Tolerance, then return immediately
     if(maxres < Newton_solver_tolerance) {LOOP_FLAG=0; count=0; continue;}
    }
   
   //If it's the block hopf solver we need to solve for both rhs's
   //simultaneously. This is because the block decomposition involves
   //solves with two different matrices and storing both at once to 
   //allow general resolves would be more expensive than necessary.
   if(dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt))
    {
     //Get the vector dresiduals/dparameter
     get_derivative_wrt_global_parameter(parameter_pt,z);
     
     // Copy rhs vector into local storage so it doesn't get overwritten
     // if the linear solver decides to initialise the solution vector, say,
     // which it's quite entitled to do!
     DoubleVector input_z(z);

     //Solve the system for the two right-hand sides.
     dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt)->
      solve_for_two_rhs(this,y,input_z,z);
    }
   //Otherwise
   else
    {
     //Solve the standard problem
     Linear_solver_pt->solve(this,y);

     //Get the vector dresiduals/dparameter
     get_derivative_wrt_global_parameter(parameter_pt,z);
     
     // Copy rhs vector into local storage so it doesn't get overwritten
     // if the linear solver decides to initialise the solution vector, say,
     // which it's quite entitled to do!
     DoubleVector input_z(z);
     
     //Now resolve the system with the new RHS
     Linear_solver_pt->resolve(input_z,z);
    }
  
   //Now we need to calculate dparam, for which we must calculate the 
   //dot product of the derivatives and y and z
   //Reset these values to zero
   uderiv_dot_y = 0.0; uderiv_dot_z=0.0;
   //Now calculate the dot products of the derivative and the solutions
   //of the linear system
   for(unsigned long l=0;l<n_dofs;l++) 
    {
     uderiv_dot_y += Dof_derivatives[l]*y[l];
     uderiv_dot_z += Dof_derivatives[l]*z[l];
    }
   uderiv_dot_y *= Theta_squared;
   uderiv_dot_z *= Theta_squared;

   //The set the change in the parameter, given by the pseudo-arclength
   //equation. Note that here we are assuming that the arc-length
   //equation is always exactly zero, 
   //which seems to work OK, and saves on some storage. 
   //In fact, it's more subtle than that. If we include this
   //proper residual then we will have to solve the eigenproblem.
   //This will make the solver more robust and *should* be done
   // ... at some point.
   double dparam = (arc_length_constraint_residual - uderiv_dot_y)
    /(Parameter_derivative - uderiv_dot_z);
   
   //Set the new value of the parameter
   *parameter_pt -= dparam;
       
   //Update the values of the other degrees of freedom
   for(unsigned long l=0;l<n_dofs;l++) 
    {
     *Dof_pt[l] -= y[l] - dparam*z[l];
    }

   //Calculate the new residuals
#ifdef OOMPH_HAS_MPI
   // Synchronise the solution on different processors (on each submesh)
   unsigned nmesh=nsub_mesh();
   if (nmesh==0)
    {
     synchronise_dofs(mesh_pt());
    }
   else
    {
     // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
     for (unsigned imesh=0; imesh<nmesh; imesh++)
      {
       synchronise_dofs(mesh_pt(imesh));
      }
     for (unsigned imesh=0; imesh<nmesh; imesh++)
      {
       synchronise_external_dofs(mesh_pt(imesh));
      }
    }
#endif

   // Do any updates that are required 
   actions_after_newton_step();
   actions_before_newton_convergence_check();

   get_residuals(y);
   
   //Get the maximum residuals
   double maxres = y.max();
                           
   //Assemble the residuals for the arc-length step
   arc_length_constraint_residual = 0.0; 
   //Add the variables
   for(unsigned long l=0;l<n_dofs;l++)
    {
     arc_length_constraint_residual +=
      Dof_derivatives[l]*(*Dof_pt[l] - Dofs_current[l]);
    }
   arc_length_constraint_residual *= Theta_squared;
   arc_length_constraint_residual +=
    Parameter_derivative*(*parameter_pt - Parameter_current) 
    - Ds_current;

   //Is it the max
   if(std::abs(arc_length_constraint_residual) > maxres)
    {
     maxres = std::abs(arc_length_constraint_residual);
    }

   //If we have converged jump straight to the test at the end of the loop
   if(maxres < Newton_solver_tolerance) {LOOP_FLAG=0; continue;} 
   
   //This section will not be reached if we have converged already
   //If the maximum number of residuals is too high or the maximum number
   //of iterations has been reached
   if((maxres > Max_residuals) || (count == Max_newton_iterations))
    {
     throw NewtonSolverError(count,maxres);
    }
   
  }
 while(LOOP_FLAG);
 
 //Now update anything that needs updating
 actions_after_newton_solve();

 //Reset the storage of the matrix on the linear solver to what it was
 //on entry to this routine
 if (enable_resolve)
  {
   Linear_solver_pt->enable_resolve();
  }
 else
  {
   Linear_solver_pt->disable_resolve();
  }

 //Return the number of Newton Steps taken
 return count;
}

//=========================================================================
/// A function to calculate the derivatives wrt the arc-length. This version
/// of the function actually does a linear solve so that the derivatives
/// are calculated "exactly" rather than using the values at the Newton
/// step just before convergence. This is only necessary in spatially adaptive
/// problems, in which the number of degrees of freedom changes and so
/// the appropriate derivatives must be calculated for the new variables.
//=========================================================================
void Problem::
calculate_continuation_derivatives(double* const &parameter_pt)
{
 //Find the number of degrees of freedom in the problem
 const unsigned long n_dofs = ndof();

 // create the distribution
 LinearAlgebraDistribution dist(Communicator_pt,n_dofs,false);

 //Assign memory for solutions of the equations
 DoubleVector z(&dist,0.0);

 //If it's the block hopf solver need to solve for both RHS
 //at once, but this would all be alleviated if we have the solve
 //for the non-residuals RHS.
 if(dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt))
  {
   //Get the vector dresiduals/dparameter
   get_derivative_wrt_global_parameter(parameter_pt,z);
   
   
   // Copy rhs vector into local storage so it doesn't get overwritten
   // if the linear solver decides to initialise the solution vector, say,
   // which it's quite entitled to do!
   DoubleVector dummy(&dist,0.0), input_z(z);

   //Solve for the two RHSs
   dynamic_cast<BlockHopfLinearSolver*>(Linear_solver_pt)->
    solve_for_two_rhs(this,dummy,input_z,z);
  }
 //Otherwise we can use the normal resolve
 else
  {
   //Save the status before entry to this routine
   bool enable_resolve = Linear_solver_pt->resolve_is_enabled();
   
   //We need to do resolves
   Linear_solver_pt->enable_resolve();

   //Solve the standard problem, we only want to make sure that
   //we factorise the matrix, if it has not been factorised. We shall
   //ignore the return value of z.
   Linear_solver_pt->solve(this,z);
   
   //Get the vector dresiduals/dparameter
   get_derivative_wrt_global_parameter(parameter_pt,z);
   
   
   // Copy rhs vector into local storage so it doesn't get overwritten
   // if the linear solver decides to initialise the solution vector, say,
   // which it's quite entitled to do!
   DoubleVector input_z(z);
   
   //Now resolve the system with the new RHS and overwrite the solution
   Linear_solver_pt->resolve(input_z,z);
   
   //Restore the storage status of the linear solver
   if (enable_resolve)
    {
     Linear_solver_pt->enable_resolve();
    }
   else 
    {
     Linear_solver_pt->disable_resolve();
    }
  }

 //Now, we can calculate the derivatives, etc
 calculate_continuation_derivatives(z);
}

//=======================================================================
/// A function to calculate the derivatives with respect to the arc-length
/// required for continuation. The arguments is the solution of the 
/// linear system,
/// Jz = dR/dparameter, that gives du/dparameter and the direction 
/// output from the newton_solve_continuation function. The derivatives
/// are stored in the ContinuationParameters namespace.
//===================================================================
void Problem::calculate_continuation_derivatives(const DoubleVector &z)
{

 //Calculate the continuation derivatives
 calculate_continuation_derivatives_helper(z);

 //Scale the value of theta if the control flag is set
 if(Scale_arc_length)
  {
   Theta_squared *= (Parameter_derivative*Parameter_derivative/
                     Desired_proportion_of_arc_length)*
    ((1.0 - Desired_proportion_of_arc_length)/
     (1.0 - Parameter_derivative*Parameter_derivative));

   //Recalculate the continuation derivatives with the new scaled values
   calculate_continuation_derivatives_helper(z);
  }
}

//=======================================================================
/// A private helper function to
/// calculate the derivatives with respect to the arc-length
/// required for continuation. The arguments is the solution of the 
/// linear system,
/// Jz = dR/dparameter, that gives du/dparameter and the direction 
/// output from the newton_solve_continuation function. The derivatives
/// are stored in the ContinuationParameters namespace.
//===================================================================
void Problem::calculate_continuation_derivatives_helper(
 const DoubleVector &z)
{
 //Find the number of degrees of freedom in the problem
 unsigned long n_dofs = ndof();

 //Work out the continuation direction
 //The idea is that (du/ds)_{old} . (du/ds)_{new} >= 0 
 //if the direction is to remain the same.
 //du/ds_{new} = [dlambda/ds; du/ds] = [dlambda/ds ; - dlambda/ds z]
 //so (du/ds)_{new} . (du/ds)_{old} 
 // = dlambda/ds [1 ; - z] . [ Parameter_derivative ; Dof_derivaives]
 // = dlambda/ds (Parameter_derivative - Dof_derivatives . z)
 Continuation_direction = Parameter_derivative;
 for(unsigned long l=0;l<n_dofs;l++) 
  {Continuation_direction -= Dof_derivatives[l]*z[l];}

 //Calculate the magnitude of the du/ds Vector

 //Note that actually, we are usually approximating by using the value at 
 //newton step just before convergence, which saves one additional 
 //Newton solve.
 
 //First calculate the magnitude of du/dparameter, chi
 double chi = 0.0; 
 for(unsigned long l=0;l<n_dofs;l++) {chi += z[l]*z[l];}
                                                   
 //Calculate the current derivative of the parameter wrt the arc-length
 Parameter_derivative = 1.0/sqrt(1.0 + Theta_squared*chi);

 //If the dot product of the current derivative wrt the Direction
 //is less than zero, switch the sign of the derivative to ensure 
 //smooth continuation
 if(Parameter_derivative*Continuation_direction < 0.0) 
  {Parameter_derivative*= -1.0;}

 //Resize the derivatives array, if necessary
 if(Dof_derivatives.size() != n_dofs) {Dof_derivatives.resize(n_dofs,0.0);}
 //Calculate the new derivatives wrt the arc-length
 for(unsigned long l=0;l<n_dofs;l++)
  {
   //This comes from the formulation J u_dot + dr/dlambda  lambda_dot = 0
   //on the curve and then it follows that.
   Dof_derivatives[l] = -Parameter_derivative*z[l];
  }
}


/// \short Virtual function that is used to symmetrise the problem so that
/// the current solution exactly satisfies any symmetries within the system.
/// Used when adpativly solving pitchfork detection problems when small 
/// asymmetries in the coarse solution can be magnified
/// leading to very inaccurate answers on the fine mesh. 
/// This is always problem-specific and must be filled in by the user
/// The default issues a warning
void Problem::symmetrise_eigenfunction_for_adaptive_pitchfork_tracking()
{
 std::ostringstream warn_message;
 warn_message
  << "Warning: This function is called after spatially adapting the\n"
  << "eigenfunction associated with a pitchfork bifurcation and should\n"
  << "ensure that the exact (anti-)symmetries of problem are enforced\n"
  << "within that eigenfunction. It is problem specific and must be\n"
  << "filled in by the user if required.\n"
  << "A sign of problems is if the slack paramter gets too large and\n"
  << "if the solution at the Pitchfork is not symmetric.\n";
 OomphLibWarning(
  warn_message.str(),
  "Problem::symmetrise_eigenfunction_for_adaptive_pitchfork_tracking()",
  OOMPH_EXCEPTION_LOCATION);
}

//====================================================================
///Return pointer to the parameter that is used in the
/// bifurcation detection. If we are not tracking a bifurcation then
/// an error will be thrown by the AssemblyHandler
//====================================================================
double* Problem::bifurcation_parameter_pt() const 
{return Assembly_handler_pt->bifurcation_parameter_pt();}

//====================================================================
/// Return the eigenfunction calculated as part of a
/// bifurcation tracking process. If we are not tracking a bifurcation
/// then an error will be thrown by the AssemblyHandler
//======================================================================
void Problem::get_bifurcation_eigenfunction(
 Vector<DoubleVector> &eigenfunction)
{
 //Simply call the appropriate assembly handler function
 Assembly_handler_pt->get_eigenfunction(eigenfunction);
}

//============================================================
/// Activate the fold tracking system by changing the assembly
/// handler and initialising it using the parameter addressed 
/// by parameter_pt.
//============================================================
void Problem::activate_fold_tracking(double* const &parameter_pt,
                                     const bool &block_solve)
{
 //Reset the assembly handler to default
 reset_assembly_handler_to_default();
 //Set the new assembly handler. Note that the constructor actually
 //solves the original problem to get some initial conditions, but
 //this is OK because the RHS is always evaluated before assignment.
 Assembly_handler_pt = new FoldHandler(this,parameter_pt);

 //If we are using a block solver, we must set the linear solver pointer
 //to the block fold solver. The present linear solver is
 //used by the block solver and so must be passed as an argument.
 //The destructor of the Fold handler returns the linear
 //solver to the original non-block version.
 if(block_solve)
  {
   Linear_solver_pt = new AugmentedBlockFoldLinearSolver(Linear_solver_pt);
  }
}


//==================================================================
/// Activate the pitchfork tracking system by changing the assembly
/// handler and initialising it using the parameter addressed 
/// by parameter_pt and a symmetry vector. The boolean flag is
/// used to specify whether a block solver is used, default is true.
//===================================================================
void Problem::activate_pitchfork_tracking(
 double* const &parameter_pt,
 const DoubleVector &symmetry_vector,const bool &block_solve)
{

 //Reset the assembly handler to default
 reset_assembly_handler_to_default();
 //Set the new assembly handler. Note that the constructor actually
 //solves the original problem to get some initial conditions, but
 //this is OK because the RHS is always evaluated before assignment.
 Assembly_handler_pt = new PitchForkHandler(this,parameter_pt,
                                            symmetry_vector);

 //If we are using a block solver, we must set the linear solver pointer
 //to the block pitchfork solver. The present linear solver is
 //used by the block solver and so must be passed as an argument.
 //The destructor of the PitchFork handler returns the linear
 //solver to the original non-block version.
 if(block_solve)
  {
   Linear_solver_pt = 
    new BlockPitchForkLinearSolver(Linear_solver_pt);
  }
}



//============================================================
/// Activate the hopf tracking system by changing the assembly
/// handler and initialising it using the parameter addressed 
/// by parameter_pt.
//============================================================
void Problem::activate_hopf_tracking(
 double* const &parameter_pt, const bool &block_solve)
{
 //Reset the assembly handler to default
 reset_assembly_handler_to_default();
 //Set the new assembly handler. Note that the constructor actually
 //solves the original problem to get some initial conditions, but
 //this is OK because the RHS is always evaluated before assignment.
 Assembly_handler_pt = new HopfHandler(this,parameter_pt);

 //If we are using a block solver, we must set the linear solver pointer
 //to the block hopf solver. The present linear solver is
 //used by the block solver and so must be passed as an argument.
 //The destructor of the Hopf handler returns the linear
 //solver to the original non-block version.
 if(block_solve)
  {
   Linear_solver_pt = new BlockHopfLinearSolver(Linear_solver_pt);
  }
}


//============================================================
/// Activate the hopf tracking system by changing the assembly
/// handler and initialising it using the parameter addressed 
/// by parameter_pt and the frequency and null vectors
/// specified.
//============================================================
void Problem::activate_hopf_tracking(
 double* const &parameter_pt, const double &omega,
 const DoubleVector &null_real, const DoubleVector &null_imag,
 const bool &block_solve)
{
 //Reset the assembly handler to default
 reset_assembly_handler_to_default();
 //Set the new assembly handler. Note that the constructor actually
 //solves the original problem to get some initial conditions, but
 //this is OK because the RHS is always evaluated before assignment.
 Assembly_handler_pt = new HopfHandler(this,parameter_pt,omega,
                                       null_real,null_imag);

 //If we are using a block solver, we must set the linear solver pointer
 //to the block hopf solver. The present linear solver is
 //used by the block solver and so must be passed as an argument.
 //The destructor of the Hopf handler returns the linear
 //solver to the original non-block version.
 if(block_solve)
  {
   Linear_solver_pt = new BlockHopfLinearSolver(Linear_solver_pt);
  }
}

 
//===============================================================
///Reset the assembly handler to default
//===============================================================
void Problem::reset_assembly_handler_to_default()
{
 //If we have a non-default handler
 if(Assembly_handler_pt != Default_assembly_handler_pt) 
  {
   //Delete the current assembly handler 
   delete Assembly_handler_pt; 
   //Reset the assembly handler
   Assembly_handler_pt = Default_assembly_handler_pt;
  }
}

//===================================================================
/// This function takes one step of length ds in pseudo-arclength.The 
/// argument parameter_pt is a pointer to the parameter (global variable) 
/// that is being traded for arc-length. The function returns the next desired
/// arc-length according to criteria based upon the desired number of Newton 
/// Iterations per solve.
//=====================================================================
double Problem::arc_length_step_solve(double* const &parameter_pt,
                                      const double &ds)
{

 //----------------------MAKE THE PROBLEM STEADY-----------------------
 //Loop over the timesteppers and make them (temporarily) steady.
 //We can only do continuation for steady problems!
 unsigned n_time_steppers = ntime_stepper();
 for(unsigned i=0;i<n_time_steppers;i++) 
  {
   time_stepper_pt(i)->make_steady();
  }
 
 //----------SAVE THE INITIAL VALUES, IN CASE THE STEP FAILS-----------
 //Find total number of dofs
 unsigned long n_dofs = ndof();
 //Safety check, set up the array of Dof_derivatives, if necessary
 if(Dof_derivatives.size() != n_dofs) {Dof_derivatives.resize(n_dofs,0.0);}
 //Save the current value of the parameter
 Parameter_current = *parameter_pt;

 //Save the current values of the degrees of freedom 
 //Safety check, set up the array of Dof_derivatives, if necessary
 if(Dofs_current.size() != n_dofs) {Dofs_current.resize(n_dofs);}
 for(unsigned long l=0;l<n_dofs;l++) {Dofs_current[l] = *Dof_pt[l];}
 //Set the value of ds_current
 Ds_current = ds;

 //----SET UP MEMORY FOR QUANTITIES THAT ARE REQUIRED OUTSIDE THE LOOP----

 //Assign memory for solutions of the equations Jz = du/dparameter
 //This is needed here (outside the loop), so that we can save on 
 //one linear solve when calculating the derivatives wrt the arc-length
 LinearAlgebraDistribution dist(Communicator_pt,n_dofs,false);
 DoubleVector z(&dist,0.0);

 //Store sign of the Jacobian, used for bifurcation detection
 //If this is the first time that we are calling the arc-length solver,
 //this should not be used.
 int previous_sign = Sign_of_jacobian;

 //Counter for the number of newton steps
 unsigned count=0;
 //Flag to indicate a successful step
 bool STEP_REJECTED=false;
 
 //Flag to indicate a sign change
 bool SIGN_CHANGE=false;

 //Loop around the step in arc-length
 do
  {
   //Check that the step has not fallen below the minimum tolerance
   if(std::abs(Ds_current) < Minimum_ds)
    {
     std::ostringstream error_message;
     error_message << "DESIRED ARC-LENGTH STEP " << Ds_current 
                   << " HAS FALLEN BELOW MINIMUM TOLERANCE, " 
                   << Minimum_ds << std::endl;
     
     throw OomphLibError(error_message.str(),
                         "Problem::arc_length_step_solve()",
                         OOMPH_EXCEPTION_LOCATION);
    }

   //Assume that we shall accept the step
   STEP_REJECTED=false;
   //Set initial value of the parameter
   *parameter_pt += Parameter_derivative*Ds_current;
   //Loop over the variables and set their initial values
   for(unsigned long l=0;l<n_dofs;l++) 
    {
     *Dof_pt[l] += Dof_derivatives[l]*Ds_current;
    }
   
   //Actually do the newton solve stage for the continuation problem
   try
    {
     count = newton_solve_continuation(parameter_pt,z);
    } 
   //Catch any exceptions thrown in the Newton solver
   catch(NewtonSolverError &error)
    {
     //Check whether it's the linear solver
     if(error.linear_solver_error)
      {
       std::ostringstream error_stream;
       error_stream << std::endl 
                    << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
       oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
       throw OomphLibError(error_stream.str(),
                           "Problem::arc_length_step_solve()",
                           OOMPH_EXCEPTION_LOCATION);
      }
     //Otherwise mark the step as having failed
     else 
      {
       oomph_info << "STEP REJECTED --- TRYING AGAIN" << std::endl;
       STEP_REJECTED=true;
       //Let's take a smaller step
       Ds_current *= (2.0/3.0);
       //Reset the dofs and parameter
       *parameter_pt = Parameter_current;
       for(unsigned long l=0;l<n_dofs;l++) {*Dof_pt[l] = Dofs_current[l];}
      }
    }
  }
 while(STEP_REJECTED); //continue until a step is accepted

 //Only recalculate the derivatives if there has been a Newton solve
 //If not, the previous values should be close enough 
 if(count>0)
  {

   //--------------------CHECK FOR POTENTIAL BIFURCATIONS-------------
   //If the sign of the jacobian is zero issue a warning
   if(Sign_of_jacobian == 0) 
    {
     std::string error_message =
      "The sign of the jacobian is zero after a linear solve\n";
     error_message +=
      "Either the matrix is singular (unlikely),\n";
     error_message +=
      "or the linear solver cannot compute the determinant of the matrix;\n";
     error_message += "e.g. an iterative linear solver.\n";
     error_message +=
      "If the latter, bifurcation detection must be via an eigensolver\n";
     OomphLibWarning(error_message,
                     "Problem::arc_length_step_solve",
                     OOMPH_EXCEPTION_LOCATION);
    }
   //If this is the first step, we cannot rely on the previous value 
   //of the jacobian so set the previous sign to the present sign
   if(!Arc_length_step_taken) {previous_sign = Sign_of_jacobian;}
   //If we have detected a sign change in the last converged Jacobian,
   //it must be a turning point or bifurcation
   if(Sign_of_jacobian != previous_sign)
    {

     //There has been, at least, one sign change
     First_jacobian_sign_change = true;

     //The sign has changed this time
     SIGN_CHANGE=true;

     //Calculate the dot product of the approximate null vector
     //of the Jacobian matrix ((badly) approximated by z)
     //and the vectors of derivatives of the residuals wrt the global parameter
     //If this is small it is a bifurcation rather than a turning point.
     //Get the derivative wrt global parameter
     DoubleVector dparam;
     get_derivative_wrt_global_parameter(parameter_pt,dparam);
     //Calculate the dot product
     double dot=0.0;
     for(unsigned long n=0;n<n_dofs;++n) {dot += dparam[n]*z[n];}

     //Write the output message
     std::ostringstream message;
     message << "-----------------------------------------------------------";
     message << std::endl << "SIGN CHANGE IN DETERMINANT OF JACOBIAN: " 
             << std::endl;
     message << "BIFURCATION OR TURNING POINT DETECTED BETWEEN "
             << Parameter_current << " AND " << *parameter_pt << std::endl;
     message << "APPROXIMATE DOT PRODUCT : " << dot << "," << std::endl;
     message << "IF CLOSE TO ZERO WE HAVE A BIFURCATION; ";
     message << "OTHERWISE A TURNING POINT" << std::endl;
     message << "-----------------------------------------------------------"
             << std::endl;

     //Write the message to standard output
     oomph_info << message.str();

     //Open the information file for appending
     std::ofstream bifurcation_info("bifurcation_info",std::ios_base::app);
     //Write the message to the file
     bifurcation_info << message.str();
     bifurcation_info.close();
    }
   
   //Calculate the derivatives required for the next stage of continuation
   //In this we pass the last value of z (i.e. approximation)
   calculate_continuation_derivatives(z); 

   //If it's the first step then the value of the next step should
   //be the change in parameter divided by the parameter derivative
   //to obtain approximately the same parameter change
   if(!Arc_length_step_taken) 
    {
     Ds_current = (*parameter_pt - Parameter_current)/Parameter_derivative;
    }
   
   //We have taken our first step
   Arc_length_step_taken = true;
  }
 //Otherwise calculate the continuation derivatives by solving the linear
 //system. We must do this to ensure that the derivatives are in sync
 //It could lead to problems near turning points when we should really be
 //solving an eigenproblem, but seems OK so far!
 else
  {

   //Save the current sign of the jacobian
   int temp_sign=Sign_of_jacobian;
   //Calculate the continuation derivatives, which includes a solve
   //of the linear system
   calculate_continuation_derivatives(parameter_pt);

   //Reset the sign of the jacobian, just in case the sign has changed when
   //solving the continuation derivatives. The sign change will be picked
   //up on the next continuation step.
   Sign_of_jacobian = temp_sign;
  }

 /*{
   Vector<double> z(n_dofs);
   //Cheeky tester
   double length=0.0;
   for(unsigned long l=0;l<n_dofs;l++)
   {
   z[l] = (*Dof_pt[l] - Dofs_current[l])/Ds_current;
   length += Theta_squared*z[l]*z[l];
   }

   double Z = (*parameter_pt - Parameter_current)/Ds_current;
   length += Z*Z;

   length = sqrt(length);
   for(unsigned long l=0;l<n_dofs;l++)
   {
   Dof_derivatives[l] = z[l]/length;
   }
  
   Parameter_derivative = Z/length;
   } */

 //If we are trying to find a bifurcation and the first sign change
 //has occured, use bisection
 if((Bifurcation_detection) && (First_jacobian_sign_change))
  {
   //If there has been a sign change we need to half the step size
   //and reverse the direction
   if(SIGN_CHANGE) {Ds_current *= -0.5;}
   //Otherwise
   else
    {
     //The size of the bracketed interval is always
     //2ds - Ds_current (this will work even if the original step failed)
     //We want our new step size to be half this
     Ds_current = ds - 0.5*Ds_current;
    }
   //Return the desired value of the step
   return Ds_current;
  }

 //If fewer than the desired number of Newton Iterations, increase the step
 if(count < Desired_newton_iterations_ds) {return Ds_current*1.5;}
 //If more than the desired number of Newton Iterations, reduce the step
 if(count > Desired_newton_iterations_ds) {return Ds_current*(2.0/3.0);}
 //Otherwise return the step just taken
 return Ds_current;
}


//=======================================================================
/// Take an explicit timestep of size dt
//======================================================================
void Problem::explicit_timestep(const double &dt, const bool &shift_values)
{
 //Firstly we shift the time values
 if(shift_values) {shift_time_values();}
 //Set the current value of dt, if we can
 if(time_pt()->ndt() > 0) {time_pt()->dt() = dt;}
 
 //Make all the implicit timestepper steady
 //Find out how many timesteppers there are
 unsigned n_time_steppers = ntime_stepper();
 
 // Vector of bools to store the is_steady status of the various
 // timesteppers when we came in here
 std::vector<bool> was_steady(n_time_steppers);
 
 //Loop over them all and make them (temporarily) static
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   was_steady[i]=time_stepper_pt(i)->is_steady();
   time_stepper_pt(i)->make_steady();
  }
 
 //Take the explicit step
 this->explicit_time_stepper_pt()->timestep(this,dt);
 
 // Reset the is_steady status of all timesteppers that
 // weren't already steady when we came in here and reset their
 // weights
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   if (!was_steady[i])
    {
     time_stepper_pt(i)->undo_make_steady();
    }
  }
}


//========================================================================
/// Do one timestep of size dt using Newton's method with the specified 
/// tolerance and linear solver defined as member data of the Problem class.
/// This will be the most commonly used version 
/// of  unsteady_newton_solve, in which the time values are always shifted
/// This does not include any kind of adaptativity. If the solution fails to 
/// converge the program will end.
//========================================================================
void Problem::unsteady_newton_solve(const double &dt)
{
 //We shift the values, so shift_values is true
 unsteady_newton_solve(dt,true);
}

//========================================================================
/// Do one timestep forward of size dt using Newton's method with the 
/// specified tolerance and linear solver defined via member data of the 
/// Problem class.
/// The boolean flag shift_values is used to control whether the time values
/// should be shifted or not. 
//========================================================================
void Problem::unsteady_newton_solve(const double &dt, const bool &shift_values)
{
 //Shift the time values and the dts, according to the control flag
 if(shift_values) {shift_time_values();}
 
 // Advance global time and set current value of dt 
 time_pt()->time()+=dt;
 time_pt()->dt()=dt;

 //Find out how many timesteppers there are
 unsigned n_time_steppers = ntime_stepper();

 //Loop over them all and set the weights
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   time_stepper_pt(i)->set_weights();
  }

 //Now update anything that needs updating before the timestep
 //This could be time-dependent boundary conditions, for example.
 actions_before_implicit_timestep();

 try
  {
   //Solve the non-linear problem for this timestep with Newton's method
   newton_solve();
  }
 //Catch any exceptions thrown in the Newton solver
 catch(NewtonSolverError &error)
  {
   oomph_info << std::endl << "USER-DEFINED ERROR IN NEWTON SOLVER " 
              << std::endl;
   //Check whether it's the linear solver
   if(error.linear_solver_error)
    {
     oomph_info << "ERROR IN THE LINEAR SOLVER" << std::endl;
    }
   //Check to see whether we have reached Max_iterations
   else if(error.iterations==Max_newton_iterations)
    {
     oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" << error.iterations 
                << ") REACHED WITHOUT CONVERGENCE " << std::endl;
    }
   //If not, it must be that we have exceeded the maximum residuals
   else
    {
     oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                << " EXCEEDS PREDEFINED MAXIMUM " << Max_residuals
                << std::endl;
    }
   //Die horribly!!
   std::ostringstream error_stream;
   error_stream << "Error occured in unsteady Newton solver. "
                << std::endl;
   throw OomphLibError(error_stream.str(),
                       "Problem::unsteady_newton_solve()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Now update anything that needs updating after the timestep
 actions_after_implicit_timestep();
}

//=======================================================================
/// Attempt to take one timestep forward using dt_desired. The error control
/// parameter, epsilon, is used to specify the desired approximate value of the
/// global error norm per timestep. The routine returns the value an estimate
/// of the next value of dt that should be taken. 
//=======================================================================
double Problem::
adaptive_unsteady_newton_solve(const double &dt_desired,
                               const double &epsilon)
{
 //We always want to shift the time values
 return adaptive_unsteady_newton_solve(dt_desired,epsilon,true);
}


//=======================================================================
/// Attempt to take  one timestep forward using the dt_desired.
/// This is the driver for a number of adaptive solvers. If the solution
/// fails to converge at a given timestep, the routine will automatically
/// halve the time step and try again, until the time step falls below the
/// specified minimum value. The routine returns the value an estimate
/// of the next value of dt that should be taken. 
//========================================================================
double Problem::
adaptive_unsteady_newton_solve(const double &dt_desired,
                               const double &epsilon, 
                               const bool &shift_values)
{
 //First, we need to backup the existing dofs, in case the timestep is 
 //rejected 
 //Find total number of dofs
 unsigned long n_dofs = ndof();
 //Now set up a Vector to hold current values
 Vector<double> dofs_current(n_dofs);
 //Load values into dofs_current
 for(unsigned i=0;i<n_dofs;i++) dofs_current[i] = *Dof_pt[i];
 //Store the time
 double time_current = time_pt()->time();

 //Flag to detect whether the timestep has been rejected or not
 unsigned REJECT_TIMESTEP=0;
 //Flag to detect whether any of the timesteppers are adaptive
 unsigned ADAPTIVE_FLAG=0;
 //The value of the actual timestep, by default the same as desired timestep
 double dt_actual=dt_desired;
 //Timestep rescaling factor, 1.0 by default
 double DTSF = 1.0;
 
 //Determine the number of timesteppers
 unsigned n_time_steppers = ntime_stepper();
 //Find out whether any of the timesteppers are adaptive
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   if(time_stepper_pt(i)->adaptive_flag())
    {
     ADAPTIVE_FLAG=1; 
     break;
    }
  }

 //Shift the time_values according to the control flag
 if(shift_values) {shift_time_values();}

 //This loop surrounds the adaptive time-stepping critera
 do
  {
   //This loop surrounds the Newton solver and will not
   //be broken until a timestep is accepted
   do
    {
     //Initially the timestep is presumed to be accepted
     REJECT_TIMESTEP=0;

     //Set the new time and value of dt
     time_pt()->time() += dt_actual;
     time_pt()->dt() = dt_actual;

     //Loop over all timesteppers and set the weights and predictor weights
     for(unsigned i=0;i<n_time_steppers;i++)
      {
       //If the time_stepper is non-adaptive, this will be zero
       time_stepper_pt(i)->set_predictor_weights();
       time_stepper_pt(i)->set_weights();
      }

     //Now calculate the predicted values for the all data and all positions
     calculate_predictions();

     //Do any updates/boundary conditions changes here
     actions_before_implicit_timestep();
   
     //Attempt to solver the non-linear system
     try
      {
       //Solve the non-linear problem at this timestep
       newton_solve();
      }
     //Catch any exceptions thrown
     catch(NewtonSolverError &error)
      {
       //If it's a solver error then die
       if(error.linear_solver_error)
        {
         std::string error_message =
          "USER-DEFINED ERROR IN NEWTON SOLVER\n";
         error_message +=  "ERROR IN THE LINEAR SOLVER\n";

         //Die
         throw OomphLibError(error_message,
                             "Problem::adaptive_unsteady_newton_solve()",
                             OOMPH_EXCEPTION_LOCATION);
        }
       else
        {
         oomph_info << "TIMESTEP REJECTED --- HALVING TIMESTEP AND TRYING AGAIN" 
                    << std::endl;
         //Reject the timestep, if we have an exception
         REJECT_TIMESTEP=1;
         //Essentially all I do here is half the next timestep
         dt_actual *= 0.5;
         //Reset the time
         time_pt()->time() = time_current;
         //Reload the dofs
         for(unsigned i=0;i<n_dofs;i++) *Dof_pt[i] = dofs_current[i];

#ifdef OOMPH_HAS_MPI
         // Synchronise the solution on different processors (on each submesh)
         unsigned nmesh=nsub_mesh();
         if (nmesh==0)
          {
           synchronise_dofs(mesh_pt());
          }
         else
          {
           // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
           for (unsigned imesh=0; imesh<nmesh; imesh++)
            {
             synchronise_dofs(mesh_pt(imesh));
            }
           for (unsigned imesh=0; imesh<nmesh; imesh++)
            {
             synchronise_external_dofs(mesh_pt(imesh));
            }
          }
#endif
         //Call all "after" actions, e.g. to handle mesh updates
         actions_after_newton_step();
         actions_before_newton_convergence_check();
         actions_after_newton_solve();
         actions_after_implicit_timestep();
         //Skip to the next iteration
         continue;
        }
      }

     //Break out of the loop if the timestep has become too small
     if(dt_actual < Minimum_dt)
      {
       std::ostringstream error_message;
       error_message 
        << "TIMESTEP (" << dt_actual 
        << ") HAS FALLEN BELOW SPECIFIED THRESHOLD: Problem::Minimum_dt=" 
        << Minimum_dt << std::endl;

       throw OomphLibError(error_message.str(),
                           "Problem::adaptive_unsteady_newton_solve()",
                           OOMPH_EXCEPTION_LOCATION);
      }

     //Update anything that needs updating after the timestep
     actions_after_implicit_timestep();
    }
   //Keep looping until we accept the timestep
   while(REJECT_TIMESTEP);
  
   //If we have an adapative timestepper
   if(ADAPTIVE_FLAG)
    {
     //Once timestep has been accepted can do fancy error processing
     //Set the error weights
     for(unsigned i=0;i<n_time_steppers;i++)
      {
       time_stepper_pt(i)->set_error_weights();
      }

     //Call a global error, at the moment I'm just going to use a square norm 
     double error = global_temporal_error_norm(); 

     //Calculate the scaling  factor
     DTSF = pow((epsilon/error),
                (1.0/(1.0+time_stepper_pt()->order())));

     oomph_info << "DTSF is  " << DTSF << std::endl;
     oomph_info << "Estimated timestepping error is " << error << std::endl;

     //Now decide what to do based upon DTSF
     //If it's small reject the timestep
     if(DTSF <= 0.8)
      {
       oomph_info << "TIMESTEP REJECTED" << std::endl;
       //Reject the timestep
       REJECT_TIMESTEP=1;
       //Modify the actual timestep
       dt_actual *= DTSF;
       //Reset the time
       time_pt()->time() = time_current;
       //Reload the dofs
       for(unsigned i=0;i<n_dofs;i++) *Dof_pt[i] = dofs_current[i];

#ifdef OOMPH_HAS_MPI
       // Synchronise the solution on different processors (on each submesh)
       unsigned nmesh=nsub_mesh();
       if (nmesh==0)
        {
         synchronise_dofs(mesh_pt());
        }
       else
        {
          // Synchronise ALL halo(ed) dofs BEFORE external halo(ed) dofs
         for (unsigned imesh=0; imesh<nmesh; imesh++)
          {
           synchronise_dofs(mesh_pt(imesh));
          }
         for (unsigned imesh=0; imesh<nmesh; imesh++)
          {
           synchronise_external_dofs(mesh_pt(imesh));
          }
        }
#endif

       //Call all "after" actions, e.g. to handle mesh updates
       actions_after_newton_step();
       actions_before_newton_convergence_check();
       actions_after_newton_solve();
       actions_after_implicit_timestep();

       continue;
      }
     //If it's large change the timestep
     if(DTSF >= 1.0)
      {
       //Restrict the increase
       if(DTSF > DTSF_max_increase)
        {
         DTSF = DTSF_max_increase;
         oomph_info << "DTSF LIMITED TO " << DTSF_max_increase << std::endl;
        }
      }

    } //End of if adaptive flag

  }
 //Keep this loop going, again until we accept the timestep
 while(REJECT_TIMESTEP);


 //Make sure timestep doesn't get too large
 if ((dt_actual*DTSF) > Maximum_dt)
  {
   oomph_info << "DTSF WOULD INCREASE TIMESTEP "
              << "ABOVE SPECIFIED THRESHOLD: Problem::Maximum_dt=" 
              <<  Maximum_dt << std::endl;
   DTSF =  Maximum_dt/dt_actual;
   oomph_info << "ADJUSTING DTSF TO " << DTSF << std::endl;
  }
 

 //Once the timestep has been accepted, return the actual timestep taken, 
 //suitably scaled, to be used the next time
 return (dt_actual*DTSF);
}



//=======================================================================
///  Unsteady "doubly" adaptive Newton solve: Does temporal
/// adaptation first, i.e. we try to do a timestep with an increment
/// of dt, and adjusting dt until the solution on the given mesh satisfies
/// the temporal error measure with tolerance epsilon. Following
/// this, we do up to max_adapt spatial adaptions (without 
/// re-examining the temporal error). If first==true, the initial conditions
/// are re-assigned after the mesh adaptations.
/// Shifting of time can be suppressed by overwriting the
/// default value of shift (true). [Shifting must be done
/// if first_timestep==true because we're constantly re-assigning
/// the initial conditions; if first_timestep==true and shift==false
/// shifting is performed anyway and a warning is issued.
//========================================================================
double Problem::doubly_adaptive_unsteady_newton_solve(const double &dt_desired,
                                                      const double &epsilon, 
                                                      const unsigned &max_adapt,
                                                      const bool &first,
                                                      const bool &shift_values)
{
 //Store the initial time
 double initial_time = time_pt()->time();

 // Take adaptive timestep, adjusting dt until tolerance is satisfied
 double new_dt=adaptive_unsteady_newton_solve(dt_desired, 
                                              epsilon,
                                              shift_values);
 double dt_taken=time_pt()->dt();
 oomph_info << "Accepted solution taken with timestep: " 
            << dt_taken << std::endl;

 // Adapt problem/mesh
 unsigned n_refined=0;
 unsigned n_unrefined=0;
 adapt(n_refined,n_unrefined); 
 
 // Re-solve the problem if the adaptation has changed anything
 if ((n_refined  !=0)||
     (n_unrefined!=0))
  {
   oomph_info << "Mesh was adapted --> we'll re-solve for current timestep." 
              << std::endl;

   // Reset time to what it was when we entered here
   // because it will be incremented again by dt_taken.
   time_pt()->time()=initial_time;

   // Shift the timesteps? No! They've been shifted already when we
   // called the solve with pure temporal adaptivity...
   bool shift=false;

   // Reset the inital condition on refined meshes
   if (first) 
    {
     //Reset the initial conditions
     oomph_info << "Re-assigning initial condition at time=" 
                << time_pt()->time()<< std::endl;
     set_initial_condition();

     // This is the first timestep so shifting
     // has to be done following the assignment of initial conditions.
     // In fact, unsteady_newton_solve(...) does that automatically.
     // We're changing the flag here to avoid warning messages.
     shift=true;
    }

   // Now take the step again on the refined mesh, using the same
   // timestep as used before.
   unsteady_newton_solve(dt_taken,max_adapt,first,shift);
  }
 else
  {
   oomph_info << "Mesh wasn't adapted --> we'll accept spatial refinement." 
              << std::endl;
  }

 return new_dt;

}




//========================================================================
/// \short Initialise the previous values of the variables for time stepping
/// corresponding to an impulsive start. Previous history for all data
/// is generated by the appropriate timesteppers. Previous nodal 
/// positions are simply copied backwards.
//========================================================================
void Problem::assign_initial_values_impulsive()
{
 //Assign the impulsive values in the "master" mesh
 Mesh_pt->assign_initial_values_impulsive();

 // Loop over global data
 unsigned Nglobal=Global_data_pt.size();
 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   Global_data_pt[iglobal]->time_stepper_pt()->
    assign_initial_values_impulsive(Global_data_pt[iglobal]);
  }
}


//=======================================================================
///Assign the values for an impulsive start and also set the initial
///values of the previous dts to both be dt
//======================================================================
void Problem::assign_initial_values_impulsive(const double &dt)
{
 //First initialise the dts and set the weights
 initialise_dt(dt);
 //Now call assign_initial_values_impulsive
 assign_initial_values_impulsive();
}

//=======================================================================
/// Return the current value of continuous time. If not Time object
/// has been assigned, then throw an error
//======================================================================
double& Problem::time()
{
 if(Time_pt==0) 
  {
   throw OomphLibError("Time object has not been set",
                       "Problem::time()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 else {return Time_pt->time();}
}

//========================================================================
/// Shift all time-dependent data along for next timestep.
//========================================================================
void Problem::shift_time_values()
{
 //Move the values of dt in the Time object
 Time_pt->shift_dt(); 

 //Only shift time values in the "master" mesh, otherwise things will
 //get shifted twice in complex problems
 Mesh_pt->shift_time_values();

 // Shift global data with their own timesteppers
 unsigned Nglobal=Global_data_pt.size();
 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   Global_data_pt[iglobal]->time_stepper_pt()->
    shift_time_values(Global_data_pt[iglobal]);
  }

}


//========================================================================
/// Calculate the predictions of all variables in problem
//========================================================================
void Problem::calculate_predictions()
{
 //Calculate all predictions in the "master" mesh
 Mesh_pt->calculate_predictions();

 //Calculate predictions for global data with their own timesteppers
 unsigned Nglobal=Global_data_pt.size();
 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   Global_data_pt[iglobal]->time_stepper_pt()->
    calculate_predicted_values(Global_data_pt[iglobal]);
  }

}



//======================================================================
///\short Enable recycling of the mass matrix in explicit timestepping
///schemes. Useful for timestepping on fixed meshes when you want
///to avoid the linear solve phase.
//=====================================================================
void Problem::enable_mass_matrix_reuse()
{
 Mass_matrix_reuse_is_enabled=true;
 Mass_matrix_has_been_computed=false;
 
 //If we have a discontinuous formulation set the elements to reuse
 //their own mass matrices
 if(Discontinuous_element_formulation)
  {
   const unsigned n_element = Problem::mesh_pt()->nelement();
   //Loop over the other elements
   for(unsigned e=0;e<n_element;e++)
    {
     //Cache the element
     DGElement* const elem_pt =
      dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));
     elem_pt->enable_mass_matrix_reuse();
    }
  }
}
 
//======================================================================
/// \short Turn off the recyling of the mass matrix in explicit 
/// time-stepping schemes
//======================================================================
void Problem::disable_mass_matrix_reuse()
{
 Mass_matrix_reuse_is_enabled=false;
 Mass_matrix_has_been_computed=false;
 
 //If we have a discontinuous formulation set the element-level
 //function
 if(Discontinuous_element_formulation)
  {
   const unsigned n_element = Problem::mesh_pt()->nelement();
   //Loop over the other elements
   for(unsigned e=0;e<n_element;e++)
    {
     //Cache the element
     DGElement* const elem_pt =
      dynamic_cast<DGElement*>(Problem::mesh_pt()->element_pt(e));
     elem_pt->disable_mass_matrix_reuse();
    }
  }
}


//=========================================================================
/// Copy Data values, nodal positions etc from specified problem.
/// Note: This is not a copy constructor. We assume that the current
/// and the "original" problem have both been created by calling
/// the same problem constructor so that all Data objects,
/// time steppers etc. in the two problems are completely independent.
/// This function copies the nodal, internal and global values
/// and the time parameters from the original problem into "this"
/// one. This functionality is required, e.g. for
/// multigrid computations.
//=========================================================================
void Problem::copy(Problem* orig_problem_pt)
{
  
 // Copy time
 //----------

 // Flag to indicate that orig problem is unsteady problem
 bool unsteady_flag=(orig_problem_pt->time_pt()!=0);

 // Copy current time and previous time increments for proper unsteady run
 if (unsteady_flag)
  {
   oomph_info << "Copying an unsteady problem." << std::endl;
   // Current time
   this->time_pt()->time()=orig_problem_pt->time_pt()->time();
   // Timesteps
   unsigned n_dt=orig_problem_pt->time_pt()->ndt();
   time_pt()->resize(n_dt);
   for (unsigned i=0;i<n_dt;i++)
    { 
     time_pt()->dt(i)=orig_problem_pt->time_pt()->dt(i);
    }

   //Find out how many timesteppers there are
   unsigned n_time_steppers = ntime_stepper();
   
   //Loop over them all and set the weights
   for(unsigned i=0;i<n_time_steppers;i++)
    {
     time_stepper_pt(i)->set_weights();
    }


  }

 // Copy nodes
 //-----------

 // Loop over submeshes:
 unsigned nmesh=nsub_mesh();
 if (nmesh==0) nmesh=1;
 for (unsigned m=0;m<nmesh;m++)
  {
   // Find number of nodes in present mesh
   unsigned long n_node = mesh_pt(m)->nnode(); 
   
   // Check # of nodes: 
   unsigned long n_node_orig=orig_problem_pt->mesh_pt(m)->nnode();
   if (n_node!=n_node_orig)
    {
     std::ostringstream error_message;
     error_message << "Number of nodes in copy " << n_node 
                   << " not equal to the number in the original "
                   << n_node_orig << std::endl;

     throw OomphLibError(error_message.str(),
                         "Problem::copy()",
                         OOMPH_EXCEPTION_LOCATION);
    }
   
   //Loop over the nodes
   for(unsigned long i=0;i<n_node;i++)
    {     
     /// Try to cast to elastic node \todo there's got to be a better way
     /// but making Problem::mesh_pt() virtual doesn't do the right thing...
     SolidNode* el_node_pt=dynamic_cast<SolidNode*>(mesh_pt(m)->node_pt(i));
     if (el_node_pt!=0)
      {
       SolidNode* el_node_orig_pt=
        dynamic_cast<SolidNode*>(orig_problem_pt->mesh_pt(m)->node_pt(i));
       el_node_pt->copy(el_node_orig_pt);
      }
     else
      {
       mesh_pt(m)->node_pt(i)->copy(orig_problem_pt->mesh_pt(m)->node_pt(i));
      }
    }
  }


 // Copy global data:
 //------------------

 // Number of global data
 unsigned n_global=Global_data_pt.size();

 // Check # of nodes in orig problem
 unsigned long n_global_orig=orig_problem_pt->nglobal_data();
 if (n_global!=n_global_orig)
  {
   std::ostringstream error_message;
   error_message << "Number of global data in copy " << n_global
                 << " not equal to the number in the original "
                 << n_global_orig << std::endl;

   throw OomphLibError(error_message.str(),
                       "Problem::copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 for (unsigned iglobal=0;iglobal<n_global;iglobal++)
  {
   Global_data_pt[iglobal]->copy(orig_problem_pt->global_data_pt(iglobal));
  }


 // Copy internal data of elements:
 //--------------------------------

 // Loop over submeshes:
 for (unsigned m=0;m<nmesh;m++)
  {
   // Loop over elements and deal with internal data
   unsigned n_element=mesh_pt(m)->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     GeneralisedElement* el_pt=mesh_pt(m)->element_pt(e);
     unsigned n_internal=el_pt->ninternal_data();
     if (n_internal>0)
      {
       // Check # of internals :
       unsigned long n_internal_orig=orig_problem_pt->
        mesh_pt(m)->element_pt(e)->ninternal_data();
       if (n_internal!=n_internal_orig)
        {
         std::ostringstream error_message;
         error_message << "Number of internal data in copy " << n_internal
                       << " not equal to the number in the original "
                       << n_internal_orig << std::endl;
         
         throw OomphLibError(error_message.str(),
                             "Problem::copy()",
                             OOMPH_EXCEPTION_LOCATION);
        }
       for (unsigned i=0;i<n_internal;i++)
        {
         el_pt->internal_data_pt(i)->copy(orig_problem_pt->
                                          mesh_pt(m)->element_pt(e)->
                                          internal_data_pt(i));
        }
      }
    }
  }

}

//=========================================================================
/// Make and return a pointer to the copy of the problem. A virtual
/// function that must be filled in by the user is they wish to perform
/// adaptive refinement in bifurcation tracking or in multigrid problems.
/// ALH: WILL NOT BE NECESSARY IN BIFURCATION TRACKING IN LONG RUN...
//=========================================================================
Problem* Problem::make_copy()
{
 std::ostringstream error_stream;
 error_stream 
  << "This function must be overloaded in your specific problem, and must\n"
  << "create an exact copy of your problem. Usually this will be achieved\n"
  << "by a call to the constructor with exactly the same arguments as used\n";

 throw OomphLibError(error_stream.str(),
                     "Problem::make_copy()",
                     OOMPH_EXCEPTION_LOCATION);
}



//=========================================================================
/// Dump refinement pattern of all refineable meshes and all  generic
/// Problem data to file for restart. 
//=========================================================================
void Problem::dump(std::ofstream& dump_file)
{

 // Number of submeshes?
 unsigned n_mesh=nsub_mesh();
 
 // Single mesh:
 //------------
 if(n_mesh==0)
  {
   // Dump single mesh refinement pattern (if mesh is refineable)
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    { 
     mmesh_pt->dump_refinement(dump_file);
    }
  }
 
 //Multiple submeshes
 //------------------
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<n_mesh;imesh++)
    {
     // Dump single mesh refinement pattern (if mesh is refineable)
     if(RefineableMeshBase* mmesh_pt =
        dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
      {
       mmesh_pt->dump_refinement(dump_file);
      }
    } // End of loop over submeshes
  } 

 // Dump time
 // ---------

 // Flag to indicate unsteady run
 bool unsteady_flag=(time_pt()!=0);
 dump_file << unsteady_flag << " # bool flag for unsteady" << std::endl;

 // Current time and previous time increments for proper unsteady run
 if (unsteady_flag)
  {
   // Current time
   dump_file << time_pt()->time() << " # Time " << std::endl;
   // Timesteps
   unsigned n_dt=time_pt()->ndt();
   dump_file << n_dt << " # Number of timesteps " << std::endl;
   for (unsigned i=0;i<n_dt;i++)
    { 
     dump_file << time_pt()->dt(i) << " # dt " << std::endl;
    }
  }
 // Dummy time and previous time increments for steady run
 else
  {
   // Current time
   dump_file << "0.0 # Dummy time from steady run " << std::endl;
   // Timesteps
   dump_file << "0 # Dummy number of timesteps from steady run" << std::endl;
  }

 // Loop over submeshes and dump their data
 unsigned nmesh=nsub_mesh();
 if (nmesh==0) nmesh=1;
 for(unsigned m=0;m<nmesh;m++) {mesh_pt(m)->dump(dump_file);}

 // Dump global data

 // Loop over global data
 unsigned Nglobal=Global_data_pt.size();
 dump_file << Nglobal << " # number of global Data items " << std::endl;
 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   Global_data_pt[iglobal]->dump(dump_file);
   dump_file << std::endl;
  }

}

//=========================================================================
/// Read refinement pattern of all refineable meshes and refine them
/// accordingly, then read all Data and nodal position info from 
/// file for restart. Return flag to indicate if the restart was from 
/// steady or unsteady solution.
//=========================================================================
void Problem::read(std::ifstream& restart_file, bool& unsteady_restart) 
{

 //Call the actions before adaptation
 actions_before_adapt();
 
 // Number of submeshes?
 unsigned n_mesh=nsub_mesh();
 
 // Single mesh:
 //------------
 if(n_mesh==0)
  {
   // Refine single mesh (if it's refineable)
   if(RefineableMeshBase* mmesh_pt =
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    { 
     // When we get in here the problem has been constructed
     // by the constructor and the mesh is its original unrefined
     // form. 
     // RefineableMeshBase::refine(...) reads the refinement pattern from the
     // specified file and performs refinements until the mesh has
     // reached the same level of refinement as the mesh that existed
     // when the problem was dumped to disk.
     mmesh_pt->refine(restart_file);
    }
  }
 
 //Multiple submeshes
 //------------------
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<n_mesh;imesh++)
    {
     // Refine single mesh (if its refineable)
     if(RefineableMeshBase* mmesh_pt
        =dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
      {
       // When we get in here the problem has been constructed
       // by the constructor and the mesh is its original unrefined
       // form. 
       // RefineableMeshBase::refine(...) reads the refinement pattern from 
       // the specified file and performs refinements until the mesh has
       // reached the same level of refinement as the mesh that existed
       // when the problem was dumped to disk.
       mmesh_pt->refine(restart_file);
      }
    } // End of loop over submeshes

   // Rebuild the global mesh
   rebuild_global_mesh();
  } 

 //Any actions after adapt
 actions_after_adapt();
    
 // Setup equation numbering scheme
 oomph_info <<"\nNumber of equations: " << assign_eqn_numbers() 
            << std::endl<< std::endl; 

 std::string input_string;
 
 // Read time
 //----------

 // Read line up to termination sign
 getline(restart_file,input_string,'#');

 // Ignore rest of line
 restart_file.ignore(80,'\n');

 // Is the restart data from an unsteady run?
 unsteady_restart=atol(input_string.c_str());

 // Read line up to termination sign
 getline(restart_file,input_string,'#');

 // Ignore rest of line
 restart_file.ignore(80,'\n');

 // Read in initial time and set
 double time=atof(input_string.c_str());
 if (unsteady_restart) time_pt()->time()=time;
 

 // Read line up to termination sign
 getline(restart_file,input_string,'#');

 // Ignore rest of line
 restart_file.ignore(80,'\n');

 // Read & set number of timesteps
 unsigned n_dt=atoi(input_string.c_str());
 if (unsteady_restart) time_pt()->resize(n_dt);
 Vector<double> dt(n_dt);


 // Read in timesteps:
 for (unsigned i=0;i<n_dt;i++)
  {
   // Read line up to termination sign
   getline(restart_file,input_string,'#');
   
   // Ignore rest of line
   restart_file.ignore(80,'\n');
   
   // Read in initial time and set
   double prev_dt=atof(input_string.c_str());
   dt[i]=prev_dt;
  }
 
 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 if (unsteady_restart) initialise_dt(dt);

 // Loop over submeshes:
 unsigned nmesh=nsub_mesh();
 if (nmesh==0) nmesh=1;
 for (unsigned m=0;m<nmesh;m++)
  {
   mesh_pt(m)->read(restart_file);
  }


 // Read global data:
 //------------------

 // Number of global data
 unsigned Nglobal=Global_data_pt.size();

 // Read line up to termination sign
 getline(restart_file,input_string,'#');

 // Ignore rest of line
 restart_file.ignore(80,'\n');

 // Check # of nodes:
 unsigned long check_nglobal=atoi(input_string.c_str());

 if (check_nglobal!=Nglobal)
  {
   std::ostringstream error_message;
   error_message << "The number of global data " << Nglobal 
                 << " is not equal to that specified in the input file " 
                 <<   check_nglobal << std::endl;

   throw OomphLibError(error_message.str(),
                       "Problem::read()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   Global_data_pt[iglobal]->read(restart_file);
  }

}

//===================================================================
/// Set all timesteps to the same value, dt, and assign 
/// weights for all timesteppers in the problem.
//===================================================================
void Problem::initialise_dt(const double& dt)
{
 // Initialise the timesteps in the Problem's time object
 Time_pt->initialise_dt(dt);
 
 //Find out how many timesteppers there are
 unsigned n_time_steppers = ntime_stepper();
 
 //Loop over them all and set the weights
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   time_stepper_pt(i)->set_weights();
  }
}

//=========================================================================
/// Set the value of the timesteps to be equal to the values passed in 
/// a vector and assign weights for all timesteppers in the problem
//========================================================================
void Problem::initialise_dt(const Vector<double>& dt)
{
 // Initialise the timesteps in the Problem's time object
 Time_pt->initialise_dt(dt);
 
 //Find out how many timesteppers there are
 unsigned n_time_steppers = ntime_stepper();
 
 //Loop over them all and set the weights
 for(unsigned i=0;i<n_time_steppers;i++)
  {
   time_stepper_pt(i)->set_weights();
  }
}

//========================================================
/// Self-test: Check meshes and global data. Return 0 for OK
//========================================================
unsigned Problem::self_test()
{ 
 // Initialise
 bool passed=true;

 // Are there any submeshes?
 unsigned Nmesh=nsub_mesh();

 // Just one mesh: Check it
 if (Nmesh==0)
  {
   if (mesh_pt()->self_test()!=0)
    {
     passed=false;
     oomph_info 
      << "\n ERROR: Failed Mesh::self_test() for single mesh in problem" 
      << std::endl;
    }
  }
 // Loop over all submeshes and check them
 else
  {
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {
     if (mesh_pt(imesh)->self_test()!=0)
      {
       passed=false;
       oomph_info << "\n ERROR: Failed Mesh::self_test() for mesh imesh" 
                  << imesh  << std::endl;
      }
    }
  }
  

 // Check global data
 unsigned Nglobal=Global_data_pt.size();
 for (unsigned iglobal=0;iglobal<Nglobal;iglobal++)
  {
   if (Global_data_pt[iglobal]->self_test()!=0)
    {
     passed=false;
     oomph_info << "\n ERROR: Failed Data::self_test() for global data iglobal" 
                << iglobal << std::endl;
    }
  }


#ifdef OOMPH_HAS_MPI

 if (Problem_has_been_distributed)
  {
   // Note: This throws an error if it fails so no return is required.
   DocInfo tmp_doc_info;
   tmp_doc_info.doc_flag()=false;
   check_halo_schemes(tmp_doc_info);
  }

#endif

 // Return verdict
 if (passed) {return 0;}
 else {return 1;}

}

//====================================================================
/// A function that is used to adapt a bifurcation-tracking
/// problem, which requires separate interpolation of the 
/// associated eigenfunction. The error measure is chosen to be
/// a suitable combination of the errors in the base flow and the
/// eigenfunction. The bifurcation type is passed as an argument
//=====================================================================
void Problem::bifurcation_adapt_helper(
 unsigned &n_refined, unsigned &n_unrefined,
 const unsigned &bifurcation_type, const bool &actually_adapt)
{
 //Storage for eigenfunction from the problem
 Vector<DoubleVector> eigenfunction;
 //Get the eigenfunction from the problem
 this->get_bifurcation_eigenfunction(eigenfunction);
 
 //Get the bifurcation parameter 
 double *parameter_pt = this->bifurcation_parameter_pt();
 
 //Get the frequency parameter if tracking a Hopf bifurcation
 double omega = 0.0;
 //If we're tracking a Hopf then also get the frequency
 if(bifurcation_type==3)
  {
   omega = dynamic_cast<HopfHandler*>(assembly_handler_pt())->omega();
  }

 //If we're tracking a Pitchfork get the slack parameter (Hack)
 double sigma = 0.0;
 if(bifurcation_type==2)
  {
   sigma = this->dof(this->ndof()-1);
  }

 //We can now deactivate the bifurcation tracking in the problem
 //to restore the degrees of freedom to the unaugmented value
 this->deactivate_bifurcation_tracking();
 
 //Next, we create copies of the present problem
 //The number of copies depends on the number of eigenfunctions 
 //One copy for each eigenfunction
 const unsigned n_copies = eigenfunction.size();
 Copy_of_problem_pt.resize(n_copies);
 
 //Loop over the number of copies
 for(unsigned c=0;c<n_copies;c++)
  {
   //If we don't already have a copy
   if(Copy_of_problem_pt[c]==0)
    {
     //Create the copy
     Copy_of_problem_pt[c] = this->make_copy();

     //Refine the copy to the same level as the current problem
     
     //Find number of submeshes
     const unsigned N_mesh = Copy_of_problem_pt[c]->nsub_mesh();
     //If there is only one mesh
     if(N_mesh==0)
      {
       //Can we refine the mesh
       if(RefineableMeshBase* mmesh_pt =
          dynamic_cast<RefineableMeshBase*>(Copy_of_problem_pt[c]->mesh_pt(0)))
        {
         //Is the adapt flag set
         if(mmesh_pt->adapt_flag())
          {
           //Now get the original problem's mesh if it's refineable
           if(RefineableMeshBase* original_mesh_pt
              = dynamic_cast<RefineableMeshBase*>(this->mesh_pt(0)))
           {
            mmesh_pt->refine_base_mesh_as_in_reference_mesh(original_mesh_pt);
           }
           else
            {
             oomph_info 
              << 
              "Info/Warning: Mesh in orginal problem is not refineable."
              << std::endl;
            }
          }
         else
          {
           oomph_info << "Info/Warning: Mesh adaptation is disabled in copy." 
                      << std::endl;
          }
        }
       else
        {
         oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                    << std::endl;
        }
      } //End of single mesh case
     //Otherwise loop over the submeshes
     else
      {
       for(unsigned m=0;m<N_mesh;m++)
        {
         //Can we refine the submesh
         if(RefineableMeshBase* mmesh_pt =
            dynamic_cast<RefineableMeshBase*>(
             Copy_of_problem_pt[c]->mesh_pt(m)))
          {
           //Is the adapt flag set
           if(mmesh_pt->adapt_flag())
            {
             //Now get the original problem's mesh
             if(RefineableMeshBase* original_mesh_pt
                = dynamic_cast<RefineableMeshBase*>(this->mesh_pt(m)))
              {
               mmesh_pt->
                refine_base_mesh_as_in_reference_mesh(original_mesh_pt);
              }
             else
              {
               oomph_info 
                << 
                "Info/Warning: Mesh in orginal problem is not refineable."
                << std::endl;
              }
            }
           else
            {
             oomph_info << 
              "Info/Warning: Mesh adaptation is disabled in copy." 
                        << std::endl;
            }
          }
         else
           {
            oomph_info << "Info/Warning: Mesh cannot be adapted in copy."
                       << std::endl;
           }
        }
       //rebuild the global mesh in the copy
       Copy_of_problem_pt[c]->rebuild_global_mesh();

      } //End of multiple mesh case
     
     //Must call actions after adapt
     Copy_of_problem_pt[c]->actions_after_adapt();
     
     //Assign the equation numbers to the copy (quietly)
     (void)Copy_of_problem_pt[c]->assign_eqn_numbers();
    }
  } //End of creation of copies
 

 //Now check some numbers
 for(unsigned c=0;c<n_copies;c++)
  {
   //Check that the dofs match for each copy
#ifdef PARANOID
   //If the problems don't match then complain
   if(Copy_of_problem_pt[c]->ndof() != this->ndof())
    {
       std::ostringstream error_stream;
       error_stream 
        << 
        "Number of unknowns in the problem copy " << c << " "
        << "not equal to number in the original:\n"
        << this->ndof() << " (original) " << Copy_of_problem_pt[c]->ndof() 
        << " (copy)\n";
       
       throw OomphLibError(error_stream.str(),
                           "Problem::bifurcation_adapt_helper()",
                         OOMPH_EXCEPTION_LOCATION);
      }
#endif    
     
   //Assign the eigenfunction(s) to the copied problems 
   Copy_of_problem_pt[c]->assign_eigenvector_to_dofs(eigenfunction[c]);
   //Set all pinned values to zero
   Copy_of_problem_pt[c]->set_pinned_values_to_zero();
  }
 
 //Symmetrise the problem if we are solving a pitchfork
 if(bifurcation_type==2)
  {Copy_of_problem_pt[0]->
    symmetrise_eigenfunction_for_adaptive_pitchfork_tracking();}
 
 //Find error estimates based on current problem and eigenproblem
 //Now we need to get the error estimates for both problems.
 Vector<Vector<double> > base_error, eigenfunction_error;
 this->get_all_error_estimates(base_error);
 //Loop over the copies
 for(unsigned c=0;c<n_copies;c++)
  {
   //Get the error estimates for the copy
   Copy_of_problem_pt[c]->get_all_error_estimates(eigenfunction_error);

   //Find the number of meshes
   unsigned n_mesh = base_error.size();

#ifdef PARANOID   
   if(n_mesh != eigenfunction_error.size())
    {
     std::ostringstream error_stream;
     error_stream << 
      "Problems do not have the same number of meshes\n"
                  << "Base : " << n_mesh 
                  << " : Eigenproblem : " 
                  << eigenfunction_error.size() << "\n";
     throw OomphLibError(error_stream.str(),
                         "Problem::bifurcation_adapt_helper()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   for(unsigned m=0;m<n_mesh;m++)
    {
     //Check the number of elements is the same
     unsigned n_element = base_error[m].size();
#ifdef PARANOID
     if(n_element != eigenfunction_error[m].size())
      {
       std::ostringstream error_stream;
       error_stream << "Mesh " << m << 
        " does not have the same number of elements in the two problems:\n"
                    << "Base: " << n_element << " :  Eigenproblem: " 
                    << eigenfunction_error[m].size() << "\n";
       throw OomphLibError(error_stream.str(),
                           "Problem::bifurcation_adapt_helper()",
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif
     //Now add all the error esimates together
     for(unsigned e=0;e<n_element;e++)
      {
       //Add the error estimates (lazy)
       base_error[m][e] += eigenfunction_error[m][e];
      }
      }
    } //End of loop over copies
 
 //Then refine all problems based on the combined measure
 //if we are actually adapting (not just estimating the errors)
 if(actually_adapt)
  {
   this->adapt_based_on_error_estimates(n_refined,n_unrefined,base_error);
   for(unsigned c=0;c<n_copies;c++)
    {
     Copy_of_problem_pt[c]->
      adapt_based_on_error_estimates(n_refined,n_unrefined,base_error);
    }
   //Symmetrise the problem (again) if we are solving for a pitchfork
   if(bifurcation_type==2)
    {Copy_of_problem_pt[0]->
      symmetrise_eigenfunction_for_adaptive_pitchfork_tracking();}
  
   //Now get the refined guess for the eigenvector
   for(unsigned c=0;c<n_copies;c++)
    {
     Copy_of_problem_pt[c]->get_dofs(eigenfunction[c]);
    }
  }

 //Reactivate the tracking
 switch(bifurcation_type)
  {
   //Fold tracking
  case 1:
   this->activate_fold_tracking(parameter_pt);
   break;
   
   //Pitchfork
  case 2:
   this->activate_pitchfork_tracking(parameter_pt,eigenfunction[0]);
   //reset the slack parameter
   this->dof(this->ndof()-1) = sigma;
   break;

   //Hopf
  case 3:
   this->activate_hopf_tracking(parameter_pt,omega,
                                eigenfunction[0],eigenfunction[1]);
   break;
   
  default:
   std::ostringstream error_stream;
   error_stream << "Bifurcation type " << bifurcation_type << " not known\n"
                << "1: Fold, 2: Pitchfork, 3: Hopf\n";
   throw OomphLibError(error_stream.str(),
                       "Problem::bifurcation_adapt_helper()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}


//====================================================================
/// A function that is used to document the errors when
/// adapting a bifurcation-tracking
/// problem, which requires separate interpolation of the 
/// associated eigenfunction. The error measure is chosen to be
/// a suitable combination of the errors in the base flow and the
/// eigenfunction. The bifurcation type is passed as an argument
//=====================================================================
void Problem::bifurcation_adapt_doc_errors(const unsigned &bifurcation_type)
{
 //Dummy arguments
 unsigned n_refined, n_unrefined;
 //Just call the bifurcation helper without actually adapting
 bifurcation_adapt_helper(n_refined,n_unrefined,bifurcation_type,false);
}



//========================================================================
/// Adapt problem:
/// Perform mesh adaptation for (all) refineable (sub)mesh(es),
/// based on their own error estimates and the target errors specified
/// in the mesh(es). Following mesh adaptation,
/// update global mesh, and re-assign equation numbers. 
/// Return # of refined/unrefined elements. On return from this
/// function, Problem can immediately be solved again.
//======================================================================
void Problem::adapt(unsigned &n_refined, unsigned &n_unrefined)
{
 //Get the bifurcation type
 int bifurcation_type = this->Assembly_handler_pt->bifurcation_type();
 //If we are tracking a bifurcation then call the bifurcation adapt function
 if(bifurcation_type!=0)
  {
   this->bifurcation_adapt_helper(n_refined,n_unrefined,bifurcation_type);
   //Return immediately
   return;
  }

 oomph_info << std::endl << std::endl;
 oomph_info << "Adapting problem:" << std::endl;
 oomph_info << "=================" << std::endl;

 //Call the actions before adaptation
 actions_before_adapt();

 // Initialise counters
 n_refined=0;
 n_unrefined=0;
 
 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();
 
 // Single mesh:
 //------------
 if(Nmesh==0)
  {
   // Refine single mesh if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    { 
     if (mmesh_pt->adapt_flag())
      {
       // Get pointer to error estimator
       ErrorEstimator* error_estimator_pt=mmesh_pt->
        spatial_error_estimator_pt();
       
#ifdef PARANOID
       if (error_estimator_pt==0)
        {
         throw OomphLibError(
          "Error estimator hasn't been set yet",
          "Problem::adapt()",
          OOMPH_EXCEPTION_LOCATION);
        }
#endif

       // Get error for all elements
       Vector<double> elemental_error(mmesh_pt->nelement());
       
       if (mmesh_pt->doc_info_pt()==0)
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                mesh_pt(0),elemental_error);
        }
       else
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                mesh_pt(0),elemental_error,
                                                *mmesh_pt->doc_info_pt());
        }
       
       // Store max./min actual error
       mmesh_pt->max_error()=
        std::abs(*std::max_element(elemental_error.begin(),
                                   elemental_error.end(),AbsCmp<double>()));
       
       mmesh_pt->min_error()=
        std::abs(*std::min_element(elemental_error.begin(),
                                   elemental_error.end(),AbsCmp<double>()));

       oomph_info << "\n Max/min error: " 
                  << mmesh_pt->max_error() << " "
                  << mmesh_pt->min_error() << std::endl;
       
       // Adapt mesh
       mmesh_pt->adapt(this->communicator_pt(),elemental_error);
        
       // Add to counters
       n_refined+=mmesh_pt->nrefined();
       n_unrefined+=mmesh_pt->nunrefined();

      }
     else
      {
       oomph_info << "Info/Warning: Mesh adaptation is disabled." << std::endl;
      }
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
    }

  }

 //Multiple submeshes
 //------------------
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {
     // Refine single mesh uniformly if possible
     if(RefineableMeshBase* mmesh_pt =
        dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
      {
       // Get pointer to error estimator
       ErrorEstimator* error_estimator_pt=mmesh_pt->
        spatial_error_estimator_pt();
        
#ifdef PARANOID
       if (error_estimator_pt==0)
        {
         throw OomphLibError(
          "Error estimator hasn't been set yet",
          "Problem::adapt()",
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
        
       if (mmesh_pt->adapt_flag())
        {
         // Get error for all elements
         Vector<double> elemental_error(mmesh_pt->nelement());
         if (mmesh_pt->doc_info_pt()==0)
          {
           error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                  mesh_pt(imesh),
                                                  elemental_error);
          }
         else
          {
           error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                  mesh_pt(imesh),
                                                  elemental_error,
                                                  *mmesh_pt->doc_info_pt());
          }
        
         // Store max./min error if the mesh has any elements
         if (mesh_pt(imesh)->nelement()>0)
          {
           mmesh_pt->max_error()=
            std::abs(*std::max_element(elemental_error.begin(),
                                       elemental_error.end(),
                                       AbsCmp<double>()));
          
           mmesh_pt->min_error()=
            std::abs(*std::min_element(elemental_error.begin(),
                                       elemental_error.end(),
                                       AbsCmp<double>()));
          }

         oomph_info << "\n Max/min error: " 
                    << mmesh_pt->max_error() << " "
                    << mmesh_pt->min_error() << std::endl;

         // Adapt mesh
         mmesh_pt->adapt(this->communicator_pt(),elemental_error); 
  
         // Add to counters
         n_refined+=mmesh_pt->nrefined();
         n_unrefined+=mmesh_pt->nunrefined();

        }
       else
        {
         oomph_info << "Info/Warning: Mesh adaptation is disabled." 
                    << std::endl;
        }
      }
     else
      {
       oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
      }
      
    } // End of loop over submeshes

   // Rebuild the global mesh
   rebuild_global_mesh();

  } 

 //Any actions after adapt
 actions_after_adapt();

 //Attach the boundary conditions to the mesh
 oomph_info <<"\nNumber of equations: " << assign_eqn_numbers() 
            << std::endl<< std::endl; 

}

//========================================================================
/// Perform mesh adaptation for (all) refineable (sub)mesh(es),
/// based on the error estimates in elemental_error 
/// and the target errors specified
/// in the mesh(es). Following mesh adaptation,
/// update global mesh, and re-assign equation numbers. 
/// Return # of refined/unrefined elements. On return from this
/// function, Problem can immediately be solved again.
//========================================================================
void Problem::adapt_based_on_error_estimates(unsigned &n_refined, 
                                             unsigned &n_unrefined,
                                             Vector<Vector<double> > 
                                             &elemental_error)
{
 oomph_info << std::endl << std::endl;
 oomph_info << "Adapting problem:" << std::endl;
 oomph_info << "=================" << std::endl;
 
 //Call the actions before adaptation
 actions_before_adapt();
 
 //Initialise counters
 n_refined = 0;
 n_unrefined = 0;
 
 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();
 
 // Single mesh:
 //------------
 if(Nmesh==0)
  {
   // Refine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(0)))
    { 
     if (mmesh_pt->adapt_flag())
      {
       // Adapt mesh
       mmesh_pt->adapt(this->communicator_pt(),elemental_error[0]);
       
       // Add to counters
       n_refined += mmesh_pt->nrefined();
       n_unrefined += mmesh_pt->nunrefined();
       
      }
     else
      {
       oomph_info << "Info/Warning: Mesh adaptation is disabled." 
                  << std::endl;
      }
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
    }
   
  }
 
 //Multiple submeshes
 //------------------
 else
  {
   // Loop over submeshes
    for (unsigned imesh=0;imesh<Nmesh;imesh++)
     {
      // Refine single mesh uniformly if possible
      if(RefineableMeshBase* mmesh_pt =
         dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(imesh)))
       {
        if (mmesh_pt->adapt_flag())
         {
          // Adapt mesh
          mmesh_pt->adapt(this->communicator_pt(),elemental_error[imesh]); 
          
          // Add to counters
          n_refined += mmesh_pt->nrefined();
          n_unrefined += mmesh_pt->nunrefined();
         }
        else
         {
          oomph_info << "Info/Warning: Mesh adaptation is disabled." 
                     << std::endl;
         }
       }
      else
       {
        oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
       }
      
     } // End of loop over submeshes
    
    // Rebuild the global mesh
    rebuild_global_mesh();
    
  }

 //Any actions after adapt
 actions_after_adapt();
 
 //Attach the boundary conditions to the mesh
 oomph_info <<"\nNumber of equations: " << assign_eqn_numbers() 
            << std::endl<< std::endl; 
 
}


//========================================================================
/// Return the error estimates computed by (all) refineable 
/// (sub)mesh(es) in the elemental_error structure, which consists of
/// a vector of elemental errors for each (sub)mesh.
//========================================================================
void Problem::get_all_error_estimates(Vector<Vector<double> > &elemental_error)
{
 // Number of submeshes?
 const unsigned Nmesh=nsub_mesh();

 // Single mesh:
 //------------
 if(Nmesh==0)
  {
   //There is only one mesh
   elemental_error.resize(1);
   // Refine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(0)))
    { 
     //If we can adapt the mesh
     if(mmesh_pt->adapt_flag())
      {
       // Get pointer to error estimator
       ErrorEstimator* error_estimator_pt=mmesh_pt->
        spatial_error_estimator_pt();
       
#ifdef PARANOID
       if (error_estimator_pt==0)
          {
           throw OomphLibError(
            "Error estimator hasn't been set yet",
            "Problem::get_all_error_estimates()",
            OOMPH_EXCEPTION_LOCATION);
          }
#endif
       
       // Get error for all elements
       elemental_error[0].resize(mmesh_pt->nelement());
       //Are we documenting the errors or not
       if(mmesh_pt->doc_info_pt()==0)
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                Problem::mesh_pt(0),
                                                elemental_error[0]);
        }
       else
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                Problem::mesh_pt(0),
                                                elemental_error[0],
                                                *mmesh_pt->doc_info_pt());
        }
       
       // Store max./min actual error
       mmesh_pt->max_error()=
        std::abs(*std::max_element(elemental_error[0].begin(),
                                   elemental_error[0].end(),AbsCmp<double>()));
       
       mmesh_pt->min_error()=
        std::abs(*std::min_element(elemental_error[0].begin(),
                                   elemental_error[0].end(),AbsCmp<double>()));
       
       oomph_info << "\n Max/min error: " 
                  << mmesh_pt->max_error() << " "
                  << mmesh_pt->min_error() << std::endl;
      }
     else
      {
       oomph_info << "Info/Warning: Mesh adaptation is disabled." 
                  << std::endl;
      }
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be adapted" << std::endl;
    }
   
  }
 
 //Multiple submeshes
 //------------------
 else
  {
   //Resize to the number of submeshes
   elemental_error.resize(Nmesh);
   
   // Loop over submeshes
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {
     // Refine single mesh uniformly if possible
     if(RefineableMeshBase* mmesh_pt =
        dynamic_cast<RefineableMeshBase*>(Problem::mesh_pt(imesh)))
      {
       // Get pointer to error estimator
       ErrorEstimator* error_estimator_pt=mmesh_pt->
        spatial_error_estimator_pt();
       
#ifdef PARANOID
       if (error_estimator_pt==0)
        {
         throw OomphLibError(
          "Error estimator hasn't been set yet",
          "Problem::get_all_error_estimates()",
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
       //If we can adapt the mesh
       if (mmesh_pt->adapt_flag())
        {
         // Get error for all elements
         elemental_error[imesh].resize(mmesh_pt->nelement());
         if (mmesh_pt->doc_info_pt()==0)
          {
           error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                  Problem::mesh_pt(imesh),
                                                  elemental_error[imesh]);
          }
         else
          {
           error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                  Problem::mesh_pt(imesh),
                                                  elemental_error[imesh],
                                                  *mmesh_pt->doc_info_pt());
          }
         
         // Store max./min error
         mmesh_pt->max_error()=
          std::abs(*std::max_element(elemental_error[imesh].begin(),
                                     elemental_error[imesh].end(),
                                     AbsCmp<double>()));
         
         mmesh_pt->min_error()=
          std::abs(*std::min_element(elemental_error[imesh].begin(),
                                     elemental_error[imesh].end(),
                                     AbsCmp<double>()));
         
         oomph_info << "\n Max/min error: " 
                    << mmesh_pt->max_error() << " "
                    << mmesh_pt->min_error() << std::endl;
        }
       else
        {
         oomph_info << "Info/Warning: Mesh adaptation is disabled." 
                    << std::endl;
        }
      }
     else
      {
       oomph_info << "Info/Warning: Mesh cannot be adapted." << std::endl;
      }
     
    } // End of loop over submeshes

  }
}

//========================================================================
/// \short Get max and min error for all elements in submeshes
//========================================================================
void Problem::doc_errors(DocInfo& doc_info)
{
 //Get the bifurcation type
 int bifurcation_type = this->Assembly_handler_pt->bifurcation_type();
 //If we are tracking a bifurcation then call the bifurcation adapt function
 if(bifurcation_type!=0)
  {
   this->bifurcation_adapt_doc_errors(bifurcation_type);
   //Return immediately
   return;
  }

 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();
 
 // Single mesh:
 //------------
 if (Nmesh==0)
  {
   // Is the single mesh refineable?
   if (RefineableMeshBase* mmesh_pt =
       dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    { 
     
     // Get pointer to error estimator
     ErrorEstimator* error_estimator_pt=mmesh_pt->
      spatial_error_estimator_pt();

#ifdef PARANOID
     if (error_estimator_pt==0)
      {
       throw OomphLibError(
        "Error estimator hasn't been set yet",
        "Problem::doc_errors()",
        OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // Get error for all elements
     Vector<double> elemental_error(mmesh_pt->nelement());
     if (!doc_info.doc_flag())
      {
       error_estimator_pt->get_element_errors(this->communicator_pt(),
                                              mesh_pt(0),
                                              elemental_error);
      }
     else
      {
       error_estimator_pt->get_element_errors(this->communicator_pt(),
                                              mesh_pt(0),
                                              elemental_error,
                                              doc_info);
      }

     // Store max./min actual error
     mmesh_pt->max_error()=
      std::abs(*std::max_element(elemental_error.begin(),
                                 elemental_error.end(),AbsCmp<double>()));
      
     mmesh_pt->min_error()=
      std::abs(*std::min_element(elemental_error.begin(),
                                 elemental_error.end(),AbsCmp<double>()));
      
     oomph_info << "\n Max/min error: " 
                << mmesh_pt->max_error() << " "
                << mmesh_pt->min_error() << std::endl;

    }

  }
  
 //Multiple submeshes
 //------------------
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {

     // Is the single mesh refineable?
     if (RefineableMeshBase* mmesh_pt=
         dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
      { 

       // Get pointer to error estimator
       ErrorEstimator* error_estimator_pt=mmesh_pt->
        spatial_error_estimator_pt();

#ifdef PARANOID
       if (error_estimator_pt==0)
        {
         throw OomphLibError(
          "Error estimator hasn't been set yet",
          "Problem::doc_errors()",
          OOMPH_EXCEPTION_LOCATION);
        }
#endif

       // Get error for all elements
       Vector<double> elemental_error(mmesh_pt->nelement());
       if (mmesh_pt->doc_info_pt()==0)
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                mesh_pt(imesh),
                                                elemental_error);
        }
       else
        {
         error_estimator_pt->get_element_errors(this->communicator_pt(),
                                                mesh_pt(imesh),
                                                elemental_error,
                                                *mmesh_pt->doc_info_pt());
        }
        
       // Store max./min error if the mesh has any elements
       if (mesh_pt(imesh)->nelement()>0)
        {
         mmesh_pt->max_error()=
          std::abs(*std::max_element(elemental_error.begin(),
                                     elemental_error.end(),AbsCmp<double>()));
        
         mmesh_pt->min_error()=
          std::abs(*std::min_element(elemental_error.begin(),
                                     elemental_error.end(),AbsCmp<double>()));
        }

       oomph_info << "\n Max/min error: " 
                  << mmesh_pt->max_error() << " "
                  << mmesh_pt->min_error() << std::endl;
      }

    } // End of loop over submeshes

  }

}

//========================================================================
/// Refine (one and only!) mesh by splitting the elements identified
/// by their numbers relative to the problems' only mesh, then rebuild 
/// the problem. 
//========================================================================
void Problem::refine_selected_elements(const Vector<unsigned>& 
                                       elements_to_be_refined)
{
 actions_before_adapt();
 
 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();

 // Single mesh:
 if (Nmesh==0)
  {
   // Refine single mesh if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    {
     mmesh_pt->refine_selected_elements(elements_to_be_refined);
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be refined " 
                << std::endl;
    }
  }
 //Multiple submeshes
 else
  {
   std::ostringstream error_message;
   error_message << "Problem::refine_selected_elements(...) only works for\n"
                 << "multiple-mesh problems if you specify the mesh\n"
                 << "number in the function argument before the Vector,\n"
                 << "or a Vector of Vectors for each submesh.\n"
                 << std::endl;
   throw OomphLibError(error_message.str(),
                       "Problem::refine_selected_elements()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Any actions after the adapatation phase
 actions_after_adapt();

 //Attach the boundary conditions to the mesh
 oomph_info <<"Number of equations: " 
            << assign_eqn_numbers() << std::endl; 
}

//========================================================================
/// Refine (one and only!) mesh by splitting the elements identified
/// by their pointers, then rebuild the problem. 
//========================================================================
void Problem::refine_selected_elements(const Vector<RefineableElement*>& 
                                       elements_to_be_refined_pt)
{
 actions_before_adapt();
 
 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();

 // Single mesh:
 if (Nmesh==0)
  {
   // Refine single mesh if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    {
     mmesh_pt->refine_selected_elements(elements_to_be_refined_pt);
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be refined " 
                << std::endl;
    }
  }
 //Multiple submeshes
 else
  {
   std::ostringstream error_message;
   error_message << "Problem::refine_selected_elements(...) only works for\n"
                 << "multiple-mesh problems if you specify the mesh\n"
                 << "number in the function argument before the Vector,\n"
                 << "or a Vector of Vectors for each submesh.\n"
                 << std::endl;
   throw OomphLibError(error_message.str(),
                       "Problem::refine_selected_elements()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Any actions after the adapatation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " << assign_eqn_numbers()
            << std::endl; 
}

//========================================================================
/// Refine specified submesh by splitting the elements identified
/// by their numbers relative to the specified mesh, then rebuild the problem. 
//========================================================================
void Problem::refine_selected_elements(const unsigned& i_mesh,
                                       const Vector<unsigned>& 
                                       elements_to_be_refined)
 {
  actions_before_adapt();
 
  // Number of submeshes?
  unsigned n_mesh=nsub_mesh();

  if (i_mesh>=n_mesh)
   {
    std::ostringstream error_message;
    error_message <<
     "Problem only has " << n_mesh << " submeshes. Cannot refine submesh " 
                         << i_mesh << std::endl;
    throw OomphLibError(error_message.str(),
                        "Problem::refine_selected_elements()",
                        OOMPH_EXCEPTION_LOCATION);   
   }

  // Refine single mesh if possible
  if (mesh_pt(i_mesh)->nelement()>0)
   {
    if(RefineableMeshBase* mmesh_pt = 
       dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
     {
      mmesh_pt->refine_selected_elements(elements_to_be_refined);
     }
    else
     {
      oomph_info << "Info/Warning: Mesh cannot be refined " 
                 << std::endl;
     }
   }

  if (n_mesh>1)
   {
    //Rebuild the global mesh
    rebuild_global_mesh();
   }

  //Any actions after the adapatation phase
  actions_after_adapt();

  //Do equation numbering
  oomph_info <<"Number of equations: " << assign_eqn_numbers()
            << std::endl; 
 }


//========================================================================
/// Refine specified submesh by splitting the elements identified
/// by their pointers, then rebuild the problem. 
//========================================================================
void Problem::refine_selected_elements(const unsigned& i_mesh,
                                       const Vector<RefineableElement*>& 
                                       elements_to_be_refined_pt)
{
 actions_before_adapt();
 
 // Number of submeshes?
 unsigned n_mesh=nsub_mesh();

 if (i_mesh>=n_mesh)
  {
   std::ostringstream error_message;
   error_message <<
    "Problem only has " << n_mesh << " submeshes. Cannot refine submesh " 
                 << i_mesh << std::endl;
   throw OomphLibError(error_message.str(),
                       "Problem::refine_selected_elements()",
                       OOMPH_EXCEPTION_LOCATION);   
  }

 // Refine single mesh if possible
 if (mesh_pt(i_mesh)->nelement()>0)
  {
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
     mmesh_pt->refine_selected_elements(elements_to_be_refined_pt);
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be refined " 
                << std::endl;
    }
  }

 if (n_mesh>1)
  {
   //Rebuild the global mesh
   rebuild_global_mesh();
  }

 //Any actions after the adapatation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " << assign_eqn_numbers()
            << std::endl; 
}

//========================================================================
/// Refine all submeshes by splitting the elements identified by their
/// numbers relative to each submesh in a Vector of Vectors, then 
/// rebuild the problem. 
//========================================================================
void Problem::refine_selected_elements(const Vector<Vector<unsigned> >&
                                       elements_to_be_refined)
 {
  actions_before_adapt();
 
  // Number of submeshes?
  unsigned n_mesh=nsub_mesh();

  // Refine all submeshes if possible
  for (unsigned i_mesh=0; i_mesh<n_mesh; i_mesh++)
   {
    if (mesh_pt(i_mesh)->nelement()>0)
     {
      if(RefineableMeshBase* mmesh_pt = 
         dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
       {
        mmesh_pt->refine_selected_elements(elements_to_be_refined[i_mesh]);
       }
      else
       {
        oomph_info << "Info/Warning: Mesh cannot be refined " 
                   << std::endl;
       }
     }
   }

  //Rebuild the global mesh
  rebuild_global_mesh();

  //Any actions after the adapatation phase
  actions_after_adapt();

  //Do equation numbering
  oomph_info <<"Number of equations: " << assign_eqn_numbers()
            << std::endl; 
 }

//========================================================================
/// Refine all submeshes by splitting the elements identified by their
/// pointers within each submesh in a Vector of Vectors, then 
/// rebuild the problem. 
//========================================================================
void Problem::refine_selected_elements(const 
                                       Vector<Vector<RefineableElement*> >&
                                       elements_to_be_refined_pt)
 {
  actions_before_adapt();
 
  // Number of submeshes?
  unsigned n_mesh=nsub_mesh();

  // Refine all submeshes if possible
  for (unsigned i_mesh=0; i_mesh<n_mesh; i_mesh++)
   {
    if (mesh_pt(i_mesh)->nelement()>0)
     {
      if(RefineableMeshBase* mmesh_pt = 
         dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
       {
        mmesh_pt->refine_selected_elements(elements_to_be_refined_pt[i_mesh]);
       }
      else
       {
        oomph_info << "Info/Warning: Mesh cannot be refined " 
                   << std::endl;
       }
     }
   }

  //Rebuild the global mesh
  rebuild_global_mesh();

  //Any actions after the adapatation phase
  actions_after_adapt();

  //Do equation numbering
  oomph_info <<"Number of equations: " << assign_eqn_numbers()
            << std::endl; 
 }


//========================================================================
/// Refine (all) refineable (sub)mesh(es) uniformly and rebuild problem
/// and doc refinement process.
//========================================================================
void Problem::refine_uniformly(DocInfo& doc_info)
{
 actions_before_adapt();

 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();
  
 // Single mesh:
 if (Nmesh==0)
  {
   // Refine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    {
     mmesh_pt->refine_uniformly(doc_info);
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be refined uniformly " 
                << std::endl;
    }
  }
 //Multiple submeshes
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {
     // If the mesh has elements...
     if (mesh_pt(imesh)->nelement()>0)
      {
       // Refine i-th submesh uniformly if possible
       if (RefineableMeshBase* mmesh_pt =
           dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
         mmesh_pt->refine_uniformly(doc_info);
        }
       else
        {
         oomph_info << "Info/Warning: Cannot refine mesh " << imesh 
                    << std::endl;
        } 
      }
    }
   //Rebuild the global mesh
   rebuild_global_mesh();
  }

 //Any actions after the adaptation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " 
            << assign_eqn_numbers() << std::endl; 

}

//========================================================================
/// Refine submesh i_mesh uniformly and rebuild problem;
/// doc refinement process.
//========================================================================
void Problem::refine_uniformly(const unsigned& i_mesh, 
                               DocInfo& doc_info)
{
 actions_before_adapt();
 
#ifdef PARANOID
 // Number of submeshes?
 if (i_mesh>=nsub_mesh())
  {
   std::ostringstream error_message;
   error_message  << "imesh " << i_mesh 
                  << " is greater than the number of sub meshes " 
                  << nsub_mesh() << std::endl;
 
   throw OomphLibError(error_message.str(),
                       "Problem::refine_uniformly()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // If the mesh has elements...
 if (mesh_pt(i_mesh)->nelement()>0)
  {
   // Refine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
     mmesh_pt->refine_uniformly(doc_info);
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be refined uniformly " 
                << std::endl;
    }
  }

 //Rebuild the global mesh
 rebuild_global_mesh();

 //Any actions after the adaptation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " 
            << assign_eqn_numbers() << std::endl; 

}
 

//========================================================================
/// Unrefine (all) refineable (sub)mesh(es) uniformly and rebuild problem.
/// Return 0 for success,
/// 1 for failure (if unrefinement has reached the coarsest permitted
/// level)
//========================================================================
unsigned Problem::unrefine_uniformly()
{
 actions_before_adapt();

 // Has unrefinement been successful?
 unsigned success_flag=0;

 // Number of submeshes?
 unsigned Nmesh=nsub_mesh();

 // Single mesh:
 if (Nmesh==0)
  {
   // Unrefine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(0)))
    {
     success_flag+=mmesh_pt->unrefine_uniformly(this->communicator_pt());
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be unrefined uniformly " 
                << std::endl;
    }
  }
 //Multiple submeshes
 else
  {
   // Loop over submeshes
   for (unsigned imesh=0;imesh<Nmesh;imesh++)
    {
     // If the mesh has elements...
     if (mesh_pt(imesh)->nelement()>0)
      {
       // Unrefine i-th submesh uniformly if possible
       if (RefineableMeshBase* mmesh_pt=
           dynamic_cast<RefineableMeshBase*>(mesh_pt(imesh)))
        {
         success_flag+=mmesh_pt->unrefine_uniformly(this->communicator_pt());
        }
       else
        {
         oomph_info << "Info/Warning: Cannot unrefine mesh " << imesh 
                    << std::endl;
        } 
      }
    }
   //Rebuild the global mesh
   rebuild_global_mesh();
  }

 //Any actions after the adaptation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " 
            << assign_eqn_numbers() << std::endl; 

 // Judge success
 if (success_flag>0)
  {
   return 1;
  }
 else
  {
   return 0;
  }

}

//========================================================================
/// Unrefine submesh i_mesh uniformly and rebuild problem.
/// Return 0 for success,
/// 1 for failure (if unrefinement has reached the coarsest permitted
/// level)
//========================================================================
unsigned Problem::unrefine_uniformly(const unsigned& i_mesh)
{
 actions_before_adapt();

 // Has unrefinement been successful?
 unsigned success_flag=0;

#ifdef PARANOID
 // Number of submeshes?
 if (i_mesh>=nsub_mesh())
  {
   std::ostringstream error_message;
   error_message  << "imesh " << i_mesh 
                  << " is greater than the number of sub meshes " 
                  << nsub_mesh() << std::endl;
 
   throw OomphLibError(error_message.str(),
                       "Problem::unrefine_uniformly()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // If the mesh has elements...
 if (mesh_pt(i_mesh)->nelement()>0)
  {
   // Unrefine single mesh uniformly if possible
   if(RefineableMeshBase* mmesh_pt = 
      dynamic_cast<RefineableMeshBase*>(mesh_pt(i_mesh)))
    {
     success_flag+=mmesh_pt->unrefine_uniformly(this->communicator_pt());
    }
   else
    {
     oomph_info << "Info/Warning: Mesh cannot be unrefined uniformly " 
                << std::endl;
    }
  }

 //Rebuild the global mesh
 rebuild_global_mesh();

 //Any actions after the adaptation phase
 actions_after_adapt();

 //Do equation numbering
 oomph_info <<"Number of equations: " 
            << assign_eqn_numbers() << std::endl; 

 // Judge success
 if (success_flag>0)
  {
   return 1;
  }
 else
  {
   return 0;
  }

}


//========================================================================
/// Do one timestep, dt, forward  using Newton's method with specified 
/// tolerance and linear solver specified via member data.
/// Keep adapting on all meshes to criteria specified in
/// these meshes (up to max_adapt adaptations are performed).
/// If first_timestep==true, re-set initial conditions after mesh adaptation.
/// Shifting of time can be suppressed by overwriting the
/// default value of shift (true). [Shifting must be done
/// if first_timestep==true because we're constantly re-assigning
/// the initial conditions; if first_timestep==true and shift==false
/// shifting is performed anyway and a warning is issued.
//========================================================================
void Problem::unsteady_newton_solve(const double &dt, 
                                    const unsigned &max_adapt, 
                                    const bool &first_timestep,
                                    const bool& shift)
{
 

 // Do shifting or not?
 bool shift_it=shift;

 // Warning:
 if (first_timestep&&(!shift))
  {
   shift_it=true;
   oomph_info 
    << "\n\n===========================================================\n";
   oomph_info << "                  ********  WARNING *********** \n";
   oomph_info 
    << "===========================================================\n";
   oomph_info << "Problem::unsteady_newton_solve() called with " << std::endl;
   oomph_info << "first_timestep: " << first_timestep << std::endl;
   oomph_info << "shift: " << shift << std::endl;
   oomph_info << "This doesn't make sense (shifting does have to be done" 
              << std::endl;
   oomph_info 
    << "since we're constantly re-assigning the initial conditions"
    << std::endl;
   oomph_info 
    << "\n===========================================================\n\n";
  }


 //Find the initial time
 double initial_time = time_pt()->time();

 // Max number of solves
 unsigned max_solve=max_adapt+1;
 
 // Adaptation loop
 //----------------
 for (unsigned isolve=0;isolve<max_solve;isolve++)
  {
   // Only adapt after the first solve has been done!
   if (isolve>0)
    {
     unsigned n_refined;
     unsigned n_unrefined;
     
     // Adapt problem 
     adapt(n_refined,n_unrefined);

#ifdef OOMPH_HAS_MPI
     // Adaptation only converges if ALL the processes have no
     // refinement or unrefinement to perform
     unsigned total_refined=0;
     unsigned total_unrefined=0;
     if (Problem_has_been_distributed)
      {
       MPI_Allreduce(&n_refined,&total_refined,1,MPI_INT,MPI_SUM,
                     this->communicator_pt()->mpi_comm());
       n_refined=total_refined;
       MPI_Allreduce(&n_unrefined,&total_unrefined,1,MPI_INT,MPI_SUM,
                     this->communicator_pt()->mpi_comm());
       n_unrefined=total_unrefined;
      }
#endif

     oomph_info << "---> " << n_refined << " elements were refined, and " 
                << n_unrefined << " were unrefined, in total." << std::endl;
     
     // Check convergence of adaptation cycle
     if ((n_refined==0)&&(n_unrefined==0))
      {
       oomph_info << "\n \n Solution is fully converged in "
                  << "Problem::unsteady_newton_solver() \n \n ";
       break;
      }       

     //Reset the time
     time_pt()->time() = initial_time;

     // Reset the inital condition on refined meshes
     if (first_timestep) 
      {
       oomph_info << "Re-setting initial condition " << std::endl;
       set_initial_condition();
      }
    }

   //Now do the actual unsteady timestep
   //If it's the first time around the loop, or the first timestep 
   //shift the timevalues, otherwise don't
   // Note: we need to shift if it's the first timestep because
   // we're constantly re-assigning the initial condition above!
   if((isolve==0) || (first_timestep))
    {
     Problem::unsteady_newton_solve(dt,shift_it);
    }
   // Subsequent solve: Have shifted already -- don't do it again.
   else
    {
     shift_it=false;
     Problem::unsteady_newton_solve(dt,shift_it);
    }

   if (isolve==max_solve-1)
    {
     oomph_info << std::endl 
                << "----------------------------------------------------------"
                << std::endl
                << "Reached max. number of adaptations in \n"
                << "Problem::unsteady_newton_solver().\n"
                << "----------------------------------------------------------"
                << std::endl
                << std::endl;
    }
   
  } // End of adaptation loop


}



//========================================================================
/// \short Adaptive Newton solver. 
/// The linear solver takes a pointer to the problem (which defines
/// the Jacobian \b J and the residual Vector \b r) and returns 
/// the solution \b x of the system 
/// \f[ {\bf J} {\bf x} = - \bf{r} \f]. 
/// Performs at most max_adapt adaptations on all meshes. 
//========================================================================
void Problem::newton_solve(const unsigned &max_adapt)
{
 // Max number of solves
 unsigned max_solve=max_adapt+1;

 // Adaptation loop
 //----------------
 for (unsigned isolve=0;isolve<max_solve;isolve++)
  {

   // Only adapt after the first solve has been done!
   if (isolve>0)
    {

     unsigned n_refined;
     unsigned n_unrefined;
       
     // Adapt problem 
     adapt(n_refined,n_unrefined);
       
#ifdef OOMPH_HAS_MPI
     // Adaptation only converges if ALL the processes have no
     // refinement or unrefinement to perform
     unsigned total_refined=0;
     unsigned total_unrefined=0;
     if (Problem_has_been_distributed)
      {
       MPI_Allreduce(&n_refined,&total_refined,1,MPI_INT,MPI_SUM,
                     this->communicator_pt()->mpi_comm());
       n_refined=total_refined;
       MPI_Allreduce(&n_unrefined,&total_unrefined,1,MPI_INT,MPI_SUM,
                     this->communicator_pt()->mpi_comm());
       n_unrefined=total_unrefined;
      }
#endif

     oomph_info << "---> " << n_refined << " elements were refined, and " 
                << n_unrefined << " were unrefined, in total." << std::endl;

     // Check convergence of adaptation cycle
     if ((n_refined==0)&&(n_unrefined==0))
      {
       oomph_info << "\n \n Solution is fully converged in "
                  << "Problem::newton_solver(). \n \n ";
       break;
      }       
    }


   // Do actual solve
   //----------------
   {
      
    //Now update anything that needs updating
    // NOT NEEDED -- IS CALLED IN newton_solve BELOW! #
    // actions_before_newton_solve();
      
    try
     {
      //Solve the non-linear problem for this timestep with Newton's method
      Problem::newton_solve();
     }
    //Catch any exceptions thrown in the Newton solver
    catch(NewtonSolverError &error)
     {
      oomph_info << std::endl
                 << "USER-DEFINED ERROR IN NEWTON SOLVER " << std::endl;
      //Check to see whether we have reached Max_iterations
      if(error.iterations==Max_newton_iterations)
       {
        oomph_info << "MAXIMUM NUMBER OF ITERATIONS (" 
                   << error.iterations 
                   << ") REACHED WITHOUT CONVERGENCE " << std::endl;
       }
      //If not, it must be that we have exceeded the maximum residuals
      else
       {
        oomph_info << "MAXIMUM RESIDUALS: " << error.maxres
                   <<"EXCEEDS PREDEFINED MAXIMUM " 
                   << Max_residuals
                   << std::endl;
       }

      //Die horribly!!
      std::ostringstream error_stream;
      error_stream << "Error occured in adaptive Newton solver. "
                   << std::endl;
      throw OomphLibError(error_stream.str(),
                          "Problem::newton_solve()",
                          OOMPH_EXCEPTION_LOCATION);
     }

    //Now update anything that needs updating
    // NOT NEEDED -- WAS CALLED IN newton_solve ABOVE
    // !actions_after_newton_solve();
      
   } //End of solve block


   if (isolve==max_solve-1)
    {
     oomph_info 
      << std::endl 
      << "----------------------------------------------------------" 
      << std::endl
      << "Reached max. number of adaptations in \n"
      << "Problem::newton_solver().\n"
      << "----------------------------------------------------------" 
      << std::endl << std::endl;
    }

  } // End of adaptation loop


}

//========================================================================
/// Flush any external storage for any submeshes
/// NB this would ordinarily take place within the adaptation procedure
/// for each submesh (See RefineableMesh::adapt_mesh(...)), but there
/// are instances where the actions_before/after_adapt routines are used
/// and no adaptive routines are called in between (e.g. when doc-ing 
/// errors at the end of an adaptive newton solver)
//========================================================================
void Problem::flush_all_external_storage()
{
 //Number of submeshes
 unsigned n_mesh=nsub_mesh();

 //External storage will only exist if there is more than one (sub)mesh
 if (n_mesh>1)
  {
   for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
    {
     mesh_pt(i_mesh)->flush_all_external_storage();
    }
  }
}


#ifdef OOMPH_HAS_MPI

//========================================================================
/// Check the halo/haloed/shared node/element schemes.
//========================================================================
void Problem::check_halo_schemes(DocInfo& doc_info)
{
 // The bulk of the stuff that was in this routine is mesh-based, and 
 // should therefore drop into the Mesh base class.  All that needs to remain
 // here is a "wrapper" which calls the function dependent upon the number
 // of (sub)meshes that may have been distributed.

 unsigned n_mesh=nsub_mesh();

 if (n_mesh==0)
  {
   oomph_info << "Checking halo schemes on single mesh" << std::endl;
   mesh_pt()->check_halo_schemes(this->communicator_pt(),doc_info,
                                 Max_permitted_error_for_halo_check);
  }
 else // there are submeshes
  {
   for (unsigned i_mesh=0; i_mesh<n_mesh; i_mesh++)
    {
     oomph_info << "Checking halo schemes on submesh " << i_mesh << std::endl;
     doc_info.number()=i_mesh;
     mesh_pt(i_mesh)->check_halo_schemes(this->communicator_pt(),doc_info,
                                         Max_permitted_error_for_halo_check);
    }
  }

}



//========================================================================
///  Synchronise the degrees of freedom by overwriting
/// the haloed values with their non-halo counterparts held
/// on other processors. This works on the (sub)mesh argument
//========================================================================
void Problem::synchronise_dofs(Mesh* &mesh_pt)
{ 
 MPI_Status status;

 // Storage for number of processors and current processor
 int n_proc=this->communicator_pt()->nproc();
 int my_rank=this->communicator_pt()->my_rank();

 // Loop over all processors whose eqn numbers are to be updated
 for (int rank=0;rank<n_proc;rank++)
  {   
   // Prepare a vector of values
   Vector<double> values_on_other_proc;
   Vector<double> internal_values_on_other_proc;
   
   // If I'm not the processor whose halo values are updated,
   // some of my nodes may be haloed: Stick their
   // values into the vector
   if (rank!=my_rank)
    {
     // How many of my nodes are haloed by the processor whose values
     // are updated?
     unsigned nnod=mesh_pt->nhaloed_node(rank);
     unsigned count=0;
     for (unsigned j=0;j<nnod;j++)
      {
       // Generalised to variable number of values per node
       Node* haloed_nod_pt=mesh_pt->haloed_node_pt(rank,j);

       // Does the node have a timestepper?  Synchronise all history values!
       unsigned n_prev=1;
       if (haloed_nod_pt->time_stepper_pt()!=0)
        {
         n_prev+=haloed_nod_pt->time_stepper_pt()->nprev_values();
        }

       unsigned nval=haloed_nod_pt->nvalue();

       for (unsigned ival=0;ival<nval;ival++)
        {
         for (unsigned t=0;t<n_prev;t++)
          {
           values_on_other_proc.push_back(haloed_nod_pt->value(t,ival));
           count++;
          }
        }

       // Synchronise positions too!
       unsigned n_dim=haloed_nod_pt->ndim();
       for (unsigned i_dim=0;i_dim<n_dim;i_dim++)
        {
         for (unsigned t=0;t<n_prev;t++)
          {
           values_on_other_proc.push_back(haloed_nod_pt->x(t,i_dim));
           count++;
          }
        }

       // Is it a solid node?
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(haloed_nod_pt);

       if (solid_nod_pt!=0)
        {
         unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
         for (unsigned ival=0; ival<nval; ival++)
          {
           for (unsigned t=0;t<n_prev;t++)
            {
             values_on_other_proc.push_back
              (solid_nod_pt->variable_position_pt()->value(t,ival));
             count++;
            }
          }

        }

      }
          
     // Since nval may vary from node to node in the most general case,
     // the best way to send/receive here is to get the size of
     // the array before the array is sent, and send that first, so that
     // the receiver knows how much data to expect.  The order will be
     // preserved since the halo/haloed nodes are already ordered correctly
     MPI_Send(&count,1,MPI_INT,rank,0,this->communicator_pt()->mpi_comm());

     // Send it across
     if (count!=0)
      {
       MPI_Send(&values_on_other_proc[0],count,MPI_DOUBLE,rank,1,
                this->communicator_pt()->mpi_comm());
      }
     
     // Now loop over haloed elements and prepare to send internal data
     Vector<FiniteElement*> haloed_elem_pt=mesh_pt->haloed_element_pt(rank);
     unsigned nelem_haloed=haloed_elem_pt.size();
     unsigned count_intern=0;

     for (unsigned e=0; e<nelem_haloed; e++)
      {
       // how many internal data values for this element?
       unsigned nintern_data = haloed_elem_pt[e]->ninternal_data();
 
       for (unsigned iintern=0; iintern<nintern_data; iintern++)
        {
         // Cache internal_data local copy
         Data* int_data_pt=haloed_elem_pt[e]->internal_data_pt(iintern);

         unsigned n_prev=1;
         if (int_data_pt->time_stepper_pt()!=0)
          {
           n_prev+=int_data_pt->time_stepper_pt()->nprev_values();
          }

         unsigned nval=int_data_pt->nvalue();

         for (unsigned ival=0;ival<nval;ival++)
          {
           for (unsigned t=0;t<n_prev;t++)
            {
             internal_values_on_other_proc.push_back(int_data_pt
                                                     ->value(t,ival));
             count_intern++;
            }
          }
        }
      }
     // send the size of the vector of internal data values to the receiver
     MPI_Send(&count_intern,1,MPI_INT,rank,2,
              this->communicator_pt()->mpi_comm());

     // now send the vector itself
     if (count_intern!=0)
      {
       MPI_Send(&internal_values_on_other_proc[0],count_intern,MPI_DOUBLE,rank,
                3,this->communicator_pt()->mpi_comm());
      }
  
     // done
    }
   // Receive the vector of values
   else
    {
     // Loop over all other processors to receive their
     // values
     for (int send_rank=0;send_rank<n_proc;send_rank++)
      {
       
       // Don't talk to yourself
       if (send_rank!=my_rank)
        {
         // How many of my nodes are halos whose non-halo counter
         // parts live on processor send_rank?
         unsigned nnod=mesh_pt->nhalo_node(send_rank);
         
         // Receive size of vector of values
         unsigned count=0;
         MPI_Recv(&count,1,MPI_INT,send_rank,0,
                  this->communicator_pt()->mpi_comm(),&status);

         if (count!=0)
          {
           // Prepare vector for receipt of values
           values_on_other_proc.resize(count);
      
           // Receive 
           MPI_Recv(&values_on_other_proc[0],count,MPI_DOUBLE,send_rank,
                    1,this->communicator_pt()->mpi_comm(),&status);

           // Copy into the values of the halo nodes
           // on the present processors
           count=0; // reset array index counter to zero
           for (unsigned j=0;j<nnod;j++)
            {
             // Generalised to variable number of values per node
             Node* halo_nod_pt=mesh_pt->halo_node_pt(send_rank,j);

             // If the node has a timestepper, synchronise all history values
             unsigned n_prev=1;
             if (halo_nod_pt->time_stepper_pt()!=0)
              {
               n_prev+=halo_nod_pt->time_stepper_pt()->nprev_values();
              }

             unsigned nval=halo_nod_pt->nvalue();
        
             for (unsigned ival=0;ival<nval;ival++)
              {
               for (unsigned t=0;t<n_prev;t++)
                {
                 halo_nod_pt->set_value(t,ival,values_on_other_proc[count]);
                 count++; // increase array index
                }
              }

             // Synchronise positions too
             unsigned n_dim=halo_nod_pt->ndim();
             for (unsigned i_dim=0;i_dim<n_dim;i_dim++)
              {
               for (unsigned t=0;t<n_prev;t++)
                {
                 halo_nod_pt->x(t,i_dim)=values_on_other_proc[count];
                 count++;
                }
              }

             // Is this a solid node?
             SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(halo_nod_pt);

             if (solid_nod_pt!=0)
              {
               unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 for (unsigned t=0;t<n_prev;t++)
                  {
                   solid_nod_pt->variable_position_pt()->set_value
                    (t,ival,values_on_other_proc[count]);;
                   count++;
                  }
                }
              }

            }
          }

         // Get number of halo elements whose non-halo is on process send_rank
         Vector<FiniteElement*> halo_elem_pt=mesh_pt->
          halo_element_pt(send_rank);
         unsigned nelem_halo=halo_elem_pt.size();
         
         // Receive size of vector of internal data values
         unsigned count_intern=0;
         MPI_Recv(&count_intern,1,MPI_INT,send_rank,2,
                  this->communicator_pt()->mpi_comm(),&status);

         // Prepare and receive vector
         if (count_intern!=0)
          {
           internal_values_on_other_proc.resize(count_intern);
           MPI_Recv(&internal_values_on_other_proc[0],count_intern,
                    MPI_DOUBLE,send_rank,3,
                    this->communicator_pt()->mpi_comm(),&status);

           // reset array counter index to zero
           count_intern=0;
           for (unsigned e=0;e<nelem_halo;e++)
            {
             unsigned nintern_data=halo_elem_pt[e]->ninternal_data();
             for (unsigned iintern=0;iintern<nintern_data;iintern++)
              {
               // Cache internal_data local copy
               Data* int_data_pt=halo_elem_pt[e]->internal_data_pt(iintern);

               // Does the data have a timestepper?
               unsigned n_prev=1;
               if (int_data_pt->time_stepper_pt()!=0)
                {
                 n_prev+=int_data_pt->time_stepper_pt()->nprev_values();
                }

               unsigned nval=int_data_pt->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 for (unsigned t=0;t<n_prev;t++)
                  {
                   int_data_pt->set_value
                    (t,ival,internal_values_on_other_proc[count_intern]);
                   count_intern++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//========================================================================
/// Synchronise the external degrees of freedom by overwriting
/// the external halo values with their external haloed counterparts held
/// on other processors and the external mesh.
//========================================================================
void Problem::synchronise_external_dofs(Mesh* &mesh_pt)
{ 
 MPI_Status status;

 // Storage for number of processors and current processor
 int n_proc=this->communicator_pt()->nproc();
 int my_rank=this->communicator_pt()->my_rank();

 // Loop over all processors whose eqn numbers are to be updated
 for (int rank=0;rank<n_proc;rank++)
  {   
   // Prepare a vector of values
   Vector<double> values_on_other_proc;
   Vector<double> internal_values_on_other_proc;
   
   // If I'm not the processor whose external halo values are updated,
   // some of my nodes may be externally haloed: Stick their
   // values into the vector
   if (rank!=my_rank)
    {
     // How many of my nodes are externally haloed by the processor whose
     // values are updated?  NB these nodes are on the external mesh.
     unsigned next_nod=mesh_pt->nexternal_haloed_node(rank);

     unsigned count=0;
     for (unsigned j=0;j<next_nod;j++)
      {
       // Generalised to variable number of values per node
       Node* ext_haloed_nod_pt=mesh_pt->external_haloed_node_pt(rank,j);

       // Does the node have a timestepper?  Synchronise all history values!
       unsigned n_prev=1;
       if (ext_haloed_nod_pt->time_stepper_pt()!=0)
        {
         n_prev+=ext_haloed_nod_pt->time_stepper_pt()->nprev_values();
        }

       // Synchronise nodal values
       unsigned nval=ext_haloed_nod_pt->nvalue();
       for (unsigned ival=0;ival<nval;ival++)
        {
         for (unsigned t=0;t<n_prev;t++)
          {
           values_on_other_proc.push_back(ext_haloed_nod_pt->value(t,ival));
           count++;
          }
        }

       // Synchronise history values for positions
       unsigned n_dim=ext_haloed_nod_pt->ndim();
       for (unsigned i_dim=0;i_dim<n_dim;i_dim++)
        {
         for (unsigned t=0;t<n_prev;t++)
          {
           values_on_other_proc.push_back(ext_haloed_nod_pt->x(t,i_dim));
           count++;
          }
        }

       // Is it a solid node? Synchronise solid node values too
       SolidNode* ext_solid_nod_pt=
        dynamic_cast<SolidNode*>(ext_haloed_nod_pt);

       if (ext_solid_nod_pt!=0)
        {
         unsigned nval=ext_solid_nod_pt->variable_position_pt()->nvalue();
         for (unsigned ival=0; ival<nval; ival++)
          {
           for (unsigned t=0;t<n_prev;t++)
            {
             values_on_other_proc.push_back
              (ext_solid_nod_pt->variable_position_pt()->value(t,ival));
             count++;
            }
          }
        }
      }
          
     // Since nval may vary from node to node in the most general case,
     // the best way to send/receive here is to get the size of
     // the array before the array is sent, and send that first, so that
     // the receiver knows how much data to expect.  The order will be
     // preserved since the halo/haloed nodes are already ordered correctly
     MPI_Send(&count,1,MPI_INT,rank,0,this->communicator_pt()->mpi_comm());

     if (count!=0)
      {
       // Send it across
       MPI_Send(&values_on_other_proc[0],count,MPI_DOUBLE,rank,1,
                this->communicator_pt()->mpi_comm());
      }
     
     // Now loop over haloed elements and prepare to send internal data
     unsigned next_elem_haloed=mesh_pt->nexternal_haloed_element(rank);
     unsigned count_intern=0;

     for (unsigned e=0; e<next_elem_haloed; e++)
      {
       // How many internal data values for this element?
       unsigned nintern_data = mesh_pt->
        external_haloed_element_pt(rank,e)->ninternal_data();
       for (unsigned iintern=0; iintern<nintern_data; iintern++)
        {
         // Cache internal_data local copy
         Data* int_data_pt=mesh_pt->
          external_haloed_element_pt(rank,e)->internal_data_pt(iintern);

         // Does the data have a timestepper?
         unsigned n_prev=1;
         if (int_data_pt->time_stepper_pt()!=0)
          {
           n_prev+=int_data_pt->time_stepper_pt()->nprev_values();
          }

         unsigned nval=int_data_pt->nvalue();
         for (unsigned ival=0;ival<nval;ival++)
          {
           for (unsigned t=0;t<n_prev;t++)
            {
             internal_values_on_other_proc.push_back(int_data_pt
                                                     ->value(t,ival));
             count_intern++;
            }
          }
        }
      }
     // Send the size of the vector of internal data values to the receiver
     MPI_Send(&count_intern,1,MPI_INT,rank,2,
              this->communicator_pt()->mpi_comm());

     if (count_intern!=0)
      {
       // Now send the vector itself
       MPI_Send(&internal_values_on_other_proc[0],count_intern,MPI_DOUBLE,rank,
                3,this->communicator_pt()->mpi_comm());
      }
    }
   // Receive the vector of values
   else
    {
     // Loop over all other processors to receive their values
     for (int send_rank=0;send_rank<n_proc;send_rank++)
      {
       // Don't talk to yourself
       if (send_rank!=my_rank)
        {
         // How many of my nodes are external halos whose external non-halo
         // counterparts live on processor send_rank? 
         unsigned next_nod=mesh_pt->nexternal_halo_node(send_rank);

         // Receive size of vector of values
         unsigned count=0;
         MPI_Recv(&count,1,MPI_INT,send_rank,0,
                  this->communicator_pt()->mpi_comm(),&status);

         if (count!=0)
          {
           // Prepare vector for receipt of values
           values_on_other_proc.resize(count);
      
           // Receive 
           MPI_Recv(&values_on_other_proc[0],count,MPI_DOUBLE,send_rank,
                    1,this->communicator_pt()->mpi_comm(),&status);

           // Copy into the values of the external halo nodes
           // on the present processors
           count=0; // reset array index counter to zero
           for (unsigned j=0;j<next_nod;j++)
            {
             // Generalised to variable number of values per node
             Node* ext_halo_nod_pt=mesh_pt->external_halo_node_pt(send_rank,j);

             // Does the node have a timestepper?  Synchronise all values!
             unsigned n_prev=1;
             if (ext_halo_nod_pt->time_stepper_pt()!=0)
              {
               n_prev+=ext_halo_nod_pt->time_stepper_pt()->nprev_values();
              }

             // Synchronise values
             unsigned nval=ext_halo_nod_pt->nvalue();        
             for (unsigned ival=0;ival<nval;ival++)
              {
               for (unsigned t=0;t<n_prev;t++)
                {
                 ext_halo_nod_pt->set_value(t,ival,
                                            values_on_other_proc[count]);
                 count++;
                }
              }

             // Synchronise history values for positions
             unsigned n_dim=ext_halo_nod_pt->ndim();
             for (unsigned i_dim=0;i_dim<n_dim;i_dim++)
              {
               for (unsigned t=0;t<n_prev;t++)
                {
                 ext_halo_nod_pt->x(t,i_dim)=values_on_other_proc[count];
                 count++;
                }
              }

             // Is this a solid node?
             SolidNode* ext_solid_nod_pt=
              dynamic_cast<SolidNode*>(ext_halo_nod_pt);
             if (ext_solid_nod_pt!=0)
              {
               unsigned nval=ext_solid_nod_pt->
                variable_position_pt()->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 for (unsigned t=0;t<n_prev;t++)
                  {
                   ext_solid_nod_pt->variable_position_pt()->
                    set_value(t,ival,values_on_other_proc[count]);
                   count++;
                  }
                }
              }
            }
          }

         // Get number of halo elements whose non-halo is on process send_rank
         unsigned next_elem_halo=mesh_pt->nexternal_halo_element(send_rank);
         
         // Receive size of vector of internal data values
         unsigned count_intern=0;
         MPI_Recv(&count_intern,1,MPI_INT,send_rank,2,
                  this->communicator_pt()->mpi_comm(),&status);

         if (count_intern!=0)
          {
           // Prepare and receive vector
           internal_values_on_other_proc.resize(count_intern);
           MPI_Recv(&internal_values_on_other_proc[0],count_intern,
                    MPI_DOUBLE,send_rank,3,
                    this->communicator_pt()->mpi_comm(),&status);

           // reset array counter index to zero
           count_intern=0;
           for (unsigned e=0;e<next_elem_halo;e++)
            {
             unsigned nintern_data=mesh_pt->
              external_halo_element_pt(send_rank,e)->ninternal_data();
             for (unsigned iintern=0;iintern<nintern_data;iintern++)
              {
               // Cache internal_data local copy
               Data* int_data_pt=mesh_pt->
                external_halo_element_pt(send_rank,e)->
                internal_data_pt(iintern);

               // Does the data have a timestepper?
               unsigned n_prev=1;
               if (int_data_pt->time_stepper_pt()!=0)
                {
                 n_prev+=int_data_pt->time_stepper_pt()->nprev_values();
                }

               unsigned nval=int_data_pt->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 for (unsigned t=0;t<n_prev;t++)
                  {
                   int_data_pt->set_value
                    (t,ival,internal_values_on_other_proc[count_intern]);
                   count_intern++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//========================================================================
///  Synchronise equation numbers and return the total
/// number of degrees of freedom in the overall problem
//========================================================================
long Problem::synchronise_eqn_numbers(const bool& assign_local_eqn_numbers)
{ 
 // number of equations on this processor
 unsigned my_n_eqn = Dof_pt.size();

 // my rank
 unsigned my_rank = Communicator_pt->my_rank();
 
 // number of processors
 unsigned nproc = Communicator_pt->nproc();

 // send my_n_eqn to with rank greater than my_rank
 unsigned n_send = nproc-my_rank-1;
 Vector<MPI_Request> send_req(n_send);
 for (unsigned p = my_rank+1; p < nproc; p++)
  {
   MPI_Isend(&my_n_eqn,1,MPI_UNSIGNED,p,0,
             Communicator_pt->mpi_comm(),&send_req[p-my_rank-1]);
  }

 // recv n_eqn from processors with rank less than my_rank
 Vector<unsigned> n_eqn_on_proc(my_rank);
 for (unsigned p = 0; p < my_rank; p++)
  {
   MPI_Recv(&n_eqn_on_proc[p],1,MPI_UNSIGNED,p,0,
            Communicator_pt->mpi_comm(),MPI_STATUS_IGNORE);
  }

 // determine the number of equation on processors with rank
 // less than my_rank
 unsigned my_eqn_num_base = 0;
 for (unsigned p = 0; p < my_rank; p++)
  {
   my_eqn_num_base += n_eqn_on_proc[p];
  }

 // Loop over all internal data (on elements) and bump up their
 // equation numbers if they exist
 unsigned nelem=mesh_pt()->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   GeneralisedElement* el_pt=mesh_pt()->element_pt(e); 

   unsigned nintern_data=el_pt->ninternal_data();
   for (unsigned iintern=0;iintern<nintern_data;iintern++)
    {
     Data* int_data_pt=el_pt->internal_data_pt(iintern);
     unsigned nval=int_data_pt->nvalue();
     for (unsigned ival=0;ival<nval;ival++)
      {
       int old_eqn_number=int_data_pt->eqn_number(ival);
       if (old_eqn_number>=0) // i.e. it's being used
        {
         // Bump up eqn number
         int new_eqn_number=old_eqn_number+my_eqn_num_base;
         int_data_pt->eqn_number(ival)=new_eqn_number;
        }
      }
    }
  }

 // Loop over all nodes on current processor and bump up their 
 // equation numbers if they're not pinned!
 unsigned nnod=mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=mesh_pt()->node_pt(j);

   // loop over ALL eqn numbers - variable number of values
   unsigned nval=nod_pt->nvalue(); 
   
   for (unsigned ival=0;ival<nval;ival++)
    {
     int old_eqn_number=nod_pt->eqn_number(ival); 
     // Include all eqn numbers
     if (old_eqn_number>=0)
      {
       // Bump up eqn number
       int new_eqn_number=old_eqn_number+my_eqn_num_base;
       nod_pt->eqn_number(ival)=new_eqn_number;
      }
    }

   // Is this a solid node? If so, need to bump up its equation number(s)
   SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(nod_pt);

   if (solid_nod_pt!=0)
    {
     // Find equation numbers
     unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
     for (unsigned ival=0;ival<nval;ival++)
      {
       int old_eqn_number=solid_nod_pt->variable_position_pt()
        ->eqn_number(ival);
       // include all eqn numbers

       if (old_eqn_number>=0)
        {
         // Bump up eqn number
         int new_eqn_number=old_eqn_number+my_eqn_num_base;
         solid_nod_pt->variable_position_pt()->eqn_number(ival)=new_eqn_number;
        }
      }
    }
  }

 // The following copies the halo(ed) eqn numbers, and therefore
 // needs to be called at a submesh level as the halo(ed) structure
 // is only visible at the submesh level
 unsigned nmesh=nsub_mesh();
 // Now copy the haloed eqn numbers across
 // This has to include the internal data equation numbers as well
 // as the solid node equation numbers
 if (nmesh==0)
  {
   copy_haloed_eqn_numbers_helper(mesh_pt());
  }
 else // nmesh!=0
  {
   for (unsigned imesh=0; imesh<nmesh; imesh++)
    {
     // Do the haloed eqn numbers for this submesh
     copy_haloed_eqn_numbers_helper(mesh_pt(imesh));
    }
  }

 // Also copy external haloed equation numbers
 for (unsigned i=0;i<nmesh;i++)
  {
   copy_external_haloed_eqn_numbers_helper(mesh_pt(i));
  }

 // Now the global equation numbers have been updated.
 //---------------------------------------------------
 // Setup the local equation numbers again.
 //----------------------------------------

 if (assign_local_eqn_numbers)
  {
   // Loop over the submeshes: Note we need to call the submeshes' own
   // assign_*_eqn_number() otherwise we miss additional functionality
   // that is implemented (e.g.) in SolidMeshes!
   unsigned n_sub_mesh=nsub_mesh();
   if (n_sub_mesh==0)
    {
     mesh_pt()->assign_local_eqn_numbers();
    }
   else
    {
     for (unsigned i=0;i<n_sub_mesh;i++)
      {
       mesh_pt(i)->assign_local_eqn_numbers();
      }
    }
  }

 // wait for the sends to complete
 if (n_send > 0)
  {
   Vector<MPI_Status> send_status(n_send);
   MPI_Waitall(n_send,&send_req[0],&send_status[0]);
  }

 // build the Dof distribution pt
 Dof_distribution_pt->rebuild(Communicator_pt,my_eqn_num_base,
                              my_n_eqn);

 // resize the sparse assemble with arrays previous allocation
 Sparse_assemble_with_arrays_previous_allocation.resize(0);
 
 // and return the total number of equations in the problem
 return (long)Dof_distribution_pt->nrow();
}


//=======================================================================
/// A private helper function to
/// copy the haloed equation numbers into the halo equation numbers.
/// The argument is the mesh/submesh which will be worked on.
//===================================================================
void Problem::copy_haloed_eqn_numbers_helper(Mesh* &mesh_pt)
{
 MPI_Status status;

 // Storage for number of processors and current processor
 int n_proc=this->communicator_pt()->nproc();
 int my_rank=this->communicator_pt()->my_rank();

 // Loop over all processors whose eqn numbers are to be updated
 for (int rank=0;rank<n_proc;rank++)
  {
   // Prepare a vector of equation numbers 
   Vector<int> eqn_numbers_on_other_proc;
   Vector<int> internal_eqn_numbers_on_other_proc;

   // If I'm not the processor whose halo eqn numbers are updated,
   // some of my nodes may be haloed: Stick their
   // eqn numbers into the vector
   if (rank!=my_rank)
    {
     // How many of my nodes are haloed by the processor whose eqn
     // numbers are updated?
     unsigned nnod=mesh_pt->nhaloed_node(rank);
     unsigned count=0;     
     for (unsigned j=0;j<nnod;j++)
      {
       // Generalise to variable number of values per node
       Node* haloed_nod_pt=mesh_pt->haloed_node_pt(rank,j);
       unsigned nval=haloed_nod_pt->nvalue();

       for (unsigned ival=0;ival<nval;ival++)
        {
         eqn_numbers_on_other_proc.push_back(haloed_nod_pt->eqn_number(ival));
         count++;
        }

       // Is it a solid node?
       SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(haloed_nod_pt);

       if (solid_nod_pt!=0)
        {
         unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
         for (unsigned ival=0; ival<nval; ival++)
          {
           eqn_numbers_on_other_proc.push_back
            (solid_nod_pt->variable_position_pt()->eqn_number(ival));
           count++;
          }
        }
      }

     // The receiving process needs to know how many values it's getting
     // since it only descends into the loop over the nodes after
     // receiving the vector
     MPI_Send(&count,1,MPI_INT,rank,0,this->communicator_pt()->mpi_comm());

     if (count!=0)
      {
       // Send it across
       MPI_Send(&eqn_numbers_on_other_proc[0],count,MPI_INT,rank,1,
                this->communicator_pt()->mpi_comm());
      }

     // now loop over haloed elements and prepare to send 
     // equation numbers for internal data
     Vector<FiniteElement*> haloed_elem_pt=mesh_pt->haloed_element_pt(rank);
     unsigned nelem_haloed=haloed_elem_pt.size();
     unsigned count_intern=0;

     for (unsigned e=0; e<nelem_haloed; e++)
      {
       // how many internal data values for this element?
       unsigned nintern_data = haloed_elem_pt[e]->ninternal_data();
 
       for (unsigned iintern=0; iintern<nintern_data; iintern++)
        {
         // Cache local copy of data
         Data* int_data_pt=haloed_elem_pt[e]->internal_data_pt(iintern); 
         unsigned nval=int_data_pt->nvalue();

         for (unsigned ival=0;ival<nval;ival++)
          {
           internal_eqn_numbers_on_other_proc.push_back
            (int_data_pt->eqn_number(ival));
           count_intern++;
          }
        }
      }
     // This vector should only be sent if count_intern is non-zero
     
     // send the size of the vector of internal data values to the receiver
     MPI_Send(&count_intern,1,MPI_INT,rank,2,
              this->communicator_pt()->mpi_comm());

     // now send the vector itself
     if (count_intern!=0)
      {
       MPI_Send(&internal_eqn_numbers_on_other_proc[0],count_intern,
                MPI_INT,rank,3,this->communicator_pt()->mpi_comm());
      }
     // done

    }
   // Receive the vector of eqn numbers
   else
    {
     // Loop over all other processors to receive their
     // eqn numbers
     for (int send_rank=0;send_rank<n_proc;send_rank++)
      {
       // Don't talk to yourself
       if (send_rank!=my_rank)
        {
         // How many of my nodes are halos whose non-halo counter
         // parts live on processor send_rank?
         unsigned nnod=mesh_pt->nhalo_node(send_rank);
         
         // Receive the size of the vector of eqn numbers
         unsigned count=0;
         MPI_Recv(&count,1,MPI_INT,send_rank,0,
                  this->communicator_pt()->mpi_comm(),&status);
         
         if (count!=0)
          {
           // Prepare vector for receipt of eqn numbers
           eqn_numbers_on_other_proc.resize(count);
         
           // Receive it
           MPI_Recv(&eqn_numbers_on_other_proc[0],count,MPI_INT,send_rank,
                    1,this->communicator_pt()->mpi_comm(),&status);

           // Copy into the equation numbers of the halo nodes
           // on the present processors
           count=0; // reset count for array index
           for (unsigned j=0;j<nnod;j++)
            {
             // Generalise to variable number of values per node
             Node* halo_nod_pt=mesh_pt->halo_node_pt(send_rank,j);
             unsigned nval=halo_nod_pt->nvalue();

             for (unsigned ival=0;ival<nval;ival++)
              {
               halo_nod_pt->eqn_number(ival)=eqn_numbers_on_other_proc[count];
               count++;
              }

             // Is this a solid node?
             SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(halo_nod_pt);

             if (solid_nod_pt!=0)
              {
               unsigned nval=solid_nod_pt->variable_position_pt()->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 solid_nod_pt->variable_position_pt()->eqn_number(ival)=
                  eqn_numbers_on_other_proc[count];
                 count++;
                }
              }
            }
          }

         // Get number of halo elements whose non-halo is on process send_rank
         Vector<FiniteElement*> halo_elem_pt=mesh_pt->
          halo_element_pt(send_rank);
         unsigned nelem_halo=halo_elem_pt.size();

         // Receive size of vector of internal data values
         unsigned count_intern=0;
         MPI_Recv(&count_intern,1,MPI_INT,send_rank,2,
                  this->communicator_pt()->mpi_comm(),&status);

         // Prepare and receive vector
         if (count_intern!=0)
          {
           internal_eqn_numbers_on_other_proc.resize(count_intern);
           MPI_Recv(&internal_eqn_numbers_on_other_proc[0],
                    count_intern,MPI_INT,send_rank,3,
                    this->communicator_pt()->mpi_comm(),&status);

           // reset array counter index to zero
           count_intern=0;
           for (unsigned e=0;e<nelem_halo;e++)
            {
             unsigned nintern_data=halo_elem_pt[e]->ninternal_data();
             for (unsigned iintern=0;iintern<nintern_data;iintern++)
              {
               Data* int_data_pt=halo_elem_pt[e]->internal_data_pt(iintern); 
               // cache internal_data_pt copy
               unsigned nval=int_data_pt->nvalue();

               for (unsigned ival=0;ival<nval;ival++)
                {
                 int_data_pt->eqn_number(ival)=
                  internal_eqn_numbers_on_other_proc[count_intern];
                 count_intern++;
                }
              }
            }
          }
         
        }
      }
    }

  }


}

//=======================================================================
/// A helper function to copy the external haloed equation 
/// numbers into the external halo equation numbers.
/// The arguments are the submesh for the external halo nodes and the
/// (external) submesh for the external haloed nodes
//===================================================================
void Problem::copy_external_haloed_eqn_numbers_helper(Mesh* &mesh_pt)
{
 MPI_Status status;

 // Storage for number of processors and current processor
 int n_proc=this->communicator_pt()->nproc();
 int my_rank=this->communicator_pt()->my_rank();

 // Loop over all processors whose eqn numbers are to be updated
 for (int rank=0;rank<n_proc;rank++)
  {
   // Prepare a vector of equation numbers 
   Vector<int> eqn_numbers_on_other_proc;
   Vector<int> internal_eqn_numbers_on_other_proc;

   // If I'm not the processor whose external halo eqn numbers are updated,
   // some of my nodes may be externally haloed: Stick their
   // eqn numbers into the vector
   if (rank!=my_rank)
    {
     // How many of my nodes are externally haloed by the processor whose
     // eqn numbers are updated?
     unsigned next_nod=mesh_pt->nexternal_haloed_node(rank);

     unsigned count=0;
     for (unsigned j=0;j<next_nod;j++)
      {
       // Variable number of values per node
       Node* ext_haloed_nod_pt=mesh_pt->external_haloed_node_pt(rank,j);
       unsigned nval=ext_haloed_nod_pt->nvalue();
       for (unsigned ival=0;ival<nval;ival++)
        {
         eqn_numbers_on_other_proc.push_back
          (ext_haloed_nod_pt->eqn_number(ival));
         count++;
        }

       // Is it a solid node?
       SolidNode* ext_solid_nod_pt=
        dynamic_cast<SolidNode*>(ext_haloed_nod_pt);

       if (ext_solid_nod_pt!=0)
        {
         unsigned nval=ext_solid_nod_pt->variable_position_pt()->nvalue();
         for (unsigned ival=0; ival<nval; ival++)
          {
           eqn_numbers_on_other_proc.push_back
            (ext_solid_nod_pt->variable_position_pt()->eqn_number(ival));
           count++;
          }
        }

      }

     // The receiving process needs to know how many values it's getting
     // since it only descends into the loop over the nodes after
     // receiving the vector
     MPI_Send(&count,1,MPI_INT,rank,0,this->communicator_pt()->mpi_comm());

     if (count!=0)
      {
       // Send it across
       MPI_Send(&eqn_numbers_on_other_proc[0],count,MPI_INT,rank,1,
                this->communicator_pt()->mpi_comm());
      }

     // now loop over external haloed elements and prepare to send 
     // equation numbers for internal data
     unsigned next_elem_haloed=mesh_pt->nexternal_haloed_element(rank);
     unsigned count_intern=0;

     for (unsigned e=0; e<next_elem_haloed; e++)
      {
       // how many internal data values for this element?
       unsigned nintern_data = mesh_pt->
        external_haloed_element_pt(rank,e)->ninternal_data();
 
       for (unsigned iintern=0; iintern<nintern_data; iintern++)
        {
         // Cache local copy of data
         Data* int_data_pt=mesh_pt->
          external_haloed_element_pt(rank,e)->internal_data_pt(iintern);
         unsigned nval=int_data_pt->nvalue();

         for (unsigned ival=0;ival<nval;ival++)
          {
           internal_eqn_numbers_on_other_proc.push_back
            (int_data_pt->eqn_number(ival));
           count_intern++;
          }
        }
      }
     
     // send the size of the vector of internal data values to the receiver
     MPI_Send(&count_intern,1,MPI_INT,rank,2,
              this->communicator_pt()->mpi_comm());

     if (count_intern!=0)
      {
       // now send the vector itself
       MPI_Send(&internal_eqn_numbers_on_other_proc[0],count_intern,
                MPI_INT,rank,3,this->communicator_pt()->mpi_comm());
      }
     // done

    }
   // Receive the vector of eqn numbers
   else
    {
     // Loop over all other processors to receive their
     // eqn numbers
     for (int send_rank=0;send_rank<n_proc;send_rank++)
      {
       // Don't talk to yourself
       if (send_rank!=my_rank)
        {
         // How many of my nodes are external halos whose external non-halo
         // counterparts live on processor send_rank?
         unsigned next_nod=mesh_pt->nexternal_halo_node(send_rank);

         // Receive the size of the vector of eqn numbers
         unsigned count=0;
         MPI_Recv(&count,1,MPI_INT,send_rank,0,
                  this->communicator_pt()->mpi_comm(),&status);
         
         if (count!=0)
          {
           // Prepare vector for receipt of eqn numbers
           eqn_numbers_on_other_proc.resize(count);
         
           // Receive it
           MPI_Recv(&eqn_numbers_on_other_proc[0],count,MPI_INT,send_rank,
                    1,this->communicator_pt()->mpi_comm(),&status);

           // Copy into the equation numbers of the external halo nodes
           // on the present processors
           count=0; // reset count for array index
           for (unsigned j=0;j<next_nod;j++)
            {
             // Generalise to variable number of values per node
             Node* ext_halo_nod_pt=mesh_pt->
              external_halo_node_pt(send_rank,j);
             unsigned nval=ext_halo_nod_pt->nvalue();
             for (unsigned ival=0;ival<nval;ival++)
              {
               ext_halo_nod_pt->eqn_number(ival)=
                eqn_numbers_on_other_proc[count];
               count++;
              }

             // Is this a solid node?
             SolidNode* ext_solid_nod_pt=
              dynamic_cast<SolidNode*>(ext_halo_nod_pt);

             if (ext_solid_nod_pt!=0)
              {
               unsigned nval=ext_solid_nod_pt->
                variable_position_pt()->nvalue();
               for (unsigned ival=0;ival<nval;ival++)
                {
                 ext_solid_nod_pt->variable_position_pt()->
                  eqn_number(ival)=eqn_numbers_on_other_proc[count];
                 count++;
                }
              }
            }
          }

         // Get number of external halo elements whose external haloed
         // counterpart is on process send_rank
         unsigned next_elem_halo=mesh_pt->nexternal_halo_element(send_rank);

         // Receive size of vector of internal data values
         unsigned count_intern=0;
         MPI_Recv(&count_intern,1,MPI_INT,send_rank,2,
                  this->communicator_pt()->mpi_comm(),&status);

         if (count_intern!=0)
          {
           // Prepare and receive vector
           internal_eqn_numbers_on_other_proc.resize(count_intern);
           MPI_Recv(&internal_eqn_numbers_on_other_proc[0],
                    count_intern,MPI_INT,send_rank,3,
                    this->communicator_pt()->mpi_comm(),&status);

           // reset array counter index to zero
           count_intern=0;
           for (unsigned e=0;e<next_elem_halo;e++)
            {
             unsigned nintern_data=mesh_pt->
              external_halo_element_pt(send_rank,e)->ninternal_data();
             for (unsigned iintern=0;iintern<nintern_data;iintern++)
              {
               Data* int_data_pt=mesh_pt->
                external_halo_element_pt(send_rank,e)->
                internal_data_pt(iintern);
               // cache internal_data_pt copy
               unsigned nval=int_data_pt->nvalue();

               for (unsigned ival=0;ival<nval;ival++)
                {
                 int_data_pt->eqn_number(ival)=
                  internal_eqn_numbers_on_other_proc[count_intern];
                 count_intern++;
                }
              }
            }
          }

        }
      }
    }

  }


}

#endif


}
