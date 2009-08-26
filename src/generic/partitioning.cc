//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#include<float.h>

#include "partitioning.h"
#include "mesh.h"
#include "refineable_mesh.h"


//#define USE_OLD_VERSION


namespace oomph
{

//====================================================================
/// Namespace for METIS graph partitioning routines
//====================================================================
namespace METIS
{



 /// \short Default function that translates spatial
 /// error into weight for METIS partitioning (unit weight regardless
 /// of input).
 void default_error_to_weight_fct(const double& spatial_error, 
                                  const double& max_error, 
                                  const double& min_error,  
                                  int& weight)
 {
  weight=1;
 }

 /// \short Function pointer to to function that translates spatial
 /// error into weight for METIS partitioning.
 ErrorToWeightFctPt Error_to_weight_fct_pt=&default_error_to_weight_fct;
 

 /// \short Partition mesh uniformly by dividing elements
 /// equally over the partitions, in the order
 /// in which they are returned by problem.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 void uniform_partition_mesh(Problem* problem_pt,
                             const unsigned& ndomain,
                             Vector<unsigned>& element_domain);
 
 
 /// \short Use METIS to assign each element to a domain.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// Partioning is based on nodal graph of mesh.
 void partition_mesh(Problem* problem_pt,
                     const unsigned& ndomain,
                     const unsigned& objective,
                     Vector<unsigned>& element_domain);

#ifdef OOMPH_HAS_MPI

 /// \short Use METIS to assign each element in an already-distributed mesh
 /// to a domain. On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which element ielem has been
 /// assigned, on every processor involved in the Problem.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// The partioning is based on the nodal graph of the mesh by taking into
 /// account which global equation numbers are affected by each element and
 /// connecting elements which affect the same global equation number.
 void partition_distributed_mesh(Problem* problem_pt,
                                 const unsigned& objective,
                                 Vector<unsigned>& element_domain);

#endif

}




//==================================================================
/// Partition mesh uniformly by dividing elements
/// equally over the partitions, in the order
/// in which they are returned by problem.
/// On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which 
/// element ielem has been assigned.
//==================================================================
void METIS::uniform_partition_mesh(Problem* problem_pt,
                                   const unsigned& ndomain,
                                   Vector<unsigned>& element_domain)
{
 
 // Number of elements
 unsigned nelem=problem_pt->mesh_pt()->nelement();
 
#ifdef PARANOID
 if (nelem!=element_domain.size())
  {
   std::ostringstream error_stream;
   error_stream << "element_domain Vector has wrong length " 
                << nelem << " " << element_domain.size() << std::endl;

   throw OomphLibError(error_stream.str(),
                       "METIS::uniform_partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 
 
 // Uniform partitioning
 unsigned nel_per_domain=int(float(nelem)/float(ndomain));
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   unsigned idomain=unsigned(float(ielem)/float(nel_per_domain));
   element_domain[ielem]=idomain;
  }
 
}


//==================================================================
/// Use METIS to assign each element to a domain.
/// On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which 
/// element ielem has been assigned.
/// - objective=0: minimise edgecut.
/// - objective=1: minimise total communications volume.
/// .
/// Partioning is based on dual graph of mesh.
//==================================================================
void METIS::partition_mesh(Problem* problem_pt, const unsigned& ndomain,
                           const unsigned& objective,
                           Vector<unsigned>& element_domain)
{
 // Communicator
 OomphCommunicator* comm_pt=problem_pt->communicator_pt();

 // Global mesh
 Mesh* mesh_pt=problem_pt->mesh_pt();

 // Number of elements
 unsigned nelem=mesh_pt->nelement();

#ifdef PARANOID
 if (nelem!=element_domain.size())
  {
   std::ostringstream error_stream;
   error_stream << "element_domain Vector has wrong length " 
                << nelem << " " << element_domain.size() << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "METIS::partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 // Setup dual graph
 //------------------

 // Start timer
 clock_t cpu_start=clock();

 // Container to collect all elements associated with given global eqn number
 std::map<unsigned,std::set<unsigned> > elements_connected_with_global_eqn;
 
 // Container for all unique global eqn numbers
 std::set<unsigned> all_global_eqns;

 // Loop over all elements
 for (unsigned e=0;e<nelem;e++)
  {
   FiniteElement* el_pt=mesh_pt->finite_element_pt(e);

   // Add all global eqn numbers
   unsigned ndof=el_pt->ndof();
   for (unsigned j=0;j<ndof;j++)
    {
     // Get global eqn number
     unsigned eqn_number=el_pt->eqn_number(j);
     elements_connected_with_global_eqn[eqn_number].insert(e);
     all_global_eqns.insert(eqn_number);
    }
  }

 // Now reverse the lookup scheme to find out all elements
 // that are connected because they share the same global eqn
 Vector<std::set<unsigned> > connected_elements(nelem);
 
 // Counter for total number of entries in connected_elements structure
 unsigned count=0; 
   
 // Loop over all global eqns
 for (std::set<unsigned>::iterator it=all_global_eqns.begin();
      it!=all_global_eqns.end();it++)
  {
   // Get set of elements connected with this data item
   std::set<unsigned> elements=elements_connected_with_global_eqn[*it];
   
   // Double loop over connnected elements: Everybody's connected to
   // everybody
   for (std::set<unsigned>::iterator it1=elements.begin();
        it1!=elements.end();it1++) 
    {
     for (std::set<unsigned>::iterator it2=elements.begin();
          it2!=elements.end();it2++) 
      {
       if ((*it1)!=(*it2))
        {
         connected_elements[(*it1)].insert(*it2);
        }
      }
    }
  }
       
 
 // Now convert into C-style packed array for interface with METIS
 int* xadj = new int[nelem+1];
 Vector<int> adjacency_vector;
 
 // Reserve (too much) space
 adjacency_vector.reserve(count); 

 // Initialise counters
 unsigned ientry=0;

 // Loop over all elements
 for (unsigned e=0;e<nelem;e++)
  {
   // First entry for current element
   xadj[e]=ientry;

   // Loop over elements that are connected to current element
   typedef std::set<unsigned>::iterator IT;
   for (IT it=connected_elements[e].begin();
        it!=connected_elements[e].end();it++)
    {
     // Copy into adjacency array
     adjacency_vector.push_back(*it);
     
     // We've just made another entry
     ientry++;
    }
   
   // Entry after last entry for current element:
   xadj[e+1]=ientry;
   
  }
 
 // End timer
 clock_t cpu_end=clock();
 
 // Doc
 double cpu0=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for setup of METIS data structures            [nelem="
  << nelem <<"]: " 
  << cpu0
  << " sec" << std::endl;
 
 
 // Call METIS graph partitioner
 //-----------------------------
 
 // Start timer
 cpu_start=clock();
 
 // Number of vertices in graph
 int nvertex=nelem;

 // No vertex weights 
 int* vwgt=0;

 // No edge weights
 int* adjwgt=0;

 // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
 // Note that wgtflag==2 requires nodal weights to be stored in vwgt. 
 int wgtflag=0;

 // Use C-style numbering (first array entry is zero)
 int numflag=0; 

 // Number of desired partitions
 int nparts=ndomain;

 // Use default options
 int* options= new int[10];
 options[0]=0;

 // Number of cut edges in graph
 int* edgecut = new int[nelem];

 // Array containing the partition information
 int* part = new int[nelem];

 // Can we get an error estimate?

 unsigned n_mesh=problem_pt->nsub_mesh();

 if (n_mesh==0)
  {

   RefineableMeshBase* mmesh_pt=dynamic_cast<RefineableMeshBase*>(mesh_pt); 
   if (mmesh_pt!=0)
    {

     // Bias distribution?
     if (Error_to_weight_fct_pt!=&default_error_to_weight_fct)
      {
       oomph_info << 
        "Biasing element distribution via spatial error estimate\n";

       // Adjust flag and provide storage for weights
       wgtflag=2;
       vwgt=new int[nelem];
     
       // Get error for all elements
       Vector<double> elemental_error(nelem);     
       mmesh_pt->spatial_error_estimator_pt()->
        get_element_errors(comm_pt,mesh_pt,elemental_error);
     
       double max_error=*(std::max_element(elemental_error.begin(),
                                           elemental_error.end()));
       double min_error=*(std::min_element(elemental_error.begin(),
                                           elemental_error.end()));

       // Bias weights
       int weight=1;
       for (unsigned e=0;e<nelem;e++)
        {
         // Translate error into weight
         Error_to_weight_fct_pt(elemental_error[e],max_error,min_error,
                                weight);
         vwgt[e]=weight;
        }
      }
    }
  }
 else // There are submeshes
  {
   // Are any of the submeshes refineable?
   bool refineable_submesh_exists=false;
   // Vector to store "start and end point" for loops in submeshes
   Vector<unsigned> loop_helper(n_mesh+1);
   loop_helper[0]=0;

   // Loop over submeshes
   for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
    {
     // Store the end of the loop
     loop_helper[i_mesh+1]=problem_pt->mesh_pt(i_mesh)->nelement()+
      loop_helper[i_mesh];

     RefineableMeshBase* mmesh_pt=
      dynamic_cast<RefineableMeshBase*>(problem_pt->mesh_pt(i_mesh));
     if (mmesh_pt!=0)
      {
       refineable_submesh_exists=true;
      }
    }

   // If a refineable submesh exists
   if (refineable_submesh_exists)
    {
     // Bias distribution?
     if (Error_to_weight_fct_pt!=&default_error_to_weight_fct)
      {
       oomph_info << 
        "Biasing element distribution via spatial error estimate\n";

       // Adjust flag and provide storage for weights
       wgtflag=2;
       vwgt=new int[nelem];

       // Loop over submeshes
       for (unsigned i_mesh=0;i_mesh<n_mesh;i_mesh++)
        {
         RefineableMeshBase* mmesh_pt=
          dynamic_cast<RefineableMeshBase*>(problem_pt->mesh_pt(i_mesh));
         if (mmesh_pt!=0)
          {
           // Get error for all elements
           unsigned nsub_elem=loop_helper[i_mesh+1]-loop_helper[i_mesh];
           Vector<double> elemental_error(nsub_elem);
           mmesh_pt->spatial_error_estimator_pt()->
            get_element_errors(comm_pt,problem_pt->mesh_pt(i_mesh),
                               elemental_error);
     
           double max_error=*(std::max_element(elemental_error.begin(),
                                               elemental_error.end()));
           double min_error=*(std::min_element(elemental_error.begin(),
                                               elemental_error.end()));

           // Bias weights
           int weight=1;
           unsigned start=loop_helper[i_mesh];
           unsigned end=loop_helper[i_mesh+1];
           for (unsigned e=start;e<end;e++)
            {
             unsigned error_index=e-start;
             // Translate error into weight
             Error_to_weight_fct_pt(elemental_error[error_index],max_error,
                                    min_error,weight);
             vwgt[e]=weight;
            }
          }
         else // This mesh is not refineable
          {
           // There's no error estimator, so use the default weight
           int weight=1;
           unsigned start=loop_helper[i_mesh];
           unsigned end=loop_helper[i_mesh+1];
           for (unsigned e=start;e<end;e++)
            {
             vwgt[e]=weight;
            }
          }
        }
      }
    }
  }

 // Call partitioner
 if (objective==0)
  {
   // Partition with the objective of minimising the edge cut
   METIS_PartGraphKway(&nvertex, xadj, &adjacency_vector[0], vwgt, adjwgt, 
                       &wgtflag, &numflag, &nparts, options, edgecut, part);
  }
 else if (objective==1)
  {
   // Partition with the objective of minimising the total communication 
   // volume  
   METIS_PartGraphVKway(&nvertex, xadj, &adjacency_vector[0], vwgt, adjwgt, 
                        &wgtflag, &numflag, &nparts, options, edgecut, part);
  }
 else
  {
   std::ostringstream error_stream;
   error_stream 
    << "Wrong objective for METIS. objective = " << objective << std::endl;

   throw OomphLibError(error_stream.str(),"METIS::partition_mesh()",
                       OOMPH_EXCEPTION_LOCATION);
  }
   
 
 // Copy across
 for (unsigned e=0;e<nelem;e++)
  {
   element_domain[e]=part[e];
  }

 // End timer
 cpu_end=clock();

 // Doc
 double cpu1=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for METIS mesh partitioning                   [nelem="
  << nelem <<"]: " 
  << cpu1
  << " sec" << std::endl;
 

 // Cleanup
 delete [] xadj;
 delete [] part;
 delete [] edgecut;
 delete [] options;

}

#ifdef OOMPH_HAS_MPI

//==================================================================
/// \short Use METIS to assign each element in an already-distributed mesh
/// to a domain. On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which element ielem has been
/// assigned, on every processor involved in the Problem.
/// - objective=0: minimise edgecut.
/// - objective=1: minimise total communications volume.
/// .
/// The partioning is based on the nodal graph of the mesh by taking into 
/// account which global equation numbers are affected by each element and
/// connecting elements which affect the same global equation number
//==================================================================
void METIS::partition_distributed_mesh
(Problem* problem_pt,const unsigned& objective,
 Vector<unsigned>& element_domain_on_this_proc)
{
 // Start timer
 clock_t cpu_start=clock();

 // Communicator
 OomphCommunicator* comm_pt=problem_pt->communicator_pt();

 // Number of processors / domains
 unsigned n_proc=comm_pt->nproc();
 unsigned my_rank=comm_pt->my_rank();

 // Global mesh
 Mesh* mesh_pt=problem_pt->mesh_pt();

 // 1) reconstruct connectivity of whole mesh

 // Storage for global eqn numbers on current processor
 std::set<unsigned> global_eqns_on_this_proc;
 std::map<unsigned,std::set<unsigned> > 
  elements_connected_with_global_eqn_on_this_proc;

 std::map<unsigned,std::set<unsigned> > global_eqns_on_each_element;

 // Loop over non-halo elements on this processor
 unsigned n_elem=mesh_pt->nelement();
 unsigned count_non_halo=0;
 for (unsigned e=0; e<n_elem; e++)
  {
   FiniteElement* el_pt=mesh_pt->finite_element_pt(e);
   if (!el_pt->is_halo())
    {
     // Get global equation numbers
     unsigned n_dof=el_pt->ndof();
     for (unsigned i=0; i<n_dof; i++)
      {
       unsigned eqn_no=el_pt->eqn_number(i);
       elements_connected_with_global_eqn_on_this_proc[eqn_no].insert(e);
       global_eqns_on_this_proc.insert(eqn_no);
       global_eqns_on_each_element[e].insert(eqn_no);
      }
     count_non_halo++;
    }
  }

 // All that's required is for each element (on the "master" process) to know
 // which equation numbers it is "connected" to, and thus which other elements
 // it is "connected" to.  In order to achieve this we need to store the
 // equation numbers associated with each element on each processor locally
 // and then do a gather/reduce of this onto the "master" process

 Vector<int> n_element_on_each_proc(n_proc,0);
 Vector<int> start_index(n_proc,0);

 // Gather info to root
 unsigned root=0;
 MPI_Allgather(&count_non_halo, 1, MPI_INT,
               &n_element_on_each_proc[0], 1, MPI_INT,
               comm_pt->mpi_comm());

 unsigned total=0; 
 for (unsigned i_proc=0; i_proc<n_proc; i_proc++)
  {
   total+=n_element_on_each_proc[i_proc];
   if (i_proc!=0)
    {
     start_index[i_proc]=total-n_element_on_each_proc[i_proc];
    }
   else
    {
     start_index[0]=0;
    }
  }

 // For all non-halo elements set up a Vector detailing equation numbers
 // associated with each element in turn (flat-packed...)
 // and an index array to detail how many eqn numbers are associated with
 // each element in turn

 Vector<unsigned> eqn_numbers_with_elements_on_this_proc(0);

 Vector<unsigned> number_of_eqn_numbers(0);

 for (unsigned e=0;e<n_elem;e++)
  {
   FiniteElement* el_pt=mesh_pt->finite_element_pt(e);   
   if (!el_pt->is_halo())
    {
     // Loop over eqn numbers for this element and add to vector
     unsigned count=0;
     for (std::set<unsigned>::iterator it=
           global_eqns_on_each_element[e].begin();
          it!=global_eqns_on_each_element[e].end(); it++)
      {
       unsigned eqn_no=*it;
       eqn_numbers_with_elements_on_this_proc.push_back(eqn_no);
       count++;
      }
     number_of_eqn_numbers.push_back(count);
    }
  }

 // 
 Vector<int> n_eqns_on_each_proc(n_proc,0);
 Vector<int> start_eqns_index(n_proc,0);

 unsigned count_eqns_on_this_proc=
  eqn_numbers_with_elements_on_this_proc.size();

 // Gather info to root
 MPI_Allgather(&count_eqns_on_this_proc, 1, MPI_INT,
               &n_eqns_on_each_proc[0], 1, MPI_INT,
               comm_pt->mpi_comm());

 unsigned total_eqns=0; 
 for (unsigned i_proc=0; i_proc<n_proc; i_proc++)
  {
   total_eqns+=n_eqns_on_each_proc[i_proc];
   if (i_proc!=0)
    {
     start_eqns_index[i_proc]=total_eqns-n_eqns_on_each_proc[i_proc];
    }
   else
    {
     start_eqns_index[0]=0;
    }
  }


 // Now communicate all to the master process

 Vector<unsigned> total_number_of_eqn_numbers(total);


 MPI_Gatherv(&number_of_eqn_numbers[0], count_non_halo, MPI_INT,
             &total_number_of_eqn_numbers[0], &n_element_on_each_proc[0],
             &start_index[0], MPI_INT, root, comm_pt->mpi_comm());

 Vector<unsigned> eqn_numbers_with_elements(total_eqns);

 MPI_Gatherv(&eqn_numbers_with_elements_on_this_proc[0], 
             count_eqns_on_this_proc, MPI_INT,
             &eqn_numbers_with_elements[0], &n_eqns_on_each_proc[0], 
             &start_eqns_index[0], MPI_INT, root, comm_pt->mpi_comm());

 // Doc
 clock_t cpu_end=clock();

 double cpu0=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for global setup of METIS data structures [nelem="
  << total <<"]: " 
  << cpu0
  << " sec" << std::endl;

 // Now use METIS to determine "partition" for non-uniformly refined mesh

 // Loop over ALL elements on master...

 Vector<unsigned> element_domain(total,0);

 if (my_rank==root)
  {
   // Start timer
   clock_t cpu_start=clock();

   std::set<unsigned> all_global_eqns_root;

   std::map<unsigned,std::set<unsigned> > 
    elements_connected_with_global_eqn_on_root;

   unsigned count_all=0;
   for (unsigned e=0; e<total; e++)
    {
     unsigned n_eqn_no=total_number_of_eqn_numbers[e];
     for (unsigned n=0; n<n_eqn_no; n++)
      {
       unsigned eqn_no=eqn_numbers_with_elements[count_all];
       count_all++;
       elements_connected_with_global_eqn_on_root[eqn_no].insert(e);
       all_global_eqns_root.insert(eqn_no);
      }
    }

   // End timer
   clock_t cpu_equations=clock();
 
   // Doc
   double cpu0a=double(cpu_equations-cpu_start)/CLOCKS_PER_SEC;
   oomph_info  
    << "CPU time for setup of equations connected with elements [nelem="
    << total <<"]: " 
    << cpu0a
    << " sec" << std::endl;
 
   // Number of domains
   unsigned ndomain=n_proc;

   // and then keep going as in partitioning.cc I guess!

   // Now reverse the lookup scheme to find out all elements
   // that are connected because they share the same global eqn
   Vector<std::set<unsigned> > connected_elements(total);
 
   // Counter for total number of entries in connected_elements structure
   unsigned count=0; 
   
   // Loop over all global eqns
   for (std::set<unsigned>::iterator it=all_global_eqns_root.begin();
        it!=all_global_eqns_root.end();it++)
    {
     // Get set of elements connected with this data item
     std::set<unsigned> elements=
      elements_connected_with_global_eqn_on_root[*it];
   
     // Double loop over connnected elements: Everybody's connected to
     // everybody
     for (std::set<unsigned>::iterator it1=elements.begin();
          it1!=elements.end();it1++) 
      {
       for (std::set<unsigned>::iterator it2=elements.begin();
            it2!=elements.end();it2++) 
        {
         if ((*it1)!=(*it2))
          {
           connected_elements[(*it1)].insert(*it2);
          }
        }
      }
    }
       
   // End timer
   clock_t cpu_connected=clock();
 
   // Doc
   double cpu0b=double(cpu_connected-cpu_equations)/CLOCKS_PER_SEC;
   oomph_info  
    << "CPU time for setup of connected elements (load balance) [nelem="
    << total <<"]: " 
    << cpu0b
    << " sec" << std::endl;
  
   // Now convert into C-style packed array for interface with METIS
   int* xadj = new int[total+1];
   Vector<int> adjacency_vector;
 
   // Reserve (too much) space
   adjacency_vector.reserve(count); 

   // Initialise counters
   unsigned ientry=0;

   // Loop over all elements
   for (unsigned e=0;e<total;e++)
    {
     // First entry for current element
     xadj[e]=ientry;

     // Loop over elements that are connected to current element
     typedef std::set<unsigned>::iterator IT;
     for (IT it=connected_elements[e].begin();
          it!=connected_elements[e].end();it++)
      {
       // Copy into adjacency array
       adjacency_vector.push_back(*it);
     
       // We've just made another entry
       ientry++;
      }
   
     // Entry after last entry for current element:
     xadj[e+1]=ientry;
   
    }
 
   // End timer
   cpu_end=clock();
 
   // Doc
   double cpu0=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   oomph_info  
    << "CPU time for setup of METIS data structures (load balance) [nelem="
    << total <<"]: " 
    << cpu0
    << " sec" << std::endl;
 
 
   // Call METIS graph partitioner
   //-----------------------------
 
   // Start timer
   cpu_start=clock();
 
   // Number of vertices in graph
   int nvertex=total;

   // No vertex weights 
   int* vwgt=0;

   // No edge weights
   int* adjwgt=0;

   // Flag indicating that graph isn't weighted: 0; vertex weights only: 2
   // Note that wgtflag==2 requires nodal weights to be stored in vwgt. 
   int wgtflag=0;

   // Use C-style numbering (first array entry is zero)
   int numflag=0; 

   // Number of desired partitions
   int nparts=ndomain;

   // Use default options
   int* options= new int[10];
   options[0]=0;

   // Number of cut edges in graph
   int* edgecut = new int[total];

   // Array containing the partition information
   int* part = new int[total];

   if (objective==0)
    {
     // Partition with the objective of minimising the edge cut
     METIS_PartGraphKway(&nvertex, xadj, &adjacency_vector[0], vwgt, adjwgt, 
                         &wgtflag, &numflag, &nparts, options, edgecut, part);
    }
   else if (objective==1)
    {
     // Partition with the objective of minimising the total communication 
     // volume  
     METIS_PartGraphVKway(&nvertex, xadj, &adjacency_vector[0], vwgt, adjwgt, 
                          &wgtflag, &numflag, &nparts, options, edgecut, part);
    }

   // Copy across
   for (unsigned e=0;e<total;e++)
    {
     element_domain[e]=part[e];
    }

   // Doc
   double cpu1=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
   oomph_info  
    << "CPU time for METIS mesh partitioning [nelem="
    << total <<"]: " 
    << cpu1
    << " sec" << std::endl;

   // Cleanup
   delete [] xadj;
   delete [] part;
   delete [] edgecut;
   delete [] options;

  }

 // Communicate element domain information to all processes
 cpu_start=clock();

 element_domain_on_this_proc.resize(count_non_halo);

 MPI_Scatterv(&element_domain[0],&n_element_on_each_proc[0],&start_index[0],
              MPI_INT,&element_domain_on_this_proc[0],count_non_halo,MPI_INT,
              root,comm_pt->mpi_comm());

 // End timer
 cpu_end=clock();

 // Doc
 double cpu2=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for communication of partition to all processors [nelem="
  << total <<"]: " 
  << cpu2
  << " sec" << std::endl;
 


}

#endif

}
