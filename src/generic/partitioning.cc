///LIC// ====================================================================
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
                     Vector<unsigned>& element_domain)
 {
  partition_mesh(problem_pt->mesh_pt(),
                 ndomain,
                 objective,
                 element_domain);
 }
 

 /// \short Use METIS to assign each element to a domain.
 /// On return, element_domain[ielem] contains the number
 /// of the domain [0,1,...,ndomain-1] to which 
 /// element ielem has been assigned.
 /// - objective=0: minimise edgecut.
 /// - objective=1: minimise total communications volume.
 /// .
 /// Partioning is based on nodal graph of mesh.
 void partition_mesh(Mesh* mesh_pt,
                     const unsigned& ndomain,
                     const unsigned& objective,
                     Vector<unsigned>& element_domain);
 
//  /// \short Use METIS to assign each element to a domain.
//  /// On return, element_domain[ielem] contains the number
//  /// of the domain [0,1,...,ndomain-1] to which 
//  /// element ielem has been assigned.
//  /// - objective=0: minimise edgecut.
//  /// - objective=1: minimise total communications volume.
//  /// .
//  /// Partioning is based on "Data" graph of mesh.
//  void partition_mesh_data(Problem* problem_pt,
//                           const unsigned& ndomain,
//                           const unsigned& objective,
//                           Vector<unsigned>& element_domain);

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


#ifdef USE_OLD_VERSION


//==================================================================
/// Use METIS to assign each element to a domain.
/// On return, element_domain[ielem] contains the number
/// of the domain [0,1,...,ndomain-1] to which 
/// element ielem has been assigned.
/// - objective=0: minimise edgecut.
/// - objective=1: minimise total communications volume.
/// .
/// Partioning is based on nodal graph of mesh.
//==================================================================
void METIS::partition_mesh(Mesh* mesh_pt, const unsigned& ndomain,
                           const unsigned& objective,
                           Vector<unsigned>& element_domain)
{

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

 // Setup nodal graph
 //------------------

 // Start timer
 clock_t cpu_start=clock();

 // Use map to associate node with number
 std::map<Node*,unsigned> node_number;

 unsigned nnode=mesh_pt->nnode();
 for (unsigned inode=0;inode<nnode;inode++)
  {
   Node* node_pt=mesh_pt->node_pt(inode);
   node_number[node_pt]=inode;
  }

 // Vector of sets which store the nodes that are connected with a given node
 Vector<std::set<unsigned> > connected_nodes(nnode);

 // Counter for total number of entries in connected_nodes structure
 unsigned count=0;
 
 // Loop over all elements
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   FiniteElement* el_pt=mesh_pt->finite_element_pt(ielem);
   unsigned nnode=el_pt->nnode();
   
   // Loop over all element-related nodes for caching
   std::map<Node*,unsigned> aux_node_number;
   
   for (unsigned inode=0;inode<nnode;inode++)
    {
     // Pointer to node
     Node* nod_pt=el_pt->node_pt(inode);
     aux_node_number[nod_pt]=node_number[nod_pt];
    }
   
   
   // Double loop over nodes to establish connectivity
   for (unsigned inode=0;inode<nnode;inode++)
    {
     // Get pointer to node
     Node* inode_pt=el_pt->node_pt(inode);
     
     // Get global node number:
     unsigned global_inode=aux_node_number[inode_pt];
     
     // Now loop over all other nodes in the element -- they are
     // connected
     for (unsigned jnode=inode;jnode<nnode;jnode++)
      {
       // Get pointer to node
       Node* jnode_pt=el_pt->node_pt(jnode);
       
       // Get global node number:
       unsigned global_jnode=aux_node_number[jnode_pt];
       
       // Node jnode is connected with node inode:
       connected_nodes[global_inode].insert(global_jnode);
       
       //...and vice versa (might not be but we need a symmetric
       // graph in METIS. Nobody (?) seems to know any better
       // ways for dealing with structurally unsymmetric matrices.
       connected_nodes[global_jnode].insert(global_inode);
       
       // We've added at most two new entries
       count+=2;
       
      }
    }
  }
 
 // Now convert into C-style packed array for interface with METIS
 int* xadj = new int[nnode+1];
 Vector<int> adjacency_vector;
 
 // Reserve (too much) space
 adjacency_vector.reserve(count);

 // Initialise counters
 unsigned ientry=0;

 // Loop over all nodes
 for (unsigned inode=0;inode<nnode;inode++)
  {
   // First entry for current node
   xadj[inode]=ientry;

   // Loop over nodes that are connected to current node
   typedef std::set<unsigned>::iterator IT;
   for (IT it=connected_nodes[inode].begin();
        it!=connected_nodes[inode].end();it++)
    {
     // Copy into adjacency array
     adjacency_vector.push_back(*it);

     // We've just made another entry
     ientry++;
    }

   // Entry after last entry for current node:
   xadj[inode+1]=ientry;
    
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
 int nvertex=nnode;

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
 int* edgecut = new int[nnode];

 // Array containing the partition information
 int* part = new int[nnode];


 // Can we get an error estimate?
 RefineableMeshBase* mmesh_pt=dynamic_cast<RefineableMeshBase*>(mesh_pt); 
 if (mmesh_pt!=0)
  {

   // Bias distribution?
   if (Error_to_weight_fct_pt!=&default_error_to_weight_fct)
    {
     oomph_info << "Biasing element distribution via spatial error estimate\n";

     // Adjust flag and provide storage for weights
     wgtflag=2;
     vwgt=new int[nnode];
     
     // Get error for all elements
     Vector<double> elemental_error(nelem);
     Vector<double> node_error(nnode,0.0);
     Vector<unsigned> node_count(nnode,0);
     
     mmesh_pt->spatial_error_estimator_pt()->
      get_element_errors(mesh_pt,elemental_error);
     
     // Now loop over all elements and add error to their nodes
     for (unsigned e=0;e<nelem;e++)
      {
       FiniteElement* el_pt=mesh_pt->finite_element_pt(e);
       unsigned nnod=el_pt->nnode();
       for (unsigned j=0;j<nnod;j++)
        {
         node_error[node_number[el_pt->node_pt(j)]]+=elemental_error[e];
         node_count[node_number[el_pt->node_pt(j)]]++;
        }
      }
     
     // Map to weight
     double max_error=-DBL_MAX;
     double min_error=DBL_MAX;
     for (unsigned j=0;j<nnode;j++)
      {
       double error=node_error[j]/double(node_count[j]);
       node_error[j]=error;
       if (error>max_error) max_error=error;
       if (error<min_error) min_error=error;
      }
     
     // Range of weights
     int weight=1;
     for (unsigned j=0;j<nnode;j++)
      {
       // Translate error into weight
       Error_to_weight_fct_pt(node_error[j],max_error,min_error,
                              weight);
       vwgt[j]=weight;
       
       // doc weight
       //mesh_pt->node_pt(j)->set_value(0,double(vwgt[j]));
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
   
 
 // End timer
 cpu_end=clock();

 // Doc
 double cpu1=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for METIS mesh partitioning                   [nelem="
  << nelem <<"]: " 
  << cpu1
  << " sec" << std::endl;

 // Now turn into mesh partitioning
 //--------------------------------

 // Start timer
 cpu_start=clock();

 // Vector holding the number of elements assigned to a given partition
 Vector<unsigned> nelements_in_partition(ndomain,0);
 
 // Boolean STL vector to check if element has been assigned to domain
 std::vector<bool> element_is_assigned(nelem,false);


 // Vector of maps holding the number of nodes in the elemnt
 //that are contained in a certain
 // partition
 Vector<std::map<unsigned,unsigned> > partition_count(nelem);
  

 // First loop over all elements to find out which partitions their
 // domains are part of
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {
   
   // Cache element pointer
   FiniteElement* el_pt=mesh_pt->finite_element_pt(ielem);

   // Loop over all element-related nodes
   unsigned nnode=el_pt->nnode();
   for (unsigned inode=0;inode<nnode;inode++)
    {
     // Get pointer to node
     Node* node_pt=el_pt->node_pt(inode);

     // This partition has just received another node
     partition_count[ielem][part[node_number[node_pt]]]++;
    }

   // Check if all nodes of the element are in the same partition
   if (partition_count[ielem].size()==1)
    {
     // Assign it
     element_domain[ielem]=partition_count[ielem].begin()->first;

     // Increment counter for elements in current partition
     nelements_in_partition[partition_count[ielem].begin()->first]++;

     // Label it as assigned
     element_is_assigned[ielem]=true;
    }
  }

 // Loop over all elements again to deal with the ones
 // that stradle two partitions
 for (unsigned ielem=0;ielem<nelem;ielem++)
  {

   if (!element_is_assigned[ielem])
    {
     // Find out which partition has most nodes in common with present element
     unsigned max_part=0;
     unsigned part_max=0;
     typedef std::map<unsigned,unsigned>::iterator IT;
     for (IT it=partition_count[ielem].begin();
          it!=partition_count[ielem].end();it++)
      {
       if (it->second>max_part)
        {
         part_max=it->first;
         max_part=it->second;
        }
      }

     // Check for dead heat:
     unsigned count_max=0;
     for (IT it=partition_count[ielem].begin();
          it!=partition_count[ielem].end();it++)
      {
       // How many partitions have the max. count?
       if (it->second==max_part)
        {
         count_max++;
        }
      }
      
     unsigned partition=0;
     if (count_max==1)
      {
       partition=part_max;
      }
     else
      {
       unsigned fewest_elements=1000000;
       for (IT it=partition_count[ielem].begin();
            it!=partition_count[ielem].end();it++)
        {
         // If this this partition is (one of) the most frequently
         // shared partitions, check if it's got the fewest elements
         // assigned
         if (it->second==max_part)
          {
           if (nelements_in_partition[it->first]<fewest_elements)
            {
             fewest_elements=nelements_in_partition[it->first];
             partition=it->first;
            }
          }
        }
      }

     // Assign to the partition that's shared by most nodes
     // and (if there are multiple ones) that has the fewest elements
     // assigned to it:
     element_domain[ielem]=partition; 

     // Increment counter for elements in current partition
     nelements_in_partition[partition]++;

    }
  }
 
 // End timer
 cpu_end=clock();

 // Doc
 double cpu2=double(cpu_end-cpu_start)/CLOCKS_PER_SEC;
 oomph_info  
  << "CPU time for conversion to oomph-lib mesh partitioning [nelem="
  << nelem <<"]: " 
  << cpu2
  << " sec" << std::endl;


  oomph_info  
  << "CPU time ratio (METIS/oomph-lib)                       [nelem="
  << nelem <<"]: " 
  << cpu1/(cpu0+cpu2)*100.0
  << "% " << std::endl;
 
 

 // Cleanup
 delete [] xadj;
 delete [] part;
 delete [] edgecut;
 delete [] options;

}


#else


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
void METIS::partition_mesh(Mesh* mesh_pt, const unsigned& ndomain,
                           const unsigned& objective,
                           Vector<unsigned>& element_domain)
{

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
 RefineableMeshBase* mmesh_pt=dynamic_cast<RefineableMeshBase*>(mesh_pt); 
 if (mmesh_pt!=0)
  {

   // Bias distribution?
   if (Error_to_weight_fct_pt!=&default_error_to_weight_fct)
    {
     oomph_info << "Biasing element distribution via spatial error estimate\n";

     // Adjust flag and provide storage for weights
     wgtflag=2;
     vwgt=new int[nelem];
     
     // Get error for all elements
     Vector<double> elemental_error(nelem);     
     mmesh_pt->spatial_error_estimator_pt()->
      get_element_errors(mesh_pt,elemental_error);
     
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
       
       // doc weight
       //mesh_pt->element_pt(j)->set_value(0,double(vwgt[j]));
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

#endif
}
