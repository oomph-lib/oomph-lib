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
#include "partitioning.h"
#include "mesh.h"


namespace oomph
{

//====================================================================
/// Namespace for METIS graph partitioning routines
//====================================================================
namespace METIS
{



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
// delete [] adjcncy;
 delete [] part;
 delete [] edgecut;
 delete [] options;

}

}
