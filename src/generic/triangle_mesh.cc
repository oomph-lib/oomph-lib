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

#include <algorithm>
#include "map_matrix.h"
#include "triangle_mesh.h"


namespace oomph
{

 //=================================================
 ///  Helper namespace for BCInfo object used
 /// in the identification of boundary elements.
 //=================================================
namespace TriangleBoundaryHelper
{

 /// Structure for Boundary Informations
 struct BCInfo
 {

  /// Face ID
  unsigned Face_id;
 
  /// Boundary ID
  unsigned Boundary;
 
  /// Pointer to bulk finite element
  FiniteElement* FE_pt;
  
 };
 
//===================================================
/// Edge class
//===================================================
 class Edge
 {
  
 public:
  
  /// Constructor: Pass in the two vertex nodes
  Edge(Node* node1_pt, Node* node2_pt)
   {
    if (node1_pt==node2_pt)
     {
#ifdef PARANOID
      std::ostringstream error_stream;
      error_stream << "Edge cannot have two identical vertex nodes\n";
      throw OomphLibError(
       error_stream.str(),
       "TriangleMeshBase::setup_boundary_element_info()",
       OOMPH_EXCEPTION_LOCATION);
#endif
     }
    
    
    // Sort lexicographically
    if (node1_pt>node2_pt)
     {
      Node1_pt=node1_pt;
      Node2_pt=node2_pt;
     }
    else
     {
      Node1_pt=node2_pt;
      Node2_pt=node1_pt;
     }
    
   }
  
  
  /// Access to the first vertex node
  Node* node1_pt() const {return Node1_pt;}
  
  /// Access to the second vertex node
  Node* node2_pt() const {return Node2_pt;}
  
  /// Comparison operator
  bool operator==(const Edge& other) const
   {
    if ((Node1_pt==other.node1_pt())&&
        (Node2_pt==other.node2_pt()))
     
     {
       return true;
     }
    else
     {
      return false;
     }
   }
  
  
  
  /// Less-than operator
  bool operator<(const Edge& other) const
   {
    if (Node1_pt<other.node1_pt())
     {
      return true;
     }
    else if (Node1_pt==other.node1_pt())
     {
      if (Node2_pt<other.node2_pt())
       {
        return true;
       }
      else
       {
        return false;
       }
     }    
    else
     {
      return false;
     }
   }
  
 private:
  
  /// First vertex node
  Node* Node1_pt;
  
  /// Second vertex node
  Node* Node2_pt;
 };
 
}

//================================================================
/// Setup lookup schemes which establish which elements are located
/// next to which boundaries (Doc to outfile if it's open).
//================================================================
void TriangleMeshBase::setup_boundary_element_info(std::ostream &outfile)
{

 //Should we document the output here
 bool doc = false;

 if(outfile) doc = true;

 // Number of boundaries
 unsigned nbound=nboundary();
 
 // Wipe/allocate storage for arrays
 Boundary_element_pt.clear();
 Face_index_at_boundary.clear();
 Boundary_element_pt.resize(nbound);
 Face_index_at_boundary.resize(nbound);
 
 // Temporary vector of sets of pointers to elements on the boundaries: 
 Vector<std::set<FiniteElement*> > set_of_boundary_element_pt;
 set_of_boundary_element_pt.resize(nbound);
 
 // Matrix map for working out the fixed face for elements on boundary
 MapMatrixMixed<unsigned,FiniteElement*, int > 
  face_identifier;
 
 // Loop over elements
 //-------------------
 unsigned nel=nelement();

      
 // Get pointer to vector of boundaries that the
 // node lives on
 Vector<std::set<unsigned>*> boundaries_pt(3,0);


 // Data needed to deal with edges through the
 // interior of the domain
 std::map<TriangleBoundaryHelper::Edge,unsigned> edge_count;
 std::map<TriangleBoundaryHelper::Edge,TriangleBoundaryHelper::BCInfo> 
  edge_bcinfo;
 std::map<TriangleBoundaryHelper::Edge,TriangleBoundaryHelper::BCInfo>
  face_info;
 MapMatrixMixed<unsigned,FiniteElement*, int > face_count;
 Vector<unsigned> bonus(nbound);
     
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   FiniteElement* fe_pt=finite_element_pt(e);
   
   if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;
   
   // Only include 2D elements! Some meshes contain interface elements too.
   if (fe_pt->dim()==2)
    {
     // Loop over the element's nodes and find out which boundaries they're on
     // ----------------------------------------------------------------------

     //We need only loop over the corner nodes

     for(unsigned i=0;i<3;i++)
      {
       fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
      }

     //Find the common boundaries of each edge
     Vector<std::set<unsigned> > edge_boundary(3);


     //Edge 0 connects points 1 and 2
     //-----------------------------

     if(boundaries_pt[1] && boundaries_pt[2])
      {
       // Create the corresponding edge
       TriangleBoundaryHelper::Edge edge0(fe_pt->node_pt(1),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;      
       info.Face_id=0;
       info.FE_pt = fe_pt;
	
       std::set_intersection(boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[0],edge_boundary[0].begin()));
       std::set<unsigned>::iterator it0=edge_boundary[0].begin(); 

       // Edge does exist:
       if ( edge_boundary[0].size() > 0 )
        {
         info.Boundary=*it0;
         // How many times this edge has been visited
         edge_count[edge0]++;

         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge0,info)); 
        }
      }

     //Edge 1 connects points 0 and 2
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[2])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[1],edge_boundary[1].begin()));

       // Create the corresponding edge
       TriangleBoundaryHelper::Edge edge1(fe_pt->node_pt(0),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=1;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it1=edge_boundary[1].begin();

       // Edge does exist:
       if ( edge_boundary[1].size() > 0)
        {
         info.Boundary=*it1;

         // How many times this edge has been visited
         edge_count[edge1]++;  

         // Update edge_bcinfo              
         edge_bcinfo.insert(std::make_pair(edge1,info));
        }
      }
 
     //Edge 2 connects points 0 and 1
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[1])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[2],edge_boundary[2].begin()));

       // Create the corresponding edge
       TriangleBoundaryHelper::Edge edge2(fe_pt->node_pt(0),fe_pt->node_pt(1));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=2;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it2=edge_boundary[2].begin();

       // Edge does exist:
       if ( edge_boundary[2].size() > 0)
        {
         info.Boundary=*it2;

         // How many times this edge has been visited
         edge_count[edge2]++;  

         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge2,info));
        }
      }


#ifdef PARANOID
     
     // Check if edge is associated with multiple boundaries

     //We now know whether any edges lay on the boundaries
     for(unsigned i=0;i<3;i++)
      {
       //How many boundaries are there
       unsigned count = 0;

       //Loop over all the members of the set and add to the count
       //and set the boundary
       for(std::set<unsigned>::iterator it=edge_boundary[i].begin();
           it!=edge_boundary[i].end();++it)
        {
         ++count;
        }

       //If we're on more than one boundary, this is weird, so die
       if(count > 1)
        {
         std::ostringstream error_stream;
         error_stream << "Edge " << i << " is located on " << 
          count << " boundaries.\n";
         error_stream << "This is rather strange, so I'm going to die\n";
         throw OomphLibError(
          error_stream.str(),
          "TriangleMeshBase::setup_boundary_element_info()",
          OOMPH_EXCEPTION_LOCATION);
        }
      }
     
#endif

     // Now we set the pointers to the boundary sets to zero
     for(unsigned i=0;i<3;i++) {boundaries_pt[i] = 0;}
     
    }
  } //end of loop over all elements




 // Loop over all edges that are located on a boundary
 typedef std::map<TriangleBoundaryHelper::Edge,
  TriangleBoundaryHelper::BCInfo>::iterator ITE;
 for (ITE it=edge_bcinfo.begin();
      it!=edge_bcinfo.end();
      it++)
  {
   TriangleBoundaryHelper::Edge current_edge = it->first;
   unsigned  bound=it->second.Boundary; 
   
   //If the edge has been visited only once
   if ((edge_count[current_edge]==1) && (bound >= 0))
    {
     // Count the edges that are on the same element and on the same boundary
     face_count(static_cast<unsigned>(bound),it->second.FE_pt)=  
      face_count(static_cast<unsigned>(bound),it->second.FE_pt) + 1;
     
     //If such edges exist, let store the corresponding element
     if( face_count(bound,it->second.FE_pt) > 1)
      {
       // Update edge's infos
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=it->second.Face_id;
       info.FE_pt = it->second.FE_pt;
       info.Boundary=it->second.Boundary;
       
       // Add it to FIinfo, that stores infos of problematic elements
       face_info.insert(std::make_pair(current_edge,info));
       
       //How many edges on which boundary have to be added
       bonus[bound]++;
      }
     else
      {
       //Add element and face to the appropriate set
       set_of_boundary_element_pt[static_cast<unsigned>(bound)].insert(
        it->second.FE_pt);
       face_identifier(static_cast<unsigned>(bound),it->second.FE_pt) = 
        it->second.Face_id;
      }
     
    }
   
  }  //End of "adding-boundaries"-loop
 
 
  
 // Now copy everything across into permanent arrays
 //-------------------------------------------------

 // Loop over boundaries
 for (unsigned i=0;i<nbound;i++)
  {
   // Number of elements on this boundary that have to be added 
   // in addition to other elements
   unsigned bonus1=bonus[i];
   
   // Number of elements on this boundary (currently stored in a set)
   unsigned nel=set_of_boundary_element_pt[i].size() + bonus1;

   // Allocate storage for the coordinate identifiers
   Face_index_at_boundary[i].resize(nel);

   unsigned e_count=0;
   typedef std::set<FiniteElement*>::iterator IT;
   for (IT it=set_of_boundary_element_pt[i].begin();
        it!=set_of_boundary_element_pt[i].end();
        it++)
    {    
     // Recover pointer to element
     FiniteElement* fe_pt=*it;

     // Add to permanent storage
     Boundary_element_pt[i].push_back(fe_pt);

     Face_index_at_boundary[i][e_count] = face_identifier(i,fe_pt);

     // Increment counter
     e_count++;

    }
   // We add the elements that have two or more edges on this boundary
   for (ITE itt= face_info.begin(); itt!= face_info.end(); itt++)
    {
     if (itt->second.Boundary==i)
      {
       // Add to permanent storage
       Boundary_element_pt[i].push_back(itt->second.FE_pt);

       Face_index_at_boundary[i][e_count] = itt->second.Face_id;

       e_count++;
      }

    }

  } //End of loop over boundaries
 


 // Doc?
 //-----
 if (doc)
  {
   // Loop over boundaries
   for (unsigned i=0;i<nbound;i++)
    {
     unsigned nel=Boundary_element_pt[i].size();
     outfile << "Boundary: " << i
             << " is adjacent to " << nel
             << " elements" << std::endl;
     
     // Loop over elements on given boundary
     for (unsigned e=0;e<nel;e++)
      {
       FiniteElement* fe_pt=Boundary_element_pt[i][e];
       outfile << "Boundary element:" <<  fe_pt
               << " Face index of boundary is " 
               <<  Face_index_at_boundary[i][e] << std::endl;
      }
    }
  }
 

 // Lookup scheme has now been setup yet
 Lookup_for_elements_next_boundary_is_setup=true;

}

}
