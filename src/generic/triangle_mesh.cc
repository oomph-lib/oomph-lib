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
     Vector<std::set<unsigned> > edge(3);
     
     //Edge 0 connnects points 1 and 2
     if(boundaries_pt[1] && boundaries_pt[2])
      {
       std::set_intersection(boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge[0],edge[0].begin()));
      }

     //Edge 1 connects points 0 and 2
     if(boundaries_pt[0] && boundaries_pt[2])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge[1],edge[1].begin()));
      }
     
     //Edge 2 connects points 0 and 1
     if(boundaries_pt[0] && boundaries_pt[1])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge[2],edge[2].begin()));
      }

     //We now know whether any edges lay on the boundaries
     for(unsigned i=0;i<3;i++)
      {
       //How many boundaries are there
       unsigned count = 0;

       //The number of the boundary
       int boundary=-1;

       //Loop over all the members of the set and add to the count
       //and set the boundary
       for(std::set<unsigned>::iterator it=edge[i].begin();
           it!=edge[i].end();++it)
        {
         ++count;
         boundary = *it;
        }


       //If we're on more than one boundary, this is weird, so die
       if(count > 1)
        {
         std::ostringstream error_stream;
         error_stream << "Edge " << i << " is on " << 
          count << " boundaries.\n";
         error_stream << "This is rather strange, so I'm going to die\n";
         throw OomphLibError(
          error_stream.str(),
          "TriangleMeshBase::setup_boundary_element_info()",
          OOMPH_EXCEPTION_LOCATION);
        }

       //If we have a boundary then add this to the appropriate set
       if(boundary >= 0)
        {
         set_of_boundary_element_pt[static_cast<unsigned>(boundary)].
          insert(fe_pt);
         //Also set the fixed face
         face_identifier(static_cast<unsigned>(boundary),fe_pt) = i;
        }
      }

     //Now we set the pointers to the boundary sets to zero
     for(unsigned i=0;i<3;i++) {boundaries_pt[i] = 0;}

    }
  }
 
 // Now copy everything across into permanent arrays
 //-------------------------------------------------

 // Loop over boundaries
 //---------------------
 for (unsigned i=0;i<nbound;i++)
  {
   // Number of elements on this boundary (currently stored in a set)
   unsigned nel=set_of_boundary_element_pt[i].size();
    
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
  }
 


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
