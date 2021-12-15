// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#include "map_matrix.h"
#include "brick_mesh.h"

namespace oomph
{
  //====================================================================
  /// Helper namespace for generation of brick from tet mesh
  //====================================================================
  namespace BrickFromTetMeshHelper
  {
    /// Tolerance for mismatch during setup of boundary coordinates
    double Face_position_tolerance = 1.0e-12;

  } // namespace BrickFromTetMeshHelper


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  void BrickMeshBase::setup_boundary_element_info(std::ostream& outfile)
  {
    // Martina's fixed and commented version of the code.

    // Define variable doc and set the initial value to False.
    bool doc = false;

    // If file declared in outfile exists,set doc=true to enable writing to
    // stream
    if (outfile) doc = true;

    // Number of boundaries. Gives the value assigned by the function
    // nboundary()
    unsigned nbound = nboundary();

    if (doc)
    {
      outfile << "The number of boundaries is " << nbound << "\n";
    }

    // Wipe/allocate storage for arrays
    // Command clear will reomve all elements of vectors
    // Command resize will adapt pointer to correct boundary size. Internal data
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);


    // Temporary vector of vectors of pointers to elements on the boundaries:
    // vector_of_boundary_element_pt[i] is the vector of pointers to all
    // elements that have nodes on boundary i.
    // This is not a set to ensure UNIQUE ordering on every processor
    // There are a number of elements with nodes on each boundary.
    // For each boundary we store the appropriate elements in a vector.
    // There is therefore a vector for each boundary, which we store in
    // another vector, i.e. a matrix (vector of vectors).
    // Then vector_of_boundary_element_pt[0] is a vector of pointers to elements
    // with nodes on boundary 0 and vector_of_boundary_element_pt[0][0] is the
    // first element with nodes on boundary 0.
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed local coord for elements on boundary
    // Calling pointer to FiniteElement and pointer to  Vector<int> defining
    // the matrix to identify each element on boundary
    MapMatrixMixed<unsigned, FiniteElement*, Vector<int>*> boundary_identifier;


    // Temporary container to store pointers to temporary vectors
    // so they can be deleted. Creating storage to store these
    // temporary vectors of vectors of pointers previously defined
    Vector<Vector<int>*> tmp_vect_pt;

    // Loop over elements
    //-------------------
    unsigned nel = nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      // and put it in local storage fe_pt.
      FiniteElement* fe_pt = finite_element_pt(e);

      // print out values of all elements to doc
      if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;

      // Loop over the element's nodes and find out which boundaries they're on
      // ----------------------------------------------------------------------
      // Return number of nodes along one edge of the current element defined by
      // fe_pt hence, loop over nodes of each element in 3D-dimension
      unsigned nnode_1d = fe_pt->nnode_1d();

      // Loop over nodes in order
      for (unsigned i0 = 0; i0 < nnode_1d; i0++)
      {
        for (unsigned i1 = 0; i1 < nnode_1d; i1++)
        {
          for (unsigned i2 = 0; i2 < nnode_1d; i2++)
          {
            // Local node number, which is an assumed ordering defined in our
            // underlying finite element scheme.
            unsigned j = i0 + i1 * nnode_1d + i2 * nnode_1d * nnode_1d;

            // Get pointer to the set of boundaries that this node lives on
            // for each node reset boundaries_pt to 0,
            // create pointer to each node in each element
            // and give boundaries_pt a value
            std::set<unsigned>* boundaries_pt = 0;
            fe_pt->node_pt(j)->get_boundaries_pt(boundaries_pt);

            // If the node lives on some boundaries....
            // If not equal to 0, node is on boundary
            // Loop through values of the current node stored in boundaries_pt
            if (boundaries_pt != 0)
            {
              // Loop over the entries in the set "it" (name we use to refer to
              // the iterator)
              for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                   it != boundaries_pt->end();
                   ++it)
              {
                // What's the boundary?
                // Define boundary_id to have the value pointed to by
                // boundaries_pt
                unsigned boundary_id = *it;

                // Add pointer to finite element to vector for the appropriate
                // boundary

                // Does the pointer already exist in the vector
                Vector<FiniteElement*>::iterator b_el_it =
                  std::find(vector_of_boundary_element_pt[*it].begin(),
                            vector_of_boundary_element_pt[*it].end(),
                            fe_pt);

                // Only insert if we have not found it (i.e. got to the end)
                if (b_el_it == vector_of_boundary_element_pt[*it].end())
                {
                  vector_of_boundary_element_pt[*it].push_back(fe_pt);
                }

                // For the current element/boundary combination, create
                // a vector that stores an indicator which element boundaries
                // the node is located (boundary_identifier=-/+1 for nodes
                // on the left/right boundary; boundary_identifier=-/+2 for
                // nodes on the lower/upper boundary. We determine these indices
                // for all corner nodes of the element and add them to a vector
                // to a vector. This allows us to decide which faces of the
                // element coincide with the boundary since each face of the
                // element must have all four corner nodes on the boundary.
                if (boundary_identifier(boundary_id, fe_pt) == 0)
                {
                  // Here we make our vector of integers and keep track of the
                  // pointer to it The vector stores information about local
                  // node position for each boundary element
                  Vector<int>* tmp_pt = new Vector<int>;

                  // Add to the local scheme that will allow us to delete it
                  // later at the very end of the function
                  tmp_vect_pt.push_back(tmp_pt);

                  // Add the pointer to the storage scheme defined by
                  // boundary_id and pointer to finite element
                  boundary_identifier(boundary_id, fe_pt) = tmp_pt;
                }

                // Are we at a corner node? (local for each boundary element)
                // i0,i1,i2 represent the current 3D coordinates- which local
                // corner node.
                if (((i0 == 0) || (i0 == nnode_1d - 1)) &&
                    ((i1 == 0) || (i1 == nnode_1d - 1)) &&
                    ((i2 == 0) || (i2 == nnode_1d - 1)))
                {
                  // Create index to represent position relative to s_0
                  // s_0,s_1,s_2 are local 3D directions
                  (*boundary_identifier(boundary_id, fe_pt))
                    .push_back(1 * (2 * i0 / (nnode_1d - 1) - 1));

                  // Create index to represent position relative to s_1
                  (*boundary_identifier(boundary_id, fe_pt))
                    .push_back(2 * (2 * i1 / (nnode_1d - 1) - 1));

                  // Create index to represent position relative to s_2
                  (*boundary_identifier(boundary_id, fe_pt))
                    .push_back(3 * (2 * i2 / (nnode_1d - 1) - 1));
                }
              }
            }
          }
        }
      }
    }


    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    //---------------------
    for (unsigned i = 0; i < nbound; i++)
    {
      // Loop over elements on given boundary
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it = vector_of_boundary_element_pt[i].begin();
           it != vector_of_boundary_element_pt[i].end();
           it++)
      {
        // Recover pointer to element for each element
        FiniteElement* fe_pt = (*it);

        // Initialise count for boundary identities (-3,-2,-1,1,2,3)
        std::map<int, unsigned> count;

        // Loop over coordinates in 3D dimension
        for (int ii = 0; ii < 3; ii++)
        {
          // Loop over upper/lower end of coordinates
          // count -/+ for each direction separately
          for (int sign = -1; sign < 3; sign += 2)
          {
            count[(ii + 1) * sign] = 0;

            // Initialise map of counts to 0 before loop
            // count boundary indicator for each element
            // and its nodes
          }
        }

        // Loop over boundary indicators for this element/boundary
        // count nodes for any one fixed boundary determined by
        // boundary indicator
        unsigned n_indicators = (*boundary_identifier(i, fe_pt)).size();
        for (unsigned j = 0; j < n_indicators; j++)
        {
          count[(*boundary_identifier(i, fe_pt))[j]]++;
        }

        // Determine the correct boundary indicator by checking that it
        // occurs four times (since four corner nodes of the element's boundary
        // need to be located on the boundary domain)
        int indicator = -10;


        // Loop over coordinates
        for (int ii = 0; ii < 3; ii++)
        {
          // Loop over upper/lower end of coordinates
          for (int sign = -1; sign < 3; sign += 2)
          {
            // If an index occurs four times then a face is on the boundary
            // But we can have multiple faces on the same boundary, so add all
            // the ones that we find!
            if (count[(ii + 1) * sign] == 4)
            {
              indicator = (ii + 1) * sign;

              // Copy into the member data structures
              Boundary_element_pt[i].push_back(*it);
              Face_index_at_boundary[i].push_back(indicator);
            }
          }
        }
      }
    }

    // Doc?
    //-----
    if (doc)
    {
      // Loop over boundaries
      for (unsigned i = 0; i < nbound; i++)
      {
        unsigned nel = Boundary_element_pt[i].size();
        outfile << "Boundary: " << i << " is adjacent to " << nel << " elements"
                << std::endl;

        // Loop over elements on given boundary
        for (unsigned e = 0; e < nel; e++)
        {
          FiniteElement* fe_pt = Boundary_element_pt[i][e];
          outfile << "Boundary element:" << fe_pt
                  << " Face index along boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }


    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;


    // Cleanup temporary vectors
    unsigned n = tmp_vect_pt.size();
    for (unsigned i = 0; i < n; i++)
    {
      delete tmp_vect_pt[i];
    }

    return;


    // Doc?
    /*bool doc=false;
    if (outfile) doc=true;

    // Number of boundaries
    unsigned nbound=nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);


    // Temporary vector of vectors of pointers to elements on the boundaries:
    // vector_of_boundary_element_pt[i] is the vector of pointers to all
    // elements that have nodes on boundary i.
    // This is not a set to ensure UNIQUE ordering on every processor
    Vector<Vector<FiniteElement*> > vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed local coord for elements on boundary
    MapMatrixMixed<unsigned,FiniteElement*,Vector<int>* >
     boundary_identifier;

    // Tmp container to store pointers to tmp vectors so they can be deleted
    Vector<Vector<int>*> tmp_vect_pt;

    // Loop over elements
    //-------------------
    unsigned nel=nelement();
    for (unsigned e=0;e<nel;e++)
     {
      // Get pointer to element
      FiniteElement* fe_pt=finite_element_pt(e);

      if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;

      // Loop over the element's nodes and find out which boundaries they're on
      // ----------------------------------------------------------------------
      unsigned nnode_1d=fe_pt->nnode_1d();

      // Loop over nodes in order
      for (unsigned i0=0;i0<nnode_1d;i0++)
       {
        for (unsigned i1=0;i1<nnode_1d;i1++)
         {
          for (unsigned i2=0;i2<nnode_1d;i2++)
           {
            // Local node number
            unsigned j=i0+i1*nnode_1d+i2*nnode_1d*nnode_1d;

            // Get pointer to set of boundaries that this
            // node lives on
            std::set<unsigned>* boundaries_pt=0;
            fe_pt->node_pt(j)->get_boundaries_pt(boundaries_pt);

            // If the node lives on some boundaries....
            if (boundaries_pt!=0)
             {
              for (std::set<unsigned>::iterator it=boundaries_pt->begin();
                   it!=boundaries_pt->end();++it)
               {

                // What's the boundary?
                unsigned boundary_id=*it;

                // Add pointer to finite element to vector for the appropriate
                // boundary

                // Does the pointer already exits in the vector
                Vector<FiniteElement*>::iterator b_el_it =
                 std::find(vector_of_boundary_element_pt[*it].begin(),
                           vector_of_boundary_element_pt[*it].end(),
                           fe_pt);

                //Only insert if we have not found it (i.e. got to the end)
                if(b_el_it == vector_of_boundary_element_pt[*it].end())
                 {
                  vector_of_boundary_element_pt[*it].push_back(fe_pt);
                 }

                // For the current element/boundary combination, create
                // a vector that stores an indicator which element boundaries
                // the node is located (boundary_identifier=-/+1 for nodes
                // on the left/right boundary; boundary_identifier=-/+2 for
    nodes
                // on the lower/upper boundary. We determine these indices
                // for all corner nodes of the element and add them to a vector
                // to a vector. This allows us to decide which face of the
    element
                // coincides with the boundary since the (brick!) element must
                // have exactly four corner nodes on the boundary.
                if (boundary_identifier(boundary_id,fe_pt)==0)
                 {
                  Vector<int>* tmp_pt=new Vector<int>;
                  tmp_vect_pt.push_back(tmp_pt);
                  boundary_identifier(boundary_id,fe_pt)=tmp_pt;
                 }

                // Are we at a corner node?
                if (((i0==0)||(i0==nnode_1d-1))&&((i1==0)||(i1==nnode_1d-1))
                    &&((i2==0)||(i2==nnode_1d-1)))
                 {
                  // Create index to represent position relative to s_0
                  (*boundary_identifier(boundary_id,fe_pt)).
                   push_back(1*(2*i0/(nnode_1d-1)-1));

                  // Create index to represent position relative to s_1
                  (*boundary_identifier(boundary_id,fe_pt)).
                   push_back(2*(2*i1/(nnode_1d-1)-1));

                  // Create index to represent position relative to s_2
                  (*boundary_identifier(boundary_id,fe_pt)).
                   push_back(3*(2*i2/(nnode_1d-1)-1));
                 }
               }
             }
           }
         }
       }
     }


    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    //---------------------
    for (unsigned i=0;i<nbound;i++)
     {
      // Loop over elements on given boundary
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it=vector_of_boundary_element_pt[i].begin();
           it!=vector_of_boundary_element_pt[i].end();
           it++)
       {
        // Recover pointer to element
        FiniteElement* fe_pt=(*it);

        // Initialise count for boundary identiers (-3,-2,-1,1,2,3)
        std::map<int,unsigned> count;

        // Loop over coordinates
        for (int ii=0;ii<3;ii++)
         {
          // Loop over upper/lower end of coordinates
          for (int sign=-1;sign<3;sign+=2)
           {
            count[(ii+1)*sign]=0;
           }
         }

        // Loop over boundary indicators for this element/boundary
        unsigned n_indicators=(*boundary_identifier(i,fe_pt)).size();
        for (unsigned j=0;j<n_indicators;j++)
         {
          count[(*boundary_identifier(i,fe_pt))[j] ]++;
         }

        // Determine the correct boundary indicator by checking that it
        // occurs four times (since four corner nodes of the element's boundary
        // need to be located on the domain boundary
        int indicator=-10;

        //Check that we're finding exactly one boundary indicator
        bool found=false;

        // Loop over coordinates
        for (int ii=0;ii<3;ii++)
         {
          // Loop over upper/lower end of coordinates
          for (int sign=-1;sign<3;sign+=2)
           {
            if (count[(ii+1)*sign]==4)
             {
              // Check that we haven't found multiple boundaries
              if (found)
               {
                throw OomphLibError(
                 "Trouble: Multiple boundary identifiers!\n",
                 OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);
               }
              found=true;
              indicator=(ii+1)*sign;
             }
           }
         }

        // Check if we've found one boundary
        if (found)
         {

          // Store element
          Boundary_element_pt[i].push_back(*it);

          // Now convert boundary indicator into information required
          // for FaceElements
          switch (indicator)
           {
           case -3:

            // s_2 is fixed at -1.0:
            Face_index_at_boundary[i].push_back(-3);
            break;

           case -2:

            // s_1 is fixed at -1.0:
            Face_index_at_boundary[i].push_back(-2);
            break;

           case -1:

            // s_0 is fixed at -1.0:
            Face_index_at_boundary[i].push_back(-1);
            break;


           case 1:

            // s_0 is fixed at 1.0:
            Face_index_at_boundary[i].push_back(1);
            break;

           case 2:

            // s_1 is fixed at 1.0:
            Face_index_at_boundary[i].push_back(2);
            break;

           case 3:

            // s_2 is fixed at 1.0:
            Face_index_at_boundary[i].push_back(3);
            break;


           default:

            throw OomphLibError("Never get here",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
           }

         }
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
                  << " Face index along boundary is "
                  <<  Face_index_at_boundary[i][e] << std::endl;
         }
       }
     }


    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup=true;


    // Cleanup temporary vectors
    unsigned n=tmp_vect_pt.size();
    for (unsigned i=0;i<n;i++)
     {
      delete tmp_vect_pt[i];
      }*/
  }


} // namespace oomph
