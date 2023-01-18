// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#include <algorithm>
#include "map_matrix.h"
#include "quad_mesh.h"


namespace oomph
{
  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  void QuadMeshBase::setup_boundary_element_info(std::ostream& outfile)
  {
    bool doc = false;
    if (outfile) doc = true;

    // Number of boundaries
    unsigned nbound = nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);

    // Temporary vector of vectors to pointers to elements on the boundaries:
    // This is not a set to ensure UNIQUE ordering
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed local coord for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, Vector<int>*> boundary_identifier;


    // Loop over elements
    //-------------------
    unsigned nel = nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;

      // Only include 2D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 2)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------
        unsigned nnode_1d = fe_pt->nnode_1d();

        // Loop over nodes in order
        for (unsigned i0 = 0; i0 < nnode_1d; i0++)
        {
          for (unsigned i1 = 0; i1 < nnode_1d; i1++)
          {
            // Local node number
            unsigned j = i0 + i1 * nnode_1d;

            // Get pointer to vector of boundaries that this
            // node lives on
            std::set<unsigned>* boundaries_pt = 0;
            fe_pt->node_pt(j)->get_boundaries_pt(boundaries_pt);

            // If the node lives on some boundaries....
            if (boundaries_pt != 0)
            {
              // Loop over boundaries
              // unsigned nbound=(*boundaries_pt).size();
              for (std::set<unsigned>::iterator it = boundaries_pt->begin();
                   it != boundaries_pt->end();
                   ++it)
              {
                // Add pointer to finite element to vector for the appropriate
                // boundary

                // Does the pointer already exits in the vector
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
                // to a vector. This allows us to decide which face of the
                // element coincides with the boundary since the (quad!) element
                // must have exactly two corner nodes on the boundary.
                if (boundary_identifier(*it, fe_pt) == 0)
                {
                  boundary_identifier(*it, fe_pt) = new Vector<int>;
                }

                // Are we at a corner node?
                if (((i0 == 0) || (i0 == nnode_1d - 1)) &&
                    ((i1 == 0) || (i1 == nnode_1d - 1)))
                {
                  // Create index to represent position relative to s_0
                  (*boundary_identifier(*it, fe_pt))
                    .push_back(1 * (2 * i0 / (nnode_1d - 1) - 1));

                  // Create index to represent position relative to s_1
                  (*boundary_identifier(*it, fe_pt))
                    .push_back(2 * (2 * i1 / (nnode_1d - 1) - 1));
                }
              }
            }
            // else
            // {
            //  oomph_info << "...does not live on any boundaries " <<
            //  std::endl;
            // }
          }
        }
      }
      // else
      //{
      // oomph_info << "Element " << e << " does not qualify" << std::endl;
      //}
    }

    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Note: vector_of_boundary_element_pt contains all elements
    // that have (at least) one corner node on a boundary -- can't copy
    // them across into Boundary_element_pt
    // yet because some of them might have only one node on the
    // the boundary in which case they don't qualify as
    // boundary elements!

    // Loop over boundaries
    //---------------------
    for (unsigned i = 0; i < nbound; i++)
    {
      // Number of elements on this boundary
      // nel is unused, so I've commented it out - RWhite.
      //   unsigned nel=vector_of_boundary_element_pt[i].size();

      // Allocate storage for the face identifiers
      // Face_index_at_boundary[i].resize(nel);

      // Loop over elements that have at least one corner node on this boundary
      //-----------------------------------------------------------------------
      // unsigned e_count=0;
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it = vector_of_boundary_element_pt[i].begin();
           it != vector_of_boundary_element_pt[i].end();
           it++)
      {
        // Recover pointer to element
        FiniteElement* fe_pt = *it;

        // Initialise count for boundary identiers (-2,-1,1,2)
        std::map<int, unsigned> count;

        // Loop over coordinates
        for (unsigned ii = 0; ii < 2; ii++)
        {
          // Loop over upper/lower end of coordinates
          for (int sign = -1; sign < 3; sign += 2)
          {
            count[(ii + 1) * sign] = 0;
          }
        }

        // Loop over boundary indicators for this element/boundary
        unsigned n_indicators = (*boundary_identifier(i, fe_pt)).size();
        for (unsigned j = 0; j < n_indicators; j++)
        {
          count[(*boundary_identifier(i, fe_pt))[j]]++;
        }
        delete boundary_identifier(i, fe_pt);

        // Determine the correct boundary indicator by checking that it
        // occurs twice (since two corner nodes of the element's boundary
        // need to be located on the domain boundary
        int indicator = -10;

        // Check that we're finding exactly one boundary indicator
        // bool found=false;

        // Loop over coordinates
        for (unsigned ii = 0; ii < 2; ii++)
        {
          // Loop over upper/lower end of coordinates
          for (int sign = -1; sign < 3; sign += 2)
          {
            // If an index occurs twice then that face is on the boundary
            // But we can have multiple faces on the same boundary id, so just
            // add all the ones that we have
            if (count[(ii + 1) * sign] == 2)
            {
              // Check that we haven't found multiple boundaries
              /*if (found)
               {
                std::string error_message=
                 "Trouble: Multiple boundary identifiers!\n";
                error_message +=
                 "Elements should only have at most 2 corner ";
                error_message +=
                 "nodes on any one boundary.\n";

                throw OomphLibError(
                 error_message,
                 OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);
               }
               found=true;*/
              indicator = (ii + 1) * sign;

              // Copy into the data structure
              Boundary_element_pt[i].push_back(*it);
              Face_index_at_boundary[i].push_back(indicator);
            }
          }
        }

        // Element has exactly two corner nodes on boundary
        /*if (found)
         {
          // Add to permanent storage
          Boundary_element_pt[i].push_back(fe_pt);

          // Now convert boundary indicator into information required
          // for FaceElements
          switch (indicator)
           {
            //South face
           case -2:

            // s_1 is fixed at -1.0:
            Face_index_at_boundary[i][e_count] = -2;
            break;

            //West face
           case -1:

            // s_0 is fixed at -1.0:
            Face_index_at_boundary[i][e_count] = -1;
            break;

            //East face
           case 1:

            // s_0 is fixed at 1.0:
            Face_index_at_boundary[i][e_count] = 1;
            break;

            //North face
           case 2:

            // s_1 is fixed at 1.0:
            Face_index_at_boundary[i][e_count] = 2;
            break;

           default:

            throw OomphLibError("Never get here",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
                                }

          // Increment counter
          e_count++;
          }*/
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
                  << " Face index on boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }


    // Lookup scheme has now been setup
    Lookup_for_elements_next_boundary_is_setup = true;
  }

} // namespace oomph
