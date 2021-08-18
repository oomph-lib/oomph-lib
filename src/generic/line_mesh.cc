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
// oomph-lib headers
#include "map_matrix.h"
#include "line_mesh.h"

namespace oomph
{
  //=======================================================================
  /// Set up lookup schemes which establish which elements are located
  /// next to which boundaries (doc to outfile if it's open)
  //=======================================================================
  void LineMeshBase::setup_boundary_element_info(std::ostream& outfile)
  {
    // Initialise documentation flag
    bool doc = false;

    // Set this to true if an open file has been passed to the function
    if (outfile)
    {
      doc = true;
    }

    // Determine number of boundaries in mesh
    const unsigned n_bound = nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(n_bound);
    Face_index_at_boundary.resize(n_bound);

    // Matrix map for working out the fixed local coord for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, Vector<int>*> boundary_identifier;

    // Determine number of elements in the mesh
    const unsigned n_element = nelement();

    // Loop over elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      // Output information to output file
      if (doc)
      {
        outfile << "Element: " << e << " " << fe_pt << std::endl;
      }

      // Only include 1D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 1)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------

        // Determine number of nodes in the element
        const unsigned n_node = fe_pt->nnode_1d();

        // Loop over nodes in order
        for (unsigned n = 0; n < n_node; n++)
        {
          // Allocate storage for pointer to set of boundaries that node lives
          // on
          std::set<unsigned>* boundaries_pt = 0;

          // Get pointer to vector of boundaries that this node lives on
          fe_pt->node_pt(n)->get_boundaries_pt(boundaries_pt);

          // If the node lives on some boundaries....
          if (boundaries_pt != 0)
          {
            // Determine number of boundaries which node lives on
            const unsigned n_bound_node_lives_on = (*boundaries_pt).size();

            // Throw an error if the node lives on more than one boundary
            if (n_bound_node_lives_on > 1)
            {
              std::string error_message =
                "In a 1D mesh a node shouldn't be able to live on more than\n";
              error_message += "one boundary, yet this node claims to.";

              throw OomphLibError(error_message,
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
            // If the node lives on just one boundary
            else if (n_bound_node_lives_on == 1)
            {
              // Determine which boundary the node lives on
              const std::set<unsigned>::iterator boundary =
                boundaries_pt->begin();

              // In 1D if an element has any nodes on a boundary then it must
              // be a boundary element. This means that (unlike in the 2D and
              // 3D cases) we can immediately add this element to permanent
              // storage.
              Boundary_element_pt[*boundary].push_back(fe_pt);

              // Record information required for FaceElements.
              // `Face_index_at_boundary' = -/+1 for nodes on the left/right
              // boundary. This allows us to decide which edge of the element
              // coincides with the boundary since the line element must have
              // precisely one vertex node on the boundary.

              // Are we at the left-hand vertex node? (left face)
              if (n == 0)
              {
                Face_index_at_boundary[*boundary].push_back(-1);
              }

              // Are we at the right-hand vertex node? (right face)
              else if (n == n_node - 1)
              {
                Face_index_at_boundary[*boundary].push_back(1);
              }
            }
          } // End of if node lives on some boundaries

        } // End of loop over nodes
      }
    } // End of loop over elements

#ifdef PARANOID

    // Check each boundary only has precisely one element next to it
    // -------------------------------------------------------------
    // Only if not distributed
#ifdef OOMPH_HAS_MPI
    if (Comm_pt == 0)
#endif
    {
      // Loop over boundaries
      for (unsigned b = 0; b < n_bound; b++)
      {
        // Evaluate number of elements adjacent to boundary b
        const unsigned n_element = Boundary_element_pt[b].size();

        switch (n_element)
        {
            // Boundary b has no adjacent elements
          case 0:
          {
            std::ostringstream error_stream;
            error_stream << "Boundary " << b
                         << " has no element adjacent to it\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
          }
          // Boundary b has one adjacent element (this is good!)
          case 1:
            break;

            // Boundary b has more than one adjacent element
          default:
          {
            std::ostringstream error_stream;
            error_stream << "Boundary " << b << " has " << n_element
                         << " elements adjacent to it.\n"
                         << "This shouldn't occur in a 1D mesh.\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
            break;
          }
        } // End of switch

        // Because each boundary should only have one element adjacent to it,
        // each `Face_index_at_boundary[b]' should be of size one.

        const unsigned face_index_at_boundary_size =
          Face_index_at_boundary[b].size();

        if (face_index_at_boundary_size != 1)
        {
          std::ostringstream error_stream;
          error_stream
            << "Face_index_at_boundary[" << b << "] has size"
            << face_index_at_boundary_size << " which does not make sense.\n"
            << "In a 1D mesh its size should always be one since only\n"
            << "one element can be adjacent to any particular boundary";
          throw OomphLibError(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      } // End of loop over boundaries
    }
#endif

    // Doc?
    // ----
    if (doc)
    {
      // Loop over boundaries
      for (unsigned b = 0; b < n_bound; b++)
      {
        const unsigned n_element = Boundary_element_pt[b].size();
        outfile << "Boundary: " << b << " is adjacent to " << n_element
                << " elements" << std::endl;

        // Loop over elements on given boundary
        for (unsigned e = 0; e < n_element; e++)
        {
          FiniteElement* fe_pt = Boundary_element_pt[b][e];
          outfile << "Boundary element:" << fe_pt
                  << " Face index on boundary is "
                  << Face_index_at_boundary[b][e] << std::endl;
        }
      }
    } // End of if(doc)


    // Lookup scheme has now been set up
    Lookup_for_elements_next_boundary_is_setup = true;

  } // End of setup_boundary_element_info()

} // namespace oomph
