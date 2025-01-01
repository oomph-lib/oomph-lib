// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
#include "unstructured_two_d_mesh_geometry_base.h"
#include "triangle_mesh.h"

namespace oomph
{
#ifdef OOMPH_HAS_TRIANGLE_LIB

  //==============================================================
  /// Dump the triangulateio structure to a dump file and
  /// record boundary coordinates of boundary nodes
  //==============================================================
  void TriangleMeshBase::dump_triangulateio(std::ostream& dump_file)
  {
    TriangleHelper::dump_triangulateio(Triangulateio, dump_file);

#ifdef OOMPH_HAS_MPI
    // If the mesh is not distributed then process what follows
    if (!this->is_mesh_distributed())
    {
#endif // #ifdef OOMPH_HAS_MPI

      // Loop over all boundary nodes and dump out boundary coordinates
      // if they exist
      Vector<double> zeta(1);
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        if (Boundary_coordinate_exists[b])
        {
          dump_file << "1 # Boundary coordinate for boundary " << b
                    << " does exist\n";
          unsigned nnod = nboundary_node(b);
          dump_file << nnod << " # Number of dumped boundary nodes\n";
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* nod_pt = boundary_node_pt(b, j);
            nod_pt->get_coordinates_on_boundary(b, zeta);
            dump_file << zeta[0] << std::endl;
          }
          dump_file << "-999 # Done boundary coords for boundary " << b << "\n";
        }
        else
        {
          dump_file << "0 # Boundary coordinate for boundary " << b
                    << " does not exist\n";
        }
      }

#ifdef OOMPH_HAS_MPI
    }
#endif // #ifdef OOMPH_HAS_MPI
  }


  //==============================================================
  /// Regenerate the mesh from a dumped triangulateio file
  /// and dumped boundary coordinates of boundary nodes
  //==============================================================
  void TriangleMeshBase::remesh_from_triangulateio(std::istream& restart_file)
  {
#ifdef PARANOID
    // Record number of boundaries
    unsigned nbound_old = nboundary();
#endif

    // Clear the existing triangulate io
    TriangleHelper::clear_triangulateio(Triangulateio);

    // Read the data into the file
    TriangleHelper::read_triangulateio(restart_file, Triangulateio);

    // Now remesh from the new data structure
    this->remesh_from_internal_triangulateio();

#ifdef OOMPH_HAS_MPI
    // If the mesh is not distributed then process what follows
    if (!this->is_mesh_distributed())
    {
#endif // #ifdef OOMPH_HAS_MPI

#ifdef PARANOID
      // Record number of boundary nodes after remesh
      unsigned nbound_new = nboundary();
      if (nbound_new != nbound_old)
      {
        std::ostringstream error_stream;
        error_stream
          << "Number of boundaries before remesh from triangulateio, "
          << nbound_new << ",\ndoesn't match number boundaries afterwards, "
          << nbound_old
          << ". Have you messed \naround with boundary nodes in the "
          << "derived mesh constructor (or after calling \nit)? If so,"
          << " the dump/restart won't work as written at the moment.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // Loop over all boundary nodes and read boundary coordinates
      // if they exist
      Vector<double> zeta(1);
      std::string input_string;
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        // Read line up to termination sign
        getline(restart_file, input_string, '#');

        // Ignore rest of line
        restart_file.ignore(80, '\n');

        // Did boundary coordinate exist?
        const unsigned bound_coord_exists = atoi(input_string.c_str());
        if (bound_coord_exists == 1)
        {
          // Remember it!
          Boundary_coordinate_exists[b] = true;

          // Read line up to termination sign
          getline(restart_file, input_string, '#');

          // Ignore rest of line
          restart_file.ignore(80, '\n');

          // How many nodes did we dump?
          const unsigned nnod_dumped = atoi(input_string.c_str());

          // Does it match?
          unsigned nnod = nboundary_node(b);
          if (nnod != nnod_dumped)
          {
            std::ostringstream error_stream;
            error_stream << "Number of dumped boundary nodes " << nnod_dumped
                         << " doesn't match number of nodes on boundary " << b
                         << ": " << nnod << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          // Loop over all nodes
          for (unsigned j = 0; j < nnod; j++)
          {
            // Read line up to termination sign
            getline(restart_file, input_string);

            // Boundary coordinate
            zeta[0] = atof(input_string.c_str());

            // Set it
            Node* nod_pt = boundary_node_pt(b, j);
            nod_pt->set_coordinates_on_boundary(b, zeta);
          }

          // Read line up to termination sign
          getline(restart_file, input_string, '#');

          // Ignore rest of line
          restart_file.ignore(80, '\n');

          // Have we reached the end?
          const int check = atoi(input_string.c_str());
          if (check != -999)
          {
            std::ostringstream error_stream;
            error_stream << "Haven't read all nodes on boundary " << b
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
        else
        {
          oomph_info << "Restart: Boundary coordinate for boundary " << b
                     << " does not exist.\n";
        }
      }
#ifdef OOMPH_HAS_MPI
    } // if (!this->is_mesh_distributed())
#endif // #ifdef OOMPH_HAS_MPI
  }


  //==============================================================
  /// Write a Triangulateio_object file of the TriangulateIO object
  /// String s is add to assign a different value for
  /// input and/or output structure.
  /// The function give the same result of the "report" function
  /// included in the tricall.c, esternal_src.
  //==============================================================
  void TriangleMeshBase::write_triangulateio(TriangulateIO& triangle,
                                             std::string& s)
  {
    std::ofstream outfile;
    char filename[100];

    sprintf(filename, "Triangulateio_object_%s.dat", s.c_str());
    outfile.open(filename);
    outfile << "# Triangulateio object values:\n\n" << std::endl;

    // Write points coordinates
    if (triangle.numberofpoints != 0)
    {
      outfile << "# Triangulateio number of points is:"
              << triangle.numberofpoints << std::endl;
    }
    if (triangle.pointlist != NULL)
    {
      outfile << "# Vertex coordinates are:" << std::endl;
      for (int k = 0; k < triangle.numberofpoints * 2; k += 2)
      {
        outfile << (k * 0.5) + 1 << " " << triangle.pointlist[k] << " "
                << triangle.pointlist[k + 1] << std::endl;
      }
    }

    // Write points attribute list
    if (triangle.numberofpointattributes != 0)
    {
      outfile << "# Triangulateio number of points attributelist is:"
              << triangle.numberofpointattributes << std::endl;
    }
    if (triangle.pointattributelist != NULL)
    {
      outfile << "# Vertex attribute are:" << std::endl;
      for (int k = 0; k < triangle.numberofpointattributes; k++)
      {
        outfile << triangle.pointattributelist[k] << std::endl;
      }
    }

    // Write point markers list
    if (triangle.pointmarkerlist != NULL)
    {
      outfile << "# Vertex Markers are:" << std::endl;
      for (int k = 0; k < triangle.numberofpoints; k++)
      {
        outfile << triangle.pointmarkerlist[k] << std::endl;
      }
    }

    // Write the 1.node file used by the showme function
    std::ofstream nodefile;
    char nodename[100];

    sprintf(nodename, "file_%s.1.node", s.c_str());
    nodefile.open(nodename);
    nodefile << triangle.numberofpoints << " 2 "
             << triangle.numberofpointattributes << " 0" << std::endl;
    for (int j = 0; j < triangle.numberofpoints * 2; j += 2)
    {
      nodefile << (j / 2) + 1 << " " << triangle.pointlist[j] << " "
               << triangle.pointlist[j + 1] << std::endl;
    }
    nodefile.close();


    // Write segments edge elements
    if (triangle.numberofsegments != 0)
    {
      outfile << "# The number of segments is:" << triangle.numberofsegments
              << std::endl;
    }
    if (triangle.segmentlist != NULL)
    {
      outfile << "# Segments are:" << std::endl;
      for (int k = 0; k < triangle.numberofsegments * 2; k += 2)
      {
        outfile << triangle.segmentlist[k] << "  "
                << triangle.segmentlist[k + 1] << std::endl;
      }
    }

    // Write segments markers list
    if (triangle.segmentmarkerlist != NULL)
    {
      outfile << "# Segments Markers are:" << std::endl;
      for (int k = 0; k < triangle.numberofsegments; k++)
      {
        outfile << triangle.segmentmarkerlist[k] << std::endl;
      }
    }

    // Write regions
    if (triangle.numberofregions != 0)
    {
      outfile << "# The number of region is:" << triangle.numberofregions
              << std::endl;
    }

    // Write holes
    if (triangle.numberofholes != 0)
    {
      outfile << "# The number of holes is:" << triangle.numberofholes
              << std::endl;
    }
    if (triangle.holelist != NULL)
    {
      outfile << "#  Holes are:" << std::endl;
      for (int k = 0; k < triangle.numberofholes * 2; k += 2)
      {
        outfile << triangle.holelist[k] << "  " << triangle.holelist[k + 1]
                << std::endl;
      }
    }

    // Write triangles
    if (triangle.numberoftriangles != 0)
    {
      outfile << "# Triangulateio number of triangles:"
              << triangle.numberoftriangles << std::endl;
    }
    if (triangle.numberofcorners != 0)
    {
      outfile << "# Triangulateio number of corners:"
              << triangle.numberofcorners << std::endl;
    }
    if (triangle.numberoftriangleattributes != 0)
    {
      outfile << "# Triangulateio number of triangles attributes:"
              << triangle.numberoftriangleattributes << std::endl;
    }
    if (triangle.trianglelist != NULL)
    {
      outfile << "# Traingles are:" << std::endl;
      for (int k = 0; k < triangle.numberoftriangles * 3; k += 3)
      {
        outfile << triangle.trianglelist[k] << " "
                << triangle.trianglelist[k + 1] << " "
                << triangle.trianglelist[k + 2] << std::endl;
      }
    }

    if (triangle.trianglearealist != NULL)
    {
      outfile << "# Triangle's areas are:" << std::endl;
      for (int k = 0; k < triangle.numberoftriangles; k++)
      {
        outfile << triangle.trianglearealist[k] << std::endl;
      }
    }

    if (triangle.trianglelist != NULL)
    {
      // Write the 1.ele file used by the showme function
      std::ofstream elefile;
      char elename[100];

      sprintf(elename, "file_%s.1.ele", s.c_str());
      elefile.open(elename);
      elefile << triangle.numberoftriangles << " 3 0" << std::endl;
      for (int j = 0; j < triangle.numberoftriangles * 3; j += 3)
      {
        elefile << (j / 3) + 1 << " " << triangle.trianglelist[j] << " "
                << triangle.trianglelist[j + 1] << " "
                << triangle.trianglelist[j + 2] << std::endl;
      }
      elefile.close();
    }

    outfile.close();
  }

#endif

  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  void TriangleMeshBase::setup_boundary_element_info(std::ostream& outfile)
  {
    // Should we document the output here
    bool doc = false;

    if (outfile) doc = true;

    // Number of boundaries
    unsigned nbound = nboundary();

    // Wipe/allocate storage for arrays
    Boundary_element_pt.clear();
    Face_index_at_boundary.clear();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);

    // Temporary vector of vectors of pointers to elements on the boundaries:
    // This is a vector to ensure that order is strictly preserved
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed face for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, int> face_identifier;

    // Loop over elements
    //-------------------
    unsigned nel = nelement();

    // Get pointer to vector of boundaries that the
    // node lives on
    Vector<std::set<unsigned>*> boundaries_pt(3, 0);

    // Data needed to deal with edges through the
    // interior of the domain
    std::map<Edge, unsigned> edge_count;
    std::map<Edge, TriangleBoundaryHelper::BCInfo> edge_bcinfo;
    std::map<Edge, TriangleBoundaryHelper::BCInfo> face_info;
    MapMatrixMixed<unsigned, FiniteElement*, int> face_count;
    Vector<unsigned> bonus(nbound);

    // When using internal boundaries, an edge can be related to more than
    // one element (because of both sides of the internal boundaries)
    std::map<Edge, Vector<TriangleBoundaryHelper::BCInfo>> edge_internal_bnd;

    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      if (doc)
      {
        outfile << "Element: " << e << " " << fe_pt << std::endl;
      }

      // Only include 2D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 2)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------

        // We need only loop over the corner nodes
        for (unsigned i = 0; i < 3; i++)
        {
          fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
        }

        // Find the common boundaries of each edge
        Vector<std::set<unsigned>> edge_boundary(3);

        // Edge 0 connects points 1 and 2
        //-----------------------------

        if (boundaries_pt[1] && boundaries_pt[2])
        {
          // Create the corresponding edge
          Edge edge0(fe_pt->node_pt(1), fe_pt->node_pt(2));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 0;
          info.FE_pt = fe_pt;

          std::set_intersection(boundaries_pt[1]->begin(),
                                boundaries_pt[1]->end(),
                                boundaries_pt[2]->begin(),
                                boundaries_pt[2]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[0], edge_boundary[0].begin()));
          std::set<unsigned>::iterator it0 = edge_boundary[0].begin();

          // Edge does exist:
          if (edge_boundary[0].size() > 0)
          {
            info.Boundary = *it0;

            // How many times this edge has been visited
            edge_count[edge0]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge0, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge0].push_back(info);
          }
        }

        // Edge 1 connects points 0 and 2
        //-----------------------------

        if (boundaries_pt[0] && boundaries_pt[2])
        {
          std::set_intersection(boundaries_pt[0]->begin(),
                                boundaries_pt[0]->end(),
                                boundaries_pt[2]->begin(),
                                boundaries_pt[2]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[1], edge_boundary[1].begin()));

          // Create the corresponding edge
          Edge edge1(fe_pt->node_pt(0), fe_pt->node_pt(2));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 1;
          info.FE_pt = fe_pt;
          std::set<unsigned>::iterator it1 = edge_boundary[1].begin();

          // Edge does exist:
          if (edge_boundary[1].size() > 0)
          {
            info.Boundary = *it1;

            // How many times this edge has been visited
            edge_count[edge1]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge1, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge1].push_back(info);
          }
        }

        // Edge 2 connects points 0 and 1
        //-----------------------------

        if (boundaries_pt[0] && boundaries_pt[1])
        {
          std::set_intersection(boundaries_pt[0]->begin(),
                                boundaries_pt[0]->end(),
                                boundaries_pt[1]->begin(),
                                boundaries_pt[1]->end(),
                                std::insert_iterator<std::set<unsigned>>(
                                  edge_boundary[2], edge_boundary[2].begin()));

          // Create the corresponding edge
          Edge edge2(fe_pt->node_pt(0), fe_pt->node_pt(1));

          // Update infos about this edge
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = 2;
          info.FE_pt = fe_pt;
          std::set<unsigned>::iterator it2 = edge_boundary[2].begin();

          // Edge does exist:
          if (edge_boundary[2].size() > 0)
          {
            info.Boundary = *it2;

            // How many times this edge has been visited
            edge_count[edge2]++;

            // Update edge_bcinfo
            edge_bcinfo.insert(std::make_pair(edge2, info));

            // ... and also update the info associated with internal bnd
            edge_internal_bnd[edge2].push_back(info);
          }
        }


#ifdef PARANOID

        // Check if edge is associated with multiple boundaries

        // We now know whether any edges lay on the boundaries
        for (unsigned i = 0; i < 3; i++)
        {
          // How many boundaries are there
          unsigned count = 0;

          // Loop over all the members of the set and add to the count
          // and set the boundary
          for (std::set<unsigned>::iterator it = edge_boundary[i].begin();
               it != edge_boundary[i].end();
               ++it)
          {
            ++count;
          }

          // If we're on more than one boundary, this is weird, so die
          if (count > 1)
          {
            std::ostringstream error_stream;
            error_stream << "Edge " << i << " is located on " << count
                         << " boundaries.\n";
            error_stream << "This is rather strange, so I'm going to die\n";
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }

#endif

        // Now we set the pointers to the boundary sets to zero
        for (unsigned i = 0; i < 3; i++)
        {
          boundaries_pt[i] = 0;
        }
      }
    } // end of loop over all elements

    // Loop over all edges that are located on a boundary
    typedef std::map<Edge, TriangleBoundaryHelper::BCInfo>::iterator ITE;
    for (ITE it = edge_bcinfo.begin(); it != edge_bcinfo.end(); it++)
    {
      Edge current_edge = it->first;
      unsigned bound = it->second.Boundary;

      // If the edge has been visited only once
      if (edge_count[current_edge] == 1)
      {
        // Count the edges that are on the same element and on the same boundary
        face_count(static_cast<unsigned>(bound), it->second.FE_pt) =
          face_count(static_cast<unsigned>(bound), it->second.FE_pt) + 1;

        // If such edges exist, let store the corresponding element
        if (face_count(bound, it->second.FE_pt) > 1)
        {
          // Update edge's infos
          TriangleBoundaryHelper::BCInfo info;
          info.Face_id = it->second.Face_id;
          info.FE_pt = it->second.FE_pt;
          info.Boundary = it->second.Boundary;

          // Add it to FIinfo, that stores infos of problematic elements
          face_info.insert(std::make_pair(current_edge, info));

          // How many edges on which boundary have to be added
          bonus[bound]++;
        }
        else
        {
          // Add element and face to the appropriate vectors
          // Does the pointer already exits in the vector
          Vector<FiniteElement*>::iterator b_el_it = std::find(
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)].begin(),
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)].end(),
            it->second.FE_pt);

          // Only insert if we have not found it (i.e. got to the end)
          if (b_el_it ==
              vector_of_boundary_element_pt[static_cast<unsigned>(bound)].end())
          {
            vector_of_boundary_element_pt[static_cast<unsigned>(bound)]
              .push_back(it->second.FE_pt);
          }

          // set_of_boundary_element_pt[static_cast<unsigned>(bound)].insert(
          // it->second.FE_pt);
          face_identifier(static_cast<unsigned>(bound), it->second.FE_pt) =
            it->second.Face_id;
        }
      }

    } // End of "adding-boundaries"-loop


    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    for (unsigned i = 0; i < nbound; i++)
    {
      // Number of elements on this boundary that have to be added
      // in addition to other elements
      unsigned bonus1 = bonus[i];

      // Number of elements on this boundary
      unsigned nel = vector_of_boundary_element_pt[i].size() + bonus1;

      // Allocate storage for the coordinate identifiers
      Face_index_at_boundary[i].resize(nel);

      unsigned e_count = 0;
      typedef Vector<FiniteElement*>::iterator IT;
      for (IT it = vector_of_boundary_element_pt[i].begin();
           it != vector_of_boundary_element_pt[i].end();
           it++)
      {
        // Recover pointer to element
        FiniteElement* fe_pt = *it;

        // Add to permanent storage
        Boundary_element_pt[i].push_back(fe_pt);

        Face_index_at_boundary[i][e_count] = face_identifier(i, fe_pt);

        // Increment counter
        e_count++;
      }
      // We add the elements that have two or more edges on this boundary
      for (ITE itt = face_info.begin(); itt != face_info.end(); itt++)
      {
        if (itt->second.Boundary == i)
        {
          // Add to permanent storage
          Boundary_element_pt[i].push_back(itt->second.FE_pt);

          Face_index_at_boundary[i][e_count] = itt->second.Face_id;

          e_count++;
        }
      }

    } // End of loop over boundaries

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
                  << " Face index of boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }

    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;
  }
} // namespace oomph
