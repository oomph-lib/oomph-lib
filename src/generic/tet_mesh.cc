// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
#include "tet_mesh.h"
#include "Telements.h"


namespace oomph
{
  //=======================================================================
  /// Constructor for a FacetedSurface created from a list of nodes
  /// and connectivity information. This is used in remeshing
  //=======================================================================
  TetMeshFacetedClosedSurfaceForRemesh::TetMeshFacetedClosedSurfaceForRemesh(
    Vector<Node*> const& vertex_node_pt,
    Vector<Vector<unsigned>> const& facet_connectivity,
    Vector<unsigned> const& facet_boundary_id)
    : TetMeshFacetedClosedSurface()
  {
    // Create the vertices
    unsigned n_vertex = vertex_node_pt.size();
    Vertex_pt.resize(n_vertex);
    for (unsigned v = 0; v < n_vertex; ++v)
    {
      Vertex_pt[v] = new TetMeshVertex(vertex_node_pt[v]);
    }

    // Create the facets
    unsigned n_facet = facet_connectivity.size();
    Facet_pt.resize(n_facet);
    for (unsigned f = 0; f < n_facet; ++f)
    {
      unsigned n_vertex_on_facet = facet_connectivity[f].size();
      Facet_pt[f] = new TetMeshFacet(n_vertex_on_facet);
      for (unsigned i = 0; i < n_vertex_on_facet; ++i)
      {
        Facet_pt[f]->set_vertex_pt(i, Vertex_pt[facet_connectivity[f][i]]);
      }
      // Add in the boundary id
      Facet_pt[f]->set_one_based_boundary_id(facet_boundary_id[f]);
    }
  }

  //=================================================================
  /// Destructor. Delete allocated memory
  //================================================================
  TetMeshFacetedClosedSurfaceForRemesh::~TetMeshFacetedClosedSurfaceForRemesh()
  {
    // Delete the facets and the vertices
    unsigned n_facet = this->nfacet();
    for (unsigned f = 0; f < n_facet; f++)
    {
      delete Facet_pt[f];
    }
    unsigned n_vertex = this->nvertex();
    for (unsigned v = 0; v < n_vertex; v++)
    {
      delete Vertex_pt[v];
    }
  }


  //================================================================
  /// Global static data that specifies the permitted
  /// error in the setup of the boundary coordinates
  //================================================================
  double TetMeshBase::Tolerance_for_boundary_finding = 1.0e-5;


  //================================================================
  /// Setup lookup schemes which establish which elements are located
  /// next to which boundaries (Doc to outfile if it's open).
  //================================================================
  void TetMeshBase::setup_boundary_element_info(std::ostream& outfile)
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
    Lookup_for_elements_next_boundary_is_setup = false;

    // Temporary vector of vectors of pointers to elements on the boundaries:
    // This is a vector to ensure UNIQUE ordering in all processors
    Vector<Vector<FiniteElement*>> vector_of_boundary_element_pt;
    vector_of_boundary_element_pt.resize(nbound);

    // Matrix map for working out the fixed face for elements on boundary
    MapMatrixMixed<unsigned, FiniteElement*, int> face_identifier;

    // Loop over elements
    //-------------------
    unsigned nel = nelement();


    // Get pointer to vector of boundaries that the
    // node lives on
    Vector<std::set<unsigned>*> boundaries_pt(4, 0);

    for (unsigned e = 0; e < nel; e++)
    {
      // Get pointer to element
      FiniteElement* fe_pt = finite_element_pt(e);

      if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;

      // Only include 3D elements! Some meshes contain interface elements too.
      if (fe_pt->dim() == 3)
      {
        // Loop over the element's nodes and find out which boundaries they're
        // on
        // ----------------------------------------------------------------------
        // We need only loop over the corner nodes
        for (unsigned i = 0; i < 4; i++)
        {
          fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
        }

        // Find the common boundaries of each face
        Vector<std::set<unsigned>> face(4);

        // NOTE: Face indices defined in Telements.h

        // Face 3 connnects points 0, 1 and 2
        if (boundaries_pt[0] && boundaries_pt[1] && boundaries_pt[2])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(face[3], face[3].begin()));
        }

        // Face 2 connects points 0, 1 and 3
        if (boundaries_pt[0] && boundaries_pt[1] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[2], face[2].begin()));
        }

        // Face 1 connects points 0, 2 and 3
        if (boundaries_pt[0] && boundaries_pt[2] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[0]->begin(),
            boundaries_pt[0]->end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[1], face[1].begin()));
        }

        // Face 0 connects points 1, 2 and 3
        if (boundaries_pt[1] && boundaries_pt[2] && boundaries_pt[3])
        {
          std::set<unsigned> aux;

          std::set_intersection(
            boundaries_pt[1]->begin(),
            boundaries_pt[1]->end(),
            boundaries_pt[2]->begin(),
            boundaries_pt[2]->end(),
            std::insert_iterator<std::set<unsigned>>(aux, aux.begin()));

          std::set_intersection(
            aux.begin(),
            aux.end(),
            boundaries_pt[3]->begin(),
            boundaries_pt[3]->end(),
            std::insert_iterator<std::set<unsigned>>(face[0], face[0].begin()));
        }


        // We now know whether any faces lay on the boundaries
        for (unsigned i = 0; i < 4; i++)
        {
          // How many boundaries are there
          unsigned count = 0;

          // The number of the boundary
          int boundary = -1;

          // Loop over all the members of the set and add to the count
          // and set the boundary
          for (std::set<unsigned>::iterator it = face[i].begin();
               it != face[i].end();
               ++it)
          {
            ++count;
            boundary = *it;
          }

          // If we're on more than one boundary, this is weird, so die
          if (count > 1)
          {
            std::ostringstream error_stream;
            fe_pt->output(error_stream);
            error_stream << "Face " << i << " is on " << count
                         << " boundaries.\n";
            error_stream << "This is rather strange.\n";
            error_stream << "Your mesh may be too coarse or your tet mesh\n";
            error_stream << "may be screwed up. I'm skipping the automated\n";
            error_stream << "setup of the elements next to the boundaries\n";
            error_stream << "lookup schemes.\n";
            OomphLibWarning(error_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
          }

          // If we have a boundary then add this to the appropriate set
          if (boundary >= 0)
          {
            // Does the pointer already exits in the vector
            Vector<FiniteElement*>::iterator b_el_it = std::find(
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .begin(),
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .end(),
              fe_pt);

            // Only insert if we have not found it (i.e. got to the end)
            if (b_el_it ==
                vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                  .end())
            {
              vector_of_boundary_element_pt[static_cast<unsigned>(boundary)]
                .push_back(fe_pt);
            }

            // Also set the fixed face
            face_identifier(static_cast<unsigned>(boundary), fe_pt) = i;
          }
        }

        // Now we set the pointers to the boundary sets to zero
        for (unsigned i = 0; i < 4; i++)
        {
          boundaries_pt[i] = 0;
        }
      }
    }

    // Now copy everything across into permanent arrays
    //-------------------------------------------------

    // Loop over boundaries
    //---------------------
    for (unsigned i = 0; i < nbound; i++)
    {
      // Number of elements on this boundary (currently stored in a set)
      unsigned nel = vector_of_boundary_element_pt[i].size();

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
                  << " Face index of boundary is "
                  << Face_index_at_boundary[i][e] << std::endl;
        }
      }
    }


    // Lookup scheme has now been setup yet
    Lookup_for_elements_next_boundary_is_setup = true;
  }


  //======================================================================
  /// Assess mesh quality: Ratio of max. edge length to min. height,
  /// so if it's very large it's BAAAAAD.
  //======================================================================
  void TetMeshBase::assess_mesh_quality(std::ofstream& some_file)
  {
    Vector<Vector<double>> edge(6);
    for (unsigned e = 0; e < 6; e++)
    {
      edge[e].resize(3);
    }
    unsigned nel = this->nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* fe_pt = this->finite_element_pt(e);
      for (unsigned i = 0; i < 3; i++)
      {
        edge[0][i] = fe_pt->node_pt(2)->x(i) - fe_pt->node_pt(1)->x(i);
        edge[1][i] = fe_pt->node_pt(0)->x(i) - fe_pt->node_pt(2)->x(i);
        edge[2][i] = fe_pt->node_pt(1)->x(i) - fe_pt->node_pt(0)->x(i);
        edge[3][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(0)->x(i);
        edge[4][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(1)->x(i);
        edge[5][i] = fe_pt->node_pt(3)->x(i) - fe_pt->node_pt(2)->x(i);
      }

      double max_length = 0.0;
      for (unsigned j = 0; j < 6; j++)
      {
        double length = 0.0;
        for (unsigned i = 0; i < 3; i++)
        {
          length += edge[j][i] * edge[j][i];
        }
        length = sqrt(length);
        if (length > max_length) max_length = length;
      }


      double min_height = DBL_MAX;
      for (unsigned j = 0; j < 4; j++)
      {
        Vector<double> normal(3);
        unsigned e0 = 0;
        unsigned e1 = 0;
        unsigned e2 = 0;
        switch (j)
        {
          case 0:
            e0 = 4;
            e1 = 5;
            e2 = 1;
            break;

          case 1:
            e0 = 1;
            e1 = 3;
            e2 = 2;
            break;

          case 2:
            e0 = 3;
            e1 = 4;
            e2 = 1;
            break;

          case 3:
            e0 = 1;
            e1 = 2;
            e2 = 3;
            break;

          default:

            oomph_info << "never get here\n";
            abort();
        }

        normal[0] = edge[e0][1] * edge[e1][2] - edge[e0][2] * edge[e1][1];
        normal[1] = edge[e0][2] * edge[e1][0] - edge[e0][0] * edge[e1][2];
        normal[2] = edge[e0][0] * edge[e1][1] - edge[e0][1] * edge[e1][0];
        double norm =
          normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
        double inv_norm = 1.0 / sqrt(norm);
        normal[0] *= inv_norm;
        normal[1] *= inv_norm;
        normal[2] *= inv_norm;

        double height = fabs(edge[e2][0] * normal[0] + edge[e2][1] * normal[1] +
                             edge[e2][2] * normal[2]);

        if (height < min_height) min_height = height;
      }

      double aspect_ratio = max_length / min_height;

      some_file << "ZONE N=4, E=1, F=FEPOINT, ET=TETRAHEDRON\n";
      for (unsigned j = 0; j < 4; j++)
      {
        for (unsigned i = 0; i < 3; i++)
        {
          some_file << fe_pt->node_pt(j)->x(i) << " ";
        }
        some_file << aspect_ratio << std::endl;
      }
      some_file << "1 2 3 4" << std::endl;
    }
  }


  //======================================================================
  /// Move the nodes on boundaries with associated Geometric Objects (if any)
  /// so that they exactly coincide with the geometric object. This requires
  /// that the boundary coordinates are set up consistently
  //======================================================================
  void TetMeshBase::snap_nodes_onto_geometric_objects()
  {
    // Backup in case elements get inverted
    std::map<Node*, Vector<double>> old_nodal_posn;
    std::map<Node*, Vector<double>> new_nodal_posn;
    unsigned nnod = nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      Node* nod_pt = node_pt(j);
      Vector<double> x(3);
      nod_pt->position(x);
      old_nodal_posn[nod_pt] = x;
    }

    // Loop over all boundaries
    unsigned n_bound = this->nboundary();
    for (unsigned b = 0; b < n_bound; b++)
    {
      bool do_it = true;

      // Accumulate reason for not snapping
      std::stringstream reason;
      reason << "Can't snap nodes on boundary " << b
             << " onto geom object because: \n";

      TetMeshFacetedSurface* faceted_surface_pt = 0;
      std::map<unsigned, TetMeshFacetedSurface*>::iterator it =
        Tet_mesh_faceted_surface_pt.find(b);
      if (it != Tet_mesh_faceted_surface_pt.end())
      {
        faceted_surface_pt = (*it).second;
      }

      // Facet associated with this boundary?
      if (faceted_surface_pt == 0)
      {
        reason << "-- no facets asssociated with boundary\n";
        do_it = false;
      }

      // Get geom object associated with facet
      GeomObject* geom_obj_pt = 0;
      if (do_it)
      {
        geom_obj_pt = faceted_surface_pt->geom_object_with_boundaries_pt();
        if (geom_obj_pt == 0)
        {
          reason << "-- no geom object associated with boundary\n";
          do_it = false;
        }
      }

      // Triangular facet?
      if (Triangular_facet_vertex_boundary_coordinate[b].size() == 0)
      {
        reason << "-- facet has to be triangular and vertex coordinates have\n"
               << "   to have been set up\n";
        do_it = false;
      }

      // We need boundary coordinates!
      if (!Boundary_coordinate_exists[b])
      {
        reason << "-- no boundary coordinates were set up\n";
        do_it = false;
      }


      // Which facet is associated with this boundary?
      unsigned facet_id_of_boundary = 0;
      TetMeshFacet* f_pt = 0;
      if (do_it)
      {
        unsigned nf = faceted_surface_pt->nfacet();
        for (unsigned f = 0; f < nf; f++)
        {
          if ((faceted_surface_pt->one_based_facet_boundary_id(f) - 1) == b)
          {
            facet_id_of_boundary = f;
            break;
          }
        }
        f_pt = faceted_surface_pt->facet_pt(facet_id_of_boundary);


        // Three vertices?
        unsigned nv = f_pt->nvertex();
        if (nv != 3)
        {
          reason << "-- number of facet vertices is " << nv
                 << " rather than 3\n";
          do_it = false;
        }

        // Have we set up zeta coordinates in geometric object?
        if ((f_pt->vertex_pt(0)->zeta_in_geom_object().size() != 2) ||
            (f_pt->vertex_pt(1)->zeta_in_geom_object().size() != 2) ||
            (f_pt->vertex_pt(2)->zeta_in_geom_object().size() != 2))
        {
          reason << "-- no boundary coordinates were set up\n";
          do_it = false;
        }
      }


      // Are we ready to go?
      if (!do_it)
      {
        const bool tell_us_why = false;
        if (tell_us_why)
        {
          oomph_info << reason.str() << std::endl;
        }
      }
      else
      {
        // Setup area coordinantes in triangular facet
        double x1 = Triangular_facet_vertex_boundary_coordinate[b][0][0];
        double y1 = Triangular_facet_vertex_boundary_coordinate[b][0][1];

        double x2 = Triangular_facet_vertex_boundary_coordinate[b][1][0];
        double y2 = Triangular_facet_vertex_boundary_coordinate[b][1][1];

        double x3 = Triangular_facet_vertex_boundary_coordinate[b][2][0];
        double y3 = Triangular_facet_vertex_boundary_coordinate[b][2][1];

        double detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);


        // Boundary coordinate (cartesian coordinates inside facet)
        Vector<double> zeta(2);

        // Loop over all nodes on that boundary
        const unsigned n_boundary_node = this->nboundary_node(b);
        for (unsigned n = 0; n < n_boundary_node; ++n)
        {
          // Get the boundary node and coordinates
          Node* const nod_pt = this->boundary_node_pt(b, n);
          nod_pt->get_coordinates_on_boundary(b, zeta);

          // Now we have zeta, the cartesian boundary coordinates
          // in the (assumed to be triangular) boundary facet; let's
          // work out the area coordinates
          // Notation as in
          // https://en.wikipedia.org/wiki/Barycentric_coordinate_system
          double s0 =
            ((y2 - y3) * (zeta[0] - x3) + (x3 - x2) * (zeta[1] - y3)) / detT;
          double s1 =
            ((y3 - y1) * (zeta[0] - x3) + (x1 - x3) * (zeta[1] - y3)) / detT;
          double s2 = 1.0 - s0 - s1;

          Vector<double> zeta_in_geom_obj(2, 0.0);
          Vector<double> position_from_geom_obj(3, 0.0);

          // Vertex zeta coordinates
          Vector<double> zeta_0(2);
          zeta_0 = f_pt->vertex_pt(0)->zeta_in_geom_object();

          Vector<double> zeta_1(2);
          zeta_1 = f_pt->vertex_pt(1)->zeta_in_geom_object();

          Vector<double> zeta_2(2);
          zeta_2 = f_pt->vertex_pt(2)->zeta_in_geom_object();


#ifdef PARANOID

          // Compute zeta values of the vertices from parametrisation of
          // boundaries
          double tol = 1.0e-12;
          Vector<double> zeta_from_boundary(2);
          Vector<double> zeta_vertex(2);
          for (unsigned v = 0; v < 3; v++)
          {
            zeta_vertex = f_pt->vertex_pt(v)->zeta_in_geom_object();
            for (unsigned alt = 0; alt < 2; alt++)
            {
              switch (v)
              {
                case 0:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta01(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta20(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  break;

                case 1:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta01(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta12(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  break;

                case 2:

                  if (alt == 0)
                  {
                    faceted_surface_pt->boundary_zeta12(
                      facet_id_of_boundary, 1.0, zeta_from_boundary);
                  }
                  else
                  {
                    faceted_surface_pt->boundary_zeta20(
                      facet_id_of_boundary, 0.0, zeta_from_boundary);
                  }
                  break;
              }

              double error =
                sqrt(pow((zeta_vertex[0] - zeta_from_boundary[0]), 2) +
                     pow((zeta_vertex[1] - zeta_from_boundary[1]), 2));
              if (error > tol)
              {
                std::ostringstream error_message;
                error_message
                  << "Error in parametrisation of boundary coordinates \n"
                  << "for vertex " << v << " [alt=" << alt << "] in facet "
                  << facet_id_of_boundary << " : \n"
                  << "zeta_vertex = [ " << zeta_vertex[0] << " "
                  << zeta_vertex[1] << " ] \n"
                  << "zeta_from_boundary      = [ " << zeta_from_boundary[0]
                  << " " << zeta_from_boundary[1] << " ] \n"
                  << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }

#endif

          // Compute zeta values of the interpolation parameters
          Vector<double> zeta_a(3, 0.0);
          Vector<double> zeta_b(3, 0.0);
          Vector<double> zeta_c(3, 0.0);
          Vector<double> zeta_d(3, 0.0);
          Vector<double> zeta_e(3, 0.0);
          Vector<double> zeta_f(3, 0.0);
          faceted_surface_pt->boundary_zeta01(facet_id_of_boundary, s1, zeta_a);
          faceted_surface_pt->boundary_zeta01(
            facet_id_of_boundary, 1.0 - s0, zeta_d);

          faceted_surface_pt->boundary_zeta12(facet_id_of_boundary, s2, zeta_c);
          faceted_surface_pt->boundary_zeta12(
            facet_id_of_boundary, 1.0 - s1, zeta_f);

          faceted_surface_pt->boundary_zeta20(
            facet_id_of_boundary, 1.0 - s2, zeta_b);
          faceted_surface_pt->boundary_zeta20(facet_id_of_boundary, s0, zeta_e);

          // Transfinite mapping
          zeta_in_geom_obj[0] = s0 * (zeta_a[0] + zeta_b[0] - zeta_0[0]) +
                                s1 * (zeta_c[0] + zeta_d[0] - zeta_1[0]) +
                                s2 * (zeta_e[0] + zeta_f[0] - zeta_2[0]);
          zeta_in_geom_obj[1] = s0 * (zeta_a[1] + zeta_b[1] - zeta_0[1]) +
                                s1 * (zeta_c[1] + zeta_d[1] - zeta_1[1]) +
                                s2 * (zeta_e[1] + zeta_f[1] - zeta_2[1]);

          unsigned n_tvalues =
            1 + nod_pt->position_time_stepper_pt()->nprev_values();
          for (unsigned t = 0; t < n_tvalues; ++t)
          {
            // Get the position according to the underlying geometric object
            geom_obj_pt->position(t, zeta_in_geom_obj, position_from_geom_obj);

            // Move the node
            for (unsigned i = 0; i < 3; i++)
            {
              nod_pt->x(t, i) = position_from_geom_obj[i];
            }
          }
        }
      }
    }

    // Check if any element is inverted
    bool some_element_is_inverted = false;
    unsigned count = 0;
    unsigned nel = nelement();
    for (unsigned e = 0; e < nel; e++)
    {
      FiniteElement* el_pt = finite_element_pt(e);
      bool passed = true;
      el_pt->check_J_eulerian_at_knots(passed);
      if (!passed)
      {
        some_element_is_inverted = true;
        char filename[100];
        std::ofstream some_file;
        sprintf(filename, "overly_distorted_element%i.dat", count);
        some_file.open(filename);
        unsigned nnod_1d = el_pt->nnode_1d();
        el_pt->output(some_file, nnod_1d);
        some_file.close();

        // Reset to old nodal position
        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          Vector<double> x_current(3);
          nod_pt->position(x_current);
          new_nodal_posn[nod_pt] = x_current;
          Vector<double> old_x(old_nodal_posn[nod_pt]);
          for (unsigned i = 0; i < 3; i++)
          {
            nod_pt->x(i) = old_x[i];
          }
        }

        // Plot
        sprintf(filename, "orig_overly_distorted_element%i.dat", count);
        some_file.open(filename);
        el_pt->output(some_file, nnod_1d);
        some_file.close();

        // Reset
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          for (unsigned i = 0; i < 3; i++)
          {
            nod_pt->x(i) = new_nodal_posn[nod_pt][i];
          }
        }

        // Bump
        count++;
      }
    }
    if (some_element_is_inverted)
    {
      std::ostringstream error_message;
      error_message
        << "A number of elements, namely: " << count
        << " are inverted after snapping. Their shapes are in "
        << " overly_distorted_element*.dat and "
           "orig_overly_distorted_element*.dat"
        << "Next person to get this error: Please implement a straightforward\n"
        << "variant of one of the functors in src/mesh_smoothing to switch\n"
        << "to harmonic mapping\n"
        << std::endl;
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      oomph_info << "No elements are inverted after snapping. Yay!"
                 << std::endl;
    }
  }


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


} // namespace oomph
