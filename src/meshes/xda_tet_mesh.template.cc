// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_XDA_TET_MESH_TEMPLATE_CC
#define OOMPH_XDA_TET_MESH_TEMPLATE_CC


#include "../generic/Telements.h"
#include "xda_tet_mesh.template.h"


namespace oomph
{
  //======================================================================
  /// Constructor: Pass name of xda file. Note that all
  /// boundary elements get their own ID -- this is required for
  /// FSI problems. The vector containing the oomph-lib
  /// boundary IDs of all oomph-lib boundaries that collectively form
  /// a given boundary as specified in the xda input file can be
  /// obtained from the access function oomph_lib_boundary_ids(...).
  /// Timestepper defaults to steady pseudo-timestepper.
  //======================================================================
  template<class ELEMENT>
  XdaTetMesh<ELEMENT>::XdaTetMesh(const std::string xda_file_name,
                                  TimeStepper* time_stepper_pt)
  {
    // Mesh can only be built with 3D Telements.
    MeshChecker::assert_geometric_element<TElementGeometricBase, ELEMENT>(3);

    // Open and process xda input file
    std::ifstream infile(xda_file_name.c_str(), std::ios_base::in);
    unsigned n_node;
    unsigned n_element;
    unsigned n_bound_face;

    if (!infile.is_open())
    {
      std::ostringstream error_stream;
      error_stream << "Failed to open " << xda_file_name << "\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Dummy storage to jump lines
    char dummy[101];

    // Ignore file format
    infile.getline(dummy, 100);

    // Get number of elements
    infile >> n_element;

    // Ignore rest of line
    infile.getline(dummy, 100);

    // Get number of nodes
    infile >> n_node;

    // Ignore rest of line
    infile.getline(dummy, 100);

    // Ignore sum of element weights (whatever that is...)
    infile.getline(dummy, 100);

    // Get number of enumerated boundary faces on which boundary conditions
    // are applied.
    infile >> n_bound_face;

    // Keep reading until "Title String"
    while (dummy[0] != 'T')
    {
      infile.getline(dummy, 100);
    }

    // Make space for nodes and elements
    Node_pt.resize(n_node);
    Element_pt.resize(n_element);

    // Read first line with node labels and count them
    std::string line;
    std::getline(infile, line);
    std::istringstream ostr(line);
    std::istream_iterator<std::string> it(ostr);
    std::istream_iterator<std::string> end;
    unsigned nnod_el = 0;
    Vector<unsigned> first_node;
    while (it != end)
    {
      first_node.push_back(atoi((*it).c_str()));
      it++;
      nnod_el++;
    }

    // Check
    if (nnod_el != 10)
    {
      std::ostringstream error_stream;
      error_stream
        << "XdaTetMesh can currently only be built with quadratic tets "
        << "with 10 nodes. The specified mesh file has " << nnod_el
        << "nodes per element.\n";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // Storage for the global node numbers listed element-by-element
    Vector<unsigned> global_node(n_element * nnod_el);

    // Read in nodes
    unsigned k = 0;

    // Copy across first nodes
    for (unsigned j = 0; j < nnod_el; j++)
    {
      global_node[k] = first_node[k];
      k++;
    }

    // Read the other ones
    for (unsigned i = 1; i < n_element; i++)
    {
      for (unsigned j = 0; j < nnod_el; j++)
      {
        infile >> global_node[k];
        k++;
      }
      infile.getline(dummy, 100);
    }


    // Create storage for coordinates
    Vector<double> x_node(n_node);
    Vector<double> y_node(n_node);
    Vector<double> z_node(n_node);

    // Get nodal coordinates
    for (unsigned i = 0; i < n_node; i++)
    {
      infile >> x_node[i];
      infile >> y_node[i];
      infile >> z_node[i];
    }


    // Read in boundaries for faces
    unsigned element_nmbr;
    unsigned bound_id;
    unsigned side_nmbr;
    unsigned max_bound = 0;

    // Make space for enumeration of sub-boundaries
    Boundary_id.resize(n_bound_face + 1);

    // Counter for number of separate boundary faces
    unsigned count = 0;
    Vector<std::set<unsigned>> boundary_id(n_node);
    for (unsigned i = 0; i < n_bound_face; i++)
    {
      // Number of the element
      infile >> element_nmbr;

      // Which side/face on the tet are we dealing with (xda enumeratation)?
      infile >> side_nmbr;

      // What's the boundary ID?
      infile >> bound_id;

      // Turn into zero-based oomph-lib mesh boundary id
      unsigned oomph_lib_bound_id = bound_id - 1;
      oomph_lib_bound_id = count;
      Boundary_id[bound_id].push_back(count);

      // Increment number of separate boundary faces
      count++;

      // Get ready for allocation of total number of boundaries
      if (oomph_lib_bound_id > max_bound) max_bound = oomph_lib_bound_id;

      // Identify the "side nodes" (i.e. the nodes on the faces of
      // the bulk tet) according to the
      // conventions in '.xda' mesh files so that orientation of the
      // faces is always the same (required for computation of
      // outer unit normals
      Vector<unsigned> side_node(6);
      switch (side_nmbr)
      {
        case 0:
          side_node[0] = global_node[nnod_el * element_nmbr + 1];
          side_node[1] = global_node[nnod_el * element_nmbr];
          side_node[2] = global_node[nnod_el * element_nmbr + 2];
          side_node[3] = global_node[nnod_el * element_nmbr + 4];
          side_node[4] = global_node[nnod_el * element_nmbr + 6];
          side_node[5] = global_node[nnod_el * element_nmbr + 5];
          break;

        case 1:
          side_node[0] = global_node[nnod_el * element_nmbr];
          side_node[1] = global_node[nnod_el * element_nmbr + 1];
          side_node[2] = global_node[nnod_el * element_nmbr + 3];
          side_node[3] = global_node[nnod_el * element_nmbr + 4];
          side_node[4] = global_node[nnod_el * element_nmbr + 8];
          side_node[5] = global_node[nnod_el * element_nmbr + 7];
          break;

        case 2:
          side_node[0] = global_node[nnod_el * element_nmbr + 1];
          side_node[1] = global_node[nnod_el * element_nmbr + 2];
          side_node[2] = global_node[nnod_el * element_nmbr + 3];
          side_node[3] = global_node[nnod_el * element_nmbr + 5];
          side_node[4] = global_node[nnod_el * element_nmbr + 9];
          side_node[5] = global_node[nnod_el * element_nmbr + 8];
          break;

        case 3:
          side_node[0] = global_node[nnod_el * element_nmbr + 2];
          side_node[1] = global_node[nnod_el * element_nmbr];
          side_node[2] = global_node[nnod_el * element_nmbr + 3];
          side_node[3] = global_node[nnod_el * element_nmbr + 6];
          side_node[4] = global_node[nnod_el * element_nmbr + 7];
          side_node[5] = global_node[nnod_el * element_nmbr + 9];
          break;

        default:
          throw OomphLibError(
            "Unexcpected side number in your '.xda' input file\n",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Associate boundaries with nodes
      for (unsigned j = 0; j < 6; j++)
      {
        boundary_id[side_node[j]].insert(oomph_lib_bound_id);
      }
    }

    // Set number of boundaries
    set_nboundary(max_bound + 1);

    // Vector of bools, will tell us if we already visited a node
    std::vector<bool> done(n_node, false);

    Vector<unsigned> translate(n_node);
    translate[0] = 0;
    translate[1] = 2;
    translate[2] = 1;
    translate[3] = 3;
    translate[4] = 5;
    translate[5] = 7;
    translate[6] = 4;
    translate[7] = 6;
    translate[8] = 8;
    translate[9] = 9;

    // Create the elements
    unsigned node_count = 0;
    for (unsigned e = 0; e < n_element; e++)
    {
      Element_pt[e] = new ELEMENT;

      // Loop over all nodes
      for (unsigned j = 0; j < nnod_el; j++)
      {
        unsigned j_global = global_node[node_count];
        if (done[j_global] == false)
        {
          if (boundary_id[j_global].size() == 0)
          {
            Node_pt[j_global] = finite_element_pt(e)->construct_node(
              translate[j], time_stepper_pt);
          }
          else
          {
            Node_pt[j_global] = finite_element_pt(e)->construct_boundary_node(
              translate[j], time_stepper_pt);
            for (std::set<unsigned>::iterator it =
                   boundary_id[j_global].begin();
                 it != boundary_id[j_global].end();
                 it++)
            {
              add_boundary_node(*it, Node_pt[j_global]);
            }
          }
          done[j_global] = true;
          Node_pt[j_global]->x(0) = x_node[j_global];
          Node_pt[j_global]->x(1) = y_node[j_global];
          Node_pt[j_global]->x(2) = z_node[j_global];
        }
        else
        {
          finite_element_pt(e)->node_pt(translate[j]) = Node_pt[j_global];
        }
        node_count++;
      }
    }

    // Figure out which elements are next to the boundaries.
    setup_boundary_element_info();


    // Setup boundary coordinates
    unsigned nb = nboundary();
    for (unsigned b = 0; b < nb; b++)
    {
      bool switch_normal = false;
      setup_boundary_coordinates(b, switch_normal);
    }
  }


  //======================================================================
  /// Setup boundary coordinate on boundary b while is
  /// temporarily flattened to simplex faces. Boundary coordinates are the
  /// x-y coordinates in the plane of that boundary with the
  /// x-axis along the line from the (lexicographically)
  /// "lower left" to the "upper right" node. The y axis
  /// is obtained by taking the cross-product of the positive
  /// x direction with the outer unit normal computed by
  /// the face elements (or its negative if switch_normal is set
  /// to true). Doc faces in output file.
  //======================================================================
  template<class ELEMENT>
  void XdaTetMesh<ELEMENT>::setup_boundary_coordinates(
    const unsigned& b, const bool& switch_normal, std::ofstream& outfile)
  {
    // Temporary storage for face elements
    Vector<FiniteElement*> face_el_pt;

    // Backup for nodal positions
    std::map<Node*, Vector<double>> backup_position;

    // Loop over all elements on boundaries
    unsigned nel = this->nboundary_element(b);
    if (nel > 0)
    {
      // Loop over the bulk elements adjacent to boundary b
      for (unsigned e = 0; e < nel; e++)
      {
        // Get pointer to the bulk element that is adjacent to boundary b
        FiniteElement* bulk_elem_pt = this->boundary_element_pt(b, e);

        // Find the index of the face of element e along boundary b
        int face_index = this->face_index_at_boundary(b, e);

        // Create new face element
        DummyFaceElement<ELEMENT>* el_pt =
          new DummyFaceElement<ELEMENT>(bulk_elem_pt, face_index);
        face_el_pt.push_back(el_pt);

        // Backup nodal position
        Vector<double> pos(3);
        for (unsigned j = 3; j < 6; j++)
        {
          if (backup_position[el_pt->node_pt(j)].size() == 0)
          {
            el_pt->node_pt(j)->position(pos);
            backup_position[el_pt->node_pt(j)] = pos;
          }
        }

        // Temporarily flatten the element to a simplex
        for (unsigned i = 0; i < 3; i++)
        {
          // Node 3 is between vertex nodes 0 and 1
          el_pt->node_pt(3)->x(i) =
            0.5 * (el_pt->node_pt(0)->x(i) + el_pt->node_pt(1)->x(i));

          // Node 4 is between vertex nodes 1 and 2
          el_pt->node_pt(4)->x(i) =
            0.5 * (el_pt->node_pt(1)->x(i) + el_pt->node_pt(2)->x(i));

          // Node 5 is between vertex nodes 2 and 0
          el_pt->node_pt(5)->x(i) =
            0.5 * (el_pt->node_pt(2)->x(i) + el_pt->node_pt(0)->x(i));
        }


        // Output faces?
        if (outfile.is_open())
        {
          face_el_pt[face_el_pt.size() - 1]->output(outfile);
        }
      }


      // Loop over all nodes to find the lower left and upper right ones
      Node* lower_left_node_pt = this->boundary_node_pt(b, 0);
      Node* upper_right_node_pt = this->boundary_node_pt(b, 0);

      // Loop over all nodes on the boundary
      unsigned nnod = this->nboundary_node(b);
      for (unsigned j = 0; j < nnod; j++)
      {
        // Get node
        Node* nod_pt = this->boundary_node_pt(b, j);

        // Primary criterion for lower left: Does it have a lower z-coordinate?
        if (nod_pt->x(2) < lower_left_node_pt->x(2))
        {
          lower_left_node_pt = nod_pt;
        }
        // ...or is it a draw?
        else if (nod_pt->x(2) == lower_left_node_pt->x(2))
        {
          // If it's a draw: Does it have a lower y-coordinate?
          if (nod_pt->x(1) < lower_left_node_pt->x(1))
          {
            lower_left_node_pt = nod_pt;
          }
          // ...or is it a draw?
          else if (nod_pt->x(1) == lower_left_node_pt->x(1))
          {
            // If it's a draw: Does it have a lower x-coordinate?
            if (nod_pt->x(0) < lower_left_node_pt->x(0))
            {
              lower_left_node_pt = nod_pt;
            }
          }
        }

        // Primary criterion for upper right: Does it have a higher
        // z-coordinate?
        if (nod_pt->x(2) > upper_right_node_pt->x(2))
        {
          upper_right_node_pt = nod_pt;
        }
        // ...or is it a draw?
        else if (nod_pt->x(2) == upper_right_node_pt->x(2))
        {
          // If it's a draw: Does it have a higher y-coordinate?
          if (nod_pt->x(1) > upper_right_node_pt->x(1))
          {
            upper_right_node_pt = nod_pt;
          }
          // ...or is it a draw?
          else if (nod_pt->x(1) == upper_right_node_pt->x(1))
          {
            // If it's a draw: Does it have a higher x-coordinate?
            if (nod_pt->x(0) > upper_right_node_pt->x(0))
            {
              upper_right_node_pt = nod_pt;
            }
          }
        }
      }

      // Prepare storage for boundary coordinates
      Vector<double> zeta(2);

      // Unit vector connecting lower left and upper right nodes
      Vector<double> b0(3);
      b0[0] = upper_right_node_pt->x(0) - lower_left_node_pt->x(0);
      b0[1] = upper_right_node_pt->x(1) - lower_left_node_pt->x(1);
      b0[2] = upper_right_node_pt->x(2) - lower_left_node_pt->x(2);

      // Normalise
      double inv_length =
        1.0 / sqrt(b0[0] * b0[0] + b0[1] * b0[1] + b0[2] * b0[2]);
      b0[0] *= inv_length;
      b0[1] *= inv_length;
      b0[2] *= inv_length;

      // Get (outer) unit normal to first face element
      Vector<double> normal(3);
      Vector<double> s(2, 0.0);
      dynamic_cast<DummyFaceElement<ELEMENT>*>(face_el_pt[0])
        ->outer_unit_normal(s, normal);

      if (switch_normal)
      {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
      }

#ifdef PARANOID

      // Check that all elements have the same normal
      for (unsigned e = 0; e < nel; e++)
      {
        // Get (outer) unit normal to face element
        Vector<double> my_normal(3);
        dynamic_cast<DummyFaceElement<ELEMENT>*>(face_el_pt[0])
          ->outer_unit_normal(s, my_normal);

        // Dot product should be one!
        double error = normal[0] * my_normal[0] + normal[1] * my_normal[1] +
                       normal[2] * my_normal[2];

        error -= 1.0;
        if (switch_normal) error += 1.0;

        if (error > Tolerance_for_boundary_finding)
        {
          std::ostringstream error_message;
          error_message
            << "Error in alingment of normals (dot product-1) " << error
            << " for element " << e << std::endl
            << "exceeds tolerance specified by the static member data\n"
            << "TetMeshBase::Tolerance_for_boundary_finding = "
            << Tolerance_for_boundary_finding << std::endl
            << "This usually means that the boundary is not planar.\n\n"
            << "You can suppress this error message by recompiling \n"
            << "recompiling without PARANOID or by changing the tolerance.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif

      // Cross-product to get second in-plane vector, normal to b0
      Vector<double> b1(3);
      b1[0] = b0[1] * normal[2] - b0[2] * normal[1];
      b1[1] = b0[2] * normal[0] - b0[0] * normal[2];
      b1[2] = b0[0] * normal[1] - b0[1] * normal[0];

      // Assign boundary coordinates: projection onto the axes
      for (unsigned j = 0; j < nnod; j++)
      {
        // Get node
        Node* nod_pt = this->boundary_node_pt(b, j);

        // Difference vector to lower left corner
        Vector<double> delta(3);
        delta[0] = nod_pt->x(0) - lower_left_node_pt->x(0);
        delta[1] = nod_pt->x(1) - lower_left_node_pt->x(1);
        delta[2] = nod_pt->x(2) - lower_left_node_pt->x(2);

        // Project
        zeta[0] = delta[0] * b0[0] + delta[1] * b0[1] + delta[2] * b0[2];
        zeta[1] = delta[0] * b1[0] + delta[1] * b1[1] + delta[2] * b1[2];

#ifdef PARANOID

        // Check:
        double error = pow(lower_left_node_pt->x(0) + zeta[0] * b0[0] +
                             zeta[1] * b1[0] - nod_pt->x(0),
                           2) +
                       pow(lower_left_node_pt->x(1) + zeta[0] * b0[1] +
                             zeta[1] * b1[1] - nod_pt->x(1),
                           2) +
                       pow(lower_left_node_pt->x(2) + zeta[0] * b0[2] +
                             zeta[1] * b1[2] - nod_pt->x(2),
                           2);

        if (sqrt(error) > Tolerance_for_boundary_finding)
        {
          std::ostringstream error_message;
          error_message
            << "Error in setup of boundary coordinate " << sqrt(error)
            << std::endl
            << "exceeds tolerance specified by the static member data\n"
            << "TetMeshBase::Tolerance_for_boundary_finding = "
            << Tolerance_for_boundary_finding << std::endl
            << "This usually means that the boundary is not planar.\n\n"
            << "You can suppress this error message by recompiling \n"
            << "recompiling without PARANOID or by changing the tolerance.\n";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Set boundary coordinate
        nod_pt->set_coordinates_on_boundary(b, zeta);
      }
    }

    // Indicate that boundary coordinate has been set up
    Boundary_coordinate_exists[b] = true;

    // Cleanup
    unsigned n = face_el_pt.size();
    for (unsigned e = 0; e < n; e++)
    {
      delete face_el_pt[e];
    }


    // Reset nodal position
    for (std::map<Node*, Vector<double>>::iterator it = backup_position.begin();
         it != backup_position.end();
         it++)
    {
      Node* nod_pt = (*it).first;
      Vector<double> pos((*it).second);
      for (unsigned i = 0; i < 3; i++)
      {
        nod_pt->x(i) = pos[i];
      }
    }
  }


} // namespace oomph


#endif
