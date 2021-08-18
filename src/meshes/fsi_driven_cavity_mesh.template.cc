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
#ifndef OOMPH_FSI_DRIVEN_CAVITY_MESH_TEMPLATE_CC
#define OOMPH_FSI_DRIVEN_CAVITY_MESH_TEMPLATE_CC

// Include the headers file for collapsible channel
#include "fsi_driven_cavity_mesh.template.h"


namespace oomph
{
  //========================================================================
  /// Constructor: Pass number of elements, lengths, pointer to GeomObject
  /// that defines the collapsible segment and pointer to TimeStepper
  /// (defaults to the default timestepper, Steady).
  //========================================================================
  template<class ELEMENT>
  FSIDrivenCavityMesh<ELEMENT>::FSIDrivenCavityMesh(
    const unsigned& nx,
    const unsigned& ny,
    const double& lx,
    const double& ly,
    const double& gap_fraction,
    GeomObject* wall_pt,
    TimeStepper* time_stepper_pt)
    : SimpleRectangularQuadMesh<ELEMENT>(nx, ny, lx, ly, time_stepper_pt),
      Nx(nx),
      Ny(ny),
      Gap_fraction(gap_fraction),
      Wall_pt(wall_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Update the boundary numbering scheme and set boundary coordinate
    //-----------------------------------------------------------------

    // (Note: The original SimpleRectangularQuadMesh had four boundaries.
    // We need to overwrite the boundary lookup scheme for the current
    // mesh so that the collapsible segment becomes identifiable).
    // While we're doing this, we're also setting up a boundary
    // coordinate for the nodes located on the collapsible segment.
    // The boundary coordinate can be used to setup FSI.

    // How many boundaries does the mesh have now?
    unsigned nbound = this->nboundary();
    for (unsigned b = 0; b < nbound; b++)
    {
      // Remove all nodes on this boundary from the mesh's lookup scheme
      // and also delete the reverse lookup scheme held by the nodes
      this->remove_boundary_nodes(b);
    }

#ifdef PARANOID
    // Sanity check
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      if (this->node_pt(j)->is_on_boundary())
      {
        std::ostringstream error_message;
        error_message << "Node " << j << "is still on boundary " << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif

    // Change the numbers of boundaries
    this->set_nboundary(6);

    // Get the number of nodes along the element edge from first element
    unsigned nnode_1d = this->finite_element_pt(0)->nnode_1d();

    // Vector of Lagrangian coordinates used as boundary coordinate
    Vector<double> zeta(1);

    // Zeta increment over elements (used for assignment of
    // boundary coordinate)
    double dzeta = lx / double(nx);

    // Manually loop over the elements near the boundaries and
    // assign nodes to boundaries. Also set up boundary coordinate
    unsigned nelem = this->nelement();
    for (unsigned e = 0; e < nelem; e++)
    {
      // Bottom row of elements
      if (e < nx)
      {
        for (unsigned i = 0; i < nnode_1d; i++)
        {
          this->add_boundary_node(0, this->finite_element_pt(e)->node_pt(i));
        }
      }
      // Collapsible bit
      if (e > ((ny - 1) * nx) - 1)
      {
        for (unsigned i = 0; i < nnode_1d; i++)
        {
          this->add_boundary_node(
            3, this->finite_element_pt(e)->node_pt(2 * nnode_1d + i));

          // What column of elements are we in?
          unsigned ix = e - (ny - 1) * nx;

          // Zeta coordinate
          zeta[0] =
            double(ix) * dzeta + double(i) * dzeta / double(nnode_1d - 1);

          // Set boundary coordinate
          this->finite_element_pt(e)
            ->node_pt(2 * nnode_1d + i)
            ->set_coordinates_on_boundary(3, zeta);
        }
      }
      // Left end
      if (e % (nx) == 0)
      {
        for (unsigned i = 0; i < nnode_1d; i++)
        {
          Node* nod_pt = this->finite_element_pt(e)->node_pt(i * nnode_1d);

          // Rigid bit?
          if (nod_pt->x(1) >= ly * Gap_fraction)
          {
            this->add_boundary_node(4, nod_pt);
          }
          // Free bit
          else
          {
            this->add_boundary_node(5, nod_pt);
          }
        }
      }
      // Right end
      if (e % nx == nx - 1)
      {
        for (unsigned i = 0; i < nnode_1d; i++)
        {
          Node* nod_pt =
            this->finite_element_pt(e)->node_pt((i + 1) * nnode_1d - 1);

          // Rigid bit?
          if (nod_pt->x(1) >= ly * Gap_fraction)
          {
            this->add_boundary_node(2, nod_pt);
          }
          // Free bit
          else
          {
            this->add_boundary_node(1, nod_pt);
          }
        }
      }
    }

    // Re-setup lookup scheme that establishes which elements are located
    // on the mesh boundaries (doesn't need to be wiped)
    this->setup_boundary_element_info();

    // We have only bothered to parametrise boundary 3
    this->Boundary_coordinate_exists[3] = true;
  }


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Perform algebraic mesh update at time level t (t=0: present;
  /// t>0: previous)
  //=================================================================
  template<class ELEMENT>
  void AlgebraicFSIDrivenCavityMesh<ELEMENT>::algebraic_node_update(
    const unsigned& t, AlgebraicNode*& node_pt)
  {
#ifdef PARANOID
    // We're updating the nodal positions (!) at time level t
    // and determine them by evaluating the wall GeomObject's
    // position at that gime level. I believe this only makes sense
    // if the t-th history value in the positional timestepper
    // actually represents previous values (rather than some
    // generalised quantity). Hence if this function is called with
    // t>nprev_values(), we issue a warning and terminate the execution.
    // It *might* be possible that the code still works correctly
    // even if this condition is violated (e.g. if the GeomObject's
    // position() function returns the appropriate "generalised"
    // position value that is required by the timestepping scheme but it's
    // probably worth flagging this up and forcing the user to manually switch
    // off this warning if he/she is 100% sure that this is kosher.
    if (t > node_pt->position_time_stepper_pt()->nprev_values())
    {
      std::string error_message =
        "Trying to update the nodal position at a time level";
      error_message += "beyond the number of previous values in the nodes'";
      error_message += "position timestepper. This seems highly suspect!";
      error_message += "If you're sure the code behaves correctly";
      error_message += "in your application, remove this warning ";
      error_message += "or recompile with PARNOID switched off.";

      std::string function_name = "AlgebraicFSIDrivenCavityMesh::";
      function_name += "algebraic_node_update()";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Extract references for update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // First reference value: Original x-position. Used as the start point
    // for the lines connecting the nodes in the vertical direction
    double x_bottom = ref_value[0];

    // Second reference value: Fractional position along
    // straight line from the bottom (at the original x position)
    // to the reference point on the upper wall
    double fract = ref_value[1];

    // Third reference value:  Reference local coordinate
    // in GeomObject that represents the upper wall (local coordinate
    // in finite element if the wall GeomObject is a finite element mesh)
    Vector<double> s(1);
    s[0] = ref_value[2];

    // Fourth reference value: zeta coordinate on the upper wall
    // If the wall is a simple GeomObject, zeta[0]=s[0]
    // but if it's a compound GeomObject (e.g. a finite element mesh)
    // zeta scales during mesh refinement, whereas s[0] and the
    // pointer to the geom object have to be re-computed.
    // double zeta=ref_value[3]; // not needed here

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to actual wall geom object (either the same as the wall object
    // or the pointer to the actual finite element)
    GeomObject* geom_obj_pt = geom_object_pt[0];

    // Get position vector to wall at previous timestep t!
    Vector<double> r_wall(2);
    geom_obj_pt->position(t, s, r_wall);

    // Assign new nodal coordinate
    node_pt->x(t, 0) = x_bottom + fract * (r_wall[0] - x_bottom);
    node_pt->x(t, 1) = fract * r_wall[1];
  }


  //=====start_setup=================================================
  /// Setup algebraic mesh update -- assumes that mesh has
  /// initially been set up with a flush upper wall.
  //=================================================================
  template<class ELEMENT>
  void AlgebraicFSIDrivenCavityMesh<ELEMENT>::setup_algebraic_node_update()
  {
    // Loop over all nodes in mesh
    unsigned nnod = this->nnode();
    for (unsigned j = 0; j < nnod; j++)
    {
      // Get pointer to node -- recall that that Mesh::node_pt(...) has been
      // overloaded in the AlgebraicMesh class to return a pointer to
      // an AlgebraicNode.
      AlgebraicNode* nod_pt = node_pt(j);

      // Get coordinates
      double x = nod_pt->x(0);
      double y = nod_pt->x(1);

      // Get zeta coordinate on the undeformed wall
      Vector<double> zeta(1);
      zeta[0] = x;

      // Get pointer to geometric (sub-)object and Lagrangian coordinate
      // on that sub-object. For a wall that is represented by
      // a single geom object, this simply returns the input.
      // If the geom object consists of sub-objects (e.g.
      // if it is a finite element mesh representing a wall,
      // then we'll obtain the pointer to the finite element
      // (in its incarnation as a GeomObject) and the
      // local coordinate in that element.
      GeomObject* geom_obj_pt;
      Vector<double> s(1);
      this->Wall_pt->locate_zeta(zeta, geom_obj_pt, s);

      // Get position vector to wall:
      Vector<double> r_wall(2);
      geom_obj_pt->position(s, r_wall);

      // Sanity check: Confirm that the wall is in its undeformed position
#ifdef PARANOID
      if ((std::fabs(r_wall[0] - x) > 1.0e-15) &&
          (std::fabs(r_wall[1] - y) > 1.0e-15))
      {
        std::ostringstream error_stream;
        error_stream << "Wall must be in its undeformed position when\n"
                     << "algebraic node update information is set up!\n "
                     << "x-discrepancy: " << std::fabs(r_wall[0] - x)
                     << std::endl
                     << "y-discrepancy: " << std::fabs(r_wall[1] - y)
                     << std::endl;

        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // One geometric object is involved in update operation
      Vector<GeomObject*> geom_object_pt(1);

      // The actual geometric object (If the wall is simple GeomObject
      // this is the same as Wall_pt; if it's a compound GeomObject
      // this points to the sub-object)
      geom_object_pt[0] = geom_obj_pt;

      // The update function requires four  parameters:
      Vector<double> ref_value(4);

      // First reference value: Original x-position
      ref_value[0] = r_wall[0];

      // Second  reference value: fractional position along
      // straight line from the bottom (at the original x position)
      // to the point on the wall)
      ref_value[1] = y / r_wall[1];

      // Third reference value: Reference local coordinate
      // in wall element (local coordinate in FE if we're dealing
      // with a wall mesh)
      ref_value[2] = s[0];

      // Fourth reference value: zeta coordinate on wall
      // If the wall is a simple GeomObject, zeta[0]=s[0]
      // but if it's a compound GeomObject (e.g. a finite element mesh)
      // zeta scales during mesh refinement, whereas s[0] and the
      // pointer to the geom object have to be re-computed.
      ref_value[3] = zeta[0];

      // Setup algebraic update for node: Pass update information
      nod_pt->add_node_update_info(this, // mesh
                                   geom_object_pt, // vector of geom objects
                                   ref_value); // vector of  ref. values
    }


  } // end of setup_algebraic_node_update


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //========start_update_node_update=================================
  /// Update the geometric references that are used
  /// to update node after mesh adaptation.
  //=================================================================
  template<class ELEMENT>
  void RefineableAlgebraicFSIDrivenCavityMesh<ELEMENT>::update_node_update(
    AlgebraicNode*& node_pt)
  {
    // Extract reference values for node update by copy construction
    Vector<double> ref_value(node_pt->vector_ref_value());

    // First reference value: Original x-position
    // double x_bottom=ref_value[0]; // not needed here

    // Second reference value: fractional position along
    // straight line from the bottom (at the original x position)
    // to the point on the wall)
    // double fract=ref_value[1]; // not needed here

    // Third reference value:  Reference local coordinate
    // in GeomObject (local coordinate in finite element if the wall
    // GeomObject is a finite element mesh)
    // Vector<double> s(1);
    // s[0]=ref_value[2]; // This needs to be re-computed!

    // Fourth reference value: intrinsic coordinate on the (possibly
    // compound) wall.
    double zeta = ref_value[3];

    // Extract geometric objects for update by copy construction
    Vector<GeomObject*> geom_object_pt(node_pt->vector_geom_object_pt());

    // Pointer to actual wall geom object (either the same as wall object
    // or the pointer to the actual finite element)
    // GeomObject* geom_obj_pt=geom_object_pt[0]; // This needs to be
    // re-computed!

    // Get zeta coordinate on wall (as vector)
    Vector<double> zeta_wall(1);
    zeta_wall[0] = zeta;

    // Get pointer to geometric (sub-)object and Lagrangian coordinate
    // on that sub-object. For a wall that is represented by
    // a single geom object, this simply returns the input.
    // If the geom object consists of sub-objects (e.g.
    // if it is a finite element mesh representing a wall,
    // then we'll obtain the pointer to the finite element
    // (in its incarnation as a GeomObject) and the
    // local coordinate in that element.
    Vector<double> s(1);
    GeomObject* geom_obj_pt;
    this->Wall_pt->locate_zeta(zeta_wall, geom_obj_pt, s);

    // Update the pointer to the (sub-)GeomObject within which the
    // reference point is located. (If the wall is simple GeomObject
    // this is the same as Wall_pt; if it's a compound GeomObject
    // this points to the sub-object)
    geom_object_pt[0] = geom_obj_pt;

    // First reference value: Original x-position
    // ref_value[0]=r_wall[0]; // unchanged

    // Second reference value: fractional position along
    // straight line from the bottom (at the original x position)
    // to the point on the wall)
    // ref_value[1]=y/r_wall[1]; // unchanged

    // Update third reference value: Reference local coordinate
    // in wall element (local coordinate in FE if we're dealing
    // with a wall mesh)
    ref_value[2] = s[0];

    // Fourth reference value: zeta coordinate on wall
    // If the wall is a simple GeomObject, zeta[0]=s[0]
    // but if it's a compound GeomObject (e.g. a finite element mesh)
    // zeta scales during mesh refinement, whereas s[0] and the
    // pointer to the geom object have to be re-computed.
    // ref_value[3]=zeta[0];     //unchanged


    // Kill the existing node update info
    node_pt->kill_node_update_info();

    // Setup algebraic update for node: Pass update information
    node_pt->add_node_update_info(this, // mesh
                                  geom_object_pt, // vector of geom objects
                                  ref_value); // vector of ref. values
  }


} // namespace oomph
#endif
