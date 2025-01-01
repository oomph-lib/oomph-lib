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
#ifndef OOMPH_CHANNEL_SPINE_MESH_TEMPLATE_CC
#define OOMPH_CHANNEL_SPINE_MESH_TEMPLATE_CC

#include "channel_spine_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{
  //===========================================================================
  /// Constructor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction, axial length and height of layer,
  /// and pointer to timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains a layer of spinified fluid elements (of type ELEMENT;
  /// e.g  SpineElement<QCrouzeixRaviartElement<2>)
  /// and a surface layer of corresponding Spine interface elements
  /// of type INTERFACE_ELEMENT, e.g.
  /// SpineLineFluidInterfaceElement<ELEMENT> for 2D planar
  /// problems.
  //===========================================================================
  template<class ELEMENT>
  ChannelSpineMesh<ELEMENT>::ChannelSpineMesh(const unsigned& nx0,
                                              const unsigned& nx1,
                                              const unsigned& nx2,
                                              const unsigned& ny,
                                              const double& lx0,
                                              const double& lx1,
                                              const double& lx2,
                                              const double& h,
                                              GeomObject* wall_pt,
                                              TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(nx0 + nx1 + nx2,
                                   ny,
                                   0.0,
                                   lx0 + lx1 + lx2,
                                   0.0,
                                   h,
                                   false,
                                   false,
                                   time_stepper_pt),
      Nx0(nx0),
      Nx1(nx1),
      Nx2(nx2),
      Lx0(lx0),
      Lx1(lx1),
      Lx2(lx2),
      Wall_pt(wall_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...

    // Build the straight line object
    Straight_wall_pt = new StraightLine(h);

    // Now build the mesh:
    build_channel_spine_mesh(time_stepper_pt);
  }


  //===========================================================================
  /// Constuctor for spine 2D mesh: Pass number of elements in x-direction,
  /// number of elements in y-direction, axial length and height of layer,
  /// a boolean flag to make the mesh periodic in the x-direction,
  /// and pointer to timestepper (defaults to Static timestepper).
  ///
  /// The mesh contains a layer of elements (of type ELEMENT)
  /// and a surface layer of corresponding Spine interface elements (of type
  /// SpineLineFluidInterfaceElement<ELEMENT>).
  //===========================================================================
  template<class ELEMENT>
  ChannelSpineMesh<ELEMENT>::ChannelSpineMesh(const unsigned& nx0,
                                              const unsigned& nx1,
                                              const unsigned& nx2,
                                              const unsigned& ny,
                                              const double& lx0,
                                              const double& lx1,
                                              const double& lx2,
                                              const double& h,
                                              GeomObject* wall_pt,
                                              const bool& periodic_in_x,
                                              TimeStepper* time_stepper_pt)
    : RectangularQuadMesh<ELEMENT>(nx0 + nx1 + nx2,
                                   ny,
                                   0.0,
                                   lx0 + lx1 + lx2,
                                   0.0,
                                   h,
                                   periodic_in_x,
                                   false,
                                   time_stepper_pt),
      Nx0(nx0),
      Nx1(nx1),
      Nx2(nx2),
      Lx0(lx0),
      Lx1(lx1),
      Lx2(lx2),
      Wall_pt(wall_pt)
  {
    // Mesh can only be built with 2D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(2);

    // Mesh can only be built with spine elements
    MeshChecker::assert_geometric_element<SpineFiniteElement, ELEMENT>(2);

    // We've called the "generic" constructor for the RectangularQuadMesh
    // which doesn't do much...

    // Build the straight line object
    Straight_wall_pt = new StraightLine(h);

    // Now build the mesh:
    build_channel_spine_mesh(time_stepper_pt);
  }


  //===========================================================================
  /// Helper function that actually builds the channel-spine mesh
  /// based on the parameters set in the various constructors
  //===========================================================================
  template<class ELEMENT>
  void ChannelSpineMesh<ELEMENT>::build_channel_spine_mesh(
    TimeStepper* time_stepper_pt)
  {
    // Build the underlying quad mesh:
    RectangularQuadMesh<ELEMENT>::build_mesh(time_stepper_pt);

    // Read out the number of elements in the x-direction and y-direction
    // and in each of the left, centre and right regions
    unsigned n_x = this->Nx;
    unsigned n_y = this->Ny;
    unsigned n_x0 = this->Nx0;
    unsigned n_x1 = this->Nx1;
    unsigned n_x2 = this->Nx2;

    // Set up the pointers to elements in the left region
    unsigned nleft = n_x0 * n_y;
    ;
    Left_element_pt.reserve(nleft);
    unsigned ncentre = n_x1 * n_y;
    ;
    Centre_element_pt.reserve(ncentre);
    unsigned nright = n_x2 * n_y;
    ;
    Right_element_pt.reserve(nright);
    for (unsigned irow = 0; irow < n_y; irow++)
    {
      for (unsigned e = 0; e < n_x0; e++)
      {
        Left_element_pt.push_back(this->finite_element_pt(irow * n_x + e));
      }
      for (unsigned e = 0; e < n_x1; e++)
      {
        Centre_element_pt.push_back(
          this->finite_element_pt(irow * n_x + (n_x0 + e)));
      }
      for (unsigned e = 0; e < n_x2; e++)
      {
        Right_element_pt.push_back(
          this->finite_element_pt(irow * n_x + (n_x0 + n_x1 + e)));
      }
    }

#ifdef PARANOID
    // Check that we have the correct number of elements
    if (nelement() != nleft + ncentre + nright)
    {
      throw OomphLibError("Incorrect number of element pointers!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Allocate memory for the spines and fractions along spines
    //---------------------------------------------------------

    // Read out number of linear points in the element
    unsigned n_p = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

    unsigned nspine;
    // Allocate store for the spines:
    if (this->Xperiodic)
    {
      nspine = (n_p - 1) * n_x;
      Spine_pt.reserve(nspine);
      // Number of spines in each region
      // NOTE that boundary spines are in both regions
      Nleft_spine = (n_p - 1) * n_x0 + 1;
      Ncentre_spine = (n_p - 1) * n_x1 + 1;
      Nright_spine = (n_p - 1) * n_x2;
    }
    else
    {
      nspine = (n_p - 1) * n_x + 1;
      Spine_pt.reserve(nspine);
      // Number of spines in each region
      // NOTE that boundary spines are in both regions
      Nleft_spine = (n_p - 1) * n_x0 + 1;
      Ncentre_spine = (n_p - 1) * n_x1 + 1;
      Nright_spine = (n_p - 1) * n_x2 + 1;
    }

    // end Allocating memory


    // set up the vectors of geometric data & objects for building spines
    Vector<double> r_wall(2), zeta(1), s_wall(1);
    GeomObject* geometric_object_pt = 0;

    // LEFT REGION
    // ===========

    // SPINES IN LEFT REGION
    // ---------------------

    // Set up zeta increments
    double zeta_lo = 0.0;
    double dzeta = Lx0 / n_x0;

    // Initialise number of elements in previous regions:
    unsigned n_prev_elements = 0;


    // FIRST SPINE
    //-----------

    // Element 0
    // Node 0
    // Assign the new spine with unit length
    Spine* new_spine_pt = new Spine(1.0);
    new_spine_pt->spine_height_pt()->pin(0);
    Spine_pt.push_back(new_spine_pt);

    // Get pointer to node
    SpineNode* nod_pt = element_node_pt(0, 0);
    // Set the pointer to the spine
    nod_pt->spine_pt() = new_spine_pt;
    // Set the fraction
    nod_pt->fraction() = 0.0;
    // Pointer to the mesh that implements the update fct
    nod_pt->spine_mesh_pt() = this;
    // Set update fct id
    nod_pt->node_update_fct_id() = 0;

    // Provide spine with additional storage for wall coordinate
    // and wall geom object:

    {
      // Get the Lagrangian coordinate in the Lower Wall
      zeta[0] = 0.0;
      // Get the geometric object and local coordinate
      Straight_wall_pt->locate_zeta(zeta, geometric_object_pt, s_wall);

      // The local coordinate is a geometric parameter
      // This needs to be set (rather than added) because the
      // same spine may be visited more than once
      Vector<double> parameters(1, s_wall[0]);
      nod_pt->spine_pt()->set_geom_parameter(parameters);

      // Get position of wall
      Straight_wall_pt->position(s_wall, r_wall);

      // Adjust spine height
      nod_pt->spine_pt()->height() = r_wall[1];

      // The sub geom object is one (and only) geom object
      // for spine:
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = geometric_object_pt;

      // Pass geom object(s) to spine
      nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
    }

    // Loop vertically along the spine
    // Loop over the elements
    for (unsigned long i = 0; i < n_y; i++)
    {
      // Loop over the vertical nodes, apart from the first
      for (unsigned l1 = 1; l1 < n_p; l1++)
      {
        // Get pointer to node
        SpineNode* nod_pt = element_node_pt(i * n_x, l1 * n_p);
        // Set the pointer to the spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() =
          (double(i) + double(l1) / double(n_p - 1)) / double(n_y);
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // Set update fct id
        nod_pt->node_update_fct_id() = 0;
      }
    } // end loop over elements


    // LOOP OVER OTHER SPINES
    //----------------------

    // Now loop over the elements horizontally
    for (unsigned long j = 0; j < n_x0; j++)
    {
      // Loop over the nodes in the elements horizontally, ignoring
      // the first column
      unsigned n_pmax = n_p;
      for (unsigned l2 = 1; l2 < n_pmax; l2++)
      {
        // Assign the new spine with unit height
        new_spine_pt = new Spine(1.0);
        new_spine_pt->spine_height_pt()->pin(0);
        Spine_pt.push_back(new_spine_pt);

        // Get the node
        SpineNode* nod_pt = element_node_pt(j, l2);
        // Set the pointer to spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // Set update fct id
        nod_pt->node_update_fct_id() = 0;

        {
          // Provide spine with additional storage for wall coordinate
          // and wall geom object:

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] =
            zeta_lo + double(j) * dzeta + double(l2) * dzeta / (n_p - 1.0);
          // Get the geometric object and local coordinate
          Straight_wall_pt->locate_zeta(zeta, geometric_object_pt, s_wall);

          // The local coordinate is a geometric parameter
          // This needs to be set (rather than added) because the
          // same spine may be visited more than once
          Vector<double> parameters(1, s_wall[0]);
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Get position of wall
          Straight_wall_pt->position(s_wall, r_wall);

          // Adjust spine height
          nod_pt->spine_pt()->height() = r_wall[1];

          // The sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = geometric_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
        }

        // Loop vertically along the spine
        // Loop over the elements
        for (unsigned long i = 0; i < n_y; i++)
        {
          // Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get the node
            SpineNode* nod_pt = element_node_pt(i * n_x + j, l1 * n_p + l2);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() =
              (double(i) + double(l1) / double(n_p - 1)) / double(n_y);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // Set update fct id
            nod_pt->node_update_fct_id() = 0;
          }
        }
      }
    }

    // Increment number of previous elements
    n_prev_elements += n_x0 * n_y;


    // CENTRE REGION
    // ===========

    zeta_lo = Lx0;
    dzeta = Lx1 / n_x1;

    // SPINES IN LEFT REGION
    // ---------------------

    // LOOP OVER OTHER SPINES
    //----------------------

    // Now loop over the elements horizontally
    for (unsigned long j = n_x0; j < n_x0 + n_x1; j++)
    {
      // Loop over the nodes in the elements horizontally, ignoring
      // the first column
      unsigned n_pmax = n_p;
      for (unsigned l2 = 1; l2 < n_pmax; l2++)
      {
        // Assign the new spine with unit height
        new_spine_pt = new Spine(1.0);
        new_spine_pt->spine_height_pt()->pin(0);
        Spine_pt.push_back(new_spine_pt);

        // Get the node
        SpineNode* nod_pt = element_node_pt(j, l2);
        // Set the pointer to spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // Set update fct id
        nod_pt->node_update_fct_id() = 1;

        {
          // Provide spine with additional storage for wall coordinate
          // and wall geom object:

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] = zeta_lo + double(j - n_x0) * dzeta +
                    double(l2) * dzeta / (n_p - 1.0);
          // Get the geometric object and local coordinate
          Wall_pt->locate_zeta(zeta, geometric_object_pt, s_wall);

          // The local coordinate is a geometric parameter
          // This needs to be set (rather than added) because the
          // same spine may be visited more than once
          Vector<double> parameters(1, s_wall[0]);
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Get position of wall
          Wall_pt->position(s_wall, r_wall);

          // Adjust spine height
          nod_pt->spine_pt()->height() = r_wall[1];

          // The sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = geometric_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
        }

        // Loop vertically along the spine
        // Loop over the elements
        for (unsigned long i = 0; i < n_y; i++)
        {
          // Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get the node
            SpineNode* nod_pt = element_node_pt(i * n_x + j, l1 * n_p + l2);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() =
              (double(i) + double(l1) / double(n_p - 1)) / double(n_y);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // Set update fct id
            nod_pt->node_update_fct_id() = 1;
          }
        }
      }
    }

    // Increment number of previous elements
    n_prev_elements += n_x1 * n_y;


    // RIGHT REGION
    // ===========

    // SPINES IN RIGHT REGION
    // ----------------------

    // Set up zeta increments
    zeta_lo = Lx0 + Lx1;
    dzeta = Lx2 / n_x2;

    // LOOP OVER OTHER SPINES
    //----------------------

    // Now loop over the elements horizontally
    for (unsigned long j = n_x0 + n_x1; j < n_x0 + n_x1 + n_x2; j++)
    {
      // Need to copy last spine in previous element over to first spine
      // in next elements, for all elements other than the first
      if (j > 0)
      {
        // For first spine, re-assign the geometric objects, since
        // we treat it as part of the right region.
        if (j == n_x0 + n_x1)
        {
          SpineNode* nod_pt = element_node_pt(j, 0);
          // Set update fct id
          nod_pt->node_update_fct_id() = 2;
          {
            // Provide spine with additional storage for wall coordinate
            // and wall geom object:

            // Get the Lagrangian coordinate in the Lower Wall
            zeta[0] = zeta_lo + double(j - n_x0 - n_x1) * dzeta;
            // Get the geometric object and local coordinate
            Straight_wall_pt->locate_zeta(zeta, geometric_object_pt, s_wall);

            // The local coordinate is a geometric parameter
            // This needs to be set (rather than added) because the
            // same spine may be visited more than once
            Vector<double> parameters(1, s_wall[0]);
            nod_pt->spine_pt()->set_geom_parameter(parameters);

            // Get position of wall
            Straight_wall_pt->position(s_wall, r_wall);

            // Adjust spine height
            nod_pt->spine_pt()->height() = r_wall[1];

            // The sub geom object is one (and only) geom object
            // for spine:
            Vector<GeomObject*> geom_object_pt(1);
            geom_object_pt[0] = geometric_object_pt;

            // Pass geom object(s) to spine
            nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
          }
        }
      }
      // Loop over the nodes in the elements horizontally, ignoring
      // the first column

      // Last spine needs special treatment in x-periodic meshes:
      unsigned n_pmax = n_p;
      if ((this->Xperiodic) && (j == n_x - 1)) n_pmax = n_p - 1;

      for (unsigned l2 = 1; l2 < n_pmax; l2++)
      {
        // Assign the new spine with unit height
        new_spine_pt = new Spine(1.0);
        new_spine_pt->spine_height_pt()->pin(0);
        Spine_pt.push_back(new_spine_pt);

        // Get the node
        SpineNode* nod_pt = element_node_pt(j, l2);
        // Set the pointer to spine
        nod_pt->spine_pt() = new_spine_pt;
        // Set the fraction
        nod_pt->fraction() = 0.0;
        // Pointer to the mesh that implements the update fct
        nod_pt->spine_mesh_pt() = this;
        // Set update fct id
        nod_pt->node_update_fct_id() = 2;

        {
          // Provide spine with additional storage for wall coordinate
          // and wall geom object:

          // Get the Lagrangian coordinate in the Lower Wall
          zeta[0] = zeta_lo + double(j - n_x0 - n_x1) * dzeta +
                    double(l2) * dzeta / (n_p - 1.0);
          // Get the geometric object and local coordinate
          Straight_wall_pt->locate_zeta(zeta, geometric_object_pt, s_wall);

          // The local coordinate is a geometric parameter
          // This needs to be set (rather than added) because the
          // same spine may be visited more than once
          Vector<double> parameters(1, s_wall[0]);
          nod_pt->spine_pt()->set_geom_parameter(parameters);

          // Get position of wall
          Straight_wall_pt->position(s_wall, r_wall);

          // Adjust spine height
          nod_pt->spine_pt()->height() = r_wall[1];

          // The sub geom object is one (and only) geom object
          // for spine:
          Vector<GeomObject*> geom_object_pt(1);
          geom_object_pt[0] = geometric_object_pt;

          // Pass geom object(s) to spine
          nod_pt->spine_pt()->set_geom_object_pt(geom_object_pt);
        }

        // Loop vertically along the spine
        // Loop over the elements
        for (unsigned long i = 0; i < n_y; i++)
        {
          // Loop over the vertical nodes, apart from the first
          for (unsigned l1 = 1; l1 < n_p; l1++)
          {
            // Get the node
            SpineNode* nod_pt = element_node_pt(i * n_x + j, l1 * n_p + l2);
            // Set the pointer to the spine
            nod_pt->spine_pt() = new_spine_pt;
            // Set the fraction
            nod_pt->fraction() =
              (double(i) + double(l1) / double(n_p - 1)) / double(n_y);
            // Pointer to the mesh that implements the update fct
            nod_pt->spine_mesh_pt() = this;
            // Set update fct id
            nod_pt->node_update_fct_id() = 2;
          }
        }
      }
    }

    // Increment number of previous elements
    n_prev_elements += n_x2 * n_y;

    // Last spine needs special treatment for periodic meshes
    // because it's the same as the first one...
    if (this->Xperiodic)
    {
      // Last spine is the same as first one...
      Spine* final_spine_pt = Spine_pt[0];

      // Get the node
      SpineNode* nod_pt = element_node_pt((n_x - 1), (n_p - 1));

      // Set the pointer for the first node
      nod_pt->spine_pt() = final_spine_pt;
      // Set the fraction to be the same as for the nodes on the first row
      nod_pt->fraction() = element_node_pt(0, 0)->fraction();
      // Pointer to the mesh that implements the update fct
      nod_pt->spine_mesh_pt() = element_node_pt(0, 0)->spine_mesh_pt();

      // Now loop vertically along the spine
      for (unsigned i = 0; i < n_y; i++)
      {
        // Loop over the vertical nodes, apart from the first
        for (unsigned l1 = 1; l1 < n_p; l1++)
        {
          // Get the node
          SpineNode* nod_pt =
            element_node_pt(i * n_x + (n_x - 1), l1 * n_p + (n_p - 1));

          // Set the pointer to the spine
          nod_pt->spine_pt() = final_spine_pt;
          // Set the fraction to be the same as in first row
          nod_pt->fraction() = element_node_pt(i * n_x, l1 * n_p)->fraction();
          // Pointer to the mesh that implements the update fct
          nod_pt->spine_mesh_pt() =
            element_node_pt(i * n_x, l1 * n_p)->spine_mesh_pt();
        }
      }
    }


  } // end of build_channel_spine_mesh


  //======================================================================
  /// Reorder the elements, so we loop over them vertically first
  /// (advantageous in "wide" domains if a frontal solver is used).
  //======================================================================
  template<class ELEMENT>
  void ChannelSpineMesh<ELEMENT>::element_reorder()
  {
    unsigned n_x = this->Nx;
    unsigned n_y = this->Ny;
    // Find out how many elements there are
    unsigned long Nelement = nelement();
    // Find out how many fluid elements there are
    unsigned long Nfluid = n_x * n_y;
    // Create a dummy array of elements
    Vector<FiniteElement*> dummy;

    // Loop over the elements in horizontal order
    for (unsigned long j = 0; j < n_x; j++)
    {
      // Loop over the elements in lower layer vertically
      for (unsigned long i = 0; i < n_y; i++)
      {
        // Push back onto the new stack
        dummy.push_back(finite_element_pt(n_x * i + j));
      }

      // Push back the line element onto the stack
      dummy.push_back(finite_element_pt(Nfluid + j));
    }

    // Now copy the reordered elements into the element_pt
    for (unsigned long e = 0; e < Nelement; e++)
    {
      Element_pt[e] = dummy[e];
    }

  } // end of element_reorder


} // namespace oomph

#endif
