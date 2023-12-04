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
#ifndef OOMPH_CHANNEL_WITH_LEAFLET_MESH_HEADER
#define OOMPH_CHANNEL_WITH_LEAFLET_MESH_HEADER

// Generic includes
#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/quad_mesh.h"

// Mesh is based on simple_rectangular_quadmesh
#include "simple_rectangular_quadmesh.template.h"
#include "simple_rectangular_quadmesh.template.cc"

// Include macro elements
#include "../generic/macro_element_node_update_element.h"

// and algebraic elements
#include "../generic/algebraic_elements.h"

// Include the headers file for domain
#include "channel_with_leaflet_domain.h"

namespace oomph
{
  //===================================================================
  /// Channel with leaflet mesh
  //===================================================================
  template<class ELEMENT>
  class ChannelWithLeafletMesh : public SimpleRectangularQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    ChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor : empty
    virtual ~ChannelWithLeafletMesh() {}

    /// Access function to domain
    ChannelWithLeafletDomain* domain_pt()
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to domain
    ChannelWithLeafletDomain* Domain_pt;

    /// Pointer to GeomObject that represents the leaflet
    GeomObject* Leaflet_pt;
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable version of ChannelWithLeafletMesh
  //===================================================================
  template<class ELEMENT>
  class RefineableChannelWithLeafletMesh
    : public virtual ChannelWithLeafletMesh<ELEMENT>,
      public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    RefineableChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }

    /// Destructor (empty)
    virtual ~RefineableChannelWithLeafletMesh() {}
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====start_of_mesh=======================================================
  /// Channel with leaflet mesh with MacroElement-based node update.
  /// The leaflet is represented by the specified geometric object.
  /// Some or all of the geometric Data in that geometric object
  /// may contain unknowns in the global Problem. The dependency
  /// on these unknowns is taken into account when setting up
  /// the Jacobian matrix of the elements. For this purpose,
  /// the element (whose type is specified by the template parameter)
  /// must inherit from MacroElementNodeUpdateElementBase.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateChannelWithLeafletMesh
    : public virtual MacroElementNodeUpdateMesh,
      public virtual ChannelWithLeafletMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    MacroElementNodeUpdateChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;
        error_message << "Base class for ELEMENT in "
                      << "MacroElementNodeUpdateChannelWithLeafletMesh needs"
                      << "to be of type MacroElementNodeUpdateElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name =
          "MacroElementNodeUpdateChannelWithLeafletMesh::\n";
        function_name += "MacroElementNodeUpdateChannelWithLeafletMesh()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      delete el_pt;
#endif

      // Setup all the information that's required for MacroElement-based
      // node update: Tell the elements that their geometry depends on the
      // wall geometric object
      unsigned n_element = this->nelement();
      for (unsigned i = 0; i < n_element; i++)
      {
        // Upcast from FiniteElement to the present element
        ELEMENT* el_pt = dynamic_cast<ELEMENT*>(this->element_pt(i));

#ifdef PARANOID
        // Check if cast is successful
        MacroElementNodeUpdateElementBase* m_el_pt =
          dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt);
        if (m_el_pt == 0)
        {
          std::ostringstream error_message;
          error_message
            << "Failed to upcast to MacroElementNodeUpdateElementBase\n";
          error_message << "Element must be derived from "
                           "MacroElementNodeUpdateElementBase\n";
          error_message << "but it is of type " << typeid(el_pt).name();

          std::string function_name =
            "MacroElementNodeUpdateChannelWithLeafletMesh::\n";
          function_name += "MacroElementNodeUpdateChannelWithLeafletMesh()";

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // There's just one GeomObject
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = this->Leaflet_pt;

        // Tell the element which geom objects its macro-element-based
        // node update depends on
        el_pt->set_node_update_info(geom_object_pt);
      }

      // Add the geometric object(s) for the wall to the mesh's storage
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = this->Leaflet_pt;
      MacroElementNodeUpdateMesh::set_geom_object_vector_pt(geom_object_pt);

      // Fill in the domain pointer to the mesh's storage in the base class
      MacroElementNodeUpdateMesh::macro_domain_pt() = this->domain_pt();

    } // end of constructor


    /// Destructor: empty
    virtual ~MacroElementNodeUpdateChannelWithLeafletMesh() {}


  }; // end of mesh


  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////


  //=====start_of_mesh=======================================================
  /// Refineable  mesh with MacroElement-based node update.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateRefineableChannelWithLeafletMesh
    : public virtual MacroElementNodeUpdateChannelWithLeafletMesh<ELEMENT>,
      public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    MacroElementNodeUpdateRefineableChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt),
        MacroElementNodeUpdateChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                                              lleft,
                                                              lright,
                                                              hleaflet,
                                                              htot,
                                                              nleft,
                                                              nright,
                                                              ny1,
                                                              ny2,
                                                              time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }


    /// Destructor: empty
    virtual ~MacroElementNodeUpdateRefineableChannelWithLeafletMesh() {}

  }; // end of mesh


  /// ////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Algebraic version of ChannelWithLeafletMesh. Leaflet is
  /// assumed to be in its undeformed (straight vertical) position
  /// when the algebraic node update is set up.
  //=================================================================
  template<class ELEMENT>
  class AlgebraicChannelWithLeafletMesh
    : public AlgebraicMesh,
      public virtual ChannelWithLeafletMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    AlgebraicChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt)
    {
      // Store origin of leaflet for fast reference
      Vector<double> zeta(1);
      zeta[0] = 0.0;
      Vector<double> r(2);
      this->Leaflet_pt->position(zeta, r);
      X_0 = r[0];

      // Store length of the leaflet for fast access (it's also available
      // through the domain, of course)
      Hleaflet = hleaflet;

      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(leaflet_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }


    /// Destructor: empty
    virtual ~AlgebraicChannelWithLeafletMesh() {}


    /// Update the geometric references that are used
    /// to update node after mesh adaptation.
    /// Empty -- no update of node update required without adaptivity
    void update_node_update(AlgebraicNode*& node_pt) {}

    /// Update nodal position at time level t (t=0: present;
    /// t>0: previous)
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt);

  protected:
    /// Function to setup the algebraic node update
    void setup_algebraic_node_update();

    /// Update function for nodes in lower left region (I)
    void node_update_I(const unsigned& t, AlgebraicNode*& node_pt);

    /// Update function for nodes in lower right region (II)
    void node_update_II(const unsigned& t, AlgebraicNode*& node_pt);

    /// Update function for nodes in upper left region (III)
    void node_update_III(const unsigned& t, AlgebraicNode*& node_pt);

    /// Update function for nodes in upper right region (IV)
    void node_update_IV(const unsigned& t, AlgebraicNode*& node_pt);

    /// Helper function
    void slanted_bound_up(const unsigned& t,
                          const Vector<double>& zeta,
                          Vector<double>& r);

    /// Origin of the wall (stored explicitly for reference in
    /// algebraic node update -- it's also stored independently in
    /// domain....)
    double X_0;

    /// Length of the leaflet (stored explicitly for reference in
    /// algebraic node update -- it's also stored independently in
    /// domain....)
    double Hleaflet;
  };

  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable version of algebraic ChannelWithLeafletMesh
  //===================================================================
  template<class ELEMENT>
  class RefineableAlgebraicChannelWithLeafletMesh
    : public RefineableQuadMesh<ELEMENT>,
      public virtual AlgebraicChannelWithLeafletMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet,
    /// the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    RefineableAlgebraicChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt),
        AlgebraicChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                                 lleft,
                                                 lright,
                                                 hleaflet,
                                                 htot,
                                                 nleft,
                                                 nright,
                                                 ny1,
                                                 ny2,
                                                 time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }

    /// Update the node update data for specified node following
    /// any mesh adapation
    void update_node_update(AlgebraicNode*& node_pt);
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //==========================================================================
  /// Channel with leaflet mesh upgraded to (pseudo-)solid mesh
  //==========================================================================
  template<class ELEMENT>
  class PseudoElasticChannelWithLeafletMesh
    : public virtual ChannelWithLeafletMesh<ELEMENT>,
      public virtual SolidMesh
  {
  public:
    /// Constructor: Pass pointer to GeomObject that represents the
    /// leaflet, the length of the domain to left and right of the leaflet, the
    /// height of the leaflet and the overall height of the channel,
    /// the number of element columns to the left and right of the leaflet,
    /// the number of rows of elements from the bottom of the channel to
    /// the end of the leaflet, the number of rows of elements above the
    /// end of the leaflet. Timestepper defaults to Steady default
    /// Timestepper defined in the Mesh base class
    PseudoElasticChannelWithLeafletMesh(
      GeomObject* leaflet_pt,
      const double& lleft,
      const double& lright,
      const double& hleaflet,
      const double& htot,
      const unsigned& nleft,
      const unsigned& nright,
      const unsigned& ny1,
      const unsigned& ny2,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : ChannelWithLeafletMesh<ELEMENT>(leaflet_pt,
                                        lleft,
                                        lright,
                                        hleaflet,
                                        htot,
                                        nleft,
                                        nright,
                                        ny1,
                                        ny2,
                                        time_stepper_pt)
    {
      // Update position of all nodes (the ones haven't been given
      // positions yet!)
      bool update_all_solid_nodes = true;
      node_update(update_all_solid_nodes);

      // Make the current configuration the undeformed one by
      // setting the nodal Lagrangian coordinates to their current
      // Eulerian ones
      set_lagrangian_nodal_coordinates();
    }

    /// Destructor : empty
    virtual ~PseudoElasticChannelWithLeafletMesh() {}
  };


} // namespace oomph

#endif
