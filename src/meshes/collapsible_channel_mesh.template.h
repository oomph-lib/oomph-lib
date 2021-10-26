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
// Include guards
#ifndef OOMPH_COLLAPSIBLE_CHANNEL_MESH_HEADER
#define OOMPH_COLLAPSIBLE_CHANNEL_MESH_HEADER

// Generic includes
#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/quad_mesh.h"

// Mesh is based on simple_rectangular_quadmesh
#include "simple_rectangular_quadmesh.template.h"
#include "simple_rectangular_quadmesh.template.cc"

// Include *.template.cc to allow building the templated member functions
#include "../generic/macro_element_node_update_element.h"

// We need to include the templated sources for algebraic meshes
// to allow the build of all templates.
#include "../generic/algebraic_elements.h"

// Include the headers file for collapsible channel
#include "collapsible_channel_domain.h"

namespace oomph
{
  //========================================================================
  /// Basic collapsible channel mesh.
  /// The mesh is derived from the \c SimpleRectangularQuadMesh
  /// so it's node and element numbering scheme is the same
  /// as in that mesh. Only the boundaries are numbered differently
  /// to allow the easy identification of the "collapsible" segment.
  /// Boundary coordinates are set up for all nodes
  /// located on boundary 3 (the collapsible segment).
  /// The curvilinear ("collapsible") segment is defined by
  /// a \c GeomObject.
  //========================================================================
  template<class ELEMENT>
  class CollapsibleChannelMesh : public SimpleRectangularQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements in upstream/collapsible/
    /// downstream segment and across the channel; lengths of upstream/
    /// collapsible/downstream segments and width of channel, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    CollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// destructor
    ~CollapsibleChannelMesh()
    {
      delete Domain_pt;
    }

    /// Access function to GeomObject representing wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }

    /// Access function to domain
    CollapsibleChannelDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Function pointer for function that squashes
    /// the mesh near the walls. Default trivial mapping (the identity)
    /// leaves vertical nodal positions unchanged. Mapping is
    /// used in underlying CollapsibleChannelDomain. Virtual
    /// so we can break it in derived classes (e.g. the Algebraic
    /// versions of this mesh where it doesn't make any sense
    /// to provide the bl_squash_fct after the mesh has been built).
    virtual CollapsibleChannelDomain::BLSquashFctPt& bl_squash_fct_pt()
    {
      return Domain_pt->bl_squash_fct_pt();
    }


    /// Function pointer for function that squashes
    /// the mesh near the walls. Default trivial mapping (the identity)
    /// leaves vertical nodal positions unchanged. Mapping is
    /// used in underlying CollapsibleChannelDomain. Const version.
    CollapsibleChannelDomain::BLSquashFctPt bl_squash_fct_pt() const
    {
      return Domain_pt->bl_squash_fct_pt();
    }


    /// Function pointer for function that redistributes the
    /// elements in the axial direction. Virtual
    /// so we can break it in derived classes (e.g. the Algebraic
    /// versions of this mesh where it doesn't make any sense
    /// to provide the bl_squash_fct after the mesh has been built).
    virtual CollapsibleChannelDomain::AxialSpacingFctPt& axial_spacing_fct_pt()
    {
      return Domain_pt->axial_spacing_fct_pt();
    }


    /// Function pointer for function that redistributes the
    /// elements in the axial direction. Const version
    virtual CollapsibleChannelDomain::AxialSpacingFctPt& axial_spacing_fct_pt()
      const
    {
      return Domain_pt->axial_spacing_fct_pt();
    }


  protected:
    /// Pointer to domain
    CollapsibleChannelDomain* Domain_pt;

    /// Number of element columns in upstream part
    unsigned Nup;

    /// Number of element columns in collapsible part
    unsigned Ncollapsible;

    /// Number of element columns in downstream part
    unsigned Ndown;

    /// Number of element rows across channel
    unsigned Ny;

    /// Pointer to geometric object that represents the moving wall
    GeomObject* Wall_pt;
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //===================================================================
  /// Refineable collapsible channel mesh.
  /// The mesh is derived from the \c SimpleRectangularQuadMesh
  /// so it's node and element numbering scheme is the same
  /// as in that mesh. Only the boundaries are numbered differently
  /// to allow the easy identification of the "collapsible" segment.
  ///  Boundary coordinates are set up for all nodes
  /// located on boundary 3 (the collapsible segment).
  /// The curvilinear ("collapsible") segment is defined by
  /// a \c GeomObject.
  //====================================================================
  template<class ELEMENT>
  class RefineableCollapsibleChannelMesh
    : public virtual CollapsibleChannelMesh<ELEMENT>,
      public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements, lengths, pointer to
    /// geometric object that describes the wall and timestepper
    RefineableCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }


    /// Destructor(empty)
    ~RefineableCollapsibleChannelMesh() {}
  };


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //=====start_of_mesh=======================================================
  /// Collapsible channel mesh with MacroElement-based node update.
  /// The collapsible segment is represented by the specified geometric object.
  /// Some or all of the geometric Data in that geometric object
  /// may contain unknowns in the global Problem. The dependency
  /// on these unknowns is taken into account when setting up
  /// the Jacobian matrix of the elements. For this purpose,
  /// the element (whose type is specified by the template parameter)
  /// must inherit from MacroElementNodeUpdateElementBase.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateCollapsibleChannelMesh
    : public virtual MacroElementNodeUpdateMesh,
      public virtual CollapsibleChannelMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass numbers of elements and dimensions of the
    /// various parts of the collapsible channel, pointer to
    /// geometric object that represents the wall and pointer to
    /// timestepper (defaults to Steady).
    MacroElementNodeUpdateCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;
        error_message << "Base class for ELEMENT in "
                      << "MacroElementNodeUpdateCollapsibleChannelMesh needs"
                      << "to be of type MacroElementNodeUpdateElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name =
          "MacroElementNodeUpdateCollapsibleChannelMesh::\n";
        function_name += "MacroElementNodeUpdateCollapsibleChannelMesh()";

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
            "MacroElementNodeUpdateCollapsibleChannelMesh::\n";
          function_name += "MacroElementNodeUpdateCollapsibleChannelMesh()";

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // There's just one GeomObject
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = this->Wall_pt;

        // Tell the element which geom objects its macro-element-based
        // node update depends on
        el_pt->set_node_update_info(geom_object_pt);
      }

      // Add the geometric object(s) for the wall to the mesh's storage
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = this->Wall_pt;
      MacroElementNodeUpdateMesh::set_geom_object_vector_pt(geom_object_pt);

      // Fill in the domain pointer to the mesh's storage in the base class
      MacroElementNodeUpdateMesh::macro_domain_pt() = this->domain_pt();

    } // end of constructor


    /// Destructor: empty
    virtual ~MacroElementNodeUpdateCollapsibleChannelMesh() {}


  }; // end of mesh


  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////


  //=====start_of_mesh=======================================================
  /// Refineable collapsible channel mesh with MacroElement-based node update.
  /// The collapsible segment is represented by the specified geometric object.
  /// Some or all of the geometric Data in that geometric object
  /// may contain unknowns in the global Problem. The dependency
  /// on these unknowns is taken into account when setting up
  /// the Jacobian matrix of the elements. For this purpose,
  /// the element (whose type is specified by the template parameter)
  /// must inherit from MacroElementNodeUpdateElementBase.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateRefineableCollapsibleChannelMesh
    : public virtual MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>,
      public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass numbers of elements and dimensions of the
    /// various parts of the collapsible channel, pointer to
    /// geometric object that represents the wall and pointer to
    /// timestepper (defaults to Steady).
    MacroElementNodeUpdateRefineableCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt),
        MacroElementNodeUpdateCollapsibleChannelMesh<ELEMENT>(nup,
                                                              ncollapsible,
                                                              ndown,
                                                              ny,
                                                              lup,
                                                              lcollapsible,
                                                              ldown,
                                                              ly,
                                                              wall_pt,
                                                              time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }

    /// Destructor: empty
    virtual ~MacroElementNodeUpdateRefineableCollapsibleChannelMesh() {}

  }; // end of mesh


  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////


  //========start_of_algebraic_collapsible_channel_mesh==============
  /// Collapsible channel mesh with algebraic node update
  //=================================================================
  template<class ELEMENT>
  class AlgebraicCollapsibleChannelMesh
    : public virtual CollapsibleChannelMesh<ELEMENT>,
      public AlgebraicMesh
  {
  public:
    /// Constructor: Pass number of elements in upstream/collapsible/
    /// downstream segment and across the channel; lengths of upstream/
    /// collapsible/downstream segments and width of channel, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    AlgebraicCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt)
    {
      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(wall_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Destructor: empty
    virtual ~AlgebraicCollapsibleChannelMesh() {}


    /// Constructor: Pass number of elements in upstream/collapsible/
    /// downstream segment and across the channel; lengths of upstream/
    /// collapsible/downstream segments and width of channel, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    AlgebraicCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      CollapsibleChannelDomain::BLSquashFctPt bl_squash_function_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt)
    {
      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(wall_pt);

      // Set boundary layer squash function
      this->Domain_pt->bl_squash_fct_pt() = bl_squash_function_pt;

      // Do MacroElement-based node update
      CollapsibleChannelMesh<ELEMENT>::node_update();

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Function pointer for function that squashes
    /// the mesh near the walls. Default trivial mapping (the identity)
    /// leaves vertical nodal positions unchanged. Mapping is
    /// used in underlying CollapsibleChannelDomain. Broken function
    /// that overloads the version in the CollapsibleChannelMesh.
    /// It does not make sense to specify the function pointer
    /// after the mesh has been set up!
    CollapsibleChannelDomain::BLSquashFctPt& bl_squash_fct_pt()
    {
      std::ostringstream error_message;
      error_message
        << "It does not make sense to set the bl_squash_fct_pt \n"
        << "outside the constructor as it's only used to set up the \n"
        << "algebraic remesh data when the algebraic mesh is first built. \n";
      std::string function_name =
        "AlgebraicCollapsibleChannelMesh::bl_squash_fct_pt()\n";

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      // Dummy return
      return Dummy_fct_pt;
    }


    /// Function pointer for function that redistributes nodes
    /// axially. Default trivial mapping (the identity)
    /// leaves vertical nodal positions unchanged. Mapping is
    /// used in underlying CollapsibleChannelDomain. Broken function
    /// that overloads the version in the CollapsibleChannelMesh.
    /// It does not make sense to specify the function pointer
    /// after the mesh has been set up!
    CollapsibleChannelDomain::BLSquashFctPt& axial_spacing_fct_pt()
    {
      std::ostringstream error_message;
      error_message
        << "It does not make sense to set the axial_spacing_fct_pt \n"
        << "outside the constructor as it's only used to set up the \n"
        << "algebraic remesh data when the algebraic mesh is first built. \n";
      std::string function_name =
        "AlgebraicCollapsibleChannelMesh::axial_spacing_fct_pt()\n";

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      // Dummy return
      return Dummy_fct_pt;
    }


    /// Update nodal position at time level t (t=0: present;
    /// t>0: previous)
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt);

    /// Update the node-udate data after mesh adaptation.
    /// Empty -- no update of node update required as this is
    /// non-refineable mesh.
    void update_node_update(AlgebraicNode*& node_pt) {}

  protected:
    /// Function to setup the algebraic node update
    void setup_algebraic_node_update();

    /// Dummy function pointer
    CollapsibleChannelDomain::BLSquashFctPt Dummy_fct_pt;
  };


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////


  //=====start_of_refineable_algebraic_collapsible_channel_mesh=====
  /// Refineable version of the CollapsibleChannel mesh with
  /// algebraic node update.
  //=================================================================
  template<class ELEMENT>
  class RefineableAlgebraicCollapsibleChannelMesh
    : public RefineableQuadMesh<ELEMENT>,
      public virtual AlgebraicCollapsibleChannelMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass number of elements in upstream/collapsible/
    /// downstream segment and across the channel; lengths of upstream/
    /// collapsible/downstream segments and width of channel, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    RefineableAlgebraicCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt),
        AlgebraicCollapsibleChannelMesh<ELEMENT>(nup,
                                                 ncollapsible,
                                                 ndown,
                                                 ny,
                                                 lup,
                                                 lcollapsible,
                                                 ldown,
                                                 ly,
                                                 wall_pt,
                                                 time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }


    /// Constructor: Pass number of elements in upstream/collapsible/
    /// downstream segment and across the channel; lengths of upstream/
    /// collapsible/downstream segments and width of channel, pointer to
    /// GeomObject that defines the collapsible segment and pointer to
    /// TimeStepper (defaults to the default timestepper, Steady).
    RefineableAlgebraicCollapsibleChannelMesh(
      const unsigned& nup,
      const unsigned& ncollapsible,
      const unsigned& ndown,
      const unsigned& ny,
      const double& lup,
      const double& lcollapsible,
      const double& ldown,
      const double& ly,
      GeomObject* wall_pt,
      CollapsibleChannelDomain::BLSquashFctPt bl_squash_function_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : CollapsibleChannelMesh<ELEMENT>(nup,
                                        ncollapsible,
                                        ndown,
                                        ny,
                                        lup,
                                        lcollapsible,
                                        ldown,
                                        ly,
                                        wall_pt,
                                        time_stepper_pt),
        AlgebraicCollapsibleChannelMesh<ELEMENT>(nup,
                                                 ncollapsible,
                                                 ndown,
                                                 ny,
                                                 lup,
                                                 lcollapsible,
                                                 ldown,
                                                 ly,
                                                 wall_pt,
                                                 bl_squash_function_pt,
                                                 time_stepper_pt)
    {
      // Build quadtree forest
      this->setup_quadtree_forest();
    }

    /// Update the node update data for specified node following
    /// any mesh adapation
    void update_node_update(AlgebraicNode*& node_pt);
  };


} // namespace oomph

#endif
