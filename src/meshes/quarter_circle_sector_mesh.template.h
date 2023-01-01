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
#ifndef OOMPH_QUARTER_CIRCLE_SECTOR_MESH_HEADER
#define OOMPH_QUARTER_CIRCLE_SECTOR_MESH_HEADER

#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/quad_mesh.h"
#include "../generic/macro_element_node_update_element.h"

// Include the headers file for domain
#include "quarter_circle_sector_domain.h"


namespace oomph
{
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  class GeomObject;


  //====================================================================
  /// 2D quarter ring mesh class.
  /// The domain is specified by the GeomObject that identifies boundary 1.
  ///
  ///
  ///  \code
  ///                       ---___
  ///                      |      ---____
  ///                      |              -   BOUNDARY 1
  ///                      |               /
  ///                      |     [2]      /  |
  ///                      |             /     |
  ///                      | N          /        |
  ///                      | |_ E      /          |
  ///       BOUNDARY 2     |-----------           |
  ///                      |          |    [1]    |
  ///                      |   [0]    |           |  ^
  ///                      |          |           | / \  direction of
  ///                      | N        |    N      |  |   Lagrangian
  ///                      | |_ E     |    |_ E   |  |   coordinate
  ///                      |__________|___________|  |   along wall GeomObject
  ///
  ///                           BOUNDARY 0
  ///
  ///  Within the elements (MacroElements), the local coordinates
  ///  are such that the (E)astern direction coincides with the positive
  ///  s_0 direction,  while the (N)orther direction coincides with the positive
  ///  s_1 direction.
  ///
  /// \endcode
  ///
  /// Domain is parametrised by three macro elements as sketched.
  /// Nodal positions are determined via macro-element-based representation
  /// of the Domain (as a QuarterCircleSectorDomain).
  //====================================================================
  template<class ELEMENT>
  class QuarterCircleSectorMesh : public virtual QuadMeshBase
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh)
    QuarterCircleSectorMesh(
      GeomObject* wall_pt,
      const double& xi_lo,
      const double& fract_mid,
      const double& xi_hi,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor:
    virtual ~QuarterCircleSectorMesh() {}

    /// Access function to GeomObject representing wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }

    /// Access function to domain
    QuarterCircleSectorDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Function pointer for function that squashes
    /// the outer two macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value (defined in the underlying Domain object)
    QuarterCircleSectorDomain::BLSquashFctPt& bl_squash_fct_pt()
    {
      return Domain_pt->bl_squash_fct_pt();
    }


  protected:
    /// Pointer to Domain
    QuarterCircleSectorDomain* Domain_pt;

    /// Pointer to the geometric object that represents the curved wall
    /// (mesh boundary 1)
    GeomObject* Wall_pt;

    /// Lower limit for the (1D) coordinates along the wall
    double Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limit for the (1D) coordinates along the wall
    double Xi_hi;
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //====================================================================
  /// 2D quarter ring mesh class inherited from RefineableQuadMesh.
  /// The domain is specified by the GeomObject that identifies boundary 1.
  ///
  ///
  ///  \code
  ///                       ---___
  ///                      |      ---____
  ///                      |              -   BOUNDARY 1
  ///                      |               /
  ///                      |     [2]      /  |
  ///                      |             /     |
  ///                      | N          /        |
  ///                      | |_ E      /          |
  ///       BOUNDARY 2     |-----------           |
  ///                      |          |    [1]    |
  ///                      |   [0]    |           |  ^
  ///                      |          |           | / \  direction of
  ///                      | N        |    N      |  |   Lagrangian
  ///                      | |_ E     |    |_ E   |  |   coordinate
  ///                      |__________|___________|  |   along wall GeomObject
  ///
  ///                           BOUNDARY 0
  ///
  ///  Within the elements (MacroElements), the local coordinates
  ///  are such that the (E)astern direction coincides with the positive
  ///  s_0 direction,  while the (N)orther direction coincides with the positive
  ///  s_1 direction.
  ///
  /// \endcode
  ///
  /// Domain is parametrised by three macro elements as sketched.
  ///
  //====================================================================
  template<class ELEMENT>
  class RefineableQuarterCircleSectorMesh
    : public QuarterCircleSectorMesh<ELEMENT>,
      public virtual RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh).
    /// Adds refinement data to elements of QuarterCircleSectorMesh.
    RefineableQuarterCircleSectorMesh(
      GeomObject* wall_pt,
      const double& xi_lo,
      const double& fract_mid,
      const double& xi_hi,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuarterCircleSectorMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, time_stepper_pt)
    {
      // Basic mesh has been built -- just need to setup the
      // adaptivity information:

      // Setup quadtree forest
      this->setup_quadtree_forest();
    }

    /// Destructor: Empty
    virtual ~RefineableQuarterCircleSectorMesh() {}
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdate-version of RefineableQuarterCircleSectorMesh
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////

  class MacroElementNodeUpdateNode;

  //========================================================================
  /// MacroElementNodeUpdate version of RefineableQuarterCircleSectorMesh
  ///
  ///  \code
  ///                       ---___
  ///                      |      ---____
  ///                      |              -   BOUNDARY 1
  ///                      |               /
  ///                      |     [2]      /  |
  ///                      |             /     |
  ///                      | N          /        |
  ///                      | |_ E      /          |
  ///       BOUNDARY 2     |-----------           |
  ///                      |          |    [1]    |
  ///                      |   [0]    |           |  ^
  ///                      |          |           | / \  direction of
  ///                      | N        |    N      |  |   Lagrangian
  ///                      | |_ E     |    |_ E   |  |   coordinate
  ///                      |__________|___________|  |   along wall GeomObject
  ///
  ///                           BOUNDARY 0
  ///
  ///  Within the elements (MacroElements), the local coordinates
  ///  are such that the (E)astern direction coincides with the positive
  ///  s_0 direction,  while the (N)orther direction coincides with the positive
  ///  s_1 direction.
  ///
  /// \endcode
  ///
  /// Domain is parametrised by three macro elements as sketched. Elements
  /// need to be derived from MacroElementNodeUpdateElementBase.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateRefineableQuarterCircleSectorMesh
    : public virtual MacroElementNodeUpdateMesh,
      public virtual RefineableQuarterCircleSectorMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object, start and
    /// end coordinates on the geometric object and the fraction along
    /// which the dividing line is to be placed when updating the nodal
    /// positions, and timestepper (defaults to (Steady) default timestepper
    /// defined in Mesh). Setup the refineable mesh (by calling the
    /// constructor for the underlying  RefineableQuarterCircleSectorMesh)
    /// and the algebraic node update functions for nodes.
    MacroElementNodeUpdateRefineableQuarterCircleSectorMesh(
      GeomObject* wall_pt,
      const double& xi_lo,
      const double& fract_mid,
      const double& xi_hi,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : MacroElementNodeUpdateMesh(),
        RefineableQuarterCircleSectorMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;
        error_message
          << "Base class for ELEMENT in "
          << "MacroElementNodeUpdateRefineableQuarterCircleSectorMesh needs"
          << "to be of type MacroElementNodeUpdateElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name =
          "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
        function_name +=
          "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      delete el_pt;
#endif

      // Setup all the information that's required for MacroElement-based
      // node update: Tell the elements that their geometry depends on the
      // fishback geometric object
      this->setup_macro_element_node_update();
    }

    /// Destructor: empty
    virtual ~MacroElementNodeUpdateRefineableQuarterCircleSectorMesh() {}

    /// Resolve mesh update: Update current nodal
    /// positions via sparse MacroElement-based update.
    /// [Doesn't make sense to use this mesh with SolidElements anyway,
    /// so we buffer the case if update_all_solid_nodes is set to
    /// true.]
    void node_update(const bool& update_all_solid_nodes = false)
    {
#ifdef PARANOID
      if (update_all_solid_nodes)
      {
        std::string error_message =
          "Doesn't make sense to use an MacroElementNodeUpdateMesh with\n";
        error_message +=
          "SolidElements so specifying update_all_solid_nodes=true\n";
        error_message += "doesn't make sense either\n";

        std::string function_name =
          "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
        function_name += "node_update()";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      MacroElementNodeUpdateMesh::node_update();
    }

  private:
    /// Setup all the information that's required for MacroElement-based
    /// node update: Tell the elements that their geometry depends on the
    /// geometric object that parametrises the wall
    void setup_macro_element_node_update()
    {
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
            "MacroElementNodeUpdateRefineableQuaterCircleSectorMesh::\n";
          function_name += "setup_macro_element_node_update()";

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
    }
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic-mesh-version of RefineableQuarterCircleSectorMesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  class AlgebraicNode;

  //========================================================================
  /// Algebraic version of RefineableQuarterCircleSectorMesh
  ///
  ///  \code
  ///                       ---___
  ///                      |      ---____
  ///                      |              -   BOUNDARY 1
  ///                      |               /
  ///                      |     [2]      /  |
  ///                      |             /     |
  ///                      | N          /        |
  ///                      | |_ E      /          |
  ///       BOUNDARY 2     |-----------           |
  ///                      |          |    [1]    |
  ///                      |   [0]    |           |  ^
  ///                      |          |           | / \  direction of
  ///                      | N        |    N      |  |   Lagrangian
  ///                      | |_ E     |    |_ E   |  |   coordinate
  ///                      |__________|___________|  |   along wall GeomObject
  ///
  ///                           BOUNDARY 0
  ///
  ///  Within the elements (MacroElements), the local coordinates
  ///  are such that the (E)astern direction coincides with the positive
  ///  s_0 direction,  while the (N)orther direction coincides with the positive
  ///  s_1 direction.
  ///
  /// \endcode
  ///
  /// Domain is parametrised by three macro elements as sketched. Elements
  /// need to be derived from AlgebraicElementBase. In addition
  /// to all the refinement procedures available for
  /// RefineableQuarterCircleSectorMesh which forms the basis for
  /// this mesh, we implement algebraic node update functions for the nodes.
  //========================================================================
  template<class ELEMENT>
  class AlgebraicRefineableQuarterCircleSectorMesh
    : public virtual AlgebraicMesh,
      public RefineableQuarterCircleSectorMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object, start and
    /// end coordinates on the geometric object and the fraction along
    /// which the dividing line is to be placed when updating the nodal
    /// positions, and timestepper (defaults to (Steady) default timestepper
    /// defined in Mesh). Setup the refineable mesh (by calling the
    /// constructor for the underlying  RefineableQuarterCircleSectorMesh)
    /// and the algebraic update functions for nodes.
    AlgebraicRefineableQuarterCircleSectorMesh(
      GeomObject* wall_pt,
      const double& xi_lo,
      const double& fract_mid,
      const double& xi_hi,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : RefineableQuarterCircleSectorMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<AlgebraicElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;

        error_message << "Base class for ELEMENT in "
                      << "AlgebraicRefineableQuarterCircleSectorMesh needs"
                      << "to be of type AlgebraicElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name =
          " AlgebraicRefineableQuarterCircleSectorMesh::\n";
        function_name += "AlgebraicRefineableQuarterCircleSectorMesh()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      delete el_pt;
#endif

      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(wall_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Run self-test for algebraic mesh -- return 0/1 for OK/failure
    unsigned self_test()
    {
      return AlgebraicMesh::self_test();
    }

    /// Resolve mesh update: Update current nodal
    /// positions via algebraic node update.
    /// [Doesn't make sense to use this mesh with SolidElements anyway,
    /// so we buffer the case if update_all_solid_nodes is set to
    /// true.]
    void node_update(const bool& update_all_solid_nodes = false)
    {
#ifdef PARANOID
      if (update_all_solid_nodes)
      {
        std::string error_message =
          "Doesn't make sense to use an AlgebraicMesh with\n";
        error_message +=
          "SolidElements so specifying update_all_solid_nodes=true\n";
        error_message += "doesn't make sense either\n";

        std::string function_name =
          " AlgebraicRefineableQuarterCircleSectorMesh::";
        function_name += "node_update()";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      AlgebraicMesh::node_update();
    }


    /// Implement the algebraic node update function for a node
    /// at time level t (t=0: present; t>0: previous): Update with
    /// the node's first (default) update function.
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt)
    {
      // Update with the update function for the node's first (default)
      // node update fct
      unsigned id = node_pt->node_update_fct_id();

      switch (id)
      {
        case Central_box:

          // Central box
          node_update_in_central_box(t, node_pt);
          break;


        case Lower_right_box:

          // Lower right box
          node_update_in_lower_right_box(t, node_pt);
          break;

        case Upper_left_box:

          // Upper left box
          node_update_in_upper_left_box(t, node_pt);
          break;

        default:

          std::ostringstream error_message;
          error_message << "The node update fct id is " << id
                        << ", but it should only be one of " << Central_box
                        << ", " << Lower_right_box << " or " << Upper_left_box
                        << std::endl;
          std::string function_name =
            " AlgebraicRefineableQuarterCircleSectorMesh::";
          function_name += "algebraic_node_update()";

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }
    }

    /// Update the node update info for specified algebraic node
    /// following any spatial mesh adaptation.
    void update_node_update(AlgebraicNode*& node_pt)
    {
      // Get all node update fct for this node  (resizes internally)
      Vector<int> id;
      node_pt->node_update_fct_id(id);

      // Loop over all update fcts
      unsigned n_update = id.size();
      for (unsigned i = 0; i < n_update; i++)
      {
        switch (id[i])
        {
          case Central_box:

            // Central box: no update
            break;


          case Lower_right_box:

            // Lower right box
            update_node_update_in_lower_right_box(node_pt);
            break;

          case Upper_left_box:

            // Upper left box
            update_node_update_in_upper_left_box(node_pt);
            break;

          default:

            // Never get here....
            std::ostringstream error_message;
            error_message << "Node update fct id is " << id[i]
                          << ", but it should only be one of" << Central_box
                          << ", " << Lower_right_box << " or " << Upper_left_box
                          << std::endl;

            std::string function_name =
              " AlgebraicRefineableQuarterCircleSectorMesh::";
            function_name += "update_node_update()";

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        }
      }
    }

  private:
    /// Remesh function ids
    enum
    {
      Central_box,
      Lower_right_box,
      Upper_left_box
    };


    /// Fractional width of central box
    double Lambda_x;

    /// Fractional height of central box
    double Lambda_y;

    /// Algebraic update function for a node that is located
    /// in the central box
    void node_update_in_central_box(const unsigned& t, AlgebraicNode*& node_pt);

    /// Algebraic update function for a node that is located
    /// in the lower right box
    void node_update_in_lower_right_box(const unsigned& t,
                                        AlgebraicNode*& node_pt);

    /// Algebraic update function for a node that is located
    /// in the upper left box
    void node_update_in_upper_left_box(const unsigned& t,
                                       AlgebraicNode*& node_pt);

    /// Setup algebraic update operation for all nodes
    void setup_algebraic_node_update();


    /// Update algebraic node update function for nodes in
    /// lower right box
    void update_node_update_in_lower_right_box(AlgebraicNode*& node_pt);

    /// Update algebraic node update function for nodes
    /// in upper left box
    void update_node_update_in_upper_left_box(AlgebraicNode*& node_pt);
  };


} // namespace oomph

#endif
