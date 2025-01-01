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
#ifndef OOMPH_QUARTER_TUBE_MESH_HEADER
#define OOMPH_QUARTER_TUBE_MESH_HEADER

// Headers
#include "../generic/refineable_brick_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/brick_mesh.h"
#include "../generic/macro_element_node_update_element.h"


// Include the headers file for domain
#include "quarter_tube_domain.h"

namespace oomph
{
  //====================================================================
  /// 3D quarter tube mesh class.
  /// The domain is specified by the GeomObject that identifies
  /// boundary 3. Non-refineable base version!
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: "Inflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$
  ///               on the geometric object that specifies the wall.
  /// - Boundary 1: Plane x=0
  /// - Boundary 2: Plane y=0
  /// - Boundary 3: The curved wall
  /// - Boundary 4: "Outflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{hi} \f$
  ///               on the geometric object that specifies the wall.
  ///
  /// IMPORTANT NOTE: The interface looks more general than it should.
  ///                 The toplogy must remain that of a quarter tube,
  ///                 or the mesh generation will break.
  //====================================================================
  template<class ELEMENT>
  class QuarterTubeMesh : public virtual BrickMeshBase
  {
  public:
    /// Constructor: Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    QuarterTubeMesh(GeomObject* wall_pt,
                    const Vector<double>& xi_lo,
                    const double& fract_mid,
                    const Vector<double>& xi_hi,
                    const unsigned& nlayer,
                    TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor: empty
    virtual ~QuarterTubeMesh()
    {
      delete Domain_pt;
    }

    /// Access function to GeomObject representing wall
    GeomObject*& wall_pt()
    {
      return Wall_pt;
    }

    /// Access function to domain
    QuarterTubeDomain* domain_pt()
    {
      return Domain_pt;
    }

    /// Function pointer for function that squashes
    /// the outer macro elements towards
    /// the wall by mapping the input value of the "radial" macro element
    /// coordinate to the return value (defined in the underlying Domain object)
    QuarterTubeDomain::BLSquashFctPt& bl_squash_fct_pt()
    {
      return Domain_pt->bl_squash_fct_pt();
    }


    /// Function pointer for function for axial spacing
    virtual QuarterTubeDomain::AxialSpacingFctPt& axial_spacing_fct_pt()
    {
      return Domain_pt->axial_spacing_fct_pt();
    }

    /// Access function to underlying domain
    QuarterTubeDomain* domain_pt() const
    {
      return Domain_pt;
    }

  protected:
    /// Pointer to domain
    QuarterTubeDomain* Domain_pt;

    /// Pointer to the geometric object that represents the curved wall
    GeomObject* Wall_pt;

    /// Lower limits for the coordinates along the wall
    Vector<double> Xi_lo;

    /// Fraction along wall where outer ring is to be divided
    double Fract_mid;

    /// Upper limits for the coordinates along the wall
    Vector<double> Xi_hi;
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //=============================================================
  /// Adaptative version of the QuarterTubeMesh base mesh.
  /// The domain is specified by the GeomObject that identifies
  /// boundary 3.
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: "Inflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$
  ///               on the geometric object that specifies the wall.
  /// - Boundary 1: Plane x=0
  /// - Boundary 2: Plane y=0
  /// - Boundary 3: The curved wall
  /// - Boundary 4: "Outflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{hi} \f$
  ///               on the geometric object that specifies the wall.
  //=============================================================
  template<class ELEMENT>
  class RefineableQuarterTubeMesh : public virtual QuarterTubeMesh<ELEMENT>,
                                    public RefineableBrickMesh<ELEMENT>

  {
  public:
    /// Constructor for adaptive deformable quarter tube mesh class.
    /// The domain is specified by the GeomObject that
    /// identifies boundary 3. Pass pointer to geometric object that
    /// specifies the wall, start and end coordinates on the
    /// geometric object, and the fraction along
    /// which the dividing line is to be placed, and the timestepper.
    /// Timestepper defaults to Steady dummy timestepper.
    RefineableQuarterTubeMesh(
      GeomObject* wall_pt,
      const Vector<double>& xi_lo,
      const double& fract_mid,
      const Vector<double>& xi_hi,
      const unsigned& nlayer,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuarterTubeMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, nlayer, time_stepper_pt)
    {
      // Loop over all elements and set macro element pointer
      for (unsigned ielem = 0; ielem < QuarterTubeMesh<ELEMENT>::nelement();
           ielem++)
      {
        dynamic_cast<RefineableQElement<3>*>(
          QuarterTubeMesh<ELEMENT>::element_pt(ielem))
          ->set_macro_elem_pt(this->Domain_pt->macro_element_pt(ielem));
      }


      // Setup Octree forest: Turn elements into individual octrees
      // and plant in forest
      Vector<TreeRoot*> trees_pt;
      for (unsigned iel = 0; iel < QuarterTubeMesh<ELEMENT>::nelement(); iel++)
      {
        FiniteElement* el_pt = QuarterTubeMesh<ELEMENT>::finite_element_pt(iel);
        ELEMENT* ref_el_pt = dynamic_cast<ELEMENT*>(el_pt);
        OcTreeRoot* octree_root_pt = new OcTreeRoot(ref_el_pt);
        trees_pt.push_back(octree_root_pt);
      }
      this->Forest_pt = new OcTreeForest(trees_pt);

#ifdef PARANOID
      // Run self test
      unsigned success_flag =
        dynamic_cast<OcTreeForest*>(this->Forest_pt)->self_test();
      if (success_flag == 0)
      {
        oomph_info << "Successfully built octree forest " << std::endl;
      }
      else
      {
        throw OomphLibError("Trouble in building octree forest ",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Destructor: empty
    virtual ~RefineableQuarterTubeMesh() {}
  };


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdate-version of RefineableQuarterTubeMesh
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////

  class MacroElementNodeUpdateNode;

  //========================================================================
  /// MacroElementNodeUpdate version of RefineableQuarterTubeMesh
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateRefineableQuarterTubeMesh
    : public virtual MacroElementNodeUpdateMesh,
      public virtual RefineableQuarterTubeMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object, start and
    /// end coordinates on the geometric object and the fraction along
    /// which the dividing line is to be placed when updating the nodal
    /// positions, and timestepper (defaults to (Steady) default timestepper
    /// defined in Mesh). Setup the refineable mesh (by calling the
    /// constructor for the underlying  RefineableQuarterTubeMesh)
    /// and the algebraic node update functions for nodes.
    MacroElementNodeUpdateRefineableQuarterTubeMesh(
      GeomObject* wall_pt,
      const Vector<double>& xi_lo,
      const double& fract_mid,
      const Vector<double>& xi_hi,
      const unsigned& nlayer,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : MacroElementNodeUpdateMesh(),
        RefineableQuarterTubeMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, nlayer, time_stepper_pt),
        QuarterTubeMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, nlayer, time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;
        error_message << "Base class for ELEMENT in "
                      << "MacroElementNodeUpdateRefineableQuarterTubeMesh needs"
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
    virtual ~MacroElementNodeUpdateRefineableQuarterTubeMesh() {}

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


  //======================================================================
  /// AlgebraicMesh version of RefineableQuarterTubeMesh
  //=====================================================================


  //====================================================================
  /// Algebraic 3D quarter tube mesh class.
  ///
  /// The mesh boundaries are numbered as follows:
  /// - Boundary 0: "Inflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{lo} \f$
  ///               on the geometric object that specifies the wall.
  /// - Boundary 1: Plane x=0
  /// - Boundary 2: Plane y=0
  /// - Boundary 3: The curved wall - specified by the GeomObject
  ///               passed to the mesh constructor.
  /// - Boundary 4: "Outflow" cross section; located along the
  ///               line parametrised by \f$ \xi_0 =  \xi_0^{hi} \f$
  ///               on the geometric object that specifies the wall.
  //====================================================================


  //========================================================================
  /// Algebraic version of RefineableQuarterTubeMesh
  ///
  /// Cross section through mesh looking along tube.........
  ///
  ///                      ---___
  ///                     |      ---____
  ///                     |              -   BOUNDARY 3
  ///                     |                /
  ///                     |  [Region 2]   /  |
  ///                     |              /     |
  ///                     | N           /        |
  ///                     | |_ E       /          |
  ///       BOUNDARY 1    |------------            |
  ///                     |            |            |
  ///                     | [Region 0] | [Region 1] |  ^
  ///                     |            |            | / \  direction of
  ///                     | N          |    N       |  |   2nd Lagrangian
  ///                     | |_ E       |    |_ E    |  |   coordinate
  ///                     |____________|____________|  |   along wall GeomObject
  ///
  ///                           BOUNDARY 2
  ///
  /// The Domain is built of slices each consisting of three
  /// MacroElements as sketched.
  /// The local coordinates are such that the (E)astern direction
  /// coincides with the positive s_0 direction, while the
  /// (N)orther direction coincides with the positive s_1 direction.
  /// The positive s_2 direction points down the tube.
  ///
  /// Elements need to be derived from AlgebraicElementBase. In
  /// addition to the refinement procedures available for the
  /// RefineableQuarterTubeMesh which forms the basis for this mesh,
  /// three algebraic node update functions are implemented for the nodes
  /// in the three regions defined by the Domain MacroElements.
  /// Note: it is assumed the cross section down the tube is
  /// uniform when setup_algebraic_node_update() is called.
  //========================================================================
  template<class ELEMENT>
  class AlgebraicRefineableQuarterTubeMesh
    : public virtual AlgebraicMesh,
      public RefineableQuarterTubeMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to geometric object, start and
    /// end coordinates of the geometric object and the fraction along
    /// the 2nd Lagrangian coordinate at which the dividing line between
    /// region 1 and region 2 is to be placed, and timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh).
    /// Sets up the refineable mesh (by calling the constructor for the
    /// underlying RefineableQuarterTubeMesh).
    AlgebraicRefineableQuarterTubeMesh(
      GeomObject* wall_pt,
      const Vector<double>& xi_lo,
      const double& fract_mid,
      const Vector<double>& xi_hi,
      const unsigned& nlayer,
      const double centre_box_size = 1.0,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : QuarterTubeMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, nlayer, time_stepper_pt),
        RefineableQuarterTubeMesh<ELEMENT>(
          wall_pt, xi_lo, fract_mid, xi_hi, nlayer, time_stepper_pt),
        Centre_box_size(centre_box_size)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<AlgebraicElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;

        error_message << "Base class for ELEMENT in "
                      << "AlgebraicRefineableQuarterTubeMesh needs"
                      << "to be of type AlgebraicElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name = " AlgebraicRefineableQuarterTubeMesh::\n";
        function_name += "AlgebraicRefineableQuarterTubeMesh()";

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

      // Ensure nodes are in their default position
      node_update();
    }

    /// Run self-test for algebraic mesh -- return 0/1 for OK/failure
    unsigned self_test()
    {
      return AlgebraicMesh::self_test();
    }

    /// Broken version of the QuarterTubeDomain function
    /// Function is broken because axial spacing isn't implemented
    /// yet for the Algebraic version of the RefineableQuarterTubeMesh.
    /// Note: this function must be used BEFORE algebraic_node_update(...)
    /// is called.
    QuarterTubeDomain::AxialSpacingFctPt& axial_spacing_fct_pt()
    {
      std::ostringstream error_message;
      error_message << "AxialSpacingFctPt has not been implemented "
                    << "for the AlgebraicRefineableQuarterTubeMesh\n";

      std::string function_name =
        " AlgebraicRefineableQuarterTubeMesh::AxialSpacingFctPt()";

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      return this->Domain_pt->axial_spacing_fct_pt();
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

        std::string function_name = " AlgebraicRefineableQuarterTubeMesh::";
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
        case Central_region:

          // Central region
          node_update_central_region(t, node_pt);
          break;

        case Lower_right_region:

          // Lower right region
          node_update_lower_right_region(t, node_pt);
          break;

        case Upper_left_region:

          // Upper left region
          node_update_upper_left_region(t, node_pt);
          break;

        default:

          std::ostringstream error_message;
          error_message << "The node update fct id is " << id
                        << ", but it should only be one of " << Central_region
                        << ", " << Lower_right_region << " or "
                        << Upper_left_region << std::endl;
          std::string function_name = " AlgebraicRefineableQuarterTubeMesh::";
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
      // Get all node update fct for this node (resizes internally)
      Vector<int> id;
      node_pt->node_update_fct_id(id);

      // Loop over all update fcts
      unsigned n_update = id.size();
      for (unsigned i = 0; i < n_update; i++)
      {
        update_node_update_in_region(node_pt, id[i]);
      }
    }

  private:
    /// Size of centre box
    double Centre_box_size;

    /// Remesh function ids
    enum
    {
      Central_region,
      Lower_right_region,
      Upper_left_region
    };

    /// Fractional width of central region
    double Lambda_x;

    /// Fractional height of central region
    double Lambda_y;

    /// Algebraic update function for a node that is located
    /// in the central region
    void node_update_central_region(const unsigned& t, AlgebraicNode*& node_pt);

    /// Algebraic update function for a node that is located
    /// in the lower-right region
    void node_update_lower_right_region(const unsigned& t,
                                        AlgebraicNode*& node_pt);

    /// Algebraic update function for a node that is located
    /// in the upper-left region
    void node_update_upper_left_region(const unsigned& t,
                                       AlgebraicNode*& node_pt);

    /// Setup algebraic update operation for all nodes
    void setup_algebraic_node_update();

    /// Update algebraic node update function for nodes in
    /// the region defined by region_id
    void update_node_update_in_region(AlgebraicNode*& node_pt, int& region_id);
  };


} // namespace oomph
#endif
