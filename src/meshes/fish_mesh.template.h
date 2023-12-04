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
#ifndef OOMPH_FISH_MESH_HEADER
#define OOMPH_FISH_MESH_HEADER

// Headers
#include "../generic/refineable_quad_mesh.h"
#include "../generic/macro_element.h"
#include "../generic/domain.h"
#include "../generic/algebraic_elements.h"
#include "../generic/quad_mesh.h"
#include "../generic/macro_element_node_update_element.h"

// Include algebraic elements
#include "../generic/algebraic_elements.h"

// Include the macro element node update elements
#include "../generic/macro_element_node_update_element.h"


// Include the headers file for domain
#include "fish_domain.h"

namespace oomph
{
  //=================================================================
  /// Fish shaped mesh. The geometry is defined by
  /// the Domain object FishDomain.
  //=================================================================
  template<class ELEMENT>
  class FishMesh : public virtual QuadMeshBase
  {
  public:
    /// Constructor: Pass pointer to timestepper
    /// (defaults to the (Steady) default timestepper defined in Mesh)
    FishMesh(TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Constructor: Pass pointer GeomObject that defines
    /// the fish's back and pointer to timestepper
    /// (defaults to the (Steady) default timestepper defined in Mesh)
    FishMesh(GeomObject* back_pt,
             TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Destructor: Kill the geom object that represents the fish's back
    /// (if necessary)
    virtual ~FishMesh()
    {
      if (Must_kill_fish_back)
      {
        delete Back_pt;
        Back_pt = 0;
      }
    }

    /// Access function to geom object that represents the fish's back
    GeomObject*& fish_back_pt()
    {
      return Back_pt;
    }


    /// Access function to FishDomain
    FishDomain*& domain_pt()
    {
      return Domain_pt;
    }

  protected:
    /// Remesh function ids
    enum
    {
      Lower_body,
      Upper_body,
      Lower_fin,
      Upper_fin
    };

    /// Build the mesh, using the geometric object identified by Back_pt
    void build_mesh(TimeStepper* time_stepper_pt);

    /// Pointer to fish back
    GeomObject* Back_pt;

    /// Pointer to domain
    FishDomain* Domain_pt;

    /// Do I need to kill the fish back geom object?
    bool Must_kill_fish_back;
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Refineable fish-shaped mesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Refineable fish shaped mesh. The geometry is defined by
  /// the Domain object FishDomain.
  //=============================start_adaptive_fish_mesh============
  template<class ELEMENT>
  class RefineableFishMesh : public virtual FishMesh<ELEMENT>,
                             public RefineableQuadMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to timestepper -- defaults to (Steady)
    /// default timestepper defined in the Mesh base class
    RefineableFishMesh(
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // FishMesh<...>. Only need to setup adaptive information.

      // Do what it says....
      setup_adaptivity();

    } // end of constructor


    /// Constructor: Pass pointer GeomObject that defines
    /// the fish's back and pointer to timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh)
    RefineableFishMesh(
      GeomObject* back_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(back_pt, time_stepper_pt)
    {
      // Nodal positions etc. were created in constructor for
      // FishMesh<...>. Only need to setup adaptive information.

      // Do what it says....
      setup_adaptivity();
    }

    /// Destructor: Empty -- all cleanup gets handled in the base
    /// classes
    virtual ~RefineableFishMesh() {}


  protected:
    /// Setup all the information that's required for spatial adaptivity:
    /// Set pointers to macro elements and build quadtree forest.
    /// (contained in separate function as this functionality is common
    /// to both constructors),
    void setup_adaptivity();

  }; // end adaptive fish mesh


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdate-version of RefineableFishMesh
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////

  // Forward declaration
  class MacroElementNodeUpdateNode;


  //========================================================================
  /// Refineable fish shaped mesh with MacroElement-based node update.
  /// The fish's back is represented by a specified geometric object.
  /// Some or all of the geometric Data in that geometric object
  /// may contain unknowns in the global Problem. The dependency
  /// on these unknowns is taken into account when setting up
  /// the Jacobian matrix of the elements. For this purpose,
  /// the element (whose type is specified by the template parameter)
  /// must inherit from MacroElementNodeUpdateElementBase.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateRefineableFishMesh
    : public virtual MacroElementNodeUpdateMesh,
      public virtual RefineableFishMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer GeomObject that defines
    /// the fish's back and pointer to timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh).
    MacroElementNodeUpdateRefineableFishMesh(
      GeomObject* back_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(back_pt, time_stepper_pt),
        RefineableFishMesh<ELEMENT>(time_stepper_pt)
    {
#ifdef PARANOID
      ELEMENT* el_pt = new ELEMENT;
      if (dynamic_cast<MacroElementNodeUpdateElementBase*>(el_pt) == 0)
      {
        std::ostringstream error_message;
        error_message << "Base class for ELEMENT in "
                      << "MacroElementNodeUpdateRefineableFishMesh needs"
                      << "to be of type MacroElementNodeUpdateElement!\n";
        error_message << "Whereas it is: typeid(el_pt).name()"
                      << typeid(el_pt).name() << std::endl;

        std::string function_name =
          "MacroElementNodeUpdateRefineableFishMesh::\n";
        function_name += "MacroElementNodeUpdateRefineableFishMesh()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      delete el_pt;
#endif


      // Setup all the information that's required for MacroElement-based
      // node update: Tell the elements that their geometry depends on the
      // fishback geometric object
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
            "MacroElementNodeUpdateRefineableFishMesh::\n";
          function_name += "MacroElementNodeUpdateRefinableFishMesh()";

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // There's just one GeomObject
        Vector<GeomObject*> geom_object_pt(1);
        geom_object_pt[0] = this->Back_pt;

        // Tell the element which geom objects its macro-element-based
        // node update depends on
        el_pt->set_node_update_info(geom_object_pt);
      }

      // Add the geometric object(s) for the wall to the mesh's storage
      Vector<GeomObject*> geom_object_pt(1);
      geom_object_pt[0] = this->Back_pt;
      MacroElementNodeUpdateMesh::set_geom_object_vector_pt(geom_object_pt);

      // Fill in the domain pointer to the mesh's storage in the base class
      MacroElementNodeUpdateMesh::macro_domain_pt() = this->domain_pt();
    }

    /// Destructor: empty
    virtual ~MacroElementNodeUpdateRefineableFishMesh() {}

    /// Resolve mesh update: NodeUpdate current nodal
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
          "MacroElementNodeUpdateRefineableFishMesh::";
        function_name += "node_update()";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      MacroElementNodeUpdateMesh::node_update();
    }
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // AlgebraicElement fish-shaped mesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Fish shaped mesh with algebraic node update function for nodes.
  //=================================================================
  template<class ELEMENT>
  class AlgebraicFishMesh : public AlgebraicMesh,
                            public virtual FishMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to timestepper.
    /// (defaults to (Steady) default timestepper defined in Mesh)
    AlgebraicFishMesh(TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(time_stepper_pt)
    {
      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Constructor: Pass pointer GeomObject that defines
    /// the fish's back and pointer to timestepper
    /// (defaults to (Steady) default timestepper defined in Mesh).
    AlgebraicFishMesh(GeomObject* back_pt,
                      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(back_pt, time_stepper_pt)
    {
      // Add the geometric object to the list associated with this AlgebraicMesh
      AlgebraicMesh::add_geom_object_list_pt(back_pt);

      // Setup algebraic node update operations
      setup_algebraic_node_update();
    }

    /// Destructor: empty
    virtual ~AlgebraicFishMesh() {}

    /// Update nodal position at time level t (t=0: present;
    /// t>0: previous)
    void algebraic_node_update(const unsigned& t, AlgebraicNode*& node_pt)
    {
      // Update with the update function for the node's first (default)
      // node update fct
      unsigned id = node_pt->node_update_fct_id();

      // Upper/lower body
      if ((id == this->Lower_body) || (id == this->Upper_body))
      {
        node_update_in_body(t, node_pt);
      }
      // Upper/lower fin
      else if ((id == this->Lower_fin) || (id == this->Upper_fin))
      {
        node_update_in_fin(t, node_pt);
      }
      else
      {
        std::ostringstream error_message;
        error_message << "The node update fct id is " << id
                      << ", but it should only be one of " << this->Lower_body
                      << ", " << this->Upper_body << ", " << this->Lower_fin
                      << " or " << this->Upper_fin << std::endl;
        std::string function_name =
          "AlgebraicFishMesh::algebraic_node_update()";

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    /// Resolve the node update function (we neither want the broken
    /// empty one in the Mesh base class nor the macro-element-based one in the
    /// RefineableQuadMesh base class but the AlgebraicElement one). [It doesn't
    /// make sense to use this mesh with SolidElements so we buffer the case if
    /// update_all_solid_nodes is set to true.]
    virtual void node_update(const bool& update_all_solid_nodes = false)
    {
#ifdef PARANOID
      if (update_all_solid_nodes)
      {
        std::string error_message =
          "Doesn't make sense to use an AlgebraicMesh with\n";
        error_message +=
          "SolidElements so specifying update_all_solid_nodes=true\n";
        error_message += "doesn't make sense either\n";

        std::string function_name = "AlgebraicFishMesh::node_update()";

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      AlgebraicMesh::node_update();
    }

    /// Update the geometric references that are used
    /// to update node after mesh adaptation.
    /// We're assuming that the GeomObject that specifies
    /// the fish back does not have sub-objects, therefore
    /// no update is required -- all reference
    /// values can simply be scaled. We simply
    /// paranoid-check that this is actually the case,
    /// by checking if locate_zeta() returns the
    /// original data.
    void update_node_update(AlgebraicNode*& node_pt)
    {
#ifdef PARANOID

      // Get the start and end Lagrangian coordinates on the
      // wall from the domain:

      /// Start coordinate on wall (near nose)
      double xi_nose = this->Domain_pt->xi_nose();

      /// End coordinate on wall (near tail)
      double xi_tail = this->Domain_pt->xi_tail();

      /// Check halfway along the object
      Vector<double> zeta(1);
      zeta[0] = 0.5 * (xi_nose + xi_tail);

      Vector<double> s(1);
      GeomObject* geom_obj_pt = 0;
      this->Back_pt->locate_zeta(zeta, geom_obj_pt, s);

      if ((geom_obj_pt != this->Back_pt) || (s[0] != zeta[0]))
      {
        std::ostringstream error_message;
        error_message << "AlgebraicFishMesh only works with GeomObjects\n"
                      << "that do not contain sub-elements (e.g. GeomObjects\n"
                      << "that represent a wall finite element mesh!\n"
                      << "Back_pt    : " << this->Back_pt << std::endl
                      << "geom_obj_pt: " << geom_obj_pt << std::endl
                      << "s[0]       : " << s[0] << std::endl
                      << "zeta[0]    : " << zeta[0] << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }


  protected:
    /// Algebraic update function for nodes in upper/lower body
    void node_update_in_body(const unsigned& t, AlgebraicNode*& node_pt);

    /// Algebraic update function for nodes in upper/lower fin
    void node_update_in_fin(const unsigned& t, AlgebraicNode*& node_pt);

    /// Setup algebraic update operation for all nodes
    /// (separate function because this task needs to be performed by
    /// both constructors)
    void setup_algebraic_node_update();
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Refineable algebraic element fish-shaped mesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Refineable fish shaped mesh with algebraic node update function.
  //=================================================================
  template<class ELEMENT>
  class AlgebraicRefineableFishMesh : public AlgebraicFishMesh<ELEMENT>,
                                      public RefineableFishMesh<ELEMENT>
  {
  public:
    /// Constructor: Pass pointer to timestepper.
    /// (defaults to (Steady) default timestepper defined in Mesh)
    // Note: FishMesh is virtual base and its constructor is automatically
    // called first! --> this is where we need to build the mesh;
    // the constructors of the derived meshes don't call the
    // base constructor again and simply add the extra functionality.
    AlgebraicRefineableFishMesh(
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(time_stepper_pt),
        AlgebraicFishMesh<ELEMENT>(time_stepper_pt),
        RefineableFishMesh<ELEMENT>(time_stepper_pt)
    {
    }


    /// Constructor: Pass pointer GeomObject that defines
    /// the fish's back and pointer to timestepper.
    /// (defaults to (Steady) default timestepper defined in Mesh)
    // Note: FishMesh is virtual base and its constructor is automatically
    // called first! --> this is where we need to build the mesh;
    // the constructors of the derived meshes don't call the
    // base constructor again and simply add the extra functionality.
    AlgebraicRefineableFishMesh(
      GeomObject* back_pt,
      TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
      : FishMesh<ELEMENT>(back_pt, time_stepper_pt),
        AlgebraicFishMesh<ELEMENT>(back_pt, time_stepper_pt),
        RefineableFishMesh<ELEMENT>(back_pt, time_stepper_pt)
    {
    }


    /// Destructor: empty
    virtual ~AlgebraicRefineableFishMesh() {}

    /// Resolve node update function: Use the one defined
    /// in the AlgebraicFishMesh (where the bool flag is explained)
    void node_update(const bool& update_all_solid_nodes = false)
    {
      AlgebraicFishMesh<ELEMENT>::node_update(update_all_solid_nodes);
    }
  };

} // namespace oomph

#endif
