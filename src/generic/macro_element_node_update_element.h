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
#ifndef OOMPH_MACRO_ELEMENT_NODE_UPDATE_ELEMENTS_HEADER
#define OOMPH_MACRO_ELEMENT_NODE_UPDATE_ELEMENTS_HEADER

#include "geom_objects.h"
#include "mesh.h"
#include "elements.h"
#include "element_with_moving_nodes.h"
#include "domain.h"

namespace oomph
{
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdate nodes
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //========================================================================
  /// MacroElementNodeUpdate nodes are nodes with a positional update
  /// function, based on their element's MacroElement representation.
  //========================================================================
  class MacroElementNodeUpdateNode : public Node
  {
  public:
    /// Constructor for steady node of spatial
    /// dimension n_dim, with n_position_type generalised coordinates
    /// and with initial_nvalue dofs.
    MacroElementNodeUpdateNode(const unsigned& n_dim,
                               const unsigned& n_position_type,
                               const unsigned& initial_nvalue)
      : Node(n_dim, n_position_type, initial_nvalue)
    {
      // By default, only the nodal position is updated and no auxiliary
      // updates of function values are performed.
    }

    /// Constructor for bog-standard node of spatial
    /// dimension n_dim, with n_position_type generalised coordinates,
    /// with initial_nvalue dofs and with time dependence.
    MacroElementNodeUpdateNode(TimeStepper* time_stepper_pt,
                               const unsigned& n_dim,
                               const unsigned& n_position_type,
                               const unsigned& initial_nvalue)
      : Node(time_stepper_pt, n_dim, n_position_type, initial_nvalue)
    {
      // By default, only the nodal position is updated and no auxiliary
      // updates of function values are performed.
    }

    /// Broken copy constructor
    MacroElementNodeUpdateNode(const MacroElementNodeUpdateNode&) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const MacroElementNodeUpdateNode&) = delete;*/

    /// Destructor (empty)
    virtual ~MacroElementNodeUpdateNode() {}

    /// Update the current nodal position. If
    /// required, perform the auxiliary update of nodal values.
    /// If update_all_time_levels_for_new_node==true, previous
    /// positions are also updated -- as indicated by the name
    /// of this flag, this should only be done for newly
    /// created nodes, when this function is called from
    /// MacroElementNodeUpdateElementBase::build_macro_element_node_update_node(...)
    void node_update(const bool& update_all_time_levels_for_new_node = false);

    ///  Pointer to finite element that performs the update by referring
    /// to its macro-element representation (Access required...)
    FiniteElement*& node_update_element_pt()
    {
      return Node_update_element_pt;
    }


    /// Vector of local coordinates of node with the finite element that
    /// performs the MacroElement-based node update operation
    Vector<double>& s_in_node_update_element()
    {
      return S_in_node_update_element;
    }

    /// Number of geometric objects involved in node update function
    unsigned ngeom_object() const
    {
      return Geom_object_pt.size();
    }

    /// Vector of (pointers to) geometric objects involved in
    /// node update function
    Vector<GeomObject*>& geom_object_pt()
    {
      return Geom_object_pt;
    }


    /// Pointer to i-th geometric object involved in
    /// node update function
    GeomObject* geom_object_pt(const unsigned& i)
    {
      return Geom_object_pt[i];
    }


    /// Return vector of geometric objects involved in
    /// node update function
    Vector<GeomObject*>& vector_geom_object_pt()
    {
      return Geom_object_pt;
    }

    /// Return all geometric objects that affect the node update
    inline GeomObject** all_geom_object_pt()
    {
      if (Geom_object_pt.size() > 0)
      {
        return &(Geom_object_pt[0]);
      }
      else
      {
        return 0;
      }
    }

    /// Set node update information for node:
    /// Pass the pointer to the element that performs the update operation,
    /// the vector containing the node's local coordinates in that
    /// element and the vector of (pointers to) the geometric objects
    /// that affect the node update.
    void set_node_update_info(FiniteElement* node_update_element_pt,
                              const Vector<double>& s_in_node_update_element,
                              const Vector<GeomObject*>& geom_object_pt)
    {
      Node_update_element_pt = node_update_element_pt;
      S_in_node_update_element = s_in_node_update_element;
      Geom_object_pt = geom_object_pt;
    }


  private:
    /// Pointer to finite element that performs the node update
    /// by referring to its macro-element representation
    FiniteElement* Node_update_element_pt;

    /// Vector containing the node's local coordinates in node update
    /// element.
    Vector<double> S_in_node_update_element;

    /// Vector of geometric objects that are involved
    /// in the node update operation
    Vector<GeomObject*> Geom_object_pt;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // MacroElementNodeUpdate elements
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Base class for elements that allow MacroElement-based node update
  //========================================================================
  class MacroElementNodeUpdateElementBase
  {
  public:
    /// Constructor (empty)
    MacroElementNodeUpdateElementBase() {}

    /// Broken copy constructor
    MacroElementNodeUpdateElementBase(
      const MacroElementNodeUpdateElementBase&) = delete;

    /// Broken assignment operator
    void operator=(const MacroElementNodeUpdateElementBase&) = delete;

    /// Virtual destructor (empty)
    virtual ~MacroElementNodeUpdateElementBase() {}

    /// Set node update information:
    /// Pass the vector of (pointers to) the geometric objects
    /// that affect the node update. This gets passed on to all nodes in
    /// the element.
    virtual void set_node_update_info(
      const Vector<GeomObject*>& geom_object_pt) = 0;

    /// Number of geometric objects involved in node update function
    inline unsigned ngeom_object()
    {
      return Geom_object_pt.size();
    }

    /// Vector of (pointers to) geometric objects involved in
    /// node update function
    Vector<GeomObject*>& geom_object_pt()
    {
      return Geom_object_pt;
    }

    /// Pointer to i-th geometric object involved in
    /// node update function
    GeomObject* geom_object_pt(const unsigned& i)
    {
      return Geom_object_pt[i];
    }


  protected:
    /// Vector of geometric objects that are involved
    /// in the node update operation
    Vector<GeomObject*> Geom_object_pt;
  };


  //========================================================================
  /// MacroElementNodeUpdate elements are elements that can not only be updated
  /// via their MacroElement representation (in princple any FiniteElement
  /// could do that...) but also allows the geometric Data contained
  /// in the GeomObjects that affect the MacroElement-based node update
  /// operations to be unknowns in the overall Problem.
  ///
  /// The element wraps around the ELEMENT specified by the template
  /// parameter and computes the derivatives of the residual vector
  /// with respect to the geometric Data (needed in the setup of
  /// the element's Jacobian matrix) by finite differencing.
  /// Otherwise the element behaves exactly like the templace element.
  //========================================================================
  template<class ELEMENT>
  class MacroElementNodeUpdateElement
    : public ElementWithSpecificMovingNodes<ELEMENT,
                                            MacroElementNodeUpdateNode>,
      public MacroElementNodeUpdateElementBase
  {
  public:
    /// Constructor: Call constructor of underlying element
    MacroElementNodeUpdateElement()
      : ElementWithSpecificMovingNodes<ELEMENT, MacroElementNodeUpdateNode>(),
        MacroElementNodeUpdateElementBase()
    {
    }

    /// Constructor used for face elements
    MacroElementNodeUpdateElement(FiniteElement* const& element_pt,
                                  const int& face_index)
      : ElementWithSpecificMovingNodes<ELEMENT, MacroElementNodeUpdateNode>(
          element_pt, face_index),
        MacroElementNodeUpdateElementBase()
    {
    }

    /// Broken copy constructor
    MacroElementNodeUpdateElement(const MacroElementNodeUpdateElement&) =
      delete;

    /// Empty destructor to clean up allocated memory
    ~MacroElementNodeUpdateElement() {}

    /// Broken assignment operator
    /*void operator=(const MacroElementNodeUpdateElement&) = delete;*/

    /// Set node update information:
    /// Pass the vector of (pointers to) the geometric objects
    /// that affect the node update. This gets passed on to all nodes in
    /// the element.
    void set_node_update_info(const Vector<GeomObject*>& geom_object_pt)
    {
      // Store local copy of geom object vector, so it can be passed on
      // to son elements (and their nodes) during refinement
      unsigned ngeom_object = geom_object_pt.size();
      Geom_object_pt.resize(ngeom_object);
      for (unsigned i = 0; i < ngeom_object; i++)
      {
        Geom_object_pt[i] = geom_object_pt[i];
      }

      // Loop over nodes in element
      unsigned n_node = this->nnode();
      for (unsigned j = 0; j < n_node; j++)
      {
        // Get local coordinate in element (Vector sets its own size)
        Vector<double> s_in_node_update_element;
        this->local_coordinate_of_node(j, s_in_node_update_element);

        // Pass the lot to the node
        static_cast<MacroElementNodeUpdateNode*>(this->node_pt(j))
          ->set_node_update_info(
            this, s_in_node_update_element, geom_object_pt);
      }
    }


    /// Rebuild after unrefinement: Reset the node update information
    /// for all nodes so that the nodes get updated by this element.
    /// If we don't do that, some nodes might still want to be updated
    /// by elements that no longer exist which leads to the most wonderful
    /// seg fault... Afterwards, call the template element's own
    /// rebuild_from_son() function (if it's a RefineableElement)
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      // First call the element's own rebuild_from_sons() function
      ELEMENT::rebuild_from_sons(mesh_pt);

      // Now loop over nodes in element
      unsigned n_node = this->nnode();
      for (unsigned j = 0; j < n_node; j++)
      {
        // Get local coordinate in element (Vector sets its own size)
        Vector<double> s_in_node_update_element;
        this->local_coordinate_of_node(j, s_in_node_update_element);

        // Pass the lot to the node
        static_cast<MacroElementNodeUpdateNode*>(this->node_pt(j))
          ->set_node_update_info(
            this, s_in_node_update_element, Geom_object_pt);
      }
    }
  };


  //========================================================================
  /// MacroElementNodeUpdateMeshes contain MacroElementNodeUpdateNodes
  /// which have their own node update functions. When the node's
  /// node_update() function is called, they also perform
  /// any auxiliary update functions, e.g. to update no-slip boundary
  /// conditions on moving domain boundaries.
  //========================================================================
  class MacroElementNodeUpdateMesh : public virtual Mesh
  {
  public:
    /// Constructor (empty)
    MacroElementNodeUpdateMesh() {}

    /// Virtual destructor (empty)
    virtual ~MacroElementNodeUpdateMesh() {}

    /// Broken copy constructor
    MacroElementNodeUpdateMesh(const MacroElementNodeUpdateMesh&) = delete;

    /// Broken assignment operator
    /*void operator=(const MacroElementNodeUpdateMesh&) = delete;*/

    /// Access to Macro_domain_pt for MacroElementNodeUpdateMesh; this
    /// must be filled in by any mesh which inherits from here
    Domain*& macro_domain_pt()
    {
      return Macro_domain_pt;
    }

    /// Update all nodal positions via sparse MacroElement-based
    /// update functions. If a Node is hanging its position is updated
    /// after updating the position of its masters first.
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

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Loop over all nodes and update their positions -- hanging nodes
      // are updated via their masters; auxiliary update function
      // is performed by node, too
      unsigned n_node = nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        MacroElementNodeUpdateNode* nod_pt =
          dynamic_cast<MacroElementNodeUpdateNode*>(node_pt(n));
#ifdef PARANOID
        if (nod_pt == 0)
        {
          std::ostringstream error_message;
          error_message << "Failed to cast to MacroElementNodeUpdateNode.\n"
                        << "Node is of type: " << typeid(node_pt(n)).name()
                        << std::endl;

          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        nod_pt->node_update();
      }

#ifdef OOMPH_HAS_MPI
      // Update positions for external halo nodes attached to this mesh
      // Loop over processors
      for (std::map<unsigned, Vector<Node*>>::iterator it =
             External_halo_node_pt.begin();
           it != External_halo_node_pt.end();
           it++)
      {
        int iproc = (*it).first;
        unsigned n_ext_halo_node = nexternal_halo_node(iproc);
        for (unsigned n = 0; n < n_ext_halo_node; n++)
        {
          MacroElementNodeUpdateNode* nod_pt =
            dynamic_cast<MacroElementNodeUpdateNode*>(
              external_halo_node_pt(iproc, n));
#ifdef PARANOID
          if (nod_pt == 0)
          {
            std::ostringstream error_message;
            error_message
              << "Failed to cast (ext. halo) to MacroElementNodeUpdateNode.\n"
              << "Node is of type: " << typeid(node_pt(n)).name() << std::endl;

            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
#endif
          nod_pt->node_update();
        }
      } // end loop over processors
#endif // (ifdef OOMPH_HAS_MPI)
    }

#ifdef OOMPH_HAS_MPI
    /// Overload the base class distribute function to deal
    /// with halo nodes on halo elements that may have pointers
    /// to macro elements that no longer exist
    void distribute(OomphCommunicator* comm_pt,
                    const Vector<unsigned>& element_domain,
                    Vector<GeneralisedElement*>& deleted_element_pt,
                    DocInfo& doc_info,
                    const bool& report_stats,
                    const bool& overrule_keep_as_halo_element_status)
    {
      // Call underlying Mesh::distribute first
      Mesh::distribute(comm_pt,
                       element_domain,
                       deleted_element_pt,
                       doc_info,
                       report_stats,
                       overrule_keep_as_halo_element_status);

      // Storage for number of processors
      int n_proc = comm_pt->nproc();

      // The original call to set_node_update_info on the
      // non-distributed problem may have set a macro element which no
      // longer exists for some halo nodes which are on halo elements
      // within the distributed Mesh; this deals with the problem by
      // recalling the set_node_update_info for every halo element
      for (int iproc = 0; iproc < n_proc; iproc++)
      {
        Vector<GeneralisedElement*> halo_el_pt = halo_element_pt(iproc);
        unsigned n_halo_el = halo_el_pt.size();
        for (unsigned e = 0; e < n_halo_el; e++)
        {
          // Cast to a MacroElementNodeUpdateElement
          MacroElementNodeUpdateElementBase* macro_el_pt =
            dynamic_cast<MacroElementNodeUpdateElementBase*>(halo_el_pt[e]);

          // The vector of GeomObjects should not change!
          Vector<GeomObject*> geom_object_pt = macro_el_pt->geom_object_pt();

          // So we can just call set_node_update_info for the element!
          macro_el_pt->set_node_update_info(geom_object_pt);
        }
      }
    }
#endif

    /// Set geometric objects associated with MacroElementNodeUpdateMesh;
    /// this must also be called from the constructor of each derived mesh
    void set_geom_object_vector_pt(Vector<GeomObject*> geom_object_vector_pt)
    {
      Geom_object_vector_pt = geom_object_vector_pt;
    }

    /// Access function to the vector of GeomObject
    Vector<GeomObject*> geom_object_vector_pt()
    {
      return Geom_object_vector_pt;
    }

  private:
    /// Vector of GeomObject associated with
    /// MacroElementNodeUpdateNodeMesh
    Vector<GeomObject*> Geom_object_vector_pt;

    /// Domain associated with MacroElementNodeUpdateNodeMesh
    Domain* Macro_domain_pt;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Explicit definition of the face geometry of
  /// MacroElementNodeUpdateElements, which is the same as the face geometry of
  /// the underlying element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<MacroElementNodeUpdateElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    /// Constructor calls the constructor of the underlying ELEMENT.
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };

} // namespace oomph

#endif
