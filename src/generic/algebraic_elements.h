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
#ifndef OOMPH_ALGEBRAIC_ELEMENTS_HEADER
#define OOMPH_ALGEBRAIC_ELEMENTS_HEADER

#include "geom_objects.h"
#include "mesh.h"
#include "elements.h"
#include "domain.h"
#include "element_with_moving_nodes.h"

namespace oomph
{
  // forward references
  class AlgebraicMesh;
  class AlgebraicElementBase;
  class DummyAlgebraicMesh;


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic nodes
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Algebraic nodes are nodes with an algebraic positional update
  /// function.
  //========================================================================
  class AlgebraicNode : public Node
  {
  public:
    /// Default Constructor
    AlgebraicNode() : Node() {}

    /// Constructor for steady algebraic node of spatial
    /// dimension n_dim, with n_position_type generalised coordinates
    /// and with initial_nvalue dofs.
    AlgebraicNode(const unsigned& n_dim,
                  const unsigned& n_position_type,
                  const unsigned& initial_nvalue)
      : Node(n_dim, n_position_type, initial_nvalue)
    {
#ifdef LEAK_CHECK
      LeakCheckNames::AlgebraicNode_build += 1;
#endif

      // Add default node update info
      add_node_update_info(Dummy_node_update_fct_id, // dummy remesh fct ID
                           Dummy_mesh_pt, // dummy mesh
                           Dummy_geom_object_pt, // dummy geom object vector
                           Dummy_ref_value, // dummy  ref vector
                           true); // flag indicating call from
                                  // constructor
    }


    /// Constructor for bog-standard algebraic node of spatial
    /// dimension n_dim, with n_position_type generalised coordinates,
    /// with initial_nvalue dofs and with time dependence.
    AlgebraicNode(TimeStepper* time_stepper_pt,
                  const unsigned& n_dim,
                  const unsigned& n_position_type,
                  const unsigned& initial_nvalue)
      : Node(time_stepper_pt, n_dim, n_position_type, initial_nvalue)
    {
#ifdef LEAK_CHECK
      LeakCheckNames::AlgebraicNode_build += 1;
#endif

      // Add default node update info
      add_node_update_info(Dummy_node_update_fct_id, // dummy remesh fct ID
                           Dummy_mesh_pt, // dummy mesh
                           Dummy_geom_object_pt, // dummy geom object vector
                           Dummy_ref_value, // dummy ref vector
                           true); // flag indicating call from
                                  // constructor
    }

    /// Destructor (empty)
    virtual ~AlgebraicNode()
    {
#ifdef LEAK_CHECK
      LeakCheckNames::AlgebraicNode_build -= 1;
#endif
    }

    /// Broken copy constructor
    AlgebraicNode(const AlgebraicNode&) = delete;

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(const AlgebraicNode&) = delete;*/


    /// Update the current nodal position, using the first
    /// (default) update function if there are multiple ones. If
    /// required perform the auxiliary update of nodal values.
    /// If update_all_time_levels_for_new_node==true, previous
    /// positions are also updated -- as indicated by the name
    /// of this flag, this should only be done for newly
    /// created nodes, when this function is called from
    /// AlgebraicElementBase::setup_algebraic_node_update(...)
    void node_update(const bool& update_all_time_levels_for_new_node = false);


    /// Number of node update fcts
    unsigned nnode_update_fcts()
    {
      // Note: We could read this information out from any one of
      // various maps that that store the information for the
      // different node update functions...
      return Mesh_pt.size();
    }


    /// Default (usually first if there are multiple ones) node update fct id
    int node_update_fct_id()
    {
      return Default_node_update_fct_id;
    }

    /// Return vector of node update fct ids (vector is
    /// resized to contain the correct number of entries). Somewhat costly
    /// to call as map needs to be copied into vector.
    void node_update_fct_id(Vector<int>& id)
    {
      // Resize vector
      id.resize(0);

      // Loop over all entries and copy them across (again, we could
      // get this information from any of the maps...)
      typedef std::map<int, AlgebraicMesh*>::iterator IT;
      for (IT it = Mesh_pt.begin(); it != Mesh_pt.end(); it++)
      {
        id.push_back(it->first);
      }
    }


    /// Default (usually first) mesh that implements update function
    AlgebraicMesh* mesh_pt()
    {
      return Default_it_mesh_pt->second;
    }


    /// Mesh that implements the id-th node update function
    AlgebraicMesh* mesh_pt(const int& id)
    {
      return Mesh_pt[id];
    }


    /// Number of geometric objects involved in id-th update function
    unsigned ngeom_object(const int& id)
    {
      return Geom_object_pt[id].size();
    }


    /// Number of geometric objects involved in default (usually first)
    /// update function
    unsigned ngeom_object() const
    {
      return Default_it_geom_object_pt->second.size();
    }


    /// Return vector of geometric objects involved in
    /// id-th update function
    Vector<GeomObject*>& vector_geom_object_pt(const int& id)
    {
      return Geom_object_pt[id];
    }


    /// Return vector of geometric objects involved in
    /// default (usually first) update function
    Vector<GeomObject*>& vector_geom_object_pt()
    {
      return Default_it_geom_object_pt->second;
    }


    /// Return the vector of all geometric objects
    GeomObject** all_geom_object_pt()
    {
      if (this->ngeom_object() == 0)
      {
        return 0;
      }
      else
      {
        return &(Default_it_geom_object_pt->second[0]);
      }
    }

    /// Return pointer to i-th geometric object involved in
    /// default (usually first) update function
    GeomObject* geom_object_pt(const unsigned& i)
    {
      return Default_it_geom_object_pt->second[i];
    }

    /// Number of reference values involved in id-th update function
    unsigned nref_value(const int& id)
    {
      return Ref_value[id].size();
    }


    /// Number of reference values involved in default
    /// (usually first) update function
    unsigned nref_value()
    {
      return Default_it_ref_value->second.size();
    }


    /// Return vector of reference values involved in
    /// default (usually first) update function
    Vector<double>& vector_ref_value()
    {
      return Default_it_ref_value->second;
    }


    /// Return vector of reference values involved in
    /// id-th update function
    Vector<double>& vector_ref_value(const int& id)
    {
      return Ref_value[id];
    }


    /// Return i-th reference value involved in
    /// default (usually first) update function
    double ref_value(const unsigned& i)
    {
      return Default_it_ref_value->second[i];
    }

    /// Add algebraic update information for node: What's the
    /// ID of the mesh update function (typically used within the mesh)
    /// Which Mesh implements the update operation? Also,
    /// pass the vector of geometric objects and
    /// the vectors of reference values that are
    /// needed for the update operation. Negative values for ID are only
    /// allowed when called from node constructor, as indicated
    /// by the final argument which defaults to false.
    void add_node_update_info(const int& id,
                              AlgebraicMesh* mesh_pt,
                              const Vector<GeomObject*>& geom_object_pt,
                              const Vector<double>& ref_value,
                              const bool& called_from_constructor = false)
    {
      // Sanity check
      if (id < 0)
      {
        if (!called_from_constructor)
        {
          std::ostringstream error_message;
          error_message << "\nNegative ID, " << id
                        << ", only allowed if called from constructor and\n"
                        << "indicated as such by optional boolean flag."
                        << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // If there's just one entry -- check if it's the default dummy one
      if (Mesh_pt.size() == 1)
      {
        if (Mesh_pt.begin()->second == Dummy_mesh_pt)
        {
          if (Default_it_mesh_pt->second == Dummy_mesh_pt)
          {
            kill_node_update_info(Dummy_node_update_fct_id);
          }
        }
      }

      // Now insert the actual info
      Mesh_pt.insert(std::make_pair(id, mesh_pt));
      Geom_object_pt.insert(std::make_pair(id, geom_object_pt));
      Ref_value.insert(std::make_pair(id, ref_value));

      // Always use the "first" update fct as default -- can be overwritten
      // outside (usually only done for self test)
      set_default_node_update(Mesh_pt.begin()->first);
    }


    /// Add algebraic update information for node:
    /// Which Mesh implements the update operation? Also,
    /// pass the vector of geometric objects and
    /// the vectors of reference values that are
    /// needed for the update operation. We're assigning a default
    /// node update fct id of 0.
    void add_node_update_info(AlgebraicMesh* mesh_pt,
                              const Vector<GeomObject*>& geom_object_pt,
                              const Vector<double>& ref_value)
    {
      // No update fct id supplied: Use a default assignment of 0.
      unsigned id = 0;

      // If there's just one entry -- check if it's the default dummy one
      if (Mesh_pt.size() == 1)
      {
        // Do we still have dummy default assignment stored as the one
        // and only entry?
        if (Mesh_pt.begin()->second == Dummy_mesh_pt)
        {
          if (Default_it_mesh_pt->second == Dummy_mesh_pt)
          {
            kill_node_update_info(Dummy_node_update_fct_id);
          }
        }
      }

      // Now insert the actual info
      Mesh_pt.insert(std::make_pair(id, mesh_pt));
      Geom_object_pt.insert(std::make_pair(id, geom_object_pt));
      Ref_value.insert(std::make_pair(id, ref_value));

      // Always use the "first" update fct as default -- can be overwritten
      // outside (usually only done for self test)
      set_default_node_update(Mesh_pt.begin()->first);
    }


    /// Erase algebraic node update information for id-th
    /// node update function. Id defaults to 0.
    void kill_node_update_info(const int& id = 0)
    {
      Mesh_pt.erase(Mesh_pt.find(id));
      Geom_object_pt.erase(Geom_object_pt.find(id));
      Ref_value.erase(Ref_value.find(id));
    }


    /// Perform self test: If the node has multiple node
    /// update functions, check that they all give the same result.
    /// Return 1/0 for failure/success. (Failure if
    /// max. difference between the nodal positions for different
    /// update functions exceeds
    /// AlgebraicNode::Max_allowed_difference_between_node_update_fcts
    unsigned self_test();


  private:
    /// Make id-th node update function the default
    void set_default_node_update(const int& id)
    {
      // Set default node update fct id
      Default_node_update_fct_id = id;


      // Set iterators for default entry

      // Iterator to default mesh:
      Default_it_mesh_pt = Mesh_pt.find(id);
#ifdef PARANOID
      if (Default_it_mesh_pt == Mesh_pt.end())
      {
        std::ostringstream error_message;
        error_message << "There is no reference mesh for node update fct id"
                      << id << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Iterator to default GeomObject vector
      Default_it_geom_object_pt = Geom_object_pt.find(id);
#ifdef PARANOID
      if (Default_it_geom_object_pt == Geom_object_pt.end())
      {
        std::ostringstream error_message;
        error_message << "There is no Geom_object_pt for node update fct id"
                      << id << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Iterator to default values vector
      Default_it_ref_value = Ref_value.find(id);
#ifdef PARANOID
      if (Default_it_ref_value == Ref_value.end())
      {
        std::ostringstream error_message;
        error_message << "There is no Ref_value for node update fct id" << id
                      << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Pointer to mesh that performs the specified node update operation
    /// (Map because this node may only use the Mesh's 116th node update fct.
    /// There's no point in wasting an entire vector for the non-existing
    /// entries)
    std::map<int, AlgebraicMesh*> Mesh_pt;

    /// Vector of geometric objects that are involved
    /// in the specified node update operation.
    /// (Map because this node may only use the Mesh's 116th node update fct.
    /// There's no point in wasting an entire vector for the non-existing
    /// entries)
    std::map<int, Vector<GeomObject*>> Geom_object_pt;

    /// Vector of reference values that are required
    /// for the specified node update operation.
    /// (Map because this node may only use the Mesh's 116th node update fct.
    /// There's no point in wasting an entire vector for the non-existing
    /// entries)
    std::map<int, Vector<double>> Ref_value;

    /// Default iterator for mesh: This mesh performs the default update
    std::map<int, AlgebraicMesh*>::iterator Default_it_mesh_pt;

    /// Default iterator for vector of geom objects. These
    /// GeomObjects are involved in the default update.
    std::map<int, Vector<GeomObject*>>::iterator Default_it_geom_object_pt;

    /// Default iterator for vector of ref values. These
    /// reference values are involved in the default update.
    std::map<int, Vector<double>>::iterator Default_it_ref_value;

    /// Default node update function ID.
    int Default_node_update_fct_id;

    /// What it says: Used in self-test to check if different
    /// node update functions produce the same result.
    static double Max_allowed_difference_between_node_update_fcts;

    /// Default (negative!) remesh fct id for nodes for which no remesh
    /// fct is defined
    static int Dummy_node_update_fct_id;

    /// Default dummy mesh to point to for nodes for which no remesh
    ///  fct is defined
    static AlgebraicMesh* Dummy_mesh_pt;

    /// Static Dummy mesh to which the pointer is addressed
    static DummyAlgebraicMesh Dummy_mesh;

    /// Default dummy vector of geom objects to point to for nodes
    /// for which no remesh fct is defined
    static Vector<GeomObject*> Dummy_geom_object_pt;

    /// Default dummy vector of reference values
    ///  to point to for nodes  for which no remesh fct is defined
    static Vector<double> Dummy_ref_value;
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic elements
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Base class for algebraic elements.
  ///
  //========================================================================
  class AlgebraicElementBase
  {
  public:
    /// Empty constructor
    AlgebraicElementBase() {}

    /// Broken copy constructor
    AlgebraicElementBase(const AlgebraicElementBase&) = delete;

    /// Broken assignment operator
    void operator=(const AlgebraicElementBase&) = delete;

    /// Set up node update info for (newly created) algebraic node:
    /// I.e. work out its node update information by interpolation from
    /// the father element. Pass pointer to father element and the
    /// newly created node's local coordinate in the father element.
    void setup_algebraic_node_update(Node*& node_pt,
                                     const Vector<double>& s_father,
                                     FiniteElement* father_el_pt) const;
  };


  //========================================================================
  /// Algebraic elements are elements that have AlgebraicNodes whose
  /// position is determined by the geometric Data in the GeomObjects
  /// that are involved in their node update functions.
  /// Algebraic Elements include the derivatives w.r.t. any unknowns
  /// that are stored in this geometric Data into the element's
  /// Jacobian matrix. Otherwise they behave exactly like the templace
  /// element.
  //========================================================================
  template<class ELEMENT>
  class AlgebraicElement
    : public ElementWithSpecificMovingNodes<ELEMENT, AlgebraicNode>,
      public AlgebraicElementBase
  {
  public:
    /// Constructor -- simply calls the constructor of the
    /// underlying ELEMENT.
    AlgebraicElement()
      : ElementWithSpecificMovingNodes<ELEMENT, AlgebraicNode>(),
        AlgebraicElementBase()
    {
    }

    /// Constructor for face elements
    AlgebraicElement(FiniteElement* const& element_pt, const int& face_index)
      : ElementWithSpecificMovingNodes<ELEMENT, AlgebraicNode>(element_pt,
                                                               face_index),
        AlgebraicElementBase()
    {
    }

    /// Broken copy constructor
    AlgebraicElement(const AlgebraicElement&) = delete;

    /// Broken assignment operator
    /*void operator=(const AlgebraicElement&) = delete;*/


    /// Empty Destructor must clean up the allocated memory
    ~AlgebraicElement() {}
  };


  //=======================================================================
  /// Explicit definition of the face geometry of algebraic elements:
  /// the same as the face geometry of the underlying element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<AlgebraicElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<ELEMENT>() {}

  protected:
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Algebraic meshes
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Algebraic meshes contain AlgebraicElements and AlgebraicNodes.
  /// They implement the node update functions that are used
  /// by the AlgebraicNodes to update their positions.
  //========================================================================
  class AlgebraicMesh : public virtual Mesh
  {
  public:
    /// Constructor: create a null zeroth entry in the Geom_object_list_pt
    /// Vector (each AlgebraicMesh's constructor should add any other
    /// geometric objects to this list)
    AlgebraicMesh()
    {
      add_geom_object_list_pt(0);
    }

    /// Broken copy constructor
    AlgebraicMesh(const AlgebraicMesh&) = delete;

    /// Broken assignment operator
    /*void operator=(const AlgebraicMesh&) = delete;*/

    /// Surely a proper destructor is required... ?
    ~AlgebraicMesh() {}

    /// Return a pointer to the n-th global AlgebraicNode
    // Can safely cast the nodes to AlgebraicNodes
    AlgebraicNode* node_pt(const unsigned long& n)
    {
#ifdef PARANOID
      if (!dynamic_cast<AlgebraicNode*>(Node_pt[n]))
      {
        std::ostringstream error_stream;
        error_stream << "Error: Node " << n << "is a "
                     << typeid(Node_pt[n]).name() << ", not an AlgebraicNode"
                     << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Return a cast to the Node_pt
      return (dynamic_cast<AlgebraicNode*>(Node_pt[n]));
    }


    /// Update the nodal position posn at time level t (t=0: present;
    /// t>0: previous). Must be implemented for every specific algebraic mesh.
    virtual void algebraic_node_update(const unsigned& t,
                                       AlgebraicNode*& node_pt) = 0;

    /// Update the node update info for given node, following
    /// mesh adaptation. Must be implemented for every specific algebraic
    /// mesh, though it may, of course, be left empty.
    virtual void update_node_update(AlgebraicNode*& node_pt) = 0;


    /// Update all nodal positions via algebraic node update functions
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

        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Initial loop over ALL nodes to setup (need to place at least
      // all master nodes before we can update the position of the
      // hanging ones)
      AlgebraicNode* alg_nod_pt = 0;
      unsigned n_node = nnode();

      // In parallel there may be no nodes on a particular process
      if (n_node > 0)
      {
        for (unsigned n = 0; n < n_node; n++)
        {
          alg_nod_pt = static_cast<AlgebraicNode*>(node_pt(n));
          alg_nod_pt->node_update();
        }

        // Figure out spatial dimension of node
        unsigned n_dim = alg_nod_pt->ndim();

        // Now loop over hanging nodes and adjust their nodal positions
        // to reflect the hanging node constraints
        for (unsigned n = 0; n < n_node; n++)
        {
          Node* nod_pt = node_pt(n);
          if (nod_pt->is_hanging())
          {
            // Initialise
            Vector<double> x(n_dim);
            for (unsigned i = 0; i < n_dim; i++)
            {
              x[i] = 0.0;
            }

            // Loop over master nodes
            unsigned nmaster = nod_pt->hanging_pt()->nmaster();
            for (unsigned imaster = 0; imaster < nmaster; imaster++)
            {
              // Loop over directions
              for (unsigned i = 0; i < n_dim; i++)
              {
                x[i] += nod_pt->hanging_pt()->master_node_pt(imaster)->x(i) *
                        nod_pt->hanging_pt()->master_weight(imaster);
              }
            }

            // Copy across
            for (unsigned i = 0; i < n_dim; i++)
            {
              nod_pt->x(i) = x[i];
            }
            nod_pt->perform_auxiliary_node_update_fct();
          }
        }
      } // end if (n_node>0)

#ifdef OOMPH_HAS_MPI
      // Update positions for external halo nodes attached to this mesh
      // Loop over processors
      for (std::map<unsigned, Vector<Node*>>::iterator it =
             External_halo_node_pt.begin();
           it != External_halo_node_pt.end();
           it++)
      {
        int iproc = (*it).first;
        AlgebraicNode* alg_nod_pt = 0;
        unsigned n_ext_halo_node = nexternal_halo_node(iproc);
        // Only act if there are any external halo nodes
        if (n_ext_halo_node > 0)
        {
          for (unsigned n = 0; n < n_ext_halo_node; n++)
          {
            alg_nod_pt =
              static_cast<AlgebraicNode*>(external_halo_node_pt(iproc, n));
            alg_nod_pt->node_update();
          }

          // Figure out spatial dimension of node
          unsigned n_dim = alg_nod_pt->ndim();

          // Now loop over hanging nodes and adjust their nodal positions
          // to reflect the hanging node constraints
          for (unsigned n = 0; n < n_ext_halo_node; n++)
          {
            Node* nod_pt = external_halo_node_pt(iproc, n);
            if (nod_pt->is_hanging())
            {
              // Initialise
              Vector<double> x(n_dim);
              for (unsigned i = 0; i < n_dim; i++)
              {
                x[i] = 0.0;
              }

              // Loop over master nodes
              unsigned nmaster = nod_pt->hanging_pt()->nmaster();
              for (unsigned imaster = 0; imaster < nmaster; imaster++)
              {
                // Loop over directions
                for (unsigned i = 0; i < n_dim; i++)
                {
                  x[i] += nod_pt->hanging_pt()->master_node_pt(imaster)->x(i) *
                          nod_pt->hanging_pt()->master_weight(imaster);
                }
              }

              // Copy across
              for (unsigned i = 0; i < n_dim; i++)
              {
                nod_pt->x(i) = x[i];
              }
            }
          }
        }

      } // end loop over processors
#endif
    }

    /// Self test: check consistentency of multiple node updates.
    unsigned self_test()
    {
      // Initialise
      bool passed = true;

      unsigned test = Mesh::self_test();
      if (test != 0)
      {
        passed = false;
      }

      // Loop over nodes
      unsigned n_node = nnode();
      for (unsigned n = 0; n < n_node; n++)
      {
        if (static_cast<AlgebraicNode*>(node_pt(n))->self_test() != 0)
        {
          passed = false;
        }
      }

      oomph_info << "Done algnode selftest in mesh" << std::endl;

      // Return verdict
      if (passed)
      {
        return 0;
      }
      else
      {
        return 1;
      }
    }

    /// Add the specified GeomObject to the list of geometric objects
    /// associated with this AlgebraicMesh; remembering that the zeroth entry
    /// is null (set in the constructor above)
    void add_geom_object_list_pt(GeomObject* geom_object_pt)
    {
      Geom_object_list_pt.push_back(geom_object_pt);
    }

    /// Return number of geometric objects associated with AlgebraicMesh
    unsigned ngeom_object_list_pt()
    {
      return Geom_object_list_pt.size();
    }

    /// Access function to the ith GeomObject
    GeomObject* geom_object_list_pt(const unsigned& i)
    {
      // Probably should be a range check in here...
      return Geom_object_list_pt[i];
    }

  private:
    /// Vector of GeomObjects associated with this AlgebraicMesh
    /// The zeroth entry is null, proper entries from the 1st index onwards...
    Vector<GeomObject*> Geom_object_list_pt;
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  // Dummy algebraic mesh
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Dummy algebraic mesh -- used for default assignements
  //========================================================================
  class DummyAlgebraicMesh : public virtual AlgebraicMesh
  {
  public:
    /// Empty constructor
    DummyAlgebraicMesh() {}

    /// Broken copy constructor
    DummyAlgebraicMesh(const DummyAlgebraicMesh&) = delete;

    /// Broken assignment operator
    /*void operator=(const DummyAlgebraicMesh&) = delete;*/

    /// Update the nodal position posn at time level t (t=0: present;
    /// t>0: previous). Do nothing
    virtual void algebraic_node_update(const unsigned& t,
                                       AlgebraicNode*& node_pt)
    {
    }


    /// Update the node update info for given node, following
    /// mesh adaptation. Must be implemented for every specific algebraic
    /// mesh, though it may, of course, be left empty which is exactly
    /// what we do here
    virtual void update_node_update(AlgebraicNode*& node_pt) {}

    /// Setup algebraic node update for specified node;
    /// do nothing in this dummy version
    virtual void setup_algebraic_node_update(AlgebraicNode*& nod_pt) {}
  };


} // namespace oomph

#endif
