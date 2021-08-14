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
// Header file for spine nodes, elements and meshes

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_SPINES_HEADER
#define OOMPH_SPINES_HEADER


#include <string>

// oomph-lib headers
#include "nodes.h"
#include "elements.h"
#include "mesh.h"
#include "geom_objects.h"
#include "element_with_moving_nodes.h"

namespace oomph
{
  //=================================================================
  /// Spines are used for algebraic node update operations in free-surface
  /// fluid problems: They form the back-bones along which nodes in a
  /// a free-surface mesh are located. Typically, the free surface is
  /// located at the "end" of the spine; the nodes in the interior
  /// of the mesh are located at fixed fractions along the spine. The
  /// key Data member of the Spine object is its "height" --
  /// usually an unknown in the problem -- which is used by the
  /// SpineNode's node update function to update the SpineNode's position.
  ///
  /// In more complex problems (such as the case where a fluid layer is
  /// deposited on an elastic body), the node update function can depend on
  /// additional information, such as the GeomObject representation of the
  /// elastic body, and additional Data objects wich
  /// specify the position on the GeomObject from which the Spine
  /// emanates. The Spine class therefore provides storage for
  /// pointers to GeomObjects and storage for any additional geometric
  /// Data that may be required during node update operations.
  //=================================================================
  class Spine
  {
  public:
    /// \short Default constructor: Create the Spine and initialise its
    /// height to zero.
    Spine()
    {
      Geom_data_pt.resize(1);
      // Create Data for height. By default it's free
      Geom_data_pt[0] = new Data(1);
    }

    /// \short Constructor: Create the Spine and initialise its
    /// height to the specified value
    Spine(const double& height)
    {
      Geom_data_pt.resize(1);
      // Create Data for height. By default it's free
      Geom_data_pt[0] = new Data(1);
      // Set value
      Geom_data_pt[0]->set_value(0, height);
    }

    /// \short Constructor: Create the Spine and initialise its
    /// height to the specified value. Store the vector of (pointers to)
    /// the additional geometric Data that is required during
    /// the node update operation for this Spine.
    Spine(const double& height, const Vector<Data*>& geom_data_pt)
    {
      // Find the number of geometric data passed
      const unsigned n_geom_data = geom_data_pt.size();
      // Now allocate enough storage for the spine height and additional
      // geometric data
      Geom_data_pt.resize(n_geom_data + 1);

      // Create Data for height. By default it's free
      Geom_data_pt[0] = new Data(1);
      // Set value
      Geom_data_pt[0]->set_value(0, height);
      // Add the additional geometric data
      for (unsigned i = 0; i < n_geom_data; i++)
      {
        Geom_data_pt[i + 1] = geom_data_pt[i];
      }
    }

    /// \short Constructor: Create the Spine and initialise its
    /// height to the specified value. Store the vector of (pointers to)
    /// the additional geometric Data that is required during
    /// the node update operation; also store vector of (pointers to)
    /// GeomObjects that is required during the node update operation
    /// for this Spine
    Spine(const double& height,
          const Vector<Data*>& geom_data_pt,
          const Vector<GeomObject*>& geom_object_pt)
      : Geom_object_pt(geom_object_pt)
    {
      // Find the number of geometric data passed
      const unsigned n_geom_data = geom_data_pt.size();
      // Now allocate enough storage for the spine height and additional
      // geometric data
      Geom_data_pt.resize(n_geom_data + 1);

      // Create Data for height. By default it's free
      Geom_data_pt[0] = new Data(1);
      // Set value
      Geom_data_pt[0]->set_value(0, height);
      // Add the additional geometric data
      for (unsigned i = 0; i < n_geom_data; i++)
      {
        Geom_data_pt[i + 1] = geom_data_pt[i];
      }
    }


    /// \short Destructor: Wipe Data object that stores the
    /// Spine height. All other objects (geometric Data and
    /// geometric objects) were created outside the Spine
    /// and must be deleted there.
    ~Spine()
    {
      // Kill spine height
      delete Geom_data_pt[0];
    }


    /// Access function to spine height
    double& height()
    {
      return *(Geom_data_pt[0]->value_pt(0));
    }

    /// Access function to Data object that stores the spine height
    Data*& spine_height_pt()
    {
      return Geom_data_pt[0];
    }

    /// \short Access function to Data object that stores the spine height
    /// (const version)
    Data* spine_height_pt() const
    {
      return Geom_data_pt[0];
    }


    /// \short Number of geometric Data that is involved in the
    /// node update operations for this Spine
    unsigned ngeom_data()
    {
      return Geom_data_pt.size();
    }

    /// \short Set vector of (pointers to) geometric Data that is
    /// involved in the node update operations for this Spine.
    /// Wipes any previously existing geometric Data.
    void set_geom_data_pt(const Vector<Data*>& geom_data_pt)
    {
      unsigned n_geom_data = geom_data_pt.size();
      Geom_data_pt.resize(n_geom_data + 1);
      for (unsigned i = 1; i < n_geom_data; i++)
      {
        Geom_data_pt[i + 1] = geom_data_pt[i];
      }
    }

    /// \short Add (pointer to) geometric Data that is
    /// involved in the node update operations for this Spine
    void add_geom_data_pt(Data* geom_data_pt)
    {
      Geom_data_pt.push_back(geom_data_pt);
    }

    /// \short Return i-th geometric Data that is involved in the
    /// node update operations for this Spine
    Data*& geom_data_pt(const unsigned& i)
    {
      return Geom_data_pt[i];
    }

    /// \short Return i-th geometric Data that is involved in the
    /// node update operations for this Spine. Const version
    Data* geom_data_pt(const unsigned& i) const
    {
      return Geom_data_pt[i];
    }

    /// \short Return the vector of geometric data
    Vector<Data*>& vector_geom_data_pt()
    {
      return Geom_data_pt;
    }

    /// \short Number of geometric objects that is involved in the
    /// node update operations for this Spine
    unsigned ngeom_object()
    {
      return Geom_object_pt.size();
    }

    /// \short Set vector of (pointers to) geometric objects that is
    /// involved in the node update operations for this Spine
    void set_geom_object_pt(const Vector<GeomObject*>& geom_object_pt)
    {
      unsigned n_geom_object = geom_object_pt.size();
      Geom_object_pt.resize(n_geom_object);
      for (unsigned i = 0; i < n_geom_object; i++)
      {
        Geom_object_pt[i] = geom_object_pt[i];
      }
    }

    /// \short Add (pointer to) geometric object that is
    /// involved in the node update  operations for this Spine
    void add_geom_object_pt(GeomObject* geom_object_pt)
    {
      Geom_object_pt.push_back(geom_object_pt);
    }

    /// \short Return i-th geometric object that is involved in the
    /// node update operations for this Spine
    GeomObject*& geom_object_pt(const unsigned& i)
    {
      return Geom_object_pt[i];
    }

    /// \short Return i-th geometric object that is involved in the
    /// node update operations for this Spine. Const version
    GeomObject* geom_object_pt(const unsigned& i) const
    {
      return Geom_object_pt[i];
    }

    /// \short Return the vector of all geometric objects that affect this
    /// spine
    Vector<GeomObject*>& vector_geom_object_pt()
    {
      return Geom_object_pt;
    }

    /// \short Number of geometric parameters that are involved in the
    /// node update operations for this Spine
    unsigned ngeom_parameter()
    {
      return Geom_parameter.size();
    }

    /// \short Set vector of geometric parameters that are
    /// involved in the node update operations for this Spine.
    /// Wipes any previously existing geometric parameters
    void set_geom_parameter(const Vector<double>& geom_parameter)
    {
      Geom_parameter = geom_parameter;
    }

    /// \short Add geometric parameter
    /// involved in the node update operations for this Spine
    void add_geom_parameter(const double& geom_parameter)
    {
      Geom_parameter.push_back(geom_parameter);
    }

    /// \short Return i-th geometric parameter that is involved in the
    /// node update operations for this Spine
    double& geom_parameter(const unsigned& i)
    {
      return Geom_parameter[i];
    }

    /// \short Return i-th geometric parameter that is involved in the
    /// node update operations for this Spine. Const version
    const double& geom_parameter(const unsigned& i) const
    {
      return Geom_parameter[i];
    }


  private:
    /// Data that stores the spine height
    // Data* Spine_height_pt;

    /// Vector that stores the pointers to additional geometric Data
    Vector<Data*> Geom_data_pt;

    /// \short Vector that stores the pointers to geometric objects that is
    /// involved in the node update operation
    Vector<GeomObject*> Geom_object_pt;

    /// \short Vector that stores doubles that are used in the geometric updates
    Vector<double> Geom_parameter;
  };


  // Forward declaration
  class SpineMesh;


  //=====================================================================
  /// Class for nodes that live on spines. The assumption is that each Node
  /// lies at a fixed fraction on a single spine (although more complex
  /// behaviour could be included by adding more variables to the spine).
  /// In general, more complex node updating should be handled by the classes
  /// implemented for algebraic node updates.
  //=====================================================================
  class SpineNode : public Node
  {
  private:
    /// Private internal data pointer to a spine
    Spine* Spine_pt;

    /// Private double that represents the fixed fraction along the spine
    double Fraction;

    /// \short Pointer to SpineMesh that this node is a part of.
    /// (The mesh implements the node update function(s))
    SpineMesh* Spine_mesh_pt;

    /// ID of node update function (within specific mesh -- useful if there
    /// are multiple node update functions, e.g. in two-layer problems.
    unsigned Node_update_fct_id;


  public:
    /// Steady Constructor, initialise pointers to zero
    SpineNode(const unsigned& n_dim,
              const unsigned& n_position_type,
              const unsigned& initial_nvalue)
      : Node(n_dim, n_position_type, initial_nvalue),
        Spine_pt(0),
        Fraction(0),
        Spine_mesh_pt(0),
        Node_update_fct_id(0)
    {
    }

    /// Unsteady Constructor, initialise pointers to zero
    SpineNode(TimeStepper* const& time_stepper_pt,
              const unsigned& n_dim,
              const unsigned& n_position_type,
              const unsigned& initial_nvalue)
      : Node(time_stepper_pt, n_dim, n_position_type, initial_nvalue),
        Spine_pt(0),
        Fraction(0),
        Spine_mesh_pt(0),
        Node_update_fct_id(0)
    {
    }

    /// Access function to spine
    Spine*& spine_pt()
    {
      return Spine_pt;
    }

    /// Set reference to fraction along spine
    double& fraction()
    {
      return Fraction;
    }

    /// Access function to ID of node update function (within specific mesh)
    unsigned& node_update_fct_id()
    {
      return Node_update_fct_id;
    }

    /// \short Access function to Pointer to SpineMesh that this node is a part
    /// of and which implements the node update function(s)
    SpineMesh*& spine_mesh_pt()
    {
      return Spine_mesh_pt;
    }

    /// Access function to  spine height
    double& h()
    {
      return Spine_pt->height();
    }

    /// Overload thet node update function, call
    /// the update function in the Node's SpineMesh
    void node_update(const bool& update_all_time_levels_for_new_node = false);

    /// \short Return the number of geometric data, zero if no spine.
    unsigned ngeom_data() const
    {
      if (Spine_pt)
      {
        return Spine_pt->ngeom_data();
      }
      else
      {
        return 0;
      }
    }

    /// Return the number of geometric objects, zero if no spine.
    unsigned ngeom_object() const
    {
      if (Spine_pt)
      {
        return Spine_pt->ngeom_object();
      }
      else
      {
        return 0;
      }
    }

    /// Return the vector of all geometric data
    Data** all_geom_data_pt()
    {
      return &(Spine_pt->geom_data_pt(0));
    }

    /// Return the vector of all geometric objects
    GeomObject** all_geom_object_pt()
    {
      return &(Spine_pt->geom_object_pt(0));
    }
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// \short A policy class that serves only to establish the interface for
  /// assigning the spine equation numbers
  //=======================================================================
  class SpineFiniteElement
  {
  public:
    /// Empty constructor
    SpineFiniteElement() {}

    /// Emtpty virtual destructor
    virtual ~SpineFiniteElement() {}
  };


  //========================================================================
  /// \short The SpineElement<ELEMENT> class takes an existing element as a
  /// template parameter and adds the necessary additional functionality to
  /// allow the element to be update using the Method of Spines.
  /// A vector of pointers to spines and storage for the local equation
  /// numbers associated with the spines are added to the element.
  //========================================================================
  template<class ELEMENT>
  class SpineElement
    : public ElementWithSpecificMovingNodes<ELEMENT, SpineNode>,
      public SpineFiniteElement
  {
  private:
    /// \short Array to hold the index of the geometric data associated with
    /// the spine height of the spine that affects the n-th node
    unsigned* Spine_geometric_index;

    /// \short Complete the setup of additional dependencies. Overloads
    /// empty virtual function in GeneralisedElement to determine the "geometric
    /// Data", i.e. the Data that affects the element's shape.
    /// This function is called (for all elements) at the very beginning of the
    /// equation numbering procedure to ensure that all dependencies
    /// are accounted for.
    void complete_setup_of_dependencies();


  public:
    /// Constructor, call the constructor of the base element
    SpineElement()
      : ElementWithSpecificMovingNodes<ELEMENT, SpineNode>(),
        SpineFiniteElement(),
        Spine_geometric_index(0)
    {
    }

    /// Constructor used for spine face elements
    SpineElement(FiniteElement* const& element_pt, const int& face_index)
      : ElementWithSpecificMovingNodes<ELEMENT, SpineNode>(element_pt,
                                                           face_index),
        SpineFiniteElement(),
        Spine_geometric_index(0)
    {
    }

    /// Destructor, clean up the storage allocated to the local equation numbers
    ~SpineElement()
    {
      if (Spine_geometric_index)
      {
        delete[] Spine_geometric_index;
      }
    }

    /// \short Return the local equation number corresponding to the height
    /// of the spine at the n-th node
    inline int spine_local_eqn(const unsigned& n)
    {
#ifdef RANGE_CHECKING
      const unsigned n_node = this->nnode();
      if (n >= n_node)
      {
        std::ostringstream error_message;
        error_message << "Range Error:  Node number " << n
                      << " is not in the range (0," << n_node - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

#ifdef PARANOID
      // If there is no spine then you can't get the local equation
      if (Spine_geometric_index[n] == this->ngeom_data())
      {
        std::ostringstream error_stream;
        error_stream << "SpineNode " << n
                     << " does not have a Spine attached,\n"
                     << "so you can't get its local equation number.\n"
                     << "Check that the Mesh is correctly associating Spines "
                        "with is Nodes\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->geometric_data_local_eqn(Spine_geometric_index[n], 0);
    }
  };

  //=======================================================================
  /// \short Explicit definition of the face geometry for spine elements:
  /// The same as the face geometry of the underlying element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<SpineElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<ELEMENT>() {}

  protected:
  };

  //=====================================================================
  /// \short Explicit definition of the face geometry for spine elements:
  /// The same as the face geometry of the underlying element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<FaceGeometry<SpineElement<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}

  protected:
  };

  //=====================================================================
  /// \short Explicit definition of the face geometry for spine elements:
  /// The same as the face geometry of the underlying element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<SpineElement<FaceGeometry<ELEMENT>>>
    : public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    /// Constructor
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}

  protected:
  };


  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////


  //========================================================================
  /// General SpineMesh class.
  ///
  /// Derived from Mesh with virtual so that
  /// spine meshes can be derived from general meshes, without
  /// multiple copies of Mesh objects.
  //========================================================================
  class SpineMesh : public virtual Mesh
  {
  protected:
    /// A Spine mesh contains a Vector of pointers to spines
    Vector<Spine*> Spine_pt;

  public:
    /// Destructor to clean up the memory allocated to the spines
    virtual ~SpineMesh();

    /// Return the i-th spine in the mesh
    Spine*& spine_pt(const unsigned long& i)
    {
      return Spine_pt[i];
    }

    /// Return the i-th spine in the mesh (const version)
    const Spine* spine_pt(const unsigned long& i) const
    {
      return Spine_pt[i];
    }

    /// Return the number of spines in the mesh
    unsigned long nspine() const
    {
      return Spine_pt.size();
    }

    /// Add a spine to the mesh
    void add_spine_pt(Spine* const& spine_pt)
    {
      Spine_pt.push_back(spine_pt);
    }

    /// Return a pointer to the n-th global SpineNode
    // Can safely cast the nodes to SpineNodes
    SpineNode* node_pt(const unsigned long& n)
    {
#ifdef PARANOID
      if (!dynamic_cast<SpineNode*>(Node_pt[n]))
      {
        std::ostringstream error_message;
        error_message << "Node " << n << "is a " << typeid(Node_pt[n]).name()
                      << ", not a SpineNode" << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Return a cast to the pointer to the node
      return (dynamic_cast<SpineNode*>(Node_pt[n]));
    }

    /// \short Return the n-th local SpineNode in element e.
    /// This is required to cast the nodes in a spine mesh to be
    /// SpineNodes and therefore allow access to the extra SpineNode data
    SpineNode* element_node_pt(const unsigned long& e, const unsigned& n)
    {
#ifdef PARANOID
      // Try to cast to FiniteElement
      FiniteElement* el_pt = dynamic_cast<FiniteElement*>(Element_pt[e]);
      if (el_pt == 0)
      {
        throw OomphLibError(
          "Can't execute element_node_pt(...) for non FiniteElements",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      if (!dynamic_cast<SpineNode*>(el_pt->node_pt(n)))
      {
        std::ostringstream error_message;
        error_message << "Node " << n << "is a "
                      << typeid(el_pt->node_pt(n)).name() << ", not a SpineNode"
                      << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Return a cast to the node pointer
      return (dynamic_cast<SpineNode*>(
        dynamic_cast<FiniteElement*>(Element_pt[e])->node_pt(n)));
    }

    /// Assign spines to Spine_pt vector of element
    // N.B.: Since SpineElement<ELEMENT>'s are templated, we need the template
    // in this function so that we can do the dynamic cast to a
    // SpineElement<ELEMENT> template<class ELEMENT> void
    // add_spine_to_element(const unsigned long &e, Spine* spine)
    // {dynamic_cast<ELEMENT*>(Element_pt[e])->add_spine(spine);}

    /// Assign equation numbers for spines
    unsigned long assign_global_spine_eqn_numbers(Vector<double*>& Dof_pt);

    /// \short Function to describe the dofs of the Spine. The ostream
    /// specifies the output stream to which the description
    /// is written; the string stores the currently
    /// assembled output that is ultimately written to the
    /// output stream by Data::describe_dofs(...); it is typically
    /// built up incrementally as we descend through the
    /// call hierarchy of this function when called from
    /// Problem::describe_dofs(...)
    void describe_spine_dofs(std::ostream& out,
                             const std::string& current_string) const;

    /// \short Overload the mesh_level timestepper function to set the
    /// timestepper data for the spines
    void set_mesh_level_time_stepper(TimeStepper* const& time_stepper_pt,
                                     const bool& preserve_existing_data)
    {
      this->set_spine_time_stepper(time_stepper_pt, preserve_existing_data);
    }

    /// \short Set the time stepper forthe spine data that is stored in
    /// the mesh.
    void set_spine_time_stepper(TimeStepper* const& time_stepper_pt,
                                const bool& preserve_existing_data);

    /// \short Set any pinned spine "history" values to be consistent for
    /// continuation problems
    void set_consistent_pinned_spine_values_for_continuation(
      ContinuationStorageScheme* const& continuation_stepper_pt);


    /// \short Check whether the pointer parameter_pt addresses data stored
    /// in the spines
    bool does_pointer_correspond_to_spine_data(double* const& parameter_pt);

    /// \short Update function to update all nodes of mesh
    /// [Doesn't make sense to use this mesh with SolidElements anyway,
    /// so we buffer the case if update_all_solid_nodes is set to
    /// true.]
    void node_update(const bool& update_all_solid_nodes = false);

    /// \short Update function for given spine node -- this must be implemented
    /// by all specific SpineMeshes.
    virtual void spine_node_update(SpineNode* spine_node_pt) = 0;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif

    /// \short Overload the dump function so that the spine data is dumped
    void dump(std::ofstream& dump_file) const;

#ifdef __clang__
#pragma clang diagnostic pop
#endif

    /// \short Overload the read function so that the spine data is read
    /// from the restart file
    void read(std::ifstream& restart_file);
  };


  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////
  // Functions for the SpineElement class
  ///////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////


  //=================================================================
  /// Construct and fill the node_update_data vector
  //=================================================================
  template<class ELEMENT>
  void SpineElement<ELEMENT>::complete_setup_of_dependencies()
  {
    // Call function of underlying element
    ElementWithSpecificMovingNodes<ELEMENT,
                                   SpineNode>::complete_setup_of_dependencies();

    // Sort out the spine index stuff
    // Find the data that correspond to spine heights
    {
      const unsigned n_node = this->nnode();
      // Allocate memory
      if (Spine_geometric_index)
      {
        delete[] Spine_geometric_index;
      }
      Spine_geometric_index = new unsigned[n_node];

      // Now loop over data and find out where it fits
      for (unsigned n = 0; n < n_node; n++)
      {
        // Find pointer to the spine
        Spine* const spine_pt =
          static_cast<SpineNode*>(this->node_pt(n))->spine_pt();

        // If there is a spine then find the pointer to the data
        if (spine_pt)
        {
          // Find the pointer to the data
          Data* spine_height_data_pt = spine_pt->spine_height_pt();

          // Now find the index of the corresponding spine
          const unsigned n_node_update_data = this->ngeom_data();
          for (unsigned i = 0; i < n_node_update_data; i++)
          {
            if (this->Geom_data_pt[i] == spine_height_data_pt)
            {
              Spine_geometric_index[n] = i;
              break;
            }
          }
        }
        // Otherwise issue a warning
        else
        {
          // Set the spine_geometric_index out of range,
          // which will cause the spine_local_eqn to return a pinned value
          Spine_geometric_index[n] = this->ngeom_data();
        }
      }
    }
  }

} // namespace oomph

#endif
