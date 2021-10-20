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

// Include guard to prevent multiple inclusions of the header
#ifndef OOMPH_FACE_MESH_PROJECT_HEADER
#define OOMPH_FACE_MESH_PROJECT_HEADER

namespace oomph
{
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////


  //======================================================================
  ///  Class that makes the finite element specified as template argument
  /// projectable -- on the assumption that all fields are interpolated
  /// by isoparametric Lagrange interpolation between the nodes.
  //======================================================================
  template<class ELEMENT>
  class GenericLagrangeInterpolatedProjectableElement
    : public virtual ProjectableElement<ELEMENT>
  {
  public:
    /// Constructor
    GenericLagrangeInterpolatedProjectableElement()
    {
      Boundary_id = UINT_MAX;
    }

    ///  Nodal value of boundary coordinate
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      // Vector in which to hold the intrinsic coordinate
      Vector<double> zeta(this->dim());

#ifdef PARANOID
      if (Boundary_id == UINT_MAX)
      {
        std::ostringstream error_message;
        error_message << "Boundary_id is (still) UINT_MAX -- please set\n"
                      << "the actual value with set_boundary_id(...)\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Get the k-th generalised boundary coordinate at node n
      this->node_pt(n)->get_coordinates_on_boundary(Boundary_id, k, zeta);

      // Return the individual coordinate
      return zeta[i];
    }

    /// Boundary id
    void set_boundary_id(const unsigned& boundary_id)
    {
      Boundary_id = boundary_id;
    }

    /// Boundary id
    unsigned boundary_id() const
    {
#ifdef PARANOID
      if (Boundary_id == UINT_MAX)
      {
        std::ostringstream error_message;
        error_message << "Boundary_id is (still) UINT_MAX -- please set\n"
                      << "the actual value with set_boundary_id(...)\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Boundary_id;
    }

    ///  Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    Vector<std::pair<Data*, unsigned>> data_values_of_field(const unsigned& fld)
    {
      // Create the vector
      unsigned nnod = this->nnode();
      Vector<std::pair<Data*, unsigned>> data_values(nnod);

      // Loop over all nodes
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated field: The node itself
        data_values[j] = std::make_pair(this->node_pt(j), fld);
      }

      // Return the vector
      return data_values;
    }

    ///  Number of fields to be projected.
    unsigned nfields_for_projection()
    {
      return this->node_pt(0)->nvalue();
    }

    ///  Number of history values to be stored for fld-th field
    /// (includes current value!). Extract from first node but assume it's
    /// the same for all.
    unsigned nhistory_values_for_projection(const unsigned& fld)
    {
      return this->node_pt(0)->ntstorage();
    }

    /// Number of positional history values (Note: count includes
    /// current value!). Extract from first node but assume it's
    /// the same for all.
    unsigned nhistory_values_for_coordinate_projection()
    {
      return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
    }


    ///  Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s.
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi)
    {
      this->shape(s, psi);
      return this->J_eulerian(s);
    }


    ///  Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s)
    {
      // Local shape function
      unsigned n_node = this->nnode();
      Shape psi(n_node);

      // Find values of shape function
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Sum over the local nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(t, l, fld) * psi[l];
      }
      return interpolated_u;
    }

    ///  Return number of values in field fld
    unsigned nvalue_of_field(const unsigned& fld)
    {
      return this->nnode();
    }


    ///  Return local equation number of value j in field fld. Assumed to
    /// be the local nodal equation.
    int local_equation(const unsigned& fld, const unsigned& j)
    {
      return this->nodal_local_eqn(j, fld);
    }


  private:
    /// Boundary id
    unsigned Boundary_id;
  };


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  //======================================================================
  /// Class that makes backup (via a deep copy) of a mesh, keeping alive
  /// enough information to allow the solution that is currently stored
  /// on the mesh to be projected onto another mesh sometime in the
  /// future (when the original mesh may already have been deleted).
  /// This is mainly useful for the projection of additional
  /// nodal values (such as Lagrange multipliers) created by FaceElements.
  /// ASSUMPTION: All fields in the element are represented by isoparametric
  ///             Lagrange interpolation between the nodal values.
  ///             Any fields that do not fall into this category will not
  ///             be copied across correctly and if you're unlucky
  ///             the code may die...).
  //======================================================================
  template<class GEOMETRIC_ELEMENT>
  class BackupMeshForProjection : public virtual Mesh
  {
  public:
    ///  Constructor: Pass existing mesh and the boundary ID (need to find
    /// the boundary coordinates that are used for the projection.
    /// Optional final argument specifies the ID of the field (i.e. the
    /// index of the relevant nodal value!) to be projected.
    /// If omitted, we project all of them.
    BackupMeshForProjection(
      Mesh* mesh_pt,
      const unsigned& boundary_id,
      const unsigned& id_of_field_to_be_projected = UINT_MAX)
      : Boundary_id(boundary_id),
        ID_of_field_to_be_projected(id_of_field_to_be_projected)
    {
      // Find unique nodes (via elements because Node_pt vector in
      // original mesh may not have been filled (typical for most
      // face element meshes)
      unsigned nel = mesh_pt->nelement();
      Element_pt.reserve(nel);
      Node_pt.reserve(mesh_pt->nnode());
      for (unsigned e = 0; e < nel; e++)
      {
        FiniteElement* el_pt = mesh_pt->finite_element_pt(e);
        if (el_pt != 0)
        {
          // Make new element
#ifdef PARANOID
          if (dynamic_cast<GEOMETRIC_ELEMENT*>(mesh_pt->element_pt(e)) == 0)
          {
            std::ostringstream error_message;
            error_message << "Element is of wrong type " << typeid(el_pt).name()
                          << " doesn't match template parameter!" << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);

            if (el_pt->ninternal_data() != 0)
            {
              std::ostringstream error_message;
              error_message
                << "Internal data will NOT be projected across!\n"
                << "If you want this functionality you'll have to \n"
                << "implement it yourself" << std::endl;
              OomphLibWarning(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
            }
          }
#endif

          // Make a new element
          GenericLagrangeInterpolatedProjectableElement<GEOMETRIC_ELEMENT>*
            new_el_pt = new GenericLagrangeInterpolatedProjectableElement<
              GEOMETRIC_ELEMENT>;

          // Set boundary ID
          new_el_pt->set_boundary_id(Boundary_id);

          // Set nodal dimension
          unsigned nodal_dim = el_pt->node_pt(0)->ndim();
          new_el_pt->set_nodal_dimension(nodal_dim);

          // Add it to mesh
          add_element_pt(new_el_pt);

          // Create new nodes if needed
          unsigned nnod = el_pt->nnode();
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* old_node_pt = el_pt->node_pt(j);
            if (New_node_pt[old_node_pt] == 0)
            {
              Node* new_nod_pt = 0;


#ifdef PARANOID
              // Check boundary node-ness
              if (!old_node_pt->is_on_boundary())
              {
                std::ostringstream error_message;
                error_message << "Node isn't on a boundary!" << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif


              // How many values are we projecting? Default: All
              unsigned nval = old_node_pt->nvalue();
              unsigned first_index_in_old_node = 0;
              if (ID_of_field_to_be_projected != UINT_MAX)
              {
                nval = dynamic_cast<BoundaryNodeBase*>(old_node_pt)
                         ->nvalue_assigned_by_face_element(
                           ID_of_field_to_be_projected);
                first_index_in_old_node =
                  dynamic_cast<BoundaryNodeBase*>(old_node_pt)
                    ->index_of_first_value_assigned_by_face_element(
                      ID_of_field_to_be_projected);
              }

              // Build new node
              new_nod_pt =
                new BoundaryNode<Node>(old_node_pt->time_stepper_pt(),
                                       old_node_pt->ndim(),
                                       old_node_pt->nposition_type(),
                                       nval);

              // Copy data across
              unsigned n_time = old_node_pt->ntstorage();
              for (unsigned t = 0; t < n_time; t++)
              {
                for (unsigned i = 0; i < nval; i++)
                {
                  new_nod_pt->set_value(
                    t, i, old_node_pt->value(t, first_index_in_old_node + i));
                }
              }

              // Copy nodal positions
              unsigned n_dim = old_node_pt->ndim();
              for (unsigned i = 0; i < n_dim; i++)
              {
                new_nod_pt->x(i) = old_node_pt->x(i);
              }

              // Add to boundary
              new_nod_pt->add_to_boundary(Boundary_id);

              // Transfer boundary coordinates
#ifdef PARANOID
              if (!old_node_pt->is_on_boundary(Boundary_id))
              {
                std::ostringstream error_message;
                error_message << "Boundary ID specified as " << Boundary_id
                              << " but node isn't actually on that boundary!"
                              << std::endl;
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
#endif

              // Get vector of coordinates on mesh boundary from old node
              unsigned n = old_node_pt->ncoordinates_on_boundary(Boundary_id);
              Vector<double> zeta(n);
              old_node_pt->get_coordinates_on_boundary(Boundary_id, zeta);

              // Set for new node
              new_nod_pt->set_coordinates_on_boundary(Boundary_id, zeta);

              // Add node
              Node_pt.push_back(new_nod_pt);

              // Setup association
              New_node_pt[old_node_pt] = new_nod_pt;
            }

            // Set node pointer from new element
            new_el_pt->node_pt(j) = New_node_pt[old_node_pt];
          }
        }
      }
    }

    ///  Project the solution that was present in the original mesh
    /// and from which this backup mesh was created onto the mesh
    /// pointed to by new_mesh_pt. Note that elements in the new mesh do
    /// not have to be projectable. The original mesh may by now have
    /// been deleted.
    void project_onto_new_mesh(Mesh* new_mesh_pt)
    {
      // Make copy of new mesh that we can project onto
      BackupMeshForProjection<GEOMETRIC_ELEMENT>* projectable_new_mesh_pt =
        new BackupMeshForProjection<GEOMETRIC_ELEMENT>(
          new_mesh_pt, Boundary_id, ID_of_field_to_be_projected);

      // Create projection problem
      ProjectionProblem<
        GenericLagrangeInterpolatedProjectableElement<GEOMETRIC_ELEMENT>>*
        proj_problem_pt = new ProjectionProblem<
          GenericLagrangeInterpolatedProjectableElement<GEOMETRIC_ELEMENT>>;

      // Set the mesh we want to project onto
      proj_problem_pt->mesh_pt() = projectable_new_mesh_pt;

      // Project from projectable copy of original mesh -- this one!
      bool dont_project_positions = true;
      proj_problem_pt->project(this, dont_project_positions);

      // Copy nodal values onto the corresponding nodes in the new mesh
      projectable_new_mesh_pt->copy_onto_original_mesh();

      // Kill!
      delete proj_problem_pt;
      delete projectable_new_mesh_pt;
    }


    ///  Copy nodal values back onto original mesh from which this
    /// mesh was built. This obviously only makes sense if the original
    /// mesh still exists!
    void copy_onto_original_mesh()
    {
      for (std::map<Node*, Node*>::iterator it = New_node_pt.begin();
           it != New_node_pt.end();
           it++)
      {
        // Get old node (in the previously existing mesh)
        Node* old_node_pt = (*it).first;

        // ...and corresponding new one (in the mesh where we did the projection
        Node* new_node_pt = (*it).second;

        // How many values are we moving across?
        unsigned nval = old_node_pt->nvalue();
        unsigned first_index_in_old_node = 0;
        if (ID_of_field_to_be_projected != UINT_MAX)
        {
          nval =
            dynamic_cast<BoundaryNodeBase*>(old_node_pt)
              ->nvalue_assigned_by_face_element(ID_of_field_to_be_projected);
          first_index_in_old_node =
            dynamic_cast<BoundaryNodeBase*>(old_node_pt)
              ->index_of_first_value_assigned_by_face_element(
                ID_of_field_to_be_projected);
        }

        // Copy data across (include offset in orig mesh)
        unsigned n_time = old_node_pt->ntstorage();
        for (unsigned t = 0; t < n_time; t++)
        {
          for (unsigned i = 0; i < nval; i++)
          {
            old_node_pt->set_value(
              t, first_index_in_old_node + i, new_node_pt->value(t, i));
          }
        }
      }
    }


  private:
    /// Map returning new node, labeled by node point in original mesh
    std::map<Node*, Node*> New_node_pt;

    /// Boundary id
    unsigned Boundary_id;

    ///  ID of field to be projected (UINT_MAX if there isn't one, in
    /// which case we're doing all of them.
    unsigned ID_of_field_to_be_projected;
  };


  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////

} // namespace oomph

#endif
