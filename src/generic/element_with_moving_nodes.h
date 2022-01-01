// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
// Header file for base class that specified common interfaces
// used in elements with moving Nodes (Spine,Algebraic,MacroElementNodeUpdate)

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_ELEMENT_WITH_MOVING_NODES
#define OOMPH_ELEMENT_WITH_MOVING_NODES


// oomph-lib headers
#include "elements.h"

namespace oomph
{
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// A policy class that serves to establish the
  /// common interfaces for elements that contain moving nodes. This
  /// class provides storage for the geometric data that affect the
  /// update of all the nodes of the element, i.e. USUALLY
  /// all data that are using during a call
  /// to the Element's node_update() function. In some cases
  /// (e.g. FluidInterfaceEdge elements), node_update() is overloaded to
  /// perform an update of the bulk element, in which case the additional
  /// bulk geometric data become external data of the element and the
  /// function GeneralisedElement::update_in_external_fd(i) is overloaded
  /// to also perform the bulk node update.
  /// The storage is populated
  /// during the assignment of the equation numbers via the
  /// complete_setup_of_dependencies() function and then local equations
  /// numbers are assigned to these data, accessible via
  /// geometric_data_local_eqn(n,i). Finally, a function is provided
  /// that calculates the terms in the jacobian matrix by due to these
  /// geometric data by finite differences.
  //=======================================================================
  class ElementWithMovingNodes : public virtual FiniteElement
  {
  public:
    /// Public enumeration to choose method for computing shape derivatives
    enum
    {
      Shape_derivs_by_chain_rule,
      Shape_derivs_by_direct_fd,
      Shape_derivs_by_fastest_method
    };

    /// Constructor
    ElementWithMovingNodes()
      : Geometric_data_local_eqn(0),
        Bypass_fill_in_jacobian_from_geometric_data(false),
        Evaluate_dresidual_dnodal_coordinates_by_fd(false),
        Method_for_shape_derivs(Shape_derivs_by_direct_fd)
    // hierher: Anything other than the fd-based method is currently broken;
    // at least for refineable elements -- this all needs to be checked
    // VERY carefully again (see instructions in commit log). Until this
    // has all been fixed/re-validated, we've frozen the default assignment
    // (by breaking the functions that change the setting). Ultimately,
    // we should use:
    // Method_for_shape_derivs(Shape_derivs_by_fastest_method)
    {
    }

    /// Broken copy constructor
    ElementWithMovingNodes(const ElementWithMovingNodes&) = delete;

    /// Broken assignment operator
    void operator=(const ElementWithMovingNodes&) = delete;

    /// Function to describe the local dofs of the element. The ostream
    /// specifies the output stream to which the description
    /// is written; the string stores the currently
    /// assembled output that is ultimately written to the
    /// output stream by Data::describe_dofs(...); it is typically
    /// built up incrementally as we descend through the
    /// call hierarchy of this function when called from
    /// Problem::describe_dofs(...)
    void describe_local_dofs(std::ostream& out,
                             const std::string& current_string) const;

    /// Virtual destructor (clean up and allocated memory)
    virtual ~ElementWithMovingNodes()
    {
      if (Geometric_data_local_eqn)
      {
        delete[] Geometric_data_local_eqn[0];
        delete[] Geometric_data_local_eqn;
      }
    }

    /// Number of geometric dofs
    unsigned ngeom_dof() const
    {
      return Ngeom_dof;
    }

    /// Return the local equation number corresponding to the i-th
    /// value at the n-th geometric data object.
    inline int geometric_data_local_eqn(const unsigned& n, const unsigned& i)
    {
#ifdef RANGE_CHECKING
      unsigned n_data = Geom_data_pt.size();
      if (n >= n_data)
      {
        std::ostringstream error_message;
        error_message << "Range Error:  Data number " << n
                      << " is not in the range (0," << n_data - 1 << ")";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        unsigned n_value = Geom_data_pt[n]->nvalue();
        if (i >= n_value)
        {
          std::ostringstream error_message;
          error_message << "Range Error: value " << i << " at data " << n
                        << " is not in the range (0," << n_value - 1 << ")";
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }
#endif
#ifdef PARANOID
      // Check that the equations have been allocated
      if (Geometric_data_local_eqn == 0)
      {
        throw OomphLibError(
          "Geometric data local equation numbers have not been allocated",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Geometric_data_local_eqn[n][i];
    }


    /// Return a set of all geometric data associated with the element
    void assemble_set_of_all_geometric_data(
      std::set<Data*>& unique_geom_data_pt);

    /// Specify Data that affects the geometry of the element
    /// by adding the element's geometric Data to the set that's passed in.
    /// (This functionality is required in FSI problems; set is used to
    /// avoid double counting).
    void identify_geometric_data(std::set<Data*>& geometric_data_pt)
    {
      // Loop over the node update data and add to the set
      const unsigned n_geom_data = Geom_data_pt.size();
      for (unsigned n = 0; n < n_geom_data; n++)
      {
        geometric_data_pt.insert(Geom_data_pt[n]);
      }
    }

    /// Return whether shape derivatives are evaluated by fd
    bool are_dresidual_dnodal_coordinates_always_evaluated_by_fd() const
    {
      return Evaluate_dresidual_dnodal_coordinates_by_fd;
    }

    /// Insist that shape derivatives are always
    /// evaluated by fd (using
    /// FiniteElement::get_dresidual_dnodal_coordinates())
    void enable_always_evaluate_dresidual_dnodal_coordinates_by_fd()
    {
      Evaluate_dresidual_dnodal_coordinates_by_fd = true;
    }

    /// Insist that shape derivatives are always
    /// evaluated using the overloaded version of this function
    /// that may have been implemented in a derived class.
    /// (The default behaviour will still be finite differences unless the
    /// function has actually been overloaded
    void disable_always_evaluate_dresidual_dnodal_coordinates_by_fd()
    {
      Evaluate_dresidual_dnodal_coordinates_by_fd = false;
    }

    /// Evaluate shape derivatives by direct finite differencing
    void evaluate_shape_derivs_by_direct_fd()
    {
      // hierher: Harmless; simply re-sets the current (and currently only
      // legal) default
      Method_for_shape_derivs = Shape_derivs_by_direct_fd;
    }

    /// Evaluate shape derivatives by chain rule. Currently disabled by
    /// default because it's broken; can re-enable use by setting optional
    /// boolean to true.
    void evaluate_shape_derivs_by_chain_rule(
      const bool& i_know_what_i_am_doing = false)
    {
      if (!i_know_what_i_am_doing)
      {
        std::ostringstream error_message;
        error_message
          << "Evaluation of shape derivatives by chain rule is currently \n"
          << "disabled because it's broken, at least for refineable \n"
          << "elements. This all needs to be checked again very carefully\n"
          << "following the instructions in the commit log\n"
          << "If you know what you're doing and want to force this "
             "methodology\n"
          << "call this function with the optional boolean set to true.\n";
        throw OomphLibError(
          error_message.str(),
          "ElementWithMovingNodes::evaluate_shape_derivs_by_chain_rule()",
          OOMPH_EXCEPTION_LOCATION);
      }
      Method_for_shape_derivs = Shape_derivs_by_chain_rule;
    }

    /// Evaluate shape derivatives by (anticipated) fastest method.
    /// Currently disabled by default because it's broken; can re-enable
    /// use by setting optional boolean to true.
    void evaluate_shape_derivs_by_fastest_method(
      const bool& i_know_what_i_am_doing = false)
    {
      if (!i_know_what_i_am_doing)
      {
        std::ostringstream error_message;
        error_message
          << "Evaluation of shape derivatives by fastest method is currently \n"
          << "disabled because it's broken, at least for refineable \n"
          << "elements. This all needs to be checked again very carefully\n"
          << "following the instructions in the commit log\n"
          << "If you know what you're doing and want to force this "
             "methodology\n"
          << "call this function with the optional boolean set to true.\n";
        throw OomphLibError(
          error_message.str(),
          "ElementWithMovingNodes::evaluate_shape_derivs_by_fastest_method()",
          OOMPH_EXCEPTION_LOCATION);
      }
      Method_for_shape_derivs = Shape_derivs_by_fastest_method;
    }

    /// Access to method (enumerated flag) for determination of shape derivs
    int& method_for_shape_derivs()
    {
      return Method_for_shape_derivs;
    }

    /// Bypass the call to fill_in_jacobian_from_geometric_data
    void enable_bypass_fill_in_jacobian_from_geometric_data()
    {
      Bypass_fill_in_jacobian_from_geometric_data = true;
    }

    /// Do not bypass the call to fill_in_jacobian_from_geometric_data
    void disable_bypass_fill_in_jacobian_from_geometric_data()
    {
      Bypass_fill_in_jacobian_from_geometric_data = false;
    }

    /// Test whether the call to fill_in_jacobian_from_geometric_data is
    /// bypassed
    bool is_fill_in_jacobian_from_geometric_data_bypassed() const
    {
      return Bypass_fill_in_jacobian_from_geometric_data;
    }

  protected:
    /// Compute derivatives of the nodal coordinates w.r.t.
    /// to the geometric dofs. Default implementation by FD can be overwritten
    /// for specific elements.
    /// dnodal_coordinates_dgeom_dofs(l,i,j) = dX_{ij} / d s_l
    virtual void get_dnodal_coordinates_dgeom_dofs(
      RankThreeTensor<double>& dnodal_coordinates_dgeom_dofs);

    /// Construct the vector of (unique) geometric data
    void complete_setup_of_dependencies();

    /// Assign local equation numbers for the geometric Data in the element
    /// If the boolean argument is true then the degrees of freedom are stored
    /// in Dof_pt
    virtual void assign_all_generic_local_eqn_numbers(
      const bool& store_local_dof_pt);

    /// Calculate the contributions to the Jacobian matrix from the
    /// geometric data. This version
    /// assumes that the (full) residuals vector has already been calculated
    /// and is passed in as the first argument -- needed in case
    /// the derivatives are computed by FD.
    void fill_in_jacobian_from_geometric_data(Vector<double>& residuals,
                                              DenseMatrix<double>& jacobian);

    /// Calculate the contributions to the Jacobian matrix from the
    /// geometric data. This version computes
    /// the residuals vector before calculating the Jacobian terms
    void fill_in_jacobian_from_geometric_data(DenseMatrix<double>& jacobian)
    {
      if (!Bypass_fill_in_jacobian_from_geometric_data)
      {
        // Allocate storage for the residuals
        const unsigned n_dof = ndof();
        Vector<double> residuals(n_dof);

        // Get the residuals for the entire element
        get_residuals(residuals);

        // Call the jacobian calculation
        fill_in_jacobian_from_geometric_data(residuals, jacobian);
      }
    }


    /// Vector that stores pointers to all Data that affect the
    /// node update operations, i.e. the variables that can affect
    /// the position of the node.
    Vector<Data*> Geom_data_pt;

  public:
    /// Return the number of geometric data upon which the shape
    /// of the element depends
    unsigned ngeom_data() const
    {
      return Geom_data_pt.size();
    }

  private:
    /// Array to hold local eqn number information for the
    /// geometric Data variables
    int** Geometric_data_local_eqn;

    /// Number of geometric dofs (computed on the fly when
    /// equation numbers are set up)
    unsigned Ngeom_dof;

    /// Set flag to true to bypass calculation of Jacobain entries
    /// resulting from geometric data.
    bool Bypass_fill_in_jacobian_from_geometric_data;

    /// Boolean to decide if shape derivatives are to be evaluated
    /// by fd (using FiniteElement::get_dresidual_dnodal_coordinates())
    /// or analytically, using the overloaded version of this function
    /// in this class.
    bool Evaluate_dresidual_dnodal_coordinates_by_fd;

    /// Choose method for evaluation of shape derivatives
    /// (this takes one of the values in the enumeration)
    int Method_for_shape_derivs;
  };


  //===============================================================
  /// Specific implementation of the class for specified element
  /// and node type.
  //==============================================================
  template<class ELEMENT, class NODE_TYPE>
  class ElementWithSpecificMovingNodes : public ELEMENT,
                                         public ElementWithMovingNodes
  {
  public:
    /// Function to describe the local dofs of the element. The ostream
    /// specifies the output stream to which the description
    /// is written; the string stores the currently
    /// assembled output that is ultimately written to the
    /// output stream by Data::describe_dofs(...); it is typically
    /// built up incrementally as we descend through the
    /// call hierarchy of this function when called from
    /// Problem::describe_dofs(...)
    void describe_local_dofs(std::ostream& out,
                             const std::string& current_string) const
    {
      ELEMENT::describe_local_dofs(out, current_string);
      ElementWithMovingNodes::describe_local_dofs(out, current_string);
    }

    /// Constructor, call the constructor of the base element
    ElementWithSpecificMovingNodes() : ELEMENT(), ElementWithMovingNodes() {}

    /// Constructor used for face elements
    ElementWithSpecificMovingNodes(FiniteElement* const& element_pt,
                                   const int& face_index)
      : ELEMENT(element_pt, face_index), ElementWithMovingNodes()
    {
    }

    /// Empty Destructor,
    ~ElementWithSpecificMovingNodes() {}

    /// Unique final overrider for describe_dofs
    void describe_local_dofs(std::ostream& out, std::string& curr_str)
    {
      ElementWithMovingNodes::describe_local_dofs(out, curr_str);
      ELEMENT::describe_local_dofs(out, curr_str);
    }


    /// Overload the node assignment routine to assign nodes of the
    /// appropriate type.
    Node* construct_node(const unsigned& n)
    {
      // Assign a node to the local node pointer
      // The dimension and number of values are taken from internal element data
      // The number of timesteps to be stored comes from the problem!
      this->node_pt(n) = new NODE_TYPE(this->nodal_dimension(),
                                       this->nnodal_position_type(),
                                       this->required_nvalue(n));
      // Now return a pointer to the node, so that the mesh can find it
      return this->node_pt(n);
    }

    /// Overloaded node allocation for unsteady problems
    Node* construct_node(const unsigned& n, TimeStepper* const& time_stepper_pt)
    {
      // Assign a node to the local node pointer
      // The dimension and number of values are taken from internal element data
      // The number of timesteps to be stored comes from the problem!
      this->node_pt(n) = new NODE_TYPE(time_stepper_pt,
                                       this->nodal_dimension(),
                                       this->nnodal_position_type(),
                                       this->required_nvalue(n));
      // Now return a pointer to the node, so that the mesh can find it
      return this->node_pt(n);
    }

    /// Overload the node assignment routine to assign boundary nodes
    Node* construct_boundary_node(const unsigned& n)
    {
      // Assign a node to the local node pointer
      // The dimension and number of values are taken from internal element data
      // The number of timesteps to be stored comes from the problem!
      this->node_pt(n) =
        new BoundaryNode<NODE_TYPE>(this->nodal_dimension(),
                                    this->nnodal_position_type(),
                                    this->required_nvalue(n));
      // Now return a pointer to the node, so that the mesh can find it
      return this->node_pt(n);
    }

    /// Overloaded boundary node allocation for unsteady problems
    Node* construct_boundary_node(const unsigned& n,
                                  TimeStepper* const& time_stepper_pt)
    {
      // Assign a node to the local node pointer
      // The dimension and number of values are taken from internal element data
      // The number of timesteps to be stored comes from the problem!
      this->node_pt(n) =
        new BoundaryNode<NODE_TYPE>(time_stepper_pt,
                                    this->nodal_dimension(),
                                    this->nnodal_position_type(),
                                    this->required_nvalue(n));
      // Now return a pointer to the node, so that the mesh can find it
      return this->node_pt(n);
    }


    /// Complete the setup of additional dependencies. Overloads
    /// empty virtual function in GeneralisedElement to determine the "geometric
    /// Data", i.e. the Data that affects the element's shape.
    /// This function is called (for all elements) at the very beginning of the
    /// equation numbering procedure to ensure that all dependencies
    /// are accounted for.
    void complete_setup_of_dependencies()
    {
      // Call function of underlying element
      ELEMENT::complete_setup_of_dependencies();

      // Call function of the element with moving nodes
      ElementWithMovingNodes::complete_setup_of_dependencies();
    }


    /// Assign local equation numbers for the underlying element, then
    /// deal with the additional geometric dofs
    void assign_all_generic_local_eqn_numbers(const bool& store_local_dof_pt)
    {
      // Call the generic local equation numbering scheme of the ELEMENT
      ELEMENT::assign_all_generic_local_eqn_numbers(store_local_dof_pt);
      ElementWithMovingNodes::assign_all_generic_local_eqn_numbers(
        store_local_dof_pt);
    }

    /// Compute the element's residuals vector and jacobian matrix
    void get_jacobian(Vector<double>& residuals, DenseMatrix<double>& jacobian)
    {
      /// Call the element's get jacobian function
      ELEMENT::get_jacobian(residuals, jacobian);

      // Now call the additional geometric Jacobian terms
      this->fill_in_jacobian_from_geometric_data(jacobian);
    }

    /// Compute the element's residuals vector and jacobian matrix
    void get_jacobian_and_mass_matrix(Vector<double>& residuals,
                                      DenseMatrix<double>& jacobian,
                                      DenseMatrix<double>& mass_matrix)
    {
      /// Call the element's get jacobian function
      ELEMENT::get_jacobian_and_mass_matrix(residuals, jacobian, mass_matrix);

      // Now call the additional geometric Jacobian terms
      this->fill_in_jacobian_from_geometric_data(jacobian);
    }


  private:
  };


} // namespace oomph
#endif
