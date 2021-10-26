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
// This header file contains elements that combine two element types in
// a generic way.

#ifndef OOMPH_PSEUDO_SOLID_REMESH_ELEMENTS_HEADER
#define OOMPH_PSEUDO_SOLID_REMESH_ELEMENTS_HEADER

#include "elements.h"

namespace oomph
{
  //===========================================================================
  /// Helper namespace for pseudo-elastic elements
  //===========================================================================
  namespace PseudoSolidHelper
  {
    /// Static variable to hold the default value for the pseudo-solid's
    /// inertia parameter Lambda^2.
    extern double Zero;

  } // namespace PseudoSolidHelper


  //==========================================================================
  /// A templated class that permits combination two different element types,
  /// for the solution of problems in deforming domains. The first template
  /// paremter BASIC is the standard element and the second SOLID solves
  /// the equations that are used to control the mesh deformation.
  //==========================================================================
  template<class BASIC, class SOLID>
  class PseudoSolidNodeUpdateElement : public virtual BASIC,
                                       public virtual SOLID

  {
    /// Boolean flag to indicate shape derivative method
    bool Shape_derivs_by_direct_fd;

  public:
    /// Constructor, call the BASIC and SOLID elements' constructors and
    /// set the "density" parameter for solid element to zero
    PseudoSolidNodeUpdateElement()
      : BASIC(), SOLID(), Shape_derivs_by_direct_fd(true)
    {
      SOLID::lambda_sq_pt() = &PseudoSolidHelper::Zero;
    }

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
      BASIC::describe_local_dofs(out, current_string);
      SOLID::describe_local_dofs(out, current_string);
    }

    /// Compute norm of solution: use the version in the BASIC
    /// class if there's any ambiguity
    void compute_norm(double& el_norm)
    {
      BASIC::compute_norm(el_norm);
    }

    /// The required number of values is the sum of the two
    unsigned required_nvalue(const unsigned& n) const
    {
      return BASIC::required_nvalue(n) + SOLID::required_nvalue(n);
    }

    /// We assume that the solid stuff is stored at the end of
    /// the nodes, i.e. its index is the number of continuously interplated
    /// values in the BASIC equations.
    int solid_p_nodal_index() const
    {
      // At the moment, we can't handle this case in generality so throw an
      // error if the solid pressure is stored at the nodes
      if (SOLID::solid_p_nodal_index() >= 0)
      {
        throw OomphLibError("Cannot handle (non-refineable) continuous solid "
                            "pressure interpolation",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      return SOLID::solid_p_nodal_index();
    }

    /// Final override for the residuals function. Contributions are
    /// added from both underlying element types
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_residuals(residuals);
      // Add the solid equations contribution
      SOLID::fill_in_contribution_to_residuals(residuals);
    }

    /// Final override for jacobian function: Contributions are
    /// included from both the underlying element types
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_jacobian(residuals, jacobian);
      // Call the solid equations
      SOLID::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);
    }

    /// Final override for mass matrix function: contributions
    /// are included from both the underlying element types
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);
      // Call the solid equations
      SOLID::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);
    }


    /// Evaluate shape derivatives by direct finite differencing
    void evaluate_shape_derivs_by_direct_fd()
    {
      Shape_derivs_by_direct_fd = true;
    }

    /// Evaluate shape derivatives by chain rule
    void evaluate_shape_derivs_by_chain_rule()
    {
      Shape_derivs_by_direct_fd = false;
    }


    /// Fill in the shape derivatives of the BASIC equations
    /// w.r.t. the solid position dofs
    void fill_in_shape_derivatives(DenseMatrix<double>& jacobian)
    {
      // Default is to use finite differences
      if (Shape_derivs_by_direct_fd)
      {
        this->fill_in_shape_derivatives_by_fd(jacobian);
      }
      // Otherwise need to do a bit more work
      else
      {
        // Calculate storage requirements
        const unsigned n_dof = this->ndof();
        const unsigned n_node = this->nnode();
        const unsigned nodal_dim = this->nodal_dimension();

        // If there are no nodes or dofs return
        if ((n_dof == 0) || (n_node == 0))
        {
          return;
        }

        // Generalised dofs have NOT been considered, shout
        if (this->nnodal_position_type() != 1)
        {
          throw OomphLibError("Shape derivatives do not (yet) allow for "
                              "generalised position dofs\n",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Storage for derivatives of residuals w.r.t. nodal coordinates
        RankThreeTensor<double> dresidual_dnodal_coordinates(
          n_dof, nodal_dim, n_node, 0.0);

        // Get the analytic derivatives for the BASIC equations
        BASIC::get_dresidual_dnodal_coordinates(dresidual_dnodal_coordinates);

        // Now add the appropriate contributions to the Jacobian
        int local_unknown = 0;

        // Loop over dofs
        //(this will include the solid dofs,
        // but all those contributions should be zero)
        for (unsigned l = 0; l < n_dof; l++)
        {
          // Loop over the nodes
          for (unsigned n = 0; n < n_node; n++)
          {
            // Loop over the position_types (only one)
            unsigned k = 0;
            // Loop over the coordinates
            for (unsigned i = 0; i < nodal_dim; i++)
            {
              // Get the equation of the local unknown
              local_unknown = this->position_local_eqn(n, k, i);

              // If not pinned, add the contribution to the Jacobian
              if (local_unknown >= 0)
              {
                jacobian(l, local_unknown) +=
                  dresidual_dnodal_coordinates(l, i, n);
              }
            }
          }
        }
      }
    }


    /// Fill in the derivatives of the BASIC equations
    /// w.r.t. the solid position dofs
    void fill_in_shape_derivatives_by_fd(DenseMatrix<double>& jacobian)
    {
      // Flag to indicate if we use first or second order FD
      // bool use_first_order_fd=false;

      // Find the number of nodes
      const unsigned n_node = this->nnode();

      // If there aren't any nodes, then return straight away
      if (n_node == 0)
      {
        return;
      }

      // Call the update function to ensure that the element is in
      // a consistent state before finite differencing starts
      this->update_before_solid_position_fd();

      // Get the number of position dofs and dimensions at the node
      const unsigned n_position_type = this->nnodal_position_type();
      const unsigned nodal_dim = this->nodal_dimension();

      // Find the number of dofs in the element
      const unsigned n_dof = this->ndof();

      // Create residual newres vectors
      Vector<double> residuals(n_dof);
      Vector<double> newres(n_dof);
      // Vector<double> newres_minus(n_dof);

      // Calculate the residuals (for the BASIC) equations
      // Need to do this using fill_in because get_residuals will
      // compute all residuals for the problem, which is
      // a little ineffecient
      for (unsigned m = 0; m < n_dof; m++)
      {
        residuals[m] = 0.0;
      }
      BASIC::fill_in_contribution_to_residuals(residuals);

      // Need to determine which degrees of freedom are solid degrees of
      // freedom
      // A vector of booleans that will be true if the dof is associated
      // with the solid equations
      std::vector<bool> dof_is_solid(n_dof, false);

      // Now set all solid positional dofs in the vector
      for (unsigned n = 0; n < n_node; n++)
      {
        for (unsigned k = 0; k < n_position_type; k++)
        {
          for (unsigned i = 0; i < nodal_dim; i++)
          {
            int local_dof = this->position_local_eqn(n, k, i);
            if (local_dof >= 0)
            {
              dof_is_solid[local_dof] = true;
            }
          }
        }
      }

      // Add the solid pressures (in solid elements without
      // solid pressure the number will be zero).
      unsigned n_solid_pres = this->npres_solid();
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        int local_dof = this->solid_p_local_eqn(l);
        if (local_dof >= 0)
        {
          dof_is_solid[local_dof] = true;
        }
      }


      // Integer storage for local unknown
      int local_unknown = 0;

      // Use default value defined in GeneralisedElement
      const double fd_step = this->Default_fd_jacobian_step;

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over position dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over dimension
          for (unsigned i = 0; i < nodal_dim; i++)
          {
            // If the variable is free
            local_unknown = this->position_local_eqn(n, k, i);
            if (local_unknown >= 0)
            {
              // Store a pointer to the (generalised) Eulerian nodal position
              double* const value_pt = &(this->node_pt(n)->x_gen(k, i));

              // Save the old value of the (generalised) Eulerian nodal position
              const double old_var = *value_pt;

              // Increment the (generalised) Eulerian nodal position
              *value_pt += fd_step;

              // Perform any auxialiary node updates
              this->node_pt(n)->perform_auxiliary_node_update_fct();

              // Calculate the new residuals
              // Need to do this using fill_in because get_residuals will
              // compute all residuals for the problem, which is
              // a little ineffecient
              for (unsigned m = 0; m < n_dof; m++)
              {
                newres[m] = 0.0;
              }
              BASIC::fill_in_contribution_to_residuals(newres);

              //          if (use_first_order_fd)
              {
                // Do forward finite differences
                for (unsigned m = 0; m < n_dof; m++)
                {
                  // Stick the entry into the Jacobian matrix
                  // but only if it's not a solid dof
                  if (dof_is_solid[m] == false)
                  {
                    jacobian(m, local_unknown) =
                      (newres[m] - residuals[m]) / fd_step;
                  }
                }
              }
              //           else
              //            {
              //             //Take backwards step for the  (generalised)
              //             Eulerian nodal
              //             // position
              //             node_pt(n)->x_gen(k,i) = old_var-fd_step;

              //             //Calculate the new residuals at backward position
              //             //BASIC::get_residuals(newres_minus);

              //             //Do central finite differences
              //             for(unsigned m=0;m<n_dof;m++)
              //              {
              //               //Stick the entry into the Jacobian matrix
              //               jacobian(m,local_unknown) =
              //                (newres[m] - newres_minus[m])/(2.0*fd_step);
              //              }
              //            }

              // Reset the (generalised) Eulerian nodal position
              *value_pt = old_var;

              // Perform any auxialiary node updates
              this->node_pt(n)->perform_auxiliary_node_update_fct();
            }
          }
        }
      }

      // End of finite difference loop
      // Final reset of any dependent data
      this->reset_after_solid_position_fd();
    }


    /// Specify Data that affects the geometry of the element
    /// by adding the position Data to the set that's passed in.
    /// (This functionality is required in FSI problems; set is used to
    /// avoid double counting).
    void identify_geometric_data(std::set<Data*>& geometric_data_pt)
    {
      // Loop over the node update data and add to the set
      const unsigned n_node = this->nnode();
      for (unsigned j = 0; j < n_node; j++)
      {
        geometric_data_pt.insert(
          dynamic_cast<SolidNode*>(this->node_pt(j))->variable_position_pt());
      }
    }


    /// Overload the output function: Call that of the basic element
    void output(std::ostream& outfile)
    {
      BASIC::output(outfile);
    }

    /// Output function: Plot at n_p plot points using the basic
    /// element's output function
    void output(std::ostream& outfile, const unsigned& n_p)
    {
      BASIC::output(outfile, n_p);
    }

    /// Overload the output function: Call that of the basic element
    void output(FILE* file_pt)
    {
      BASIC::output(file_pt);
    }

    /// Output function is just the same as the basic equations
    void output(FILE* file_pt, const unsigned& n_p)
    {
      BASIC::output(file_pt, n_p);
    }

    /// Number of 'flux' terms for Z2 error estimation: Error estimation
    /// is based on error in BASIC element
    unsigned num_Z2_flux_terms()
    {
      return BASIC::num_Z2_flux_terms();
    }


    /// Plot the error when compared against a given exact flux.
    /// Also calculates the norm of the error and that of the exact flux.
    /// Use version in BASIC element
    void compute_exact_Z2_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
      double& error,
      double& norm)
    {
      BASIC::compute_exact_Z2_error(outfile, exact_flux_pt, error, norm);
    }

    /// 'Flux' vector for Z2 error estimation: Error estimation
    /// is based on error in BASIC element
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      BASIC::get_Z2_flux(s, flux);
    }

    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return BASIC::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return BASIC::vertex_node_pt(j);
    }

    /// Order of recovery shape functions for Z2 error estimation: Done
    /// for BASIC element since it determines the refinement
    unsigned nrecovery_order()
    {
      return BASIC::nrecovery_order();
    }


    /// The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into.
    unsigned ndof_types() const
    {
      return BASIC::ndof_types() + SOLID::ndof_types();
    }

    /// return the number of DOF types associated with the BASIC
    /// elements in this combined element
    unsigned nbasic_dof_types() const
    {
      return BASIC::ndof_types();
    }

    /// return the number of DOF types associated with the SOLID
    /// elements in this combined element
    unsigned nsolid_dof_types() const
    {
      return SOLID::ndof_types();
    }

    /// Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF type" that this unknown is associated with.
    /// This method combines the get_dof_numbers_for_unknowns(...)
    /// method for the BASIC and SOLID elements. The basic elements
    /// retain their DOF type numbering and the SOLID elements
    /// DOF type numbers are incremented by nbasic_dof_types().
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
    {
      // get the solid list
      std::list<std::pair<unsigned long, unsigned>> solid_list;
      SOLID::get_dof_numbers_for_unknowns(solid_list);

      // get the basic list
      BASIC::get_dof_numbers_for_unknowns(dof_lookup_list);

      // get the number of basic dof types
      unsigned nbasic_dof_types = BASIC::ndof_types();

      // add the solid lookup list to the basic lookup list
      // incrementing the solid dof numbers by nbasic_dof_types
      typedef std::list<std::pair<unsigned long, unsigned>>::iterator IT;
      for (IT it = solid_list.begin(); it != solid_list.end(); it++)
      {
        std::pair<unsigned long, unsigned> new_pair;
        new_pair.first = it->first;
        new_pair.second = it->second + nbasic_dof_types;
        dof_lookup_list.push_front(new_pair);
      }
    }
  };

  /// Explicit definition of the face geometry of these elements
  template<class BASIC, class SOLID>
  class FaceGeometry<PseudoSolidNodeUpdateElement<BASIC, SOLID>>
    : public virtual FaceGeometry<SOLID>
  {
  public:
    /// Constuctor calls the constructor of the SolidQElement
    /// (Only the Intel compiler seems to need this!)
    FaceGeometry() : FaceGeometry<SOLID>() {}
  };

  /// Explicit definition of the face geometry of these elements
  template<class BASIC, class SOLID>
  class FaceGeometry<FaceGeometry<PseudoSolidNodeUpdateElement<BASIC, SOLID>>>
    : public virtual FaceGeometry<FaceGeometry<SOLID>>
  {
  public:
    /// Constuctor calls the constructor of the SolidQElement
    /// (Only the Intel compiler seems to need this!)
    FaceGeometry() : FaceGeometry<FaceGeometry<SOLID>>() {}

  protected:
  };


  //===================================================================
  /// Refineable version of the PseudoSolidNodeUpdateELement
  //===================================================================
  template<class BASIC, class SOLID>
  class RefineablePseudoSolidNodeUpdateElement : public virtual BASIC,
                                                 public virtual SOLID
  {
  public:
    /// Constructor, call the BASIC and SOLID elements' constructors and
    /// set the "density" parameter for solid element to zero
    RefineablePseudoSolidNodeUpdateElement()
      : RefineableElement(), BASIC(), SOLID()
    {
      SOLID::lambda_sq_pt() = &PseudoSolidHelper::Zero;
    }

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
      BASIC::describe_local_dofs(out, current_string);
      SOLID::describe_local_dofs(out, current_string);
    }

    /// The required number of values is the sum of the two
    unsigned required_nvalue(const unsigned& n) const
    {
      return BASIC::required_nvalue(n) + SOLID::required_nvalue(n);
    }

    /// The number of continuously interpolated values is the
    /// sum of the SOLID and BASIC values
    unsigned ncont_interpolated_values() const
    {
      return BASIC::ncont_interpolated_values() +
             SOLID::ncont_interpolated_values();
    }

    /// We assume that the solid stuff is stored at the end of
    /// the nodes, i.e. its index is the number of continuously interplated
    /// values in the BASIC equations.
    int solid_p_nodal_index() const
    {
      // Find the index in the solid
      int solid_p_index = SOLID::solid_p_nodal_index();
      // If there is a solid pressure at the nodes, return the
      // index after all the BASIC stuff
      if (solid_p_index >= 0)
      {
        return BASIC::ncont_interpolated_values() +
               SOLID::solid_p_nodal_index();
      }
      else
      {
        return solid_p_index;
      }
    }

    /// Final override for residuals function: adds contributions
    /// from both underlying element types
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_residuals(residuals);
      // Call the solid equations
      SOLID::fill_in_contribution_to_residuals(residuals);
    }

    /// Final override for jacobian function: Calls get_jacobian() for
    /// both of the underlying element types
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Call the solid equations
      SOLID::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives_by_fd(jacobian);
    }

    /// Final override for mass matrix function: contributions
    /// are included from both the underlying element types
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call the basic equations first
      BASIC::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);
      // Call the solid equations
      SOLID::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives_by_fd(jacobian);
    }


    /// Fill in the derivatives of the BASIC equations
    /// w.r.t. to the solid position dofs, taking hanging nodes
    /// into account
    void fill_in_shape_derivatives_by_fd(DenseMatrix<double>& jacobian)
    {
      // Find the number of nodes
      const unsigned n_node = this->nnode();

      // If there are no nodes, return straight away
      if (n_node == 0)
      {
        return;
      }

      // Call the update function to ensure that the element is in
      // a consistent state before finite differencing starts
      this->update_before_solid_position_fd();

      //  bool use_first_order_fd=false;

      // Find the number of positional dofs and nodal dimension
      const unsigned n_position_type = this->nnodal_position_type();
      const unsigned nodal_dim = this->nodal_dimension();

      // Find the number of dofs in the element
      const unsigned n_dof = this->ndof();

      // Create residual newres vectors
      Vector<double> residuals(n_dof);
      Vector<double> newres(n_dof);
      // Vector<double> newres_minus(n_dof);

      // Calculate the residuals (for the BASIC) equations
      // Need to do this using fill_in because get_residuals will
      // compute all residuals for the problem, which is
      // a little ineffecient
      for (unsigned m = 0; m < n_dof; m++)
      {
        residuals[m] = 0.0;
      }
      BASIC::fill_in_contribution_to_residuals(residuals);

      // Need to determine which degrees of freedom are solid degrees of
      // freedom
      // A vector of booleans that will be true if the dof is associated
      // with the solid equations
      std::vector<bool> dof_is_solid(n_dof, false);

      // Now set all solid positional dofs in the vector
      // This is a bit more involved because we need to take account of
      // any hanging nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Get pointer to the local node
        Node* const local_node_pt = this->node_pt(n);

        // If the node is not a hanging node
        if (local_node_pt->is_hanging() == false)
        {
          for (unsigned k = 0; k < n_position_type; k++)
          {
            for (unsigned i = 0; i < nodal_dim; i++)
            {
              int local_dof = this->position_local_eqn(n, k, i);
              if (local_dof >= 0)
              {
                dof_is_solid[local_dof] = true;
              }
            }
          }
        }
        // Otherwise the node is hanging
        else
        {
          // Find the local hanging object
          HangInfo* hang_info_pt = local_node_pt->hanging_pt();
          // Loop over the master nodes
          const unsigned n_master = hang_info_pt->nmaster();
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the local equation numbers for the master node
            DenseMatrix<int> Position_local_eqn_at_node =
              this->local_position_hang_eqn(hang_info_pt->master_node_pt(m));

            // Loop over position dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
              // Loop over dimension
              for (unsigned i = 0; i < nodal_dim; i++)
              {
                int local_dof = Position_local_eqn_at_node(k, i);
                if (local_dof >= 0)
                {
                  dof_is_solid[local_dof] = true;
                }
              }
            }
          }
        }
      } // End of loop over nodes

      // Add the solid pressures (in solid elements without
      // solid pressure the number will be zero).
      unsigned n_solid_pres = this->npres_solid();
      // Now is the solid pressure hanging
      const int solid_p_index = this->solid_p_nodal_index();
      // Find out whether the solid node is hanging
      std::vector<bool> solid_p_is_hanging(n_solid_pres);
      // If we have nodal solid pressures then read out the hanging status
      if (solid_p_index >= 0)
      {
        // Loop over the solid dofs
        for (unsigned l = 0; l < n_solid_pres; l++)
        {
          solid_p_is_hanging[l] =
            this->solid_pressure_node_pt(l)->is_hanging(solid_p_index);
        }
      }
      // Otherwise the pressure is not nodal, so cannot hang
      else
      {
        for (unsigned l = 0; l < n_solid_pres; l++)
        {
          solid_p_is_hanging[l] = false;
        }
      }

      // Now we can loop of the dofs again to actually set that the appropriate
      // dofs are solid
      for (unsigned l = 0; l < n_solid_pres; l++)
      {
        // If the solid pressure is not hanging
        // we just read out the local equation numbers directly
        if (solid_p_is_hanging[l] == false)
        {
          int local_dof = this->solid_p_local_eqn(l);
          if (local_dof >= 0)
          {
            dof_is_solid[local_dof] = true;
          }
        }
        // Otherwise solid pressure is hanging and we need to take
        // care of the master nodes
        else
        {
          // Find the local hanging object
          HangInfo* hang_info_pt =
            this->solid_pressure_node_pt(l)->hanging_pt(solid_p_index);
          // Loop over the master nodes
          const unsigned n_master = hang_info_pt->nmaster();
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the local dof
            int local_dof = this->local_hang_eqn(
              hang_info_pt->master_node_pt(m), solid_p_index);

            if (local_dof >= 0)
            {
              dof_is_solid[local_dof] = true;
            }
          }
        }
      } // end of loop over solid pressure dofs


      // Used default value defined in GeneralisedElement
      const double fd_step = this->Default_fd_jacobian_step;

      // Integer storage for local unknowns
      int local_unknown = 0;

      // Loop over the nodes
      for (unsigned l = 0; l < n_node; l++)
      {
        // Get the pointer to the node
        Node* const local_node_pt = this->node_pt(l);

        // If the node is not a hanging node
        if (local_node_pt->is_hanging() == false)
        {
          // Loop over position dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Loop over dimension
            for (unsigned i = 0; i < nodal_dim; i++)
            {
              local_unknown = this->position_local_eqn(l, k, i);
              // If the variable is free
              if (local_unknown >= 0)
              {
                // Store a pointer to the (generalised) Eulerian nodal position
                double* const value_pt = &(local_node_pt->x_gen(k, i));

                // Save the old value of the (generalised) Eulerian nodal
                // position
                const double old_var = *value_pt;

                // Increment the  (generalised) Eulerian nodal position
                *value_pt += fd_step;

                // Perform any auxialiary node updates
                local_node_pt->perform_auxiliary_node_update_fct();

                // Calculate the new residuals
                // Need to do this using fill_in because get_residuals will
                // compute all residuals for the problem, which is
                // a little ineffecient
                for (unsigned m = 0; m < n_dof; m++)
                {
                  newres[m] = 0.0;
                }
                BASIC::fill_in_contribution_to_residuals(newres);


                //           if (use_first_order_fd)
                {
                  // Do forward finite differences
                  for (unsigned m = 0; m < n_dof; m++)
                  {
                    // Stick the entry into the Jacobian matrix
                    // But only if it's not a solid dof
                    if (dof_is_solid[m] == false)
                    {
                      jacobian(m, local_unknown) =
                        (newres[m] - residuals[m]) / fd_step;
                    }
                  }
                }
                //             else
                //              {
                //               //Take backwards step for the  (generalised)
                //               Eulerian nodal
                //               // position
                //               node_pt(l)->x_gen(k,i) = old_var-fd_step;

                //               //Calculate the new residuals at backward
                //               position BASIC::get_residuals(newres_minus);

                //               //Do central finite differences
                //               for(unsigned m=0;m<n_dof;m++)
                //                {
                //                 //Stick the entry into the Jacobian matrix
                //                 jacobian(m,local_unknown) =
                //                  (newres[m] - newres_minus[m])/(2.0*fd_step);
                //                }
                //              }

                // Reset the (generalised) Eulerian nodal position
                *value_pt = old_var;

                // Perform any auxialiary node updates
                local_node_pt->perform_auxiliary_node_update_fct();
              }
            }
          }
        }
        // Otherwise it's a hanging node
        else
        {
          // Find the local hanging object
          HangInfo* hang_info_pt = local_node_pt->hanging_pt();
          // Loop over the master nodes
          const unsigned n_master = hang_info_pt->nmaster();
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the pointer to the master node
            Node* const master_node_pt = hang_info_pt->master_node_pt(m);

            // Get the local equation numbers for the master node
            DenseMatrix<int> Position_local_eqn_at_node =
              this->local_position_hang_eqn(master_node_pt);

            // Loop over position dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
              // Loop over dimension
              for (unsigned i = 0; i < nodal_dim; i++)
              {
                local_unknown = Position_local_eqn_at_node(k, i);
                // If the variable is free
                if (local_unknown >= 0)
                {
                  // Store a pointer to the (generalised) Eulerian nodal
                  // position
                  double* const value_pt = &(master_node_pt->x_gen(k, i));

                  // Save the old value of the (generalised) Eulerian nodal
                  // position
                  const double old_var = *value_pt;

                  // Increment the  (generalised) Eulerian nodal position
                  *value_pt += fd_step;

                  // Perform any auxialiary node updates
                  master_node_pt->perform_auxiliary_node_update_fct();

                  // Calculate the new residuals
                  // Need to do this using fill_in because get_residuals will
                  // compute all residuals for the problem, which is
                  // a little ineffecient
                  for (unsigned m = 0; m < n_dof; m++)
                  {
                    newres[m] = 0.0;
                  }
                  BASIC::fill_in_contribution_to_residuals(newres);

                  //            if (use_first_order_fd)
                  {
                    // Do forward finite differences
                    for (unsigned m = 0; m < n_dof; m++)
                    {
                      // Stick the entry into the Jacobian matrix
                      // But only if it's not a solid dof
                      if (dof_is_solid[m] == false)
                      {
                        jacobian(m, local_unknown) =
                          (newres[m] - residuals[m]) / fd_step;
                      }
                    }
                  }
                  //               else
                  //                {
                  //                 //Take backwards step for the (generalised)
                  //                 Eulerian nodal
                  //                 // position
                  //                 master_node_pt->x_gen(k,i) =
                  //                 old_var-fd_step;

                  //                 //Calculate the new residuals at backward
                  //                 position
                  //                 BASIC::get_residuals(newres_minus);

                  //                 //Do central finite differences
                  //                 for(unsigned m=0;m<n_dof;m++)
                  //                  {
                  //                   //Stick the entry into the Jacobian
                  //                   matrix jacobian(m,local_unknown) =
                  //                    (newres[m] -
                  //                    newres_minus[m])/(2.0*fd_step);
                  //                  }
                  //                }

                  // Reset the (generalised) Eulerian nodal position
                  *value_pt = old_var;

                  // Perform any auxialiary node updates
                  master_node_pt->perform_auxiliary_node_update_fct();
                }
              }
            }
          }
        } // End of hanging node case

      } // End of loop over nodes

      // End of finite difference loop

      // Final reset of any dependent data
      this->reset_after_solid_position_fd();
    }


    /// Specify Data that affects the geometry of the element
    /// by adding the position Data to the set that's passed in.
    /// (This functionality is required in FSI problems; set is used to
    /// avoid double counting). Refineable version includes hanging nodes
    void identify_geometric_data(std::set<Data*>& geometric_data_pt)
    {
      // Loop over the node update data and add to the set
      const unsigned n_node = this->nnode();
      for (unsigned j = 0; j < n_node; j++)
      {
        // If the node is a hanging node
        if (this->node_pt(j)->is_hanging())
        {
          // Find the local hang info object
          HangInfo* hang_info_pt = this->node_pt(j)->hanging_pt();

          // Find the number of master nodes
          unsigned n_master = hang_info_pt->nmaster();

          // Loop over the master nodes
          for (unsigned m = 0; m < n_master; m++)
          {
            // Get the m-th master node
            Node* Master_node_pt = hang_info_pt->master_node_pt(m);

            // Add to set
            geometric_data_pt.insert(
              dynamic_cast<SolidNode*>(Master_node_pt)->variable_position_pt());
          }
        }
        // Not hanging
        else
        {
          // Add node itself to set
          geometric_data_pt.insert(
            dynamic_cast<SolidNode*>(this->node_pt(j))->variable_position_pt());
        }
      }
    }


    /// Final override for the assign__additional_local_eqn_numbers():
    ///  Call the version for the BASIC element
    void assign_additional_local_eqn_numbers()
    {
      BASIC::assign_additional_local_eqn_numbers();
      SOLID::assign_additional_local_eqn_numbers();
    }

    /// Call rebuild_from_sons() for both of the underlying element types
    void rebuild_from_sons(Mesh*& mesh_pt)
    {
      BASIC::rebuild_from_sons(mesh_pt);
      SOLID::rebuild_from_sons(mesh_pt);
    }

    /// Call get_interpolated_values(...) for both of the underlying element
    /// types
    void get_interpolated_values(const unsigned& t,
                                 const Vector<double>& s,
                                 Vector<double>& values)
    {
      Vector<double> basic_values;
      BASIC::get_interpolated_values(t, s, basic_values);
      Vector<double> solid_values;
      SOLID::get_interpolated_values(t, s, solid_values);

      // Now add the basic value first
      for (Vector<double>::iterator it = basic_values.begin();
           it != basic_values.end();
           ++it)
      {
        values.push_back(*it);
      }
      // Then the solid
      for (Vector<double>::iterator it = solid_values.begin();
           it != solid_values.end();
           ++it)
      {
        values.push_back(*it);
      }
    }


    /// Call get_interpolated_values(...) for both of the underlying element
    /// types
    void get_interpolated_values(const Vector<double>& s,
                                 Vector<double>& values)
    {
      Vector<double> basic_values;
      BASIC::get_interpolated_values(s, basic_values);
      Vector<double> solid_values;
      SOLID::get_interpolated_values(s, solid_values);

      // Now add the basic value first
      for (Vector<double>::iterator it = basic_values.begin();
           it != basic_values.end();
           ++it)
      {
        values.push_back(*it);
      }
      // Then the solid
      for (Vector<double>::iterator it = solid_values.begin();
           it != solid_values.end();
           ++it)
      {
        values.push_back(*it);
      }
    }

    /// We must compose the underlying interpolating nodes from
    /// the BASIC and SOLID equations, the BASIC ones are first
    Node* interpolating_node_pt(const unsigned& n, const int& value_id)
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::interpolating_node_pt(n, value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::interpolating_node_pt(n, (value_id - n_basic_values));
      }
    }

    /// The pressure nodes are the corner nodes, so when value_id==0,
    /// the fraction is the same as the 1d node number, 0 or 1.
    double local_one_d_fraction_of_interpolating_node(const unsigned& n1d,
                                                      const unsigned& i,
                                                      const int& value_id)
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::local_one_d_fraction_of_interpolating_node(
          n1d, i, value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::local_one_d_fraction_of_interpolating_node(
          n1d, i, (value_id - n_basic_values));
      }
    }


    /// The velocity nodes are the same as the geometric nodes. The
    /// pressure nodes must be calculated by using the same methods as
    /// the geometric nodes, but by recalling that there are only two pressure
    /// nodes per edge.
    Node* get_interpolating_node_at_local_coordinate(const Vector<double>& s,
                                                     const int& value_id)
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::get_interpolating_node_at_local_coordinate(s, value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::get_interpolating_node_at_local_coordinate(
          s, (value_id - n_basic_values));
      }
    }

    /// The number of 1d pressure nodes is 2, otherwise we have
    /// the positional nodes
    unsigned ninterpolating_node_1d(const int& value_id)
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::ninterpolating_node_1d(value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::ninterpolating_node_1d((value_id - n_basic_values));
      }
    }

    /// The number of pressure nodes is 2^DIM. The number of
    /// velocity nodes is the same as the number of geometric nodes.
    unsigned ninterpolating_node(const int& value_id)
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::ninterpolating_node(value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::ninterpolating_node((value_id - n_basic_values));
      }
    }

    /// The basis interpolating the pressure is given by pshape().
    /// / The basis interpolating the velocity is shape().
    void interpolating_basis(const Vector<double>& s,
                             Shape& psi,
                             const int& value_id) const
    {
      // Find the number of interpolated values in the BASIC equations
      int n_basic_values = BASIC::ncont_interpolated_values();
      // If the id is below this number, we call the BASIC functon
      if (value_id < n_basic_values)
      {
        return BASIC::interpolating_basis(s, psi, value_id);
      }
      // Otherwise it's the solid and its value_id is the the current
      // it minus n_basic_values
      else
      {
        return SOLID::interpolating_basis(s, psi, (value_id - n_basic_values));
      }
    }


    /// Number of 'flux' terms for Z2 error estimation: Error estimation
    /// is based on error in BASIC element
    unsigned num_Z2_flux_terms()
    {
      return BASIC::num_Z2_flux_terms();
    }

    /// 'Flux' vector for Z2 error estimation: Error estimation
    /// is based on error in BASIC element
    void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
    {
      BASIC::get_Z2_flux(s, flux);
    }

    /// Perform additional hanging node procedures for variables
    /// that are not interpolated by all nodes. Done for both of the
    /// underlying element types.
    void further_setup_hanging_nodes()
    {
      BASIC::further_setup_hanging_nodes();
      SOLID::further_setup_hanging_nodes();
    }


    /// Build function: Call the one for the SOLID element since it
    /// calls the one basic build function automatically.
    void build(Mesh*& mesh_pt,
               Vector<Node*>& new_node_pt,
               bool& was_already_built,
               std::ofstream& new_nodes_file)
    {
      SOLID::build(mesh_pt, new_node_pt, was_already_built, new_nodes_file);
    }


    /// Build function: Call the one for the SOLID element since it
    /// calls the one basic build function automatically.
    void build(Mesh*& mesh_pt,
               Vector<Node*>& new_node_pt,
               bool& was_already_built)
    {
      SOLID::build(mesh_pt, new_node_pt, was_already_built);
    }

    ///  Further build: Done for both of the
    /// underlying element types.
    void further_build()
    {
      BASIC::further_build();
      SOLID::further_build();
    }


    /// Number of vertex nodes in the element
    unsigned nvertex_node() const
    {
      return BASIC::nvertex_node();
    }

    /// Pointer to the j-th vertex node in the element
    Node* vertex_node_pt(const unsigned& j) const
    {
      return BASIC::vertex_node_pt(j);
    }

    /// Compute norm of solution. Use version in BASIC element.
    void compute_norm(double& norm)
    {
      BASIC::compute_norm(norm);
    }

    /// Plot the error when compared against a given exact flux.
    /// Also calculates the norm of the error and that of the exact flux.
    /// Use version in BASIC element
    void compute_exact_Z2_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_flux_pt,
      double& error,
      double& norm)
    {
      BASIC::compute_exact_Z2_error(outfile, exact_flux_pt, error, norm);
    }

    /// Order of recovery shape functions for Z2 error estimation: Done
    /// for BASIC element since it determines the refinement
    unsigned nrecovery_order()
    {
      return BASIC::nrecovery_order();
    }

    /// Overload the output function: Use that of the BASIC element
    void output(std::ostream& outfile)
    {
      BASIC::output(outfile);
    }

    /// Output function, plotting at n_p points: Use that of the BASIC element
    void output(std::ostream& outfile, const unsigned& n_p)
    {
      BASIC::output(outfile, n_p);
    }

    /// Overload the output function: Use that of the BASIC element
    void output(FILE* file_pt)
    {
      BASIC::output(file_pt);
    }

    /// Output function: Use that of the BASIC element
    void output(FILE* file_pt, const unsigned& n_p)
    {
      BASIC::output(file_pt, n_p);
    }

    /// The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into.
    unsigned ndof_types() const
    {
      return BASIC::ndof_types() + SOLID::ndof_types();
    }

    /// return the number of DOF types associated with the BASIC
    /// elements in this combined element
    unsigned nbasic_dof_types() const
    {
      return BASIC::ndof_types();
    }

    /// return the number of DOF types associated with the SOLID
    /// elements in this combined element
    unsigned nsolid_dof_types() const
    {
      return SOLID::ndof_types();
    }

    /// Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF type" that this unknown is associated with.
    /// This method combines the get_dof_numbers_for_unknowns(...)
    /// method for the BASIC and SOLID elements. The basic elements
    /// retain their DOF type numbering and the SOLID elements
    /// DOF type numbers are incremented by nbasic_dof_types().
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
    {
      // get the solid list
      std::list<std::pair<unsigned long, unsigned>> solid_list;
      SOLID::get_dof_numbers_for_unknowns(solid_list);

      // get the basic list
      BASIC::get_dof_numbers_for_unknowns(dof_lookup_list);

      // get the number of basic dof types
      unsigned nbasic_dof_types = BASIC::ndof_types();

      // add the solid lookup list to the basic lookup list
      // incrementing the solid dof numbers by nbasic_dof_types
      typedef std::list<std::pair<unsigned long, unsigned>>::iterator IT;
      for (IT it = solid_list.begin(); it != solid_list.end(); it++)
      {
        std::pair<unsigned long, unsigned> new_pair;
        new_pair.first = it->first;
        new_pair.second = it->second + nbasic_dof_types;
        dof_lookup_list.push_front(new_pair);
      }
    }
  };


  /// Explicit definition of the face geometry of these elements
  template<class BASIC, class SOLID>
  class FaceGeometry<RefineablePseudoSolidNodeUpdateElement<BASIC, SOLID>>
    : public virtual FaceGeometry<SOLID>
  {
  public:
    /// Constructor calls the constructor of the SolidQElement
    /// (Only the Intel compiler seems to need this!)
    FaceGeometry() : FaceGeometry<SOLID>() {}

  protected:
  };

  /// Explicit definition of the face geometry of these elements
  template<class BASIC, class SOLID>
  class FaceGeometry<
    FaceGeometry<RefineablePseudoSolidNodeUpdateElement<BASIC, SOLID>>>
    : public virtual FaceGeometry<FaceGeometry<SOLID>>
  {
  public:
    /// Constuctor calls the constructor of the SolidQElement
    /// (Only the Intel compiler seems to need this!)
    FaceGeometry() : FaceGeometry<FaceGeometry<SOLID>>() {}

  protected:
  };


} // namespace oomph

#endif
