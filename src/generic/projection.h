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
// Header file for a classes used to represent projectable elements

#ifndef OOMPH_PROJECTION_HEADER
#define OOMPH_PROJECTION_HEADER


#include "mesh.h"
#include "problem.h"
#include "multi_domain.h"
#include "shape.h"
#include "element_with_external_element.h"
#include "linear_solver.h"

// Using CG to solve the projection problem
#ifdef OOMPH_HAS_TRILINOS
#include "trilinos_solver.h"
#endif // #ifdef OOMPH_HAS_TRILINOS
#include "iterative_linear_solver.h"

// Use a preconditioner for the iterative solver
#include "preconditioner.h"
#include "general_purpose_preconditioners.h"

namespace oomph
{
  //==================================================================
  /// Template-free Base class for projectable elements
  //==================================================================
  class ProjectableElementBase
  {
  protected:
    /// Enumerated collection to specify which projection problem
    /// is to be solved.
    enum Projection_Type
    {
      Coordinate,
      Lagrangian,
      Value
    };

    /// Field that is currently being projected
    unsigned Projected_field;

    /// Time level we are projecting  (0=current values; >0: history values)
    unsigned Time_level_for_projection;

    /// When projecting the history values of the nodal coordinates,
    /// this is the coordinate we're projecting
    unsigned Projected_coordinate;

    /// When projecting the Lagrangain coordinates indicate which
    /// coordiante is to be projected
    unsigned Projected_lagrangian;

    /// Variable to indicate if we're projecting the history values of
    /// the nodal coordinates (Coordinate) the values themselves (Value), or the
    /// Lagrangian coordinates in Solid Mechanics problems (Lagrangian)
    Projection_Type Projection_type;

    /// Bool to know if we do projection or not. If false (the default)
    /// we solve the element's "real" equations rather than the projection
    /// equations
    bool Do_projection;


    /// Store number of "external" interactions that were assigned to
    /// the element before doing the projection.
    unsigned Backup_ninteraction;

    /// Remember if the element includes external geometric data
    /// when used in  non-projection mode (this is temporarily disabled during
    /// the projection)
    bool Backup_external_geometric_data;


    /// Remember if the element includes external data when used in
    /// non-projection mode (this is temporarily disabled during the
    /// projection)
    bool Backup_external_interaction_data;


  public:
    /// Constructor: Initialise data so that we don't project but solve
    /// the "normal" equations associated with the element.
    ProjectableElementBase()
      : Projected_field(0),
        Time_level_for_projection(0),
        Projected_coordinate(0),
        Projected_lagrangian(0),
        Projection_type(Value),
        Do_projection(false),
        Backup_ninteraction(0),
        Backup_external_geometric_data(false)
    {
    }

    /// Virtual destructor
    virtual ~ProjectableElementBase() {}

    /// Pure virtual function in which the element writer
    /// must specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    /// E.g. in Taylor Hood elements the fld-th velocities are stored
    /// at the fld-th value of the nodes; the pressures (the DIM-th
    /// field) are the DIM-th values at the vertex nodes etc.
    virtual Vector<std::pair<Data*, unsigned>> data_values_of_field(
      const unsigned& fld) = 0;

    /// Number of fields of the problem, so e.g. for 2D Navier Stokes
    /// this would be 3 (for the two velocities and one pressure)
    virtual unsigned nfields_for_projection() = 0;

    /// Number of history values to be stored for fld-th field
    /// (includes current value!)
    virtual unsigned nhistory_values_for_projection(const unsigned& fld) = 0;

    /// Number of history values to be stored when projecting
    /// the history values of the nodal coordinates (includes current value!)
    virtual unsigned nhistory_values_for_coordinate_projection() = 0;

    /// Return number of values (pinned or not) associated with
    /// field fld within the element. This must correspond to the
    /// number of shape functions returned in jacobian_and_shape_of_field(...).
    virtual unsigned nvalue_of_field(const unsigned& fld) = 0;

    /// Return local equation numbers associated with value ivalue
    /// of field fld within the element.
    virtual int local_equation(const unsigned& fld, const unsigned& ivalue) = 0;

    /// Return Jacobian of mapping and the shape functions associated
    /// with field fld. The number of shape functions must match the
    /// number of values specified in nvalue_of_field(...). For
    /// Lagrange-type interpolations the shape functinos are simply
    /// the "normal" nodal shape functions; if the element contains
    /// internal Data that is not associated with shape functions,
    /// simply set the corresonding shape function to 1.
    virtual double jacobian_and_shape_of_field(const unsigned& fld,
                                               const Vector<double>& s,
                                               Shape& psi) = 0;

    /// Return the fld-th field at local coordinates s
    /// at time-level time (time=0: current value; time>0: history values)
    virtual double get_field(const unsigned& time,
                             const unsigned& fld,
                             const Vector<double>& s) = 0;
  };


  //=====================================================================
  /// Wrapper class for projectable elements. Adds "projectability"
  /// to the underlying ELEMENT.
  //=====================================================================
  template<class ELEMENT>
  class ProjectableElement : public virtual ELEMENT,
                             public virtual ProjectableElementBase,
                             public virtual ElementWithExternalElement
  {
  protected:
    /// Overloaded version of fill_in_contribution_to_residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Do projection
      if (Do_projection)
      {
        this->residual_for_projection(
          residuals, GeneralisedElement::Dummy_matrix, 0);
      }
      // solve problem normally
      else
      {
        ELEMENT::fill_in_contribution_to_residuals(residuals);
      }
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
      ElementWithExternalElement::describe_local_dofs(out, current_string);
      ELEMENT::describe_local_dofs(out, current_string);
    }

    /// Overloaded version of fill_in_contribution_to_jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Do projection
      if (Do_projection)
      {
        this->residual_for_projection(residuals, jacobian, 1);
      }
      else
      {
        ELEMENT::fill_in_contribution_to_jacobian(residuals, jacobian);
      }
    }


  public:
    /// Constructor [this was only required explicitly
    /// from gcc 4.5.2 onwards...]
    ProjectableElement() {}

    /// Residual for the projection step. Flag indicates if we
    /// want the Jacobian (1) or not (0). Virtual so it can be
    /// overloaded if necessary
    virtual void residual_for_projection(Vector<double>& residuals,
                                         DenseMatrix<double>& jacobian,
                                         const unsigned& flag)
    {
      // Am I a solid element
      SolidFiniteElement* solid_el_pt = dynamic_cast<SolidFiniteElement*>(this);

      unsigned n_dim = dim();

      // Allocate storage for local coordinates
      Vector<double> s(n_dim);

      // Current field
      unsigned fld = Projected_field;

      // Number of nodes
      const unsigned n_node = this->nnode();
      // Number of positional dofs
      const unsigned n_position_type = this->nnodal_position_type();

      // Number of dof for current field
      const unsigned n_value = nvalue_of_field(fld);

      // Set the value of n_intpt
      const unsigned n_intpt = integral_pt()->nweight();

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the local coordinates of Gauss point
        for (unsigned i = 0; i < n_dim; i++) s[i] = integral_pt()->knot(ipt, i);

        // Get the integral weight
        double w = integral_pt()->weight(ipt);

        // Find same point in base mesh using external storage
        FiniteElement* other_el_pt = 0;
        other_el_pt = this->external_element_pt(0, ipt);
        Vector<double> other_s(n_dim);
        other_s = this->external_element_local_coord(0, ipt);

        ProjectableElement<ELEMENT>* cast_el_pt =
          dynamic_cast<ProjectableElement<ELEMENT>*>(other_el_pt);

        // Storage for the local equation and local unknown
        int local_eqn = 0, local_unknown = 0;

        // Now set up the three different projection problems
        switch (Projection_type)
        {
          case Lagrangian:
          {
            // If we don't have a solid element there is a problem
            if (solid_el_pt == 0)
            {
              throw OomphLibError("Trying to project Lagrangian coordinates in "
                                  "non-SolidElement\n",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }

            // Find the position shape functions
            Shape psi(n_node, n_position_type);
            // Get the position shape functions
            this->shape(s, psi);
            // Get the jacobian
            double J = this->J_eulerian(s);

            // Premultiply the weights and the Jacobian
            double W = w * J;

            // Get the value of the current position of the 0th coordinate
            // in the current element
            // at the current time level (which is the unkown)
            double interpolated_xi_proj = this->interpolated_x(s, 0);

            // Get the Lagrangian position in the other element
            double interpolated_xi_bar =
              dynamic_cast<SolidFiniteElement*>(cast_el_pt)
                ->interpolated_xi(other_s, Projected_lagrangian);

            // Now loop over the nodes and position dofs
            for (unsigned l = 0; l < n_node; ++l)
            {
              // Loop over position unknowns
              for (unsigned k = 0; k < n_position_type; ++k)
              {
                // The local equation is going to be the positional local
                // equation
                local_eqn = solid_el_pt->position_local_eqn(l, k, 0);

                // Now assemble residuals
                if (local_eqn >= 0)
                {
                  // calculate residuals
                  residuals[local_eqn] +=
                    (interpolated_xi_proj - interpolated_xi_bar) * psi(l, k) *
                    W;

                  // Calculate the jacobian
                  if (flag == 1)
                  {
                    for (unsigned l2 = 0; l2 < n_node; l2++)
                    {
                      // Loop over position dofs
                      for (unsigned k2 = 0; k2 < n_position_type; k2++)
                      {
                        local_unknown =
                          solid_el_pt->position_local_eqn(l2, k2, 0);
                        if (local_unknown >= 0)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            psi(l2, k2) * psi(l, k) * W;
                        }
                      }
                    }
                  } // end of jacobian
                }
              }
            }
          } // End of Lagrangian coordinate case

          break;

          // Now the coordinate history case
          case Coordinate:
          {
            // Find the position shape functions
            Shape psi(n_node, n_position_type);
            // Get the position shape functions
            this->shape(s, psi);
            // Get the jacobian
            double J = this->J_eulerian(s);

            // Premultiply the weights and the Jacobian
            double W = w * J;

            // Get the value of the current position in the current element
            // at the current time level (which is the unkown)
            double interpolated_x_proj = 0.0;
            // If we are a solid element read it out directly from the data
            if (solid_el_pt != 0)
            {
              interpolated_x_proj =
                this->interpolated_x(s, Projected_coordinate);
            }
            // Otherwise it's the field value at the current time
            else
            {
              interpolated_x_proj = this->get_field(0, fld, s);
            }

            // Get the position in the other element
            double interpolated_x_bar = cast_el_pt->interpolated_x(
              Time_level_for_projection, other_s, Projected_coordinate);

            // Now loop over the nodes and position dofs
            for (unsigned l = 0; l < n_node; ++l)
            {
              // Loop over position unknowns
              for (unsigned k = 0; k < n_position_type; ++k)
              {
                // If I'm a solid element
                if (solid_el_pt != 0)
                {
                  // The local equation is going to be the positional local
                  // equation
                  local_eqn =
                    solid_el_pt->position_local_eqn(l, k, Projected_coordinate);
                }
                // Otherwise just pick the local unknown of a field to
                // project into
                else
                {
                  // Complain if using generalised position types
                  // but this is not a solid element, because the storage
                  // is then not clear
                  if (n_position_type != 1)
                  {
                    throw OomphLibError("Trying to project generalised "
                                        "positions not in SolidElement\n",
                                        OOMPH_CURRENT_FUNCTION,
                                        OOMPH_EXCEPTION_LOCATION);
                  }
                  local_eqn = local_equation(fld, l);
                }

                // Now assemble residuals
                if (local_eqn >= 0)
                {
                  // calculate residuals
                  residuals[local_eqn] +=
                    (interpolated_x_proj - interpolated_x_bar) * psi(l, k) * W;

                  // Calculate the jacobian
                  if (flag == 1)
                  {
                    for (unsigned l2 = 0; l2 < n_node; l2++)
                    {
                      // Loop over position dofs
                      for (unsigned k2 = 0; k2 < n_position_type; k2++)
                      {
                        // If I'm a solid element
                        if (solid_el_pt != 0)
                        {
                          local_unknown = solid_el_pt->position_local_eqn(
                            l2, k2, Projected_coordinate);
                        }
                        else
                        {
                          local_unknown = local_equation(fld, l2);
                        }

                        if (local_unknown >= 0)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            psi(l2, k2) * psi(l, k) * W;
                        }
                      }
                    }
                  } // end of jacobian
                }
              }
            }
          } // End of coordinate case
          break;

          // Now the value case
          case Value:
          {
            // Field shape function
            Shape psi(n_value);

            // Calculate jacobian and shape functions for that field
            double J = jacobian_and_shape_of_field(fld, s, psi);

            // Premultiply the weights and the Jacobian
            double W = w * J;

            // Value of field in current element at current time level
            //(the unknown)
            unsigned now = 0;
            double interpolated_value_proj = this->get_field(now, fld, s);

            // Value of the interpolation of element located in base mesh
            double interpolated_value_bar =
              cast_el_pt->get_field(Time_level_for_projection, fld, other_s);

            // Loop over dofs of field fld
            for (unsigned l = 0; l < n_value; l++)
            {
              local_eqn = local_equation(fld, l);
              if (local_eqn >= 0)
              {
                // calculate residuals
                residuals[local_eqn] +=
                  (interpolated_value_proj - interpolated_value_bar) * psi[l] *
                  W;

                // Calculate the jacobian
                if (flag == 1)
                {
                  for (unsigned l2 = 0; l2 < n_value; l2++)
                  {
                    local_unknown = local_equation(fld, l2);
                    if (local_unknown >= 0)
                    {
                      jacobian(local_eqn, local_unknown) +=
                        psi[l2] * psi[l] * W;
                    }
                  }
                } // end of jacobian
              }
            }
          }
          break;

          default:
            throw OomphLibError("Wrong type specificied in Projection_type. "
                                "This should never happen\n",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
        } // End of the switch statement

      } // End of loop over ipt

    } // End of residual_for_projection function


    /// Use Eulerian coordinates for matching in locate_zeta
    /// when doing projection
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      if (Do_projection)
      {
        return nodal_position_gen(n, k, i);
      }
      else
      {
        return ELEMENT::zeta_nodal(n, k, i);
      }
    }


    /// Backup the element's state and
    /// switch it to projection mode.
    void enable_projection()
    {
      // Backup number of interaction
      Backup_ninteraction = ninteraction();

      // Backup flag for inclusion of geometric data
      if (add_external_geometric_data())
      {
        Backup_external_geometric_data = true;
      }
      else
      {
        Backup_external_geometric_data = false;
      }

      // Backup flag for inclusion of interaction data
      if (add_external_interaction_data())
      {
        Backup_external_interaction_data = true;
      }
      else
      {
        Backup_external_interaction_data = false;
      }

      // Actions to enable projection
      Do_projection = true;
      ignore_external_geometric_data();
      ignore_external_interaction_data();
      set_ninteraction(1);
    }

    /// Helper function to restore the element to the state
    /// it was in before we entered the projection mode and switch off
    /// projection mode.
    void disable_projection()
    {
      // Restore number of interaction
      set_ninteraction(Backup_ninteraction);

      // Restore geometric data
      if (Backup_external_geometric_data)
      {
        include_external_geometric_data();
      }
      else
      {
        ignore_external_geometric_data();
      }

      // Restore interaction data
      if (Backup_external_interaction_data)
      {
        include_external_interaction_data();
      }
      else
      {
        ignore_external_interaction_data();
      }

      Do_projection = false;
    }


    /// Project (history values of) coordintes
    void set_project_coordinates()
    {
      Projection_type = Coordinate;
    }

    /// Project (history values of) values
    void set_project_values()
    {
      Projection_type = Value;
    }

    /// Project (current and only values of) lagrangian coordinates
    void set_project_lagrangian()
    {
      Projection_type = Lagrangian;
    }


    /// Field that is currently being projected
    unsigned& projected_field()
    {
      return Projected_field;
    }

    /// Which history value are we projecting?
    unsigned& time_level_for_projection()
    {
      return Time_level_for_projection;
    }

    /// When projecting the history values of the nodal coordinates,
    /// this is the coordinate we're projecting
    unsigned& projected_coordinate()
    {
      return Projected_coordinate;
    }

    /// When projecting the Lagrangian coordinates this is
    /// the coordinate that is being projected
    unsigned& projected_lagrangian_coordinate()
    {
      return Projected_lagrangian;
    }


  }; // End of class


  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectableElement<ELEMENT>>
    : public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };

  // Forward definition of the friends of the class

  // The RefineableTriangleMesh
  // template<class FRIEND_PROJECTABLE_ELEMENT>
  // class RefineableTriangleMesh;

  // The RefineableTetgenMesh
  // template<class FRIEND_PROJECTABLE_ELEMENT>
  // class RefineableTetgenMesh;

  // The BackupMeshForProjection
  // template<class FRIEND_PROJECTABLE_ELEMENT>
  // class BackupMeshForProjection;

  //=======================================================================
  /// Projection problem. This is created during the adaptation
  /// of unstructured meshes and it is assumed that no boundary conditions
  /// have been set. If they have, they will be unset during the projection
  /// and must be reset afterwards.
  //=======================================================================
  template<class PROJECTABLE_ELEMENT>
  class ProjectionProblem : public virtual Problem
  {
    // The classes are friend whether the templated element of the
    // friend class is the same or not as the templated element of the
    // ProjectionProblem class
    template<class FRIEND_PROJECTABLE_ELEMENT>
    friend class RefineableTriangleMesh;
    template<class FRIEND_PROJECTABLE_ELEMENT>
    friend class RefineableTetgenMesh;
    template<class FRIEND_PROJECTABLE_ELEMENT>
    friend class BackupMeshForProjection;
    template<class FRIEND_PROJECTABLE_ELEMENT>
    friend class RefineableGmshTetMesh;

    // The classes are friend only when the templated element of the
    // friend class matches the templated element of the
    // ProjectionProblem class
    //  friend class RefineableTriangleMesh<class PROJECTABLE_ELEMENT>;
    //  friend class RefineableTetgenMesh<class PROJECTABLE_ELEMENT>;
    //  friend class BackupMeshForProjection<class PROJECTABLE_ELEMENT>;

  public:
    /// Suppress all output during projection phases
    void enable_suppress_output_during_projection()
    {
      Output_during_projection_suppressed = true;
    }

    /// Undo suppression of all output during projection phases
    void disable_suppress_output_during_projection()
    {
      Output_during_projection_suppressed = false;
    }

    /// Return the value of the flag about using an iterative solver for
    /// projection
    bool use_iterative_solver_for_projection()
    {
      return Use_iterative_solver_for_projection;
    }

    /// Enables the use of an iterative solver for projection
    void enable_use_iterative_solver_for_projection()
    {
      Use_iterative_solver_for_projection = true;
    }

    /// Disbales the use of an iterative solver for projection
    void disable_use_iterative_solver_for_projection()
    {
      Use_iterative_solver_for_projection = false;
    }

    /// Project from base into the problem's own mesh.
    void project(Mesh* base_mesh_pt, const bool& dont_project_positions = false)
    {
      // Use an iterative solver?
      if (Use_iterative_solver_for_projection)
      {
        // If oomph-lib has Trilinos installed then use the CG version
        // from Trilinos, otherwise use oomph-lib's own CG implementation
#ifdef OOMPH_HAS_TRILINOS
        // Check whether the problem is distributed?
        if (MPI_Helpers::mpi_has_been_initialised())
        {
          // Create a Trilinos Solver
          Iterative_solver_projection_pt = new TrilinosAztecOOSolver;
          // Select CG as the linear solver
          dynamic_cast<TrilinosAztecOOSolver*>(Iterative_solver_projection_pt)
            ->solver_type() = TrilinosAztecOOSolver::CG;
        }
        else
        {
          // Use CG to solve the projection problem
          Iterative_solver_projection_pt = new CG<CRDoubleMatrix>;
        }

        // Create the preconditioner
        Preconditioner_projection_pt = new MatrixBasedDiagPreconditioner();
        // Set the preconditioner for the solver
        Iterative_solver_projection_pt->preconditioner_pt() =
          Preconditioner_projection_pt;

        // Set CG as the linear solver
        Problem::linear_solver_pt() = Iterative_solver_projection_pt;

#else
        // Check whether the problem is distributed?
        if (!(MPI_Helpers::mpi_has_been_initialised()))
        {
          // If we did not installed Trilinos and the problem is not
          // distributed then we can use a (serial) preconditioned
          // iterative solver, otherwise, if we did not installed Trilinos
          // but the problem is distributed then we cannot use a
          // preconditioned iterative solver. Matrix multiplication in a
          // distributed environment is only performed by Trilinos. We
          // then use a direct solver for the projection problem.

          // Use CG to solve the projection problem
          Iterative_solver_projection_pt = new CG<CRDoubleMatrix>;

          // Create the preconditioner
          Preconditioner_projection_pt = new MatrixBasedDiagPreconditioner();
          // Set the preconditioner for the solver
          Iterative_solver_projection_pt->preconditioner_pt() =
            Preconditioner_projection_pt;

          // Set CG as the linear solver
          Problem::linear_solver_pt() = Iterative_solver_projection_pt;
        }
        else
        {
          // Use a direct solver. Do nothing
        }

#endif

      } // if (Use_iterative_solver_for_projection)

      // Backup verbosity in Newton solve status
      bool shut_up_in_newton_solve_backup = Shut_up_in_newton_solve;

      // Disable documentation of solve times
      bool backed_up_doc_time_enabled =
        linear_solver_pt()->is_doc_time_enabled();
      if (Output_during_projection_suppressed)
      {
        linear_solver_pt()->disable_doc_time();
      }

      // Display stats
      unsigned n_element = Problem::mesh_pt()->nelement();
      unsigned n_element1 = base_mesh_pt->nelement();
      unsigned n_node = Problem::mesh_pt()->nnode();
      unsigned n_node1 = base_mesh_pt->nnode();
      if (!Output_during_projection_suppressed)
      {
        oomph_info << "\n=============================\n";
        oomph_info << "Base mesh has " << n_element1 << " elements\n";
        oomph_info << "Target mesh has " << n_element << " elements\n";
        oomph_info << "Base mesh has " << n_node1 << " nodes\n";
        oomph_info << "Target mesh has " << n_node << " nodes\n";
        oomph_info << "=============================\n\n";
      }
      else
      {
        // Make Newton solver shut up too
        disable_info_in_newton_solve();
      }


      if (n_element == 0)
      {
        oomph_info << "Very odd -- no elements in target mesh; "
                   << " not doing anything in ProjectionProblem::project()\n";
        return;
      }

#ifdef PARANOID
      unsigned nnod = Problem::mesh_pt()->nnode();
      if (nnod == 0)
      {
        std::ostringstream error_stream;
        error_stream
          << "Mesh has no nodes! Please populate the Node_pt vector\n"
          << "otherwise projection won't work!\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // How many fields do we have to project?
      unsigned n_fields =
        dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(0))
          ->nfields_for_projection();

      // Spatial dimension of the problem
      unsigned n_dim = Problem::mesh_pt()->node_pt(0)->ndim();

      // Default number of history values
      unsigned n_history_values = 0;

      // Set things up for coordinate projection
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Switch to projection
        el_pt->enable_projection();
      }

      // Switch elements in base mesh to projection mode (required
      // to switch to use of Eulerian coordinates when identifying
      // corresponding points in the two meshes)
      for (unsigned e = 0; e < n_element1; e++)
      {
        PROJECTABLE_ELEMENT* el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(base_mesh_pt->element_pt(e));

        // Switch to projection
        el_pt->enable_projection();
      }


      // Set up multi domain interactions so we can locate the
      // values in the base mesh.
      // Note that it's important to switch elements to projection
      // mode first to ensure that matching is done based on Eulerian
      // rather than Lagrangian coordinates if pseudo-solid elements
      // are used.
      double t_start = TimingHelpers::timer();
      Multi_domain_functions::setup_multi_domain_interaction<
        PROJECTABLE_ELEMENT>(this, Problem::mesh_pt(), base_mesh_pt);
      if (!Output_during_projection_suppressed)
      {
        oomph_info
          << "CPU for setup of multi-domain interaction for projection: "
          << TimingHelpers::timer() - t_start << std::endl;
      }
      t_start = TimingHelpers::timer();


      // Let us first pin every degree of freedom
      // We shall unpin selected dofs for each different projection problem
      this->pin_all();

      if (!dont_project_positions)
      {
        //------------------Project coordinates first------------------------
        // If we have a solid element then we should also project Lagrangian
        // coordinates, but we can use the storage that MUST be provided for
        // the unknown positions for this.
        // If we can cast the first element of the mesh to a solid element,
        // then assume that we have solid elements
        if (dynamic_cast<SolidFiniteElement*>(
              Problem::mesh_pt()->element_pt(0)))
        {
          // Store current positions
          this->store_positions();

          // There are no history values for the Lagrangian coordinates
          // Set coordinate 0 for projection
          this->set_coordinate_for_projection(0);
          this->unpin_dofs_of_coordinate(0);

          // Loop over the Lagrangian coordinate
          const unsigned n_lagrangian =
            dynamic_cast<SolidNode*>(Problem::mesh_pt()->node_pt(0))
              ->nlagrangian();

          // Now loop over the lagrangian coordinates
          for (unsigned i = 0; i < n_lagrangian; ++i)
          {
            if (!Output_during_projection_suppressed)
            {
              oomph_info
                << "\n\n=============================================\n";
              oomph_info << "Projecting values for Lagrangian coordinate " << i
                         << std::endl;
              oomph_info << "=============================================\n\n";
            }

            // Set the coordinate for projection
            this->set_lagrangian_coordinate_for_projection(i);

            // Assign equation number
            unsigned ndof_tmp = assign_eqn_numbers();
            if (!Output_during_projection_suppressed)
            {
              oomph_info << "Number of equations for projection of Lagrangian "
                            "coordinate "
                         << " : " << ndof_tmp << std::endl
                         << std::endl;
            }


            if (Problem_is_nonlinear)
            {
              std::ostringstream error_stream;
              error_stream
                << "Solid mechanics problems will break if "
                   "Problem_is_nonlinear is\n"
                << "set to true in the projection problem because we're using "
                   "the\n "
                << "actual nodal positions to store the values of the "
                   "Lagrangian\n"
                << "coords. There shouldn't be any need for \n"
                << "Problem_is_nonlinear=true anyway, apart from debugging in "
                   "\n"
                << "which case you now know why this case will break!\n";
              OomphLibWarning(error_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
            }


            // Projection and interpolation
            Problem::newton_solve();

            // Move values back into Lagrangian coordinate for all nodes
            unsigned n_node = Problem::mesh_pt()->nnode();
            for (unsigned n = 0; n < n_node; ++n)
            {
              // Cast it to a solid node
              SolidNode* solid_node_pt =
                dynamic_cast<SolidNode*>(Problem::mesh_pt()->node_pt(n));
              // Now find number of types
              const unsigned n_position_type = solid_node_pt->nposition_type();
              // Find number of lagrangian types
              const unsigned n_lagrangian_type =
                solid_node_pt->nlagrangian_type();

              // If these are not the same, throw an error
              if (n_position_type != n_lagrangian_type)
              {
                std::ostringstream error_stream;
                error_stream
                  << "The number of generalised position dofs "
                  << n_position_type
                  << "\n not the same as the number of generalised lagrangian "
                     "dofs "
                  << n_lagrangian_type << ".\n"
                  << "Projection cannot be done at the moment, sorry.\n";

                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }

              // Now transfer the information across
              // from the first coordinate which was used during the projection
              for (unsigned k = 0; k < n_position_type; ++k)
              {
                solid_node_pt->xi_gen(k, i) = solid_node_pt->x_gen(k, 0);
                // Reset real values so that the Jacobians are correctly
                // computed next time around
                solid_node_pt->x_gen(k, 0) = Solid_backup[n](k, 0);
              }
            }
          } // End of loop over lagrangian coordinates

          // Now repin the dofs
          this->pin_dofs_of_coordinate(0);

          // Now project the position histories

          // Check number of history values for coordinates
          n_history_values = dynamic_cast<PROJECTABLE_ELEMENT*>(
                               Problem::mesh_pt()->element_pt(0))
                               ->nhistory_values_for_coordinate_projection();

          // Projection the coordinates only if there are history values
          if (n_history_values > 1)
          {
            for (unsigned i = 0; i < n_dim; i++)
            {
              if (!Output_during_projection_suppressed)
              {
                oomph_info
                  << "\n\n=============================================\n";
                oomph_info << "Projecting history values for coordinate " << i
                           << std::endl;
                oomph_info
                  << "=============================================\n\n";
              }

              // Set the coordinate for projection
              this->set_coordinate_for_projection(i);
              this->unpin_dofs_of_coordinate(i);

              // Loop over number of history values, beginning with the latest
              // one. Don't deal with current time.
              for (unsigned h_tim = n_history_values; h_tim > 1; h_tim--)
              {
                unsigned time_level = h_tim - 1;

                // Set time_level we are dealing with
                this->set_time_level_for_projection(time_level);

                // Assign equation number
                unsigned ndof_tmp = assign_eqn_numbers();
                if (!Output_during_projection_suppressed)
                {
                  oomph_info
                    << "Number of equations for projection of coordinate " << i
                    << " at time level " << time_level << " : " << ndof_tmp
                    << std::endl
                    << std::endl;
                }

                // Projection and interpolation
                Problem::newton_solve();

                // Move values back into history value of coordinate
                unsigned n_node = Problem::mesh_pt()->nnode();
                for (unsigned n = 0; n < n_node; ++n)
                {
                  // Cache the pointer to the node
                  Node* nod_pt = Problem::mesh_pt()->node_pt(n);
                  // Find the number of generalised dofs
                  const unsigned n_position_type = nod_pt->nposition_type();
                  // Now copy all back
                  for (unsigned k = 0; k < n_position_type; ++k)
                  {
                    nod_pt->x_gen(time_level, k, i) = nod_pt->x_gen(0, k, i);
                    // Reset real values so that the Jacobians are computed
                    // correctly next time around
                    nod_pt->x_gen(0, k, i) = Solid_backup[n](k, i);
                  }
                }
              }
              // Repin the positions
              this->pin_dofs_of_coordinate(i);
            }
          } // End of history value projection
        } // End of SolidElement case

        // Now for non solid elements, we are going to hijack the
        // first value as storage to be used for the projection of the history
        // values
        else
        {
          // Prepare for projection in value 0
          this->set_current_field_for_projection(0);
          this->unpin_dofs_of_field(0);

#ifdef PARANOID

          // The machinery used below  assumes that field 0 can actually
          // hold the coordinates i.e. that the field is interpolated
          // isoparametrically! The method will fail if there are no values
          // stored at the nodes.  Currently there are no examples of that --
          // the problem would only arise for elements that have all their
          // fields represented by internal data. Will fix this if/when such a
          // case arises...
          unsigned n_element = Problem::mesh_pt()->nelement();
          for (unsigned e = 0; e < n_element; e++)
          {
            FiniteElement* el_pt = Problem::mesh_pt()->finite_element_pt(e);
            unsigned nnod = el_pt->nnode();
            for (unsigned j = 0; j < nnod; j++)
            {
              Node* nod_pt = el_pt->node_pt(j);
              if (nod_pt->nvalue() == 0)
              {
                std::ostringstream error_stream;
                error_stream << "Node at  ";
                unsigned n = nod_pt->ndim();
                for (unsigned i = 0; i < n; i++)
                {
                  error_stream << nod_pt->x(i) << " ";
                }
                error_stream
                  << "\nhas no values. Projection (of coordinates) doesn't "
                     "work\n"
                  << "for such cases (at the moment), sorry! Send us the code\n"
                  << "where the problem arises and we'll try to implement "
                     "this\n"
                  << "(up to now unnecessary) capability...\n";
                throw OomphLibError(error_stream.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
            }
          }

#endif

          // Check number of history values for coordinates
          n_history_values = dynamic_cast<PROJECTABLE_ELEMENT*>(
                               Problem::mesh_pt()->element_pt(0))
                               ->nhistory_values_for_coordinate_projection();

          // Projection the coordinates only if there are history values
          if (n_history_values > 1)
          {
            for (unsigned i = 0; i < n_dim; i++)
            {
              if (!Output_during_projection_suppressed)
              {
                oomph_info
                  << "\n\n=============================================\n";
                oomph_info << "Projecting history values for coordinate " << i
                           << std::endl;
                oomph_info
                  << "=============================================\n\n";
              }

              // Set the coordinate for projection
              this->set_coordinate_for_projection(i);

              // Loop over number of history values, beginning with the latest
              // one. Don't deal with current time.
              for (unsigned h_tim = n_history_values; h_tim > 1; h_tim--)
              {
                unsigned time_level = h_tim - 1;

                // Set time_level we are dealing with
                this->set_time_level_for_projection(time_level);

                // Assign equation number
                unsigned ndof_tmp = assign_eqn_numbers();
                if (!Output_during_projection_suppressed)
                {
                  oomph_info
                    << "Number of equations for projection of coordinate " << i
                    << " at time level " << time_level << " : " << ndof_tmp
                    << std::endl
                    << std::endl;
                }

                // Projection and interpolation
                Problem::newton_solve();

                // Move values back into history value of coordinate
                unsigned n_element = Problem::mesh_pt()->nelement();
                for (unsigned e = 0; e < n_element; e++)
                {
                  PROJECTABLE_ELEMENT* new_el_pt =
                    dynamic_cast<PROJECTABLE_ELEMENT*>(
                      Problem::mesh_pt()->element_pt(e));

                  Vector<std::pair<Data*, unsigned>> data =
                    new_el_pt->data_values_of_field(0);

                  unsigned d_size = data.size();
                  for (unsigned d = 0; d < d_size; d++)
                  {
                    // Replace as coordinates
                    double coord = data[d].first->value(0, 0);
                    dynamic_cast<Node*>(data[d].first)->x(time_level, i) =
                      coord;
                  }
                }
              }
            }
          } // End of history value projection

          // Repin the dofs for field 0
          this->pin_dofs_of_field(0);

        } // End of non-SolidElement case


      } // end if for projection of coordinates

      // Disable projection of coordinates
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        el_pt->set_project_values();
      }

      // Loop over fields
      for (unsigned fld = 0; fld < n_fields; fld++)
      {
        // Let us first pin every degree of freedom
        // We shall unpin selected dofs for each different projection problem
        this->pin_all();

        // Do actions for this field
        this->set_current_field_for_projection(fld);
        this->unpin_dofs_of_field(fld);

        // Check number of history values
        n_history_values =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(0))
            ->nhistory_values_for_projection(fld);

        // Loop over number of history values
        // Beginning with the latest one
        for (unsigned h_tim = n_history_values; h_tim > 0; h_tim--)
        {
          unsigned time_level = h_tim - 1;
          if (!Output_during_projection_suppressed)
          {
            oomph_info << "\n=========================================\n";
            oomph_info << "Projecting field " << fld << " at time level "
                       << time_level << std::endl;
            oomph_info << "========================================\n";
          }

          // Set time_level we are dealing with
          this->set_time_level_for_projection(time_level);

          // Assign equation number
          unsigned ndof_tmp = assign_eqn_numbers();
          if (!Output_during_projection_suppressed)
          {
            oomph_info << "Number of equations for projection of field " << fld
                       << " at time level " << time_level << " : " << ndof_tmp
                       << std::endl
                       << std::endl;
          }

          // Projection and interpolation
          Problem::newton_solve();

          // Move computed values into the required time-level (not needed
          // for  current values which are done last -- they simply
          // stay where they are)
          if (time_level != 0)
          {
            for (unsigned e = 0; e < n_element; e++)
            {
              PROJECTABLE_ELEMENT* new_el_pt =
                dynamic_cast<PROJECTABLE_ELEMENT*>(
                  Problem::mesh_pt()->element_pt(e));

              Vector<std::pair<Data*, unsigned>> data =
                new_el_pt->data_values_of_field(fld);

              unsigned d_size = data.size();
              for (unsigned d = 0; d < d_size; d++)
              {
                // Move into time level
                double c_value = data[d].first->value(0, data[d].second);
                data[d].first->set_value(time_level, data[d].second, c_value);
              }
            }
          }
        } // End of loop over time levels

      } // End of loop over fields


      // Reset parameters of external storage and interactions
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Flush information associated with the external elements
        new_el_pt->flush_all_external_element_storage();

        new_el_pt->disable_projection();
      }

      for (unsigned e = 0; e < n_element1; e++)
      {
        PROJECTABLE_ELEMENT* el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(base_mesh_pt->element_pt(e));

        // Flush information associated with the external elements
        el_pt->flush_all_external_element_storage();

        // Disable  projection
        el_pt->disable_projection();
      }

      // Now unpin everything to restore the problem to its virgin state
      this->unpin_all();

      // Now cleanup the storage
      Solid_backup.clear();

      /* unsigned ndof_tmp=this->assign_eqn_numbers(); */
      if (!Output_during_projection_suppressed)
      {
        /* oomph_info << "Number of unknowns after project: "  */
        /*            << ndof_tmp << std::endl; */
        // std::ostringstream warn_message;
        // warn_message
        // << "WARNING: Ensure to assign equations numbers in the new mesh,\n"
        // << "this is done by calling the assign_eqn_numbers() method from\n"
        // << "the original Problem object that has an instance of the mesh.\n"
        // << "This may be performed in actions_after_adapt() if the
        // projection\n"
        // << "was performed as part of the mesh adaptation process\n\n";
        // OomphLibWarning(warn_message.str(),
        //                OOMPH_CURRENT_FUNCTION,
        //                OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Reset verbosity in Newton solver
        Shut_up_in_newton_solve = shut_up_in_newton_solve_backup;

        /// Disable documentation of solve times
        if (backed_up_doc_time_enabled)
        {
          linear_solver_pt()->enable_doc_time();
        }
      }

    } // End of function Projection

  private:
    /// Default constructor (made this private so only the friend class
    /// can call the constructor)
    ProjectionProblem()
    {
      // This is a linear problem so avoid checking the residual
      // after the solve
      Problem_is_nonlinear = false; // DO NOT CHANGE THIS -- EVER -- IN
                                    // SOLID MECHANICS PROBLEMS

      // By default suppress output during projection
      Output_during_projection_suppressed = true;

      // By default we use an iterative solver for projection
      Use_iterative_solver_for_projection = true;

      // Initialise the pointer to the solver and the preconditioner
      Iterative_solver_projection_pt = 0;
      Preconditioner_projection_pt = 0;
    }

    // Destructor
    ~ProjectionProblem()
    {
      if (Iterative_solver_projection_pt != 0)
      {
        // Clean up memory
        delete Iterative_solver_projection_pt;
        Iterative_solver_projection_pt = 0;
      }

      if (Preconditioner_projection_pt != 0)
      {
        delete Preconditioner_projection_pt;
        Preconditioner_projection_pt = 0;
      }
    }


    /// Helper function to store positions (the only things that
    /// have been set before doing projection
    void store_positions()
    {
      // No need to do anything if there are no elements (in fact, we
      // probably never get here...)
      if (Problem::mesh_pt()->nelement() == 0) return;

      // Deal with positional dofs if (pseudo-)solid element
      // If we can cast the first element to a SolidFiniteElement then
      // assume that we have a solid mesh
      SolidFiniteElement* solid_el_pt =
        dynamic_cast<SolidFiniteElement*>(Problem::mesh_pt()->element_pt(0));
      if (solid_el_pt != 0)
      {
        const unsigned n_node = this->mesh_pt()->nnode();
        Solid_backup.resize(n_node);
        // Read dimension and number of position values from the first node
        const unsigned n_dim = this->mesh_pt()->node_pt(0)->ndim();
        const unsigned n_position_type =
          this->mesh_pt()->node_pt(0)->nposition_type();
        // Loop over the nodes
        for (unsigned n = 0; n < n_node; n++)
        {
          // Cache a pointer to a solid node
          SolidNode* const solid_nod_pt =
            dynamic_cast<SolidNode*>(this->mesh_pt()->node_pt(n));
          // Now resize the appropriate storage
          Solid_backup[n].resize(n_position_type, n_dim);

          for (unsigned i = 0; i < n_dim; i++)
          {
            for (unsigned k = 0; k < n_position_type; k++)
            {
              Solid_backup[n](k, i) = solid_nod_pt->x_gen(k, i);
            }
          }
        }
      }
    }

    /// Helper function to restore positions (the only things that
    /// have been set before doing projection
    void restore_positions()
    {
      // No need to do anything if there are no elements (in fact, we
      // probably never get here...)
      if (Problem::mesh_pt()->nelement() == 0) return;

      // Deal with positional dofs if (pseudo-)solid element
      // If we can cast the first element to a SolidFiniteElement then
      // assume that we have a solid mesh
      SolidFiniteElement* solid_el_pt =
        dynamic_cast<SolidFiniteElement*>(Problem::mesh_pt()->element_pt(0));
      if (solid_el_pt != 0)
      {
        const unsigned n_node = this->mesh_pt()->nnode();
        // Read dimension and number of position values from the first node
        const unsigned n_dim = this->mesh_pt()->node_pt(0)->ndim();
        const unsigned n_position_type =
          this->mesh_pt()->node_pt(0)->nposition_type();
        // Loop over the nodes
        for (unsigned n = 0; n < n_node; n++)
        {
          // Cache a pointer to a solid node
          SolidNode* const solid_nod_pt =
            dynamic_cast<SolidNode*>(this->mesh_pt()->node_pt(n));

          for (unsigned i = 0; i < n_dim; i++)
          {
            for (unsigned k = 0; k < n_position_type; k++)
            {
              solid_nod_pt->x_gen(k, i) = Solid_backup[n](k, i);
            }
          }
        }
      }
    }

    /// Pin all the field values and position unknowns (bit inefficient)
    void pin_all()
    {
      // No need to do anything if there are no elements (in fact, we
      // probably never get here...)
      if (Problem::mesh_pt()->nelement() == 0) return;

      // Loop over all the elements
      const unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; ++e)
      {
        FiniteElement* el_pt = Problem::mesh_pt()->finite_element_pt(e);
        unsigned nint = el_pt->ninternal_data();
        for (unsigned j = 0; j < nint; j++)
        {
          Data* int_data_pt = el_pt->internal_data_pt(j);
          unsigned nval = int_data_pt->nvalue();
          for (unsigned i = 0; i < nval; i++)
          {
            int_data_pt->pin(i);
          }
        }

        unsigned nnod = el_pt->nnode();
        for (unsigned j = 0; j < nnod; j++)
        {
          Node* nod_pt = el_pt->node_pt(j);
          unsigned nval = nod_pt->nvalue();
          for (unsigned i = 0; i < nval; i++)
          {
            nod_pt->pin(i);
          }
        }
      }

      /// Do we have a solid mesh?
      SolidFiniteElement* solid_el_pt =
        dynamic_cast<SolidFiniteElement*>(Problem::mesh_pt()->element_pt(0));
      if (solid_el_pt != 0)
      {
        // Find number of nodes
        const unsigned n_node = this->mesh_pt()->nnode();
        // If no nodes then return
        if (n_node == 0)
        {
          return;
        }

        // Read dimension and number of position values from the first node
        const unsigned n_dim = this->mesh_pt()->node_pt(0)->ndim();
        const unsigned n_position_type =
          this->mesh_pt()->node_pt(0)->nposition_type();

        // Loop over the nodes
        for (unsigned n = 0; n < n_node; n++)
        {
          SolidNode* solid_nod_pt =
            dynamic_cast<SolidNode*>(this->mesh_pt()->node_pt(n));
          for (unsigned i = 0; i < n_dim; i++)
          {
            for (unsigned k = 0; k < n_position_type; k++)
            {
              solid_nod_pt->pin_position(k, i);
            }
          }
        }
      }
    }


    /// Unpin all the field values and position unknowns (bit inefficient)
    void unpin_all()
    {
      // No need to do anything if there are no elements (in fact, we
      // probably never get here...)
      if (Problem::mesh_pt()->nelement() == 0) return;

      // Loop over all the elements
      const unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; ++e)
      {
        // Cast the first element
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));
        // Find the number of fields
        unsigned n_fields = new_el_pt->nfields_for_projection();

        // Now loop over all fields
        for (unsigned f = 0; f < n_fields; f++)
        {
          // Extract the data and value for the field
          Vector<std::pair<Data*, unsigned>> data =
            new_el_pt->data_values_of_field(f);
          // Now loop over all the data and unpin the appropriate values
          unsigned d_size = data.size();
          for (unsigned d = 0; d < d_size; d++)
          {
            data[d].first->unpin(data[d].second);
          }
        }
      }

      /// Do we have a solid mesh?
      SolidFiniteElement* solid_el_pt =
        dynamic_cast<SolidFiniteElement*>(Problem::mesh_pt()->element_pt(0));
      if (solid_el_pt != 0)
      {
        // Find number of nodes
        const unsigned n_node = this->mesh_pt()->nnode();
        // If no nodes then return
        if (n_node == 0)
        {
          return;
        }

        // Read dimension and number of position values from the first node
        const unsigned n_dim = this->mesh_pt()->node_pt(0)->ndim();
        const unsigned n_position_type =
          this->mesh_pt()->node_pt(0)->nposition_type();

        // Loop over the nodes
        for (unsigned n = 0; n < n_node; n++)
        {
          SolidNode* solid_nod_pt =
            dynamic_cast<SolidNode*>(this->mesh_pt()->node_pt(n));
          for (unsigned i = 0; i < n_dim; i++)
          {
            for (unsigned k = 0; k < n_position_type; k++)
            {
              solid_nod_pt->unpin_position(k, i);
            }
          }
        }
      }
    }

    /// Helper function to unpin the values of coordinate coord
    void unpin_dofs_of_coordinate(const unsigned& coord)
    {
      // Loop over the nodes
      const unsigned n_node = Problem::mesh_pt()->nnode();
      // If there are no nodes return immediately
      if (n_node == 0)
      {
        return;
      }

      // Find the number of position values (should be the same for all nodes)
      const unsigned n_position_type =
        Problem::mesh_pt()->node_pt(0)->nposition_type();

      for (unsigned n = 0; n < n_node; ++n)
      {
        // Cache access to the Node as a solid node
        SolidNode* solid_nod_pt =
          static_cast<SolidNode*>(Problem::mesh_pt()->node_pt(n));
        // Unpin all position types associated with the given coordinate
        for (unsigned k = 0; k < n_position_type; ++k)
        {
          solid_nod_pt->unpin_position(k, coord);
        }
      }
    }

    /// Helper function to unpin the values of coordinate coord
    void pin_dofs_of_coordinate(const unsigned& coord)
    {
      // Loop over the nodes
      const unsigned n_node = Problem::mesh_pt()->nnode();
      // If there are no nodes return immediately
      if (n_node == 0)
      {
        return;
      }

      // Find the number of position values (should be the same for all nodes)
      const unsigned n_position_type =
        Problem::mesh_pt()->node_pt(0)->nposition_type();

      for (unsigned n = 0; n < n_node; ++n)
      {
        // Cache access to the Node as a solid node
        SolidNode* solid_nod_pt =
          static_cast<SolidNode*>(Problem::mesh_pt()->node_pt(n));
        // Unpin all position types associated with the given coordinate
        for (unsigned k = 0; k < n_position_type; ++k)
        {
          solid_nod_pt->pin_position(k, coord);
        }
      }
    }


    /// Helper function to unpin dofs of fld-th field
    void unpin_dofs_of_field(const unsigned& fld)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        Vector<std::pair<Data*, unsigned>> data =
          new_el_pt->data_values_of_field(fld);

        unsigned d_size = data.size();
        for (unsigned d = 0; d < d_size; d++)
        {
          data[d].first->unpin(data[d].second);
        }
      }
    }

    /// Helper function to unpin dofs of fld-th field
    void pin_dofs_of_field(const unsigned& fld)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        Vector<std::pair<Data*, unsigned>> data =
          new_el_pt->data_values_of_field(fld);

        unsigned d_size = data.size();
        for (unsigned d = 0; d < d_size; d++)
        {
          data[d].first->pin(data[d].second);
        }
      }
    }

    /// Helper function to set time level for projection
    void set_time_level_for_projection(const unsigned& time_level)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Set what time we are dealing with
        el_pt->time_level_for_projection() = time_level;
      }
    }

    /// Set the coordinate for projection
    void set_coordinate_for_projection(const unsigned& i)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Set that we are solving a projected coordinate problem
        new_el_pt->set_project_coordinates();

        // Set the projected coordinate
        new_el_pt->projected_coordinate() = i;
      }
    }

    /// Set the Lagrangian coordinate for projection
    void set_lagrangian_coordinate_for_projection(const unsigned& i)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Set that we are solving a projected Lagrangian coordinate problem
        new_el_pt->set_project_lagrangian();

        // Set the projected coordinate
        new_el_pt->projected_lagrangian_coordinate() = i;
      }
    }

    /// Set current field for projection
    void set_current_field_for_projection(const unsigned& fld)
    {
      unsigned n_element = Problem::mesh_pt()->nelement();
      for (unsigned e = 0; e < n_element; e++)
      {
        PROJECTABLE_ELEMENT* new_el_pt =
          dynamic_cast<PROJECTABLE_ELEMENT*>(Problem::mesh_pt()->element_pt(e));

        // Set current field
        new_el_pt->projected_field() = fld;
      }
    }

  private:
    /// Backup for pin status of solid node's position Data
    Vector<DenseMatrix<double>> Solid_backup;

    /// Flag to suppress output during projection
    bool Output_during_projection_suppressed;

    // Use an iterative solver for solving the system of equations
    bool Use_iterative_solver_for_projection;

    // The iterative solver to solve the projection problem
    IterativeLinearSolver* Iterative_solver_projection_pt;

    // The preconditioner for the solver
    Preconditioner* Preconditioner_projection_pt;
  };


} // namespace oomph

#endif
