// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// Header file for elements that are used to apply surface loads to
// the equations of axisymmetric linear elasticity

#ifndef OOMPH_AXISYMMETRIC_LINEAR_ELASTICITY_TRACTION_ELEMENTS_HEADER
#define OOMPH_AXISYMMETRIC_LINEAR_ELASTICITY_TRACTION_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
//#include "../generic/Qelements.h"
#include "src/generic/Qelements.h"
#include "src/generic/element_with_external_element.h"


namespace oomph
{
  //=======================================================================
  /// Namespace containing the zero traction function for
  /// axisymmetric linear elasticity traction elements
  //=======================================================================
  namespace AxisymmetricLinearElasticityTractionElementHelper
  {
    //=======================================================================
    /// Default load function (zero traction)
    //=======================================================================
    void Zero_traction_fct(const double& time,
                           const Vector<double>& x,
                           const Vector<double>& N,
                           Vector<double>& load)
    {
      unsigned n_dim = load.size();
      for (unsigned i = 0; i < n_dim; i++)
      {
        load[i] = 0.0;
      }
    }

  } // namespace AxisymmetricLinearElasticityTractionElementHelper


  //======================================================================
  /// A class for elements that allow the imposition of an applied traction
  /// in the equations of axisymmetric linear elasticity.
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELEMENT>
  class AxisymmetricLinearElasticityTractionElement
    : public virtual FaceGeometry<ELEMENT>,
      public virtual FaceElement
  {
  protected:
    /// Index at which the i-th displacement component is stored
    Vector<unsigned> U_index_axisymmetric_linear_elasticity_traction;

    /// Pointer to an imposed traction function. Arguments:
    /// Eulerian coordinate; outer unit normal;
    /// applied traction. (Not all of the input arguments will be
    /// required for all specific load functions but the list should
    /// cover all cases)
    void (*Traction_fct_pt)(const double& time,
                            const Vector<double>& x,
                            const Vector<double>& n,
                            Vector<double>& result);


    /// Get the traction vector: Pass number of integration point
    /// (dummy), Eulerian coordinate and normal vector and return the load
    /// vector (not all of the input arguments will be required for all specific
    /// load functions but the list should cover all cases). This function is
    /// virtual so it can be overloaded for FSI.
    virtual void get_traction(const double& time,
                              const unsigned& intpt,
                              const Vector<double>& x,
                              const Vector<double>& n,
                              Vector<double>& traction)
    {
      Traction_fct_pt(time, x, n, traction);
    }


    /// Helper function that actually calculates the residuals
    // This small level of indirection is required to avoid calling
    // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
    // which causes all kinds of pain if overloading later on
    void fill_in_contribution_to_residuals_axisymmetric_linear_elasticity_traction(
      Vector<double>& residuals);


  public:
    /// Constructor, which takes a "bulk" element and the
    /// value of the index and its limit
    AxisymmetricLinearElasticityTractionElement(
      FiniteElement* const& element_pt, const int& face_index)
      : FaceGeometry<ELEMENT>(), FaceElement()
    {
      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Find the dimension of the problem
      unsigned n_dim = element_pt->nodal_dimension();

      // Find the index at which the displacement unknowns are stored
      ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
      this->U_index_axisymmetric_linear_elasticity_traction.resize(n_dim + 1);
      for (unsigned i = 0; i < n_dim + 1; i++)
      {
        this->U_index_axisymmetric_linear_elasticity_traction[i] =
          cast_element_pt->u_index_axisymmetric_linear_elasticity(i);
      }

      // Zero traction
      Traction_fct_pt =
        &AxisymmetricLinearElasticityTractionElementHelper::Zero_traction_fct;
    }


    /// Reference to the traction function pointer
    void (*&traction_fct_pt())(const double& time,
                               const Vector<double>& x,
                               const Vector<double>& n,
                               Vector<double>& traction)
    {
      return Traction_fct_pt;
    }


    /// Return the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_residuals_axisymmetric_linear_elasticity_traction(
        residuals);
    }


    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the residuals
      fill_in_contribution_to_residuals_axisymmetric_linear_elasticity_traction(
        residuals);
    }

    /// Specify the value of nodal zeta from the face geometry
    /// The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Number of dimensions
      unsigned n_dim = this->nodal_dimension();

      // Find out how many nodes there are
      const unsigned n_node = nnode();

      // Get continuous time from timestepper of first node
      double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

      // Set up memory for the shape functions
      Shape psi(n_node);

      // Local and global coordinates
      Vector<double> s(n_dim - 1);
      Vector<double> interpolated_x(n_dim);

      // Tecplot header info
      outfile << this->tecplot_zone_string(n_plot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);

        // Outer unit normal
        Vector<double> unit_normal(n_dim);
        outer_unit_normal(s, unit_normal);

        // Find the shape functions
        shape(s, psi);

        // Initialise to zero
        for (unsigned i = 0; i < n_dim; i++)
        {
          interpolated_x[i] = 0.0;
        }

        // Calculate stuff
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += nodal_position(l, i) * psi[l];
          }
        }

        // Get the imposed flux
        Vector<double> traction(3);

        // Dummy integration point
        unsigned ipt = 0;

        get_traction(time, ipt, interpolated_x, unit_normal, traction);

        // Output the x,y,..
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_x[i] << " ";
        }

        // Output the traction components
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          outfile << traction[i] << " ";
        }

        // Output normal
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << unit_normal[i] << " ";
        }
        outfile << std::endl;
      }
    }

    /// C_style output function
    void output(FILE* file_pt)
    {
      FiniteElement::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FiniteElement::output(file_pt, n_plot);
    }


    /// Compute traction vector at specified local coordinate
    /// Should only be used for post-processing; ignores dependence
    /// on integration point!
    void traction(const double& time,
                  const Vector<double>& s,
                  Vector<double>& traction);
  };

  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Compute traction vector at specified local coordinate
  /// Should only be used for post-processing; ignores dependence
  /// on integration point!
  //=====================================================================
  template<class ELEMENT>
  void AxisymmetricLinearElasticityTractionElement<ELEMENT>::traction(
    const double& time, const Vector<double>& s, Vector<double>& traction)
  {
    unsigned n_dim = this->nodal_dimension();

    // Position vector
    Vector<double> x(n_dim);
    interpolated_x(s, x);

    // Outer unit normal (only in r and z direction!)
    Vector<double> unit_normal(n_dim);
    outer_unit_normal(s, unit_normal);

    // Dummy
    unsigned ipt = 0;

    // Traction vector
    get_traction(time, ipt, x, unit_normal, traction);
  }


  //=====================================================================
  /// Return the residuals for the
  /// AxisymmetricLinearElasticityTractionElement equations
  //=====================================================================
  template<class ELEMENT>
  void AxisymmetricLinearElasticityTractionElement<ELEMENT>::
    fill_in_contribution_to_residuals_axisymmetric_linear_elasticity_traction(
      Vector<double>& residuals)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Get continuous time from timestepper of first node
    double time = node_pt(0)->time_stepper_pt()->time_pt()->time();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();
    if (n_position_type != 1)
    {
      throw OomphLibError("AxisymmetricLinearElasticity is not yet implemented "
                          "for more than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Cache the nodal indices at which the displacement components are stored
    unsigned u_nodal_index[n_dim + 1];
    for (unsigned i = 0; i < n_dim + 1; i++)
    {
      u_nodal_index[i] =
        this->U_index_axisymmetric_linear_elasticity_traction[i];
    }

    // Integer to hold the local equation number
    int local_eqn = 0;

    // Set up memory for the shape functions
    // Note that in this case, the number of lagrangian coordinates is always
    // equal to the dimension of the nodes
    Shape psi(n_node);
    DShape dpsids(n_node, n_dim - 1);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Only need to call the local derivatives
      dshape_local_at_knot(ipt, psi, dpsids);

      // Calculate the Eulerian and Lagrangian coordinates
      Vector<double> interpolated_x(n_dim, 0.0);

      // Also calculate the surface Vectors (derivatives wrt local coordinates)
      DenseMatrix<double> interpolated_A(n_dim - 1, n_dim, 0.0);

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Calculate the Eulerian coords
          const double x_local = nodal_position(l, i);
          interpolated_x[i] += x_local * psi(l);

          // Loop over LOCAL derivative directions, to calculate the tangent(s)
          for (unsigned j = 0; j < n_dim - 1; j++)
          {
            interpolated_A(j, i) += x_local * dpsids(l, j);
          }
        }
      }

      // Now find the local metric tensor from the tangent Vectors
      DenseMatrix<double> A(n_dim - 1);
      for (unsigned i = 0; i < n_dim - 1; i++)
      {
        for (unsigned j = 0; j < n_dim - 1; j++)
        {
          // Initialise surface metric tensor to zero
          A(i, j) = 0.0;

          // Take the dot product
          for (unsigned k = 0; k < n_dim; k++)
          {
            A(i, j) += interpolated_A(i, k) * interpolated_A(j, k);
          }
        }
      }

      // Get the outer unit normal
      Vector<double> interpolated_normal(n_dim);
      outer_unit_normal(ipt, interpolated_normal);

      // Find the determinant of the metric tensor
      double Adet = 0.0;
      switch (n_dim)
      {
        case 2:
          Adet = A(0, 0);
          break;
        case 3:
          Adet = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
          break;
        default:
          throw OomphLibError(
            "Wrong dimension in AxisymmetricLinearElasticityTractionElement",
            "AxisymmetricLinearElasticityTractionElement::fill_in_contribution_"
            "to_residuals()",
            OOMPH_EXCEPTION_LOCATION);
      }

      // Premultiply the weights and the square-root of the determinant of
      // the metric tensor
      double W = w * sqrt(Adet);

      // Now calculate the load
      Vector<double> traction(n_dim + 1);
      get_traction(time, ipt, interpolated_x, interpolated_normal, traction);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          // Equation number
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[i]);
          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Add the loading terms to the residuals
            residuals[local_eqn] -=
              traction[i] * psi(l) * interpolated_x[0] * W;
          }
        }
      } // End of loop over shape functions
    } // End of loop over integration points
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //======================================================================
  /// A class for elements that allow the imposition of an applied traction
  /// in the equations of axisymmetric linear elasticity from an adjacent
  /// axisymmetric Navier Stokes element in a linearised FSI problem.
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELASTICITY_BULK_ELEMENT, class NAVIER_STOKES_BULK_ELEMENT>
  class FSIAxisymmetricLinearElasticityTractionElement
    : public virtual FaceGeometry<ELASTICITY_BULK_ELEMENT>,
      public virtual FaceElement,
      public virtual ElementWithExternalElement
  {
  protected:
    /// Pointer to the ratio, \f$ Q \f$ , of the stress used to
    /// non-dimensionalise the fluid stresses to the stress used to
    /// non-dimensionalise the solid stresses.
    double* Q_pt;

    /// Static default value for the ratio of stress scales
    /// used in the fluid and solid equations (default is 1.0)
    static double Default_Q_Value;

    /// Index at which the i-th displacement component is stored
    Vector<unsigned> U_index_axisym_fsi_traction;

    /// Helper function that actually calculates the residuals
    // This small level of indirection is required to avoid calling
    // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
    // which causes all kinds of pain if overloading later on
    void fill_in_contribution_to_residuals_axisym_fsi_traction(
      Vector<double>& residuals);

  public:
    /// Constructor, which takes a "bulk" element and the
    /// value of the index and its limit
    FSIAxisymmetricLinearElasticityTractionElement(
      FiniteElement* const& element_pt, const int& face_index)
      : FaceGeometry<ELASTICITY_BULK_ELEMENT>(),
        FaceElement(),
        Q_pt(&Default_Q_Value)
    {
      // Set source element storage: one interaction with an external
      // element -- the Navier Stokes bulk element that provides the traction
      this->set_ninteraction(1);

      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Find the dimension of the problem
      unsigned n_dim = element_pt->nodal_dimension();

      // Find the index at which the displacement unknowns are stored
      ELASTICITY_BULK_ELEMENT* cast_element_pt =
        dynamic_cast<ELASTICITY_BULK_ELEMENT*>(element_pt);
      this->U_index_axisym_fsi_traction.resize(n_dim + 1);
      for (unsigned i = 0; i < n_dim + 1; i++)
      {
        this->U_index_axisym_fsi_traction[i] =
          cast_element_pt->u_index_axisymmetric_linear_elasticity(i);
      }
    }


    /// Return the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_residuals_axisym_fsi_traction(residuals);
    }


    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the residuals
      fill_in_contribution_to_residuals_axisym_fsi_traction(residuals);

      // Derivatives w.r.t. external data
      fill_in_jacobian_from_external_interaction_by_fd(residuals, jacobian);
    }

    /// Return the ratio of the stress scales used to non-dimensionalise
    /// the fluid and elasticity equations. E.g.
    /// \f$ Q = (\omega a)^2 \rho/E \f$, i.e. the
    /// ratio between the inertial fluid stress and the solid's elastic
    /// modulus E.
    const double& q() const
    {
      return *Q_pt;
    }

    /// Return a pointer the ratio of stress scales used to
    /// non-dimensionalise the fluid and solid equations.
    double*& q_pt()
    {
      return Q_pt;
    }


    /// Output function
    void output(std::ostream& outfile)
    {
      /// Dummy
      unsigned nplot = 0;
      output(outfile, nplot);
    }

    /// Output function: Plot traction etc at Gauss points
    /// nplot is ignored.
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Dimension
      unsigned n_dim = this->nodal_dimension();

      // Get FSI parameter
      const double q_value = q();

      // Storage for traction (includes swirl!)
      Vector<double> traction(n_dim + 1);

      outfile << "ZONE\n";

      // Set the value of n_intpt
      unsigned n_intpt = integral_pt()->nweight();

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        Vector<double> s_int(n_dim - 1);
        s_int[0] = integral_pt()->knot(ipt, 0);

        // Get the outer unit normal
        Vector<double> interpolated_normal(n_dim);
        outer_unit_normal(ipt, interpolated_normal);

        // Boundary coordinate
        Vector<double> zeta(1);
        interpolated_zeta(s_int, zeta);

        // Get bulk element for traction
        NAVIER_STOKES_BULK_ELEMENT* ext_el_pt =
          dynamic_cast<NAVIER_STOKES_BULK_ELEMENT*>(
            this->external_element_pt(0, ipt));
        Vector<double> s_ext(this->external_element_local_coord(0, ipt));

        // Get traction from bulk element (on fluid scale)
        ext_el_pt->traction(s_ext, interpolated_normal, traction);

        outfile << ext_el_pt->interpolated_x(s_ext, 0) << " "
                << ext_el_pt->interpolated_x(s_ext, 1) << " "
                << q_value * traction[0] << " " << q_value * traction[1] << " "
                << q_value * traction[2] << " " << interpolated_normal[0] << " "
                << interpolated_normal[1] << " " << zeta[0] << std::endl;
      }
    }

    /// C_style output function
    void output(FILE* file_pt)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FaceGeometry<ELASTICITY_BULK_ELEMENT>::output(file_pt, n_plot);
    }
  };

  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //=================================================================
  /// Static default value for the ratio of stress scales
  /// used in the fluid and solid equations (default is 1.0)
  //=================================================================
  template<class ELASTICITY_BULK_ELEMENT, class NAVIER_STOKES_BULK_ELEMENT>
  double FSIAxisymmetricLinearElasticityTractionElement<
    ELASTICITY_BULK_ELEMENT,
    NAVIER_STOKES_BULK_ELEMENT>::Default_Q_Value = 1.0;


  //=====================================================================
  /// Return the residuals
  //=====================================================================
  template<class ELASTICITY_BULK_ELEMENT, class NAVIER_STOKES_BULK_ELEMENT>
  void FSIAxisymmetricLinearElasticityTractionElement<
    ELASTICITY_BULK_ELEMENT,
    NAVIER_STOKES_BULK_ELEMENT>::
    fill_in_contribution_to_residuals_axisym_fsi_traction(
      Vector<double>& residuals)
  {
    // Find out how many nodes there are
    unsigned n_node = nnode();

#ifdef PARANOID
    // Find out how many positional dofs there are
    unsigned n_position_type = this->nnodal_position_type();
    if (n_position_type != 1)
    {
      throw OomphLibError("LinearElasticity is not yet implemented for more "
                          "than one position type.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Cache the nodal indices at which the displacement components are stored
    unsigned u_nodal_index[n_dim + 1];
    for (unsigned i = 0; i < n_dim + 1; i++)
    {
      u_nodal_index[i] = this->U_index_axisym_fsi_traction[i];
    }

    // Integer to hold the local equation number
    int local_eqn = 0;

    // Set up memory for the shape functions
    Shape psi(n_node);
    DShape dpsids(n_node, n_dim - 1);

    // Get FSI parameter
    const double q_value = q();

    // Storage for traction
    Vector<double> traction(3);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Only need to call the local derivatives
      dshape_local_at_knot(ipt, psi, dpsids);

      // Calculate the coordinates
      Vector<double> interpolated_x(n_dim, 0.0);

      // Also calculate the surface tangent vectors
      DenseMatrix<double> interpolated_A(n_dim - 1, n_dim, 0.0);

      // Calculate displacements and derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over directions
        for (unsigned i = 0; i < n_dim; i++)
        {
          // Calculate the Eulerian coords
          const double x_local = nodal_position(l, i);
          interpolated_x[i] += x_local * psi(l);

          // Loop over LOCAL derivative directions, to calculate the tangent(s)
          for (unsigned j = 0; j < n_dim - 1; j++)
          {
            interpolated_A(j, i) += x_local * dpsids(l, j);
          }
        }
      }

      // Now find the local metric tensor from the tangent Vectors
      DenseMatrix<double> A(n_dim - 1);
      for (unsigned i = 0; i < n_dim - 1; i++)
      {
        for (unsigned j = 0; j < n_dim - 1; j++)
        {
          // Initialise surface metric tensor to zero
          A(i, j) = 0.0;

          // Take the dot product
          for (unsigned k = 0; k < n_dim; k++)
          {
            A(i, j) += interpolated_A(i, k) * interpolated_A(j, k);
          }
        }
      }

      // Get the outer unit normal
      Vector<double> interpolated_normal(n_dim);
      outer_unit_normal(ipt, interpolated_normal);

      // Find the determinant of the metric tensor
      double Adet = 0.0;
      switch (n_dim)
      {
        case 2:
          Adet = A(0, 0);
          break;
        case 3:
          Adet = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
          break;
        default:
          throw OomphLibError(
            "Wrong dimension in TimeHarmonicLinElastLoadedByPressureElement",
            "TimeHarmonicLinElastLoadedByPressureElement::fill_in_contribution_"
            "to_residuals()",
            OOMPH_EXCEPTION_LOCATION);
      }

      // Premultiply the weights and the square-root of the determinant of
      // the metric tensor
      double W = w * sqrt(Adet);

      // Get bulk element for traction
      NAVIER_STOKES_BULK_ELEMENT* ext_el_pt =
        dynamic_cast<NAVIER_STOKES_BULK_ELEMENT*>(
          this->external_element_pt(0, ipt));
      Vector<double> s_ext(this->external_element_local_coord(0, ipt));

      // Get traction from bulk element
      ext_el_pt->traction(s_ext, interpolated_normal, traction);

      // Loop over the test functions, nodes of the element
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over the displacement components
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          // Equation number
          local_eqn = this->nodal_local_eqn(l, u_nodal_index[i]);
          /*IF it's not a boundary condition*/
          if (local_eqn >= 0)
          {
            // Add the loading terms (multiplied by fsi scaling factor)
            // to the residuals
            residuals[local_eqn] -=
              q_value * traction[i] * psi(l) * interpolated_x[0] * W;
          }
        }
      } // End of loop over shape functions

    } // End of loop over integration points
  }


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


} // namespace oomph

#endif
