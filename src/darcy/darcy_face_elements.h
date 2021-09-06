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
// Header file for elements that are used to apply surface loads to
// the Darcy equations

#ifndef OOMPH_DARCY_FACE_ELEMENTS_HEADER
#define OOMPH_DARCY_FACE_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// OOMPH-LIB headers
#include "generic/Qelements.h"

namespace oomph
{
  //=======================================================================
  /// Namespace containing the zero pressure function for Darcy pressure
  /// elements
  //=======================================================================
  namespace DarcyFaceElementHelper
  {
    //=======================================================================
    /// Default load function (zero pressure)
    //=======================================================================
    void Zero_pressure_fct(const double& time,
                           const Vector<double>& x,
                           const Vector<double>& N,
                           double& load)
    {
      load = 0.0;
    }

  } // namespace DarcyFaceElementHelper


  //======================================================================
  /// A class for elements that allow the imposition of an applied pressure
  /// in the Darcy equations.
  /// The geometrical information can be read from the FaceGeometry<ELEMENT>
  /// class and and thus, we can be generic enough without the need to have
  /// a separate equations class.
  //======================================================================
  template<class ELEMENT>
  class DarcyFaceElement : public virtual FaceGeometry<ELEMENT>,
                           public virtual FaceElement
  {
  protected:
    /// \short Pointer to an imposed pressure function. Arguments:
    /// Eulerian coordinate; outer unit normal; applied pressure.
    /// (Not all of the input arguments will be required for all specific load
    /// functions but the list should cover all cases)
    void (*Pressure_fct_pt)(const double& time,
                            const Vector<double>& x,
                            const Vector<double>& n,
                            double& result);


    /// \short Get the pressure value: Pass number of integration point (dummy),
    /// Eulerlian coordinate and normal vector and return the pressure
    /// (not all of the input arguments will be required for all specific load
    /// functions but the list should cover all cases). This function is virtual
    /// so it can be overloaded for FSI.
    virtual void get_pressure(const double& time,
                              const unsigned& intpt,
                              const Vector<double>& x,
                              const Vector<double>& n,
                              double& pressure)
    {
      Pressure_fct_pt(time, x, n, pressure);
    }


    /// \short Helper function that actually calculates the residuals
    // This small level of indirection is required to avoid calling
    // fill_in_contribution_to_residuals in fill_in_contribution_to_jacobian
    // which causes all kinds of pain if overloading later on
    void fill_in_contribution_to_residuals_darcy_face(
      Vector<double>& residuals);


  public:
    /// \short Constructor, which takes a "bulk" element and the value of the
    /// index and its limit
    DarcyFaceElement(FiniteElement* const& element_pt, const int& face_index)
      : FaceGeometry<ELEMENT>(), FaceElement()
    {
#ifdef PARANOID
      {
        // Check that the element is not a refineable 3d element
        ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(element_pt);
        // If it's three-d
        if (elem_pt->dim() == 3)
        {
          // Is it refineable
          RefineableElement* ref_el_pt =
            dynamic_cast<RefineableElement*>(elem_pt);
          if (ref_el_pt != 0)
          {
            if (this->has_hanging_nodes())
            {
              throw OomphLibError("This flux element will not work correctly "
                                  "if nodes are hanging\n",
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
      }
#endif

      // Attach the geometrical information to the element. N.B. This function
      // also assigns nbulk_value from the required_nvalue of the bulk element
      element_pt->build_face_element(face_index, this);

      // Zero pressure
      Pressure_fct_pt = &DarcyFaceElementHelper::Zero_pressure_fct;
    }


    /// Reference to the pressure function pointer
    void (*&pressure_fct_pt())(const double& time,
                               const Vector<double>& x,
                               const Vector<double>& n,
                               double& pressure)
    {
      return Pressure_fct_pt;
    }


    /// Return the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_contribution_to_residuals_darcy_face(residuals);
    }


    /// Fill in contribution from Jacobian
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the residuals
      fill_in_contribution_to_residuals_darcy_face(residuals);
    }

    /// Specify the value of nodal zeta from the face geometry
    /// \short The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    /// \short Output function
    void output(std::ostream& outfile)
    {
      FaceGeometry<ELEMENT>::output(outfile);
    }

    /// \short Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FaceGeometry<ELEMENT>::output(outfile, n_plot);
    }

    /// \short C_style output function
    void output(FILE* file_pt)
    {
      FaceGeometry<ELEMENT>::output(file_pt);
    }

    /// \short C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FaceGeometry<ELEMENT>::output(file_pt, n_plot);
    }


    /// \short Compute pressure value at specified local coordinate
    /// Should only be used for post-processing; ignores dependence
    /// on integration point!
    void pressure(const double& time,
                  const Vector<double>& s,
                  double& pressure);
  };

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=====================================================================
  /// Compute pressure value at specified local coordinate
  /// Should only be used for post-processing; ignores dependence
  /// on integration point!
  //=====================================================================
  template<class ELEMENT>
  void DarcyFaceElement<ELEMENT>::pressure(const double& time,
                                           const Vector<double>& s,
                                           double& pressure)
  {
    unsigned n_dim = this->nodal_dimension();

    // Position vector
    Vector<double> x(n_dim);
    interpolated_x(s, x);

    // Outer unit normal
    Vector<double> unit_normal(n_dim);
    outer_unit_normal(s, unit_normal);

    // Dummy
    unsigned ipt = 0;

    // Pressure value
    get_pressure(time, ipt, x, unit_normal, pressure);
  }


  //=====================================================================
  /// Return the residuals for the DarcyFaceElement equations
  //=====================================================================
  template<class ELEMENT>
  void DarcyFaceElement<ELEMENT>::fill_in_contribution_to_residuals_darcy_face(
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
      throw OomphLibError("Darcy equations are not yet implemented for more "
                          "than one position type",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find out the dimension of the node
    unsigned n_dim = this->nodal_dimension();

    // Set the pointer to the bulk element
    ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

    unsigned n_q_basis = bulk_el_pt->nq_basis();
    unsigned n_q_basis_edge = bulk_el_pt->nq_basis_edge();

    // Integer to hold the local equation number
    int local_eqn = 0;

    // Set up memory for the shape functions
    // Note that in this case, the number of lagrangian coordinates is always
    // equal to the dimension of the nodes
    Shape psi(n_node);
    DShape dpsids(n_node, n_dim - 1);

    Shape q_basis(n_q_basis, n_dim);

    // Set the value of n_intpt
    unsigned n_intpt = integral_pt()->nweight();

    // Storage for the local coordinates
    Vector<double> s_face(n_dim - 1, 0.0), s_bulk(n_dim, 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Only need to call the local derivatives
      dshape_local_at_knot(ipt, psi, dpsids);

      // Assign values of s in FaceElement and local coordinates in bulk element
      for (unsigned i = 0; i < n_dim - 1; i++)
      {
        s_face[i] = integral_pt()->knot(ipt, i);
      }

      s_bulk = local_coordinate_in_bulk(s_face);

      // Get the q basis at bulk local coordinate s_bulk, corresponding to
      // face local coordinate s_face
      bulk_el_pt->get_q_basis(s_bulk, q_basis);

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

      // Now find the local metric tensor from the tangent vectors
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
          throw OomphLibError("Wrong dimension in DarcyFaceElement",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }

      // Premultiply the weights and the square-root of the determinant of
      // the metric tensor
      double W = w * sqrt(Adet);

      // Now calculate the load
      double pressure;
      get_pressure(time, ipt, interpolated_x, interpolated_normal, pressure);

      // Loop over the q edge test functions only (internal basis functions
      // have zero normal component on the boundary)
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        local_eqn = this->nodal_local_eqn(1, bulk_el_pt->q_edge_index(l));

        /*IF it's not a boundary condition*/
        if (local_eqn >= 0)
        {
          // Loop over the displacement components
          for (unsigned i = 0; i < n_dim; i++)
          {
            // Add the loading terms to the residuals
            residuals[local_eqn] +=
              pressure * q_basis(l, i) * interpolated_normal[i] * W;
          }
        } // End of if not boundary condition
      } // End of loop over shape functions
    } // End of loop over integration points
  }

} // namespace oomph

#endif
