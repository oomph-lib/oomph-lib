#ifndef LINEARISED_CONTACT_ANGLE_ELEMENTS_HEADER
#define LINEARISED_CONTACT_ANGLE_ELEMENTS_HEADER

#include "linearised_elastic_axisym_fluid_interface_element.h"

namespace oomph
{
  template<class FREE_SURFACE_ELEMENT>
  class LinearisedContactAngleElement : public virtual FaceElement,
                                        public virtual PointElement
  {
  private:
    /// Pointer to the desired value of the capillary number
    double* Ca_pt;

    /// Pointer to the desired value of the contact angle
    double* Contact_angle_pt;

    /// Pointer to azimuthal mode number k in e^ik(theta) decomposition
    int* Azimuthal_Mode_Number_pt;

    /// Function pointer to a wall unit normal function. Returns the
    /// unit normal on the wall, at the specified Eulerian coordinate.
    typedef void (*WallUnitNormalFctPt)(const Vector<double>& x,
                                        Vector<double>& unit_normal);

    /// Pointer to a wall normal function that returns
    /// the wall unit normal as a function of position in global
    /// Eulerian coordinates.
    WallUnitNormalFctPt Wall_unit_normal_fct_pt;

    Vector<double> Base_tangent;
    Vector<double> Base_normal;

    /// Index at which the i-th velocity component is stored in the
    /// element's nodes
    Vector<unsigned> U_index_velocity_decomposed;

    Vector<Vector<unsigned>> U_index_displacement_decomposed;

    int Face_index;

  public:
    LinearisedContactAngleElement(FiniteElement* const& element_pt,
                                  const int& face_index)
      : Ca_pt(0),
        Contact_angle_pt(0),
        Azimuthal_Mode_Number_pt(0),
        Wall_unit_normal_fct_pt(0),
        Face_index(face_index)
    {
      // Attach the geometrical information to the element, by
      // making the face element from the bulk element
      element_pt->build_face_element(face_index, this);

      FREE_SURFACE_ELEMENT* cast_bulk_element_pt =
        dynamic_cast<FREE_SURFACE_ELEMENT*>(element_pt);

      U_index_velocity_decomposed.resize(6);
      for (int n = 0; n < 6; n++)
      {
        U_index_velocity_decomposed[n] =
          cast_bulk_element_pt->U_index_interface[n];
      }

      U_index_displacement_decomposed.resize(2);
      for (int i = 0; i < 2; i++)
      {
        U_index_displacement_decomposed[i].resize(2);
        for (int j = 0; j < 2; j++)
        {
          U_index_displacement_decomposed[i][j] =
            cast_bulk_element_pt->u_index_linear_elasticity(
              this->bulk_node_number(0), i, j);
        }
      }

      Base_normal = cast_bulk_element_pt->normal(Face_index);
      Base_normal.resize(3);
      Base_tangent = cast_bulk_element_pt->tangent(Face_index);
      Base_tangent.resize(3);
    }

    /// Access function to the pointer specifying the capillary number
    double*& ca_pt()
    {
      return Ca_pt;
    }

    /// The value of the Capillary number
    const double& ca() const
    {
#ifdef PARANOID
      if (Ca_pt != 0)
      {
#endif
        return *Ca_pt;
#ifdef PARANOID
      }
      else
      {
        throw OomphLibError("Capillary number has not been set",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Access function to the pointer specifying the prescribed contact angle
    double*& contact_angle_pt()
    {
      return Contact_angle_pt;
    }

    /// Return value of the contact angle
    double& contact_angle()
    {
#ifdef PARANOID
      if (Contact_angle_pt == 0)
      {
        std::string error_message = "Contact angle not set\n";
        error_message +=
          "Please use FluidInterfaceBoundingElement::set_contact_angle()\n";
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Contact_angle_pt;
    }

    /// Azimuthal mode number k in e^ik(theta) decomposition
    const int& azimuthal_mode_number() const
    {
      return *Azimuthal_Mode_Number_pt;
    }

    /// Pointer to azimuthal mode number k in e^ik(theta) decomposition
    int*& azimuthal_mode_number_pt()
    {
      return Azimuthal_Mode_Number_pt;
    }

    /// Function that returns the unit normal of the bounding wall
    /// directed out of the fluid
    void wall_unit_normal(const Vector<double>& x, Vector<double>& normal)
    {
#ifdef PARANOID
      if (Wall_unit_normal_fct_pt)
      {
#endif
        (*Wall_unit_normal_fct_pt)(x, normal);
#ifdef PARANOID
      }
      else
      {
        throw OomphLibError("Wall unit normal fct has not been set",
                            "FluidInterfaceBoundingElement::wall_unit_normal()",
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }

    /// Access function: Pointer to wall unit normal function
    WallUnitNormalFctPt& wall_unit_normal_fct_pt()
    {
      return Wall_unit_normal_fct_pt;
    }

    /// Calculate the residuals
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Add the residual contributions using a dummy matrix
      fill_in_generic_residual_contribution_linear_contact_angle(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    }

    /// Calculate the generic residuals contribution
    void fill_in_generic_residual_contribution_linear_contact_angle(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag)
    {
      // Working in 3 dimensions
      const unsigned n_dim = 3;

      // The number of dimensions in the reduced system is 2
      const unsigned n_reduced_dim = 2;

      // Get the capillary number
      const double ca_local = ca();

      // Get the contact angle
      const double contact_angle_local = contact_angle();

      // Dummy local coordinate, of size zero
      Vector<double> s_local(0);

      // Storage for the coordinate
      Vector<double> x(n_dim, 0.0);

      // Get the x coordinate
      this->interpolated_x(s_local, x);

      // Outer unit normal to the wall
      Vector<double> wall_normal(n_dim);

      // Get the unit normal to the wall
      wall_unit_normal(x, wall_normal);

      // Need to find the current outer normal from the surface
      // which does not necessarily correspond to an imposed angle.
      // It is whatever it is...
      Vector<double> m(n_dim, 0.0);

      // Get the shape function at the single integration point
      Shape psif(nnode());
      unsigned ipt = 0;
      this->shape_at_knot(ipt, psif);
      // There is only one shape function
      unsigned l = 0;
      const double psif_ = psif(l);

      // Find the dot product of the two vectors
      Vector<double> dot(n_dim, 0.0);
      // Cosine = 0/Sine = 1
      for (unsigned j = 0; j < n_dim; j++)
      {
        // m = 1/J * d\v{R}/ds
        m = dynamic_cast<FREE_SURFACE_ELEMENT*>(
              this->bulk_element_pt())
              ->displacement_gradient(Face_index, j);
        for (unsigned i = 0; i < n_reduced_dim; i++)
        {
          dot[j] += m[i] * Base_normal[i];
        }
      }

      Vector<double> e_theta(3, 0.0);
      e_theta[2] = 1.0;

      // Just add the appropriate contribution to the momentum equations
      // This will, of course, not be added if the equation is pinned
      // (no slip)
      // cosine/sine
      for (unsigned i = 0; i < 2; i++)
      {
        // Compute the inner expression
        Vector<double> expression(n_dim, 0.0);
        for (unsigned j = 0; j < n_dim; j++)
        {
          expression[j] -= dot[i] * Base_tangent[j];

          if (i == 0)
          {
            expression[j] -= azimuthal_mode_number() * dot[1 - i] * e_theta[j];
          }
          else
          {
            expression[j] += azimuthal_mode_number() * dot[1 - i] * e_theta[j];
          }
        }
        Vector<double> projected_expression = project_onto_wall(x, expression);

        // Compute A
        Vector<double> vector_A = project_onto_wall(x, Base_normal);
        double mag_A = 0;
        for (unsigned d = 0; d < n_dim; d++)
        {
          vector_A[d] = vector_A[d] * cos(contact_angle_local) +
                        wall_normal[d] * sin(contact_angle_local);
          mag_A += vector_A[d] * vector_A[d];
        }
        mag_A = pow(mag_A, 0.5);

        // Horizontal contribution
        for (unsigned d = 0; d < n_dim; d++)
        {
          int local_eqn =
            nodal_local_eqn(0, this->U_index_velocity_decomposed[i + 2 * d]);
          if (local_eqn >= 0)
          {
            residuals[local_eqn] -= 1.0 / (ca_local * mag_A) *
                                    projected_expression[d] *
                                    cos(contact_angle_local) * psif_;
          }
        }
      }
    }

    Vector<double> project_onto_wall(Vector<double> x,
                                     Vector<double> expression)
    {
      // Working in 3 dimensions
      const unsigned n_dim = 3;

      // Outer unit normal to the wall
      Vector<double> wall_normal(n_dim);

      // Get the unit normal to the wall
      wall_unit_normal(x, wall_normal);

      // Create a n_dim x n_dim matrix of zeros
      DenseDoubleMatrix A(n_dim, n_dim, 0.0);
      for (unsigned i = 0; i < n_dim; i++)
      {
        A(i, i) = 1.0;
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        for (unsigned j = 0; j < n_dim; j++)
        {
          A(i, j) -= wall_normal[i] * wall_normal[j];
        }
      }

      Vector<double> result(n_dim, 0.0);
      for (unsigned i = 0; i < n_dim; i++)
      {
        for (unsigned j = 0; j < n_dim; j++)
        {
          result[j] += A(i, j) * expression[i];
        }
      }

      return result;
    }

    void output(std::ostream& outfile)
    {
      // Working in 3 dimensions
      const unsigned n_dim = 3;

      // The number of dimensions in the reduced system is 2
      const unsigned n_reduced_dim = 2;

      // Dummy local coordinate, of size zero
      Vector<double> s_local(0);

      // Storage for the coordinate
      Vector<double> x(n_dim, 0.0);

      // Get the x coordinate
      this->interpolated_x(s_local, x);

      // Outer unit normal to the wall
      Vector<double> wall_normal(n_dim);

      // Get the unit normal to the wall
      wall_unit_normal(x, wall_normal);

      // Need to find the current outer normal from the surface
      // which does not necessarily correspond to an imposed angle.
      // It is whatever it is...
      Vector<double> m(n_dim, 0.0);

      for (unsigned i = 0; i < n_reduced_dim; i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << Base_normal[i] << " ";
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << Base_tangent[i] << " ";
      }
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << wall_normal[i] << " ";
      }
      for (unsigned i = 0; i < 2; i++)
      {
        m = dynamic_cast<FREE_SURFACE_ELEMENT*>(
              this->bulk_element_pt())
              ->displacement_gradient(Face_index, i);

        for (unsigned j = 0; j < 2; j++)
        {
          outfile << m[j] << " ";
        }
      }
      outfile << ca() << " ";
      outfile << contact_angle() << " ";
      outfile << azimuthal_mode_number();
      outfile << endl;
    }
  };
} // namespace oomph
#endif
