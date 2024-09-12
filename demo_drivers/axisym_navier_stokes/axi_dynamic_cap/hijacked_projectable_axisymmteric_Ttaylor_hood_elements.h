#ifndef HIJACKED_PROJECTABLE_AXISYMMTERIC_TTAYLOR_HOOD_ELEMENTS_HEADER
#define HIJACKED_PROJECTABLE_AXISYMMTERIC_TTAYLOR_HOOD_ELEMENTS_HEADER

#include "solid/solid_elements.h"
#include "navier_stokes.h"
#include "fluid_interface/constrained_volume_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  template<class ELEMENT>
  class DebugElasticAxisymmetricVolumeConstraintBoundingElement
    : public virtual ElasticAxisymmetricVolumeConstraintBoundingElement<
        ELEMENT>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    DebugElasticAxisymmetricVolumeConstraintBoundingElement(
      FiniteElement* const& element_pt, const int& face_index)
      : ElasticAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(element_pt,
                                                                    face_index)
    {
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Compute the element's residual Vector and the jacobian matrix
    /// Virtual function can be overloaded by hanging-node version
    void fill_in_contribution_to_djacobian_dparameter(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam)
    {
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }
  };

  template<class ELEMENT>
  class ImposeImpenetrabilityElementWithOutput
    : public ImposeImpenetrabilityElement<ELEMENT>
  {
  public:
    // Use base class constructor
    using ImposeImpenetrabilityElement<ELEMENT>::ImposeImpenetrabilityElement;

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Compute the element's residual Vector and the jacobian matrix
    /// Virtual function can be overloaded by hanging-node version
    void fill_in_contribution_to_djacobian_dparameter(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam)
    {
    }

    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      // Number of dimensions
      unsigned n_dim = this->nodal_dimension();

      // Find out how many nodes there are
      const unsigned n_node = this->nnode();

      // Set up memory for the shape functions
      Shape psi(n_node);

      // Local and global coordinates
      Vector<double> s(n_dim - 1);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);

        // Find the shape functions
        this->shape(s, psi);

        // Initialise to zero
        Vector<double> interpolated_x(n_dim);
        Vector<double> interpolated_u(n_dim + 1);
        double interpolated_lambda = 0.0;

        // Calculate stuff
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over directions
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psi[l];
          }

          for (unsigned i = 0; i < n_dim + 1; i++)
          {
            interpolated_u[i] +=
              this->nodal_value(l, this->u_index_nst(l, i)) * psi[l];
          }
          interpolated_lambda +=
            this->nodal_value(l, this->lagrange_multiplier_index(l)) * psi[l];
        }

        // Output x,y,..
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << interpolated_x[i] << " ";
        }

        // Output u,v,w
        for (unsigned i = 0; i < n_dim + 1; i++)
        {
          outfile << interpolated_u[i] << " ";
        }

        // Output the lagrange multiplier
        outfile << interpolated_lambda << std::endl;
      }
    }
  };

  template<class ELEMENT>
  class DebugImposeImpenetrabilityElement
    : public virtual ImposeImpenetrabilityElementWithOutput<ELEMENT>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    DebugImposeImpenetrabilityElement(FiniteElement* const& element_pt,
                                      const int& face_index,
                                      const unsigned& id = 0)
      : ImposeImpenetrabilityElementWithOutput<ELEMENT>(
          element_pt, face_index, id)
    {
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

    /// Compute the element's residual Vector and the jacobian matrix
    /// Virtual function can be overloaded by hanging-node version
    void fill_in_contribution_to_djacobian_dparameter(
      double* const& parameter_pt,
      Vector<double>& dres_dparam,
      DenseMatrix<double>& djac_dparam)
    {
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
    }
  };

  class AxisymmetricTTaylorHoodPVDElement
    : public virtual PseudoSolidNodeUpdateElement<
        AxisymmetricTTaylorHoodElement,
        TPVDElement<2, 3>>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      AxisymmetricNavierStokesEquations::fill_in_contribution_to_jacobian(
        residuals, jacobian);

      PVDEquations<2>::fill_in_contribution_to_jacobian(residuals, jacobian);

      // Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);

      // Fill in the external data entries by finite difference.
      fill_in_jacobian_from_external_by_fd(residuals, jacobian, true);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      AxisymmetricNavierStokesEquations::
        fill_in_contribution_to_jacobian_and_mass_matrix(
          residuals, jacobian, mass_matrix);

      PVDEquations<2>::fill_in_contribution_to_jacobian(residuals, jacobian);

      //   Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);
    }
  };


  class HijackedProjectableAxisymmetricTTaylorHoodPVDElement
    : public virtual // Hijacked<
      ProjectableAxisymmetricTaylorHoodElement<
        AxisymmetricTTaylorHoodPVDElement> //>
  {
  private:
    double Error;

  public:
    HijackedProjectableAxisymmetricTTaylorHoodPVDElement() : Error(0.0) {}

    // Set error value for post-processing
    void set_error(const double& error)
    {
      Error = error;
    }

    void get_error(double& error)
    {
      error = Error;
    }

    void pin()
    {
      for (unsigned i = 0; i < 6; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          this->node_pt(i)->pin(j);
        }
      }
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 1; j++)
        {
          this->node_pt(i)->pin(3 + j);
        }
      }
    }

    void pin_pressure()
    {
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 1; j++)
        {
          this->node_pt(i)->pin(3 + j);
        }
      }
    }

    // Overload output function
    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(2);
      Vector<double> xi(2);

      // Tecplot header info
      outfile << this->tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        this->get_s_plot(iplot, nplot, s);
        this->interpolated_xi(s, xi);

        std::streamsize old_precision = outfile.precision();
        outfile.precision(16);
        // Coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_x(s, i) << " ";
        }
        outfile.precision(old_precision);

        // Velocities
        for (unsigned i = 0; i < 3; i++)
        {
          outfile << this->interpolated_u_axi_nst(s, i) << " ";
        }

        // Pressure
        outfile << this->interpolated_p_axi_nst(s) << " ";

        // Lagrange coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_x(s, i) - xi[i] << " ";
        }

        // Error
        outfile << Error << " ";

        // Size
        outfile << this->size() << " ";

        outfile << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      this->write_tecplot_zone_footer(outfile, nplot);
    }

    // Overload output function
    void output(std::ostream& outfile)
    {
      unsigned nplot = 3;
      output(outfile, nplot);
    }
  }; // namespace >

  // Ensure that the FaceGeometry of the new element has been set up
  template<>
  class FaceGeometry<HijackedProjectableAxisymmetricTTaylorHoodPVDElement>
    : public virtual SolidTElement<1, 3>
  {
  public:
    FaceGeometry() : SolidTElement<1, 3>() {}
  };

  // Ensure that the FaceGeometry of the new element has been set up
  template<>
  class FaceGeometry<
    FaceGeometry<HijackedProjectableAxisymmetricTTaylorHoodPVDElement>>
    : public virtual SolidPointElement
  {
  public:
    FaceGeometry() : SolidPointElement() {}
  };
} // namespace oomph

#endif
