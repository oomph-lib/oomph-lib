#ifndef SINGULAR_AXISYM_FLUID_TRACTION_ELEMENTS_HEADER
#define SINGULAR_AXISYM_FLUID_TRACTION_ELEMENTS_HEADER

#include "axisym_navier_stokes.h"
#include "singular_navier_stokes_solution_elements.h"

namespace oomph
{

  template<class ELEMENT>
  class SingularAxisymNavierStokesTractionElement
    : public virtual AxisymmetricNavierStokesTractionElement<ELEMENT>
  {
  public:
    SingularAxisymNavierStokesTractionElement(
      FiniteElement* const& element_pt,
      const int& face_index,
      Data* const& singular_scaling_data_pt)
      : AxisymmetricNavierStokesTractionElement<ELEMENT>(element_pt, face_index)
    {
      this->add_external_data(singular_scaling_data_pt);
    }

    void fill_in_contribution_to_dresiduals_dparameter(
      double* const& parameter_pt, Vector<double>& dres_dparam)
    {
    }

  protected:
    void get_traction(const double& time,
                      const Vector<double>& x,
                      const Vector<double>& n,
                      Vector<double>& result)
    {
      // Dummy integration point
      unsigned ipt = 0;
      AxisymmetricNavierStokesTractionElement<ELEMENT>::get_traction(
        time, ipt, x, n, result);
      const unsigned N = result.size();
      for (unsigned i = 0; i < N; i++)
      {
        result[i] *= this->external_data_pt(0)->value(0);
      }
    }

  public:
    /// This function returns the residuals and the jacobian
    inline void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                                 DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      AxisymmetricNavierStokesTractionElement<ELEMENT>::
        fill_in_contribution_to_residuals_axisymmetric_nst_traction(residuals);

      this->fill_in_jacobian_from_external_by_fd(residuals, jacobian, false);
    }

    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      fill_in_contribution_to_jacobian(residuals, jacobian);
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
      // Get continuous time from timestepper of first node
      const double local_time =
        this->node_pt(0)->time_stepper_pt()->time_pt()->time();

      // Loop over plot points
      unsigned num_plot_points = this->nplot_points(n_plot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get the nodal dimension
        unsigned n_dim = this->nodal_dimension();

        // Initialise local coordinates
        Vector<double> s(n_dim - 1);

        // Get local coordinates of plot point
        this->get_s_plot(iplot, n_plot, s);

        // Declare output variables
        Vector<double> x(n_dim);
        this->interpolated_x(s, x);

        Vector<double> n(n_dim);
        this->outer_unit_normal(s, n);

        Vector<double> traction(n_dim);
        get_traction(local_time, x, n, traction);

        // Output the time
        outfile << local_time << " ";
        // Output the position
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << x[i] << " ";
        }
        // Output the normal
        for (unsigned i = 0; i < n_dim; i++)
        {
          outfile << n[i] << " ";
        }
        // Output traction, skipping the comma at the end of line
        for (unsigned i = 0; i < n_dim - 1; i++)
        {
          outfile << traction[i] << " ";
        }
        outfile << traction[n_dim - 1];
        // End of line
        outfile << std::endl;
      }
    }
  };
} // namespace oomph
#endif
