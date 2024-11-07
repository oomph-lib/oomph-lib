#ifndef PROJECTABLE_TTAYLOR_HOOD_ELEMENTS_HEADER
#define PROJECTABLE_TTAYLOR_HOOD_ELEMENTS_HEADER

#include "solid/solid_elements.h"
#include "navier_stokes.h"
#include "fluid_interface/constrained_volume_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  class TTaylorHoodPVDElement
    : public virtual PseudoSolidNodeUpdateElement<TTaylorHoodElement,
                                                  TPVDElement<2, 3>>,
      public virtual DebugJacobianSolidFiniteElement
  {
  public:
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      NavierStokesEquations::fill_in_contribution_to_jacobian(residuals,
                                                              jacobian);

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
      NavierStokesEquations::fill_in_contribution_to_jacobian_and_mass_matrix(
        residuals, jacobian, mass_matrix);

      PVDEquations<2>::fill_in_contribution_to_jacobian(residuals, jacobian);

      //   Now fill in the off-diagonal entries (the shape derivatives),
      fill_in_shape_derivatives(jacobian);
    }
  };


  class ProjectableTTaylorHoodPVDElement
    : public virtual ProjectableTaylorHoodElement<TTaylorHoodPVDElement>,
      public virtual DebugJacobianSolidFiniteElement
  {
  private:
    double Error;

  public:
    ProjectableTTaylorHoodPVDElement() : Error(0.0) {}

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
          outfile << this->interpolated_u_nst(s, i) << " ";
        }

        // Pressure
        outfile << this->interpolated_p_nst(s) << " ";

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
  class FaceGeometry<ProjectableTTaylorHoodPVDElement>
    : public virtual SolidTElement<1, 3>
  {
  public:
    FaceGeometry() : SolidTElement<1, 3>() {}
  };

  // Ensure that the FaceGeometry of the new element has been set up
  template<>
  class FaceGeometry<FaceGeometry<ProjectableTTaylorHoodPVDElement>>
    : public virtual SolidPointElement
  {
  public:
    FaceGeometry() : SolidPointElement() {}
  };
} // namespace oomph

#endif
