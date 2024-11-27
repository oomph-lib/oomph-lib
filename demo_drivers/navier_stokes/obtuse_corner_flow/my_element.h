#ifndef my_element_header
#define my_element_header

#include "navier_stokes.h"
#include "store_error_element.h"

namespace oomph
{
  class MyElement : public StoreErrorElement<TTaylorHoodElement<2>>
  {
  private:
    unsigned Region_id;

  public:
    void set_region_id(const unsigned& value)
    {
      Region_id = value;
    }

    const unsigned get_region_id()
    {
      return Region_id;
    }

    void output(std::ostream& outfile, const unsigned& nplot)
    {
      // Vector of local coordinates
      Vector<double> s(2);

      // Tecplot header info
      outfile << tecplot_zone_string(nplot);

      // Loop over plot points
      unsigned num_plot_points = nplot_points(nplot);
      for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
      {
        // Get local coordinates of plot point
        get_s_plot(iplot, nplot, s);

        // Coordinates
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << interpolated_x(s, i) << " ";
        }

        // Velocities
        for (unsigned i = 0; i < 2; i++)
        {
          outfile << this->interpolated_u_nst(s, i) << " ";
        }

        // Pressure
        outfile << interpolated_p_nst(s) << " ";

        // Output the stored error
        outfile << this->get_error() << " ";

        // Output the element size
        outfile << this->size() << " ";

        // Output the continuity residual
        const double cont_res =
          interpolated_dudx_nst(s, 0, 0) + interpolated_dudx_nst(s, 1, 1);
        outfile << cont_res << " ";

        outfile << std::endl;
      }
      outfile << std::endl;

      // Write tecplot footer (e.g. FE connectivity lists)
      write_tecplot_zone_footer(outfile, nplot);
    }
  };


  //=======================================================================
  /// Face geometry of the 2D Taylor_Hood elements
  //=======================================================================
  template<>
  class FaceGeometry<MyElement> : public virtual TElement<1, 3>
  {
  public:
    /// Constructor: Call constructor of base
    FaceGeometry() : TElement<1, 3>() {}
  };
} // namespace oomph

#endif
