#ifndef my_element_header
#define my_element_header

#include "axisym_navier_stokes.h"
#include "store_error_element.h"

namespace oomph
{
  class MyElement : public StoreErrorElement<AxisymmetricTTaylorHoodElement>
  {
  public:
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
        for (unsigned i = 0; i < 3; i++)
        {
          outfile << interpolated_u_axi_nst(s, i) << " ";
        }

        // Pressure
        outfile << interpolated_p_axi_nst(s) << " ";

        // Output the stored error
        outfile << this->get_error() << " ";

        // Output the element size
        outfile << this->size() << " ";

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
