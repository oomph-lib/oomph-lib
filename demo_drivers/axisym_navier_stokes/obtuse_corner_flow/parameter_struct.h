#ifndef PARAMETER_STRUCT_HEADER
#define PARAMETER_STRUCT_HEADER

#include <iostream>
#include <fstream>
#include <string>

namespace oomph
{
  struct Params
  {
    double slip_length = 1.0;
    double geometric_base = 1.2;
    double sector_radius = 1;
    double sector_angle = 90;
    double reynolds_number = 0;
    double strouhal_reynolds_number = 0;
    unsigned n_radial = 10;
    unsigned n_azimuthal = 10;
    double inner_radius = 0;
  };

  Params create_parameters_from_file(const std::string& filename);
}; // namespace oomph

#endif
