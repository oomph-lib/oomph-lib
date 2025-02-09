#include "parameter_struct.h"

using namespace oomph;

int main()
{
  Params params = create_parameters_from_file("parameters.dat");

  if (params.slip_length != 0.1)
  {
    return 1;
  }
  if (params.n_radial != 20)
  {
    return 1;
  }
  if (params.n_azimuthal != 15)
  {
    return 1;
  }
  if (params.geometric_base != 1.05)
  {
    return 1;
  }
  if (params.sector_radius != 0.9)
  {
    return 1;
  }
  if (params.sector_angle != 135)
  {
    return 1;
  }
  if (params.reynolds_number != 0)
  {
    return 1;
  }
  if (params.strouhal_reynolds_number != 0)
  {
    return 1;
  }
  if (params.inner_radius != 0.5)
  {
    return 1;
  }

  return 0;
};
