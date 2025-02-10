#include "parameter_struct.h"

namespace oomph
{
  Params create_parameters_from_file(const std::string& filename)
  {
    Params params;

    std::ifstream parameter_filestream(filename);
    std::string input_string;

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.slip_length = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.n_radial = std::stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.n_azimuthal = std::stoi(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.geometric_base = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.sector_radius = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.sector_angle = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.reynolds_number = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.strouhal_reynolds_number = std::stod(input_string);

    getline(parameter_filestream, input_string, '#');
    parameter_filestream.ignore(80, '\n');
    params.inner_radius = std::stod(input_string);

    return params;
  }
}; // namespace oomph
