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

#endif
