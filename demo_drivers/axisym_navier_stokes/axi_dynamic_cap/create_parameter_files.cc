#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "generic.h"

#include "parameters.h"

using namespace std;
using namespace oomph;

enum
{
  FAIL = -1,
  SUCCESS = 0,
};

// function declarations
vector<double> range_string_to_numbers(string range_str);
void create_parameter_files_over_range_for_variable_pt(
  const std::string& range_str,
  const std::string& folder_base,
  oomph::Parameters& parameters,
  const bool& script_is_creating_folders,
  std::ofstream& filestream,
  double* const variable_pt);

// function definitions
void create_parameter_files_over_range_for_variable_pt(
  const std::string& range_str,
  const std::string& folder_base,
  oomph::Parameters& parameters,
  const bool& script_is_creating_folders,
  std::ofstream& filestream,
  double* const variable_pt)
{
  vector<double> values = range_string_to_numbers(range_str);
  double lower = values[0];
  double upper = values[1];
  int N = int(values[2]);
  if (N < 0) std::cout << "Negative number of points in range" << std::endl;

  for (int n = 0; n < N; n++)
  {
    parameters.dir_name = folder_base + "/" + std::to_string(n);
    // if (script_is_creating_folders)
    //{
    const int success_flag = mkdir(parameters.dir_name.c_str(), 0777);
    if (success_flag == FAIL)
      std::cout << "Failed to create directory: " << parameters.dir_name << std::endl;
    //}

    if (CommandLineArgs::command_line_flag_has_been_set("--log"))
    {
      *variable_pt =
        lower * exp(double(n) / (double(N) - 1.0) * log(upper / lower));
    }
    else
    {
      *variable_pt = lower + double(n) / double(N - 1) * (upper - lower);
    }

    filestream.open(parameters.dir_name + "/parameters.dat");
    parameters.doc(filestream);
    filestream.close();
  }
}

void create_parameter_files_over_range_for_angle_pt(
  const std::string& range_str,
  const std::string& folder_base,
  oomph::Parameters& parameters,
  const bool& script_is_creating_folders,
  std::ofstream& filestream,
  double* const variable_pt)
{
  vector<double> values = range_string_to_numbers(range_str);
  double lower = values[0];
  double upper = values[1];
  int N = int(values[2]);
  if (N < 0) std::cout << "Negative number of points in range" << std::endl;

  for (int n = 0; n < N; n++)
  {
    parameters.dir_name = folder_base + "/" + std::to_string(n);

    const int success_flag = mkdir(parameters.dir_name.c_str(), 0777);
    if (success_flag == FAIL)
      std::cout << "Failed to create directory: " << parameters.dir_name << std::endl;

    if (CommandLineArgs::command_line_flag_has_been_set("--log"))
    {
      *variable_pt = lower *
                     exp(double(n) / (double(N) - 1.0) * log(upper / lower)) *
                     MathematicalConstants::Pi / 180.0;
    }
    else
    {
      *variable_pt = (lower + double(n) / double(N - 1) * (upper - lower)) *
                     MathematicalConstants::Pi / 180.0;
    }

    filestream.open(parameters.dir_name + "/parameters.dat");
    parameters.doc(filestream);
    filestream.close();
  }
}


vector<double> range_string_to_numbers(string range_str)
{
  std::stringstream ss(range_str);
  vector<double> result;
  while (ss.good())
  {
    std::string substr;
    getline(ss, substr, ',');
    result.push_back(stod(substr));
  }
  return result;
}


int main(int argc, char** argv)
{
  std::cout << "Create parameter files." << std::endl;
  CommandLineArgs::setup(argc, argv);

  std::string folder_base = "";
  CommandLineArgs::specify_command_line_flag(
    "--folder", &folder_base, "Optional: Folder base");

  CommandLineArgs::specify_command_line_flag(
    "--overwrite",
    "Optional: Overwrite values. Assumes folders have been created.");

  std::string parameters_filename = "";
  CommandLineArgs::specify_command_line_flag(
    "--parameters", &parameters_filename, "Optional: Base parameter file");

  // Range arguments
  std::string bo_range = "";
  CommandLineArgs::specify_command_line_flag(
    "--bo",
    &bo_range,
    "Optional: Bond number range. ie Bo0,Bo1,N = \"0,10,11\"");

  std::string sl_range = "";
  CommandLineArgs::specify_command_line_flag(
    "--slip_length",
    &sl_range,
    "Optional: Slip length range. ie sl0,sl1,N = \"0.1,0.01,11\"");

  std::string ca_range = "";
  CommandLineArgs::specify_command_line_flag(
    "--angle",
    &ca_range,
    "Optional: Contact angle range. ie lower,upper,N = \"90,60,11\"");

  std::string min_element_range = "";
  CommandLineArgs::specify_command_line_flag(
    "--min_el",
    &min_element_range,
    "Optional: Min element length range. ie l0,l1,N = \"1e-1,1e-4,11\"");

  std::string flux_range = "";
  CommandLineArgs::specify_command_line_flag(
    "--flux", &flux_range, "Optional: Flux range. ie Q0,Q1,N = \"0,10,11\"");

  CommandLineArgs::specify_command_line_flag(
    "--log", "Optional: Log rather than linear");

  // Parse and assign command line arguments
  bool has_unrecognised_arg = false;
  CommandLineArgs::parse_and_assign(argc, argv, has_unrecognised_arg);
  if (has_unrecognised_arg)
  {
    std::cout << "Unrecognised args." << std::endl;
    return -1;
  }

  // If overwrite has been set then we don't need to create folders.
  bool script_is_creating_folders = true;
  if (CommandLineArgs::command_line_flag_has_been_set("--overwrite"))
  {
    script_is_creating_folders = false;
  }

  Parameters parameters;
  if (parameters_filename != "")
  {
    read_parameters_from_file(parameters_filename, parameters);
  }

  //    oomph::Global_Physical_Parameters::Ca = 1.0;
  //    oomph::Global_Physical_Parameters::Equilibrium_contact_angle =
  //      60.0 * oomph::MathematicalConstants::Pi / 180.0;

  //    parameters.max_adapt = 5;

  //    oomph::Mesh_Control_Parameters::min_element_length = 1e-3;
  //    oomph::Mesh_Control_Parameters::Min_element_size = 1e-16;
  //    // oomph::Mesh_Control_Parameters::Max_free_surface_polyline_length =
  //    0.002; oomph::Slip_Parameters::slip_length = 0.01;
  //    oomph::Slip_Parameters::wall_velocity = 0.01;

  std::ofstream filestream;
  if (script_is_creating_folders)
  {
    std::cout << "Creating super directory." << std::endl;
    const int success_flag = mkdir(folder_base.c_str(), 0777);
    if (success_flag == FAIL)
    {
      std::cout << "Failed to create directory: " << folder_base << std::endl;

      return 1;
    }
  }
  if (script_is_creating_folders)
  {
    std::cout << "Creating sub directories and parameter files." << std::endl;
  }
  else
  {
    std::cout << "Overwriting parameter files in folders, creating the folders if "
            "they don't exist."
         << std::endl;
  }

  if (bo_range != "")
  {
    create_parameter_files_over_range_for_variable_pt(
      bo_range,
      folder_base,
      parameters,
      script_is_creating_folders,
      filestream,
      &oomph::Global_Physical_Parameters::Bo);
  }
  else if (sl_range != "")
  {
    create_parameter_files_over_range_for_variable_pt(
      sl_range,
      folder_base,
      parameters,
      script_is_creating_folders,
      filestream,
      &oomph::Slip_Parameters::slip_length);
  }
  else if (min_element_range != "")
  {
    create_parameter_files_over_range_for_variable_pt(
      min_element_range,
      folder_base,
      parameters,
      script_is_creating_folders,
      filestream,
      &oomph::Mesh_Control_Parameters::min_element_length);
  }
  else if (flux_range != "")
  {
    create_parameter_files_over_range_for_variable_pt(
      flux_range,
      folder_base,
      parameters,
      script_is_creating_folders,
      filestream,
      &oomph::Slip_Parameters::wall_velocity);
  }
  else if (ca_range != "")
  {
    create_parameter_files_over_range_for_angle_pt(
      ca_range,
      folder_base,
      parameters,
      script_is_creating_folders,
      filestream,
      &oomph::Global_Physical_Parameters::Equilibrium_contact_angle);
  }
  else
  {
    parameters.dir_name = folder_base;

    const int success_flag = mkdir(parameters.dir_name.c_str(), 0777);
    if (success_flag == FAIL)
    {
      std::cout << "Failed to create directory: " << parameters.dir_name << std::endl;

      // return 1;
    }

    filestream.open(folder_base + "/parameters.dat");
    parameters.doc(filestream);
    filestream.close();
  }

  std::cout << "End of creation of parameter files." << std::endl;

  return 0;
}
