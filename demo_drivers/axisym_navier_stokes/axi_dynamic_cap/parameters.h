#ifndef PARAMETERS_HEADER
#define PARAMETERS_HEADER

#include <iostream>
#include <string>
#include <sys/stat.h>
#include <cmath>

#include "generic.h"

namespace oomph
{
  namespace Flux_Parameters
  {
    double Withdraw_speed = 0e0;
    double Duration = 1e0;
    double Ramp_up_time = 0.1;

    // Flux function
    void flux_fct(const double& t, double& flux)
    {
      // Divide by two to get the flux multipler due to the parabolic flux
      // shape and internal factor of 2*pi
      double flux_multiplier = Withdraw_speed / 2;

      flux = flux_multiplier;

      // flux = sin(20 * (t - MathematicalConstants::Pi / 4)) * std::exp(-10 *
      // t);

      // double Ramp_down_time = Ramp_up_time;
      // if (t < Ramp_up_time)
      //{
      //  flux = t * flux_multiplier / Ramp_up_time;
      //}
      // else if (t >= Ramp_up_time && t < Ramp_up_time + Duration)
      //{
      //  flux = flux_multiplier;
      //}
      // else if (t >= Ramp_up_time + Duration &&
      //         t < Ramp_up_time + Duration + Ramp_down_time)
      //{
      //  flux = (Ramp_down_time - (t - (Ramp_up_time + Duration))) *
      //         flux_multiplier / Ramp_down_time;
      //}
      // else
      //{
      //  flux = 0.0;
      //}
    }
  }; // namespace Flux_Parameters

  namespace Mesh_Control_Parameters
  {
    double Max_permitted_mesh_residual = 1e0;
    double Min_permitted_mesh_residual = 5e-2;
    double min_element_length = 2e-4;
    double inner_min_element_length = 2e-4;
    double element_length_ratio = 1.5;
    unsigned interval_between_adapts = 5;
    unsigned initial_number_of_free_surface_points = 32;
    double Max_free_surface_polyline_length = 5e-2;
    double Max_slip_polyline_length = 1e-1;
    double Max_element_size = 0.5 * std::pow(5e-1, 2.0);
    double Min_element_size = 0.5 * std::pow(min_element_length, 2.0);
    double Min_permitted_angle = 15;
    double Max_permitted_z2_error = 1e-3;
    double Min_permitted_z2_error = 1e-7;

    double Polyline_refinement_tolerence = 4e-3; // 8e-3
    double Polyline_unrefinement_tolerence = 2e-3; // 4e-3

    double Uniform_element_area = 0.5 * std::pow(5e-1, 2.0);
    unsigned Max_number_of_adapts_for_refinement = 0;
    unsigned Error_estimator_flag = 1;

    // Pseudo-solid Poisson ratio
    double Nu = 0.1;

    // Absolute, not scaled properly
    bool Use_adaptive_timestepping = false;
    double Temporal_tolerance = 1e0;
    double Max_timestep = 1e-1;
  }; // namespace Mesh_Control_Parameters

  namespace Plot_Parameters
  {
    unsigned Bulk_element_number_of_plot_points = 3;
    unsigned Surface_element_number_of_plot_points = 3;
  } // namespace Plot_Parameters

  namespace Newton_Solver_Parameters
  {
    unsigned Max_newton_iterations = 40;
    double Max_residual = 1e3;
    double Newton_solver_tolerance = 1e-8; // default: 1-e8
  }; // namespace Newton_Solver_Parameters

  namespace Doc_Parameters
  {
    bool Doc_bulk = true;
    bool Doc_free_surface = true;
    bool Doc_slip = true;
    bool Doc_no_penetration = true;
    bool Doc_volume = true;
    bool Doc_flux = true;
    bool Doc_contact_angle = true;
    bool Doc_trace = true;
    unsigned interval_between_doc = 1;
  }; // namespace Doc_Parameters

  namespace Restart_Parameters
  {
    int interval_between_dump = 10;
  };

  namespace Slip_Parameters
  {
    // Slip length
    double slip_length = 5e-3;
    double wall_velocity = 0;
    double* wall_velocity_pt = &wall_velocity;

    // Slip length near the on the outer wall near the contact line
    void prescribed_slip_fct(const double& t,
                             const oomph::Vector<double>& x,
                             const oomph::Vector<double>& n,
                             oomph::Vector<double>& slip)
    {
      const double local_slip = slip_length;

      // Assign solution
      slip[0] = -1.0; // local_slip;
      slip[1] = local_slip;
      slip[2] = local_slip;
    }

    // Slip length near the on the outer wall near the contact line
    void prescribed_wall_velocity_fct(const double& t,
                                      const oomph::Vector<double>& x,
                                      const oomph::Vector<double>& n,
                                      oomph::Vector<double>& slip)
    {
      // Assign solution
      slip[0] = 0;
      slip[1] = *wall_velocity_pt;
      slip[2] = 0;
    }
  }; // namespace Slip_Parameters

  namespace Global_Physical_Parameters
  {
    double Bo = 1.0e0;
    double Ca = 1e-2;
    double Re = 0e0;

    double Initial_fluid_height = 2.5 + 1.0;

    // Contact angle variables
    double Equilibrium_contact_angle =
      90.0 * oomph::MathematicalConstants::Pi / 180.0;

    // Used in dynamic_contact_angle function otherwise ignored
    double Advancing_contact_angle =
      90.0 * oomph::MathematicalConstants::Pi / 180.0;
    double Receding_contact_angle =
      90.0 * oomph::MathematicalConstants::Pi / 180.0;
    double dynamic_contact_angle_dependence = -0.0;

    bool Use_strong_imposition = false;

    // Direction of gravity, always vertically downward
    oomph::Vector<double> G{0.0, -1.0, 0.0};

    // Direction of the wall normal vector
    oomph::Vector<double> Wall_normal{-1.0, 0.0};

    // Function that specifies the wall unit normal
    void wall_unit_normal_fct(const oomph::Vector<double>& x,
                              oomph::Vector<double>& normal)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        normal[i] = Wall_normal[i];
      }
    }

    // Function pointer to a contact angle function. Returns the
    // contact angle between the fluid and the wall, for a given contact line
    // velocity
    void contact_angle_fct_pt(const oomph::Vector<double>& local_u,
                              double& contact_angle)
    {
      double u = local_u[1];
      if (u > 0)
      {
        contact_angle =
          Receding_contact_angle +
          dynamic_contact_angle_dependence * std::pow(u, 1.0 / 3.0);
      }
      // Not quite. If u is zero then the contact angle can be anywhere from the
      // receding to the advancing
      else if (u == 0)
      {
        contact_angle = Equilibrium_contact_angle;
      }
      else
      {
        contact_angle =
          Advancing_contact_angle -
          dynamic_contact_angle_dependence * std::pow(-u, 1.0 / 3.0);
      }
    }

  } // namespace Global_Physical_Parameters

  class Parameters
  {
  public:
    int max_adapt;
    std::string dir_name;
    double ft;
    double dt;
    std::string restart_filename;
    unsigned azimuthal_mode_number;

    Parameters()
      : max_adapt(0),
        dir_name(""),
        ft(0.0),
        dt(0.0),
        restart_filename(""),
        azimuthal_mode_number(0)
    {
    }

    void set_from_input_stream(std::ifstream& parameter_filestream)
    {
      std::string input_string;

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Global_Physical_Parameters::Bo = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Global_Physical_Parameters::Ca = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Global_Physical_Parameters::Re = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Global_Physical_Parameters::Equilibrium_contact_angle =
        stod(input_string) * oomph::MathematicalConstants::Pi / 180.0;

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      this->max_adapt = stoi(input_string);

      // Output directory name
      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      if (input_string.size() > 0) input_string.resize(input_string.size() - 1);
      this->dir_name = input_string;

      check_for_directory();

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Slip_Parameters::slip_length = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::min_element_length = stod(input_string);
      Mesh_Control_Parameters::inner_min_element_length = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::element_length_ratio = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Flux_Parameters::Withdraw_speed = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Flux_Parameters::Duration = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Flux_Parameters::Ramp_up_time = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_element_size = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Min_element_size = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_free_surface_polyline_length =
        stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_slip_polyline_length = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::interval_between_adapts = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Error_estimator_flag = stoi(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Nu = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Temporal_tolerance = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_timestep = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Use_adaptive_timestepping = stoi(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_permitted_z2_error = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Min_permitted_z2_error = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      this->ft = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      this->dt = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      if (input_string.size() > 0) input_string.resize(input_string.size() - 1);
      this->restart_filename = input_string;

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Slip_Parameters::wall_velocity = stod(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      this->azimuthal_mode_number = stoi(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Mesh_Control_Parameters::Max_number_of_adapts_for_refinement =
        stoi(input_string);

      getline(parameter_filestream, input_string, '#');
      parameter_filestream.ignore(80, '\n');
      Global_Physical_Parameters::Use_strong_imposition = stoi(input_string);
    }

    void doc(std::ofstream& output_stream)
    {
      output_stream.precision(16);
      output_stream << Global_Physical_Parameters::Bo << " # Bond number"
                    << std::endl;
      output_stream << Global_Physical_Parameters::Ca << " # Capillary number"
                    << std::endl;
      output_stream << Global_Physical_Parameters::Re << " # Reynolds number"
                    << std::endl;
      output_stream << Global_Physical_Parameters::Equilibrium_contact_angle *
                         180 / oomph::MathematicalConstants::Pi
                    << " # Contact angle" << std::endl;
      output_stream << this->max_adapt << " # Max number of adapt steps"
                    << std::endl;
      output_stream << this->dir_name << " # Output directory" << std::endl;
      output_stream << Slip_Parameters::slip_length << " # Slip length"
                    << std::endl;
      output_stream << Mesh_Control_Parameters::min_element_length
                    << " # Mininum element length" << std::endl;
      output_stream << Mesh_Control_Parameters::element_length_ratio
                    << " # Element length ratio" << std::endl;
      output_stream << Flux_Parameters::Withdraw_speed << " # Withdraw speed "
                    << std::endl;
      output_stream << Flux_Parameters::Duration << " # Flux duration"
                    << std::endl;
      output_stream << Flux_Parameters::Ramp_up_time << " # Flux ramp up time"
                    << std::endl;
      output_stream << Mesh_Control_Parameters::Max_element_size
                    << " # Max element area" << std::endl;
      output_stream << Mesh_Control_Parameters::Min_element_size
                    << " # Min element area" << std::endl;
      output_stream << Mesh_Control_Parameters::Max_free_surface_polyline_length
                    << " # Max_free_surface_polyline_length " << std::endl;
      output_stream << Mesh_Control_Parameters::Max_slip_polyline_length
                    << " # Max_slip_polyline_length" << std::endl;
      output_stream << Mesh_Control_Parameters::interval_between_adapts
                    << " # interval_between_adapts" << std::endl;
      output_stream << Mesh_Control_Parameters::Error_estimator_flag
                    << " # Error estimator flag, 0 ContactLine, 1 Z2, 2 Corner"
                    << std::endl;
      output_stream << Mesh_Control_Parameters::Nu
                    << " # Pseudo-solid Poisson ratio (Nu)" << std::endl;
      output_stream << Mesh_Control_Parameters::Temporal_tolerance
                    << " # Temporal tolerance" << std::endl;
      output_stream << Mesh_Control_Parameters::Max_timestep
                    << " # Max timestep" << std::endl;
      output_stream << Mesh_Control_Parameters::Use_adaptive_timestepping
                    << " # Use adaptive timestepping" << std::endl;
      output_stream << Mesh_Control_Parameters::Max_permitted_z2_error
                    << " # Max permitted Z2 error" << std::endl;
      output_stream << Mesh_Control_Parameters::Min_permitted_z2_error
                    << " # Min permitted Z2 error" << std::endl;
      output_stream << this->ft << " # Target final time" << std::endl;
      output_stream << this->dt << " # Time step" << std::endl;
      output_stream << this->restart_filename << " # Restart filename"
                    << std::endl;
      output_stream << Slip_Parameters::wall_velocity << " # Wall velocity"
                    << std::endl;
      output_stream << this->azimuthal_mode_number << " # Azimuthal mode number"
                    << std::endl;
      output_stream
        << Mesh_Control_Parameters::Max_number_of_adapts_for_refinement
        << " # Max number of adapts for initial mesh refinement" << std::endl;
      output_stream << Global_Physical_Parameters::Use_strong_imposition
                    << " # Use strong contact angle" << std::endl;
    }

    void check_for_directory()
    {
      struct stat stat_buffer;
      if (stat(this->dir_name.c_str(), &stat_buffer))
      {
        oomph_info << "WARNING: Directory doesn't exist." << std::endl;
      }
    }
  };

  // Read parameters from the passed filename and set them in the passed
  // parameters object.
  void read_parameters_from_file(const std::string& parameters_filename,
                                 Parameters& parameters)
  {
    try
    {
      // Open and read file
      std::ifstream parameter_filestream(parameters_filename);
      parameters.set_from_input_stream(parameter_filestream);
      parameter_filestream.close();
    }
    // if there is an error
    catch (std::exception& e)
    {
      // Warn and throw an exception
      std::cout << "Parameters' filename can't be set, or opened, or read."
                << std::endl;
      std::cout << "File: " << parameters_filename << std::endl;
      throw("Error loading parameters.");
    }
  }


} // namespace oomph
#endif
