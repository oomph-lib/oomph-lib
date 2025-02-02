#ifndef BASE_PROBLEM_HEADER
#define BASE_PROBLEM_HEADER

#include <string>
#include <iostream>
#include "generic.h"
#include "parameter_struct.h"

namespace oomph
{
  // Base problem class
  class BaseProblem : public Problem
  {
  private:
    DocInfo Doc_info;
    Params Parameters;
    bool Is_trace_file_created;

  public:
    BaseProblem(std::string parameter_filename = "parameters.dat")
      : Parameters(create_parameters_from_file(parameter_filename)),
        Is_trace_file_created(false)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      Doc_info.number() = 0;
    }

    BaseProblem(Params& params) : Parameters(params)
    {
      // Create and add the timestepper
      add_time_stepper_pt(new BDF<2>);

      // Assign doc info pointer
      Doc_info.number() = 0;
    }

    virtual void setup() = 0;

    DocInfo* doc_info_pt()
    {
      return &Doc_info;
    }

    void set_output_directory(const std::string& directory)
    {
      Doc_info.set_directory(directory);
    }

    void set_output_number(const unsigned& number)
    {
      Doc_info.number() = number;
    }

    Params& parameters()
    {
      return Parameters;
    }

    virtual void doc_solution()
    {
      doc_trace();
      Doc_info.number()++;
    }

    virtual void doc_trace()
    {
      // Open the trace file
      std::string trace_filename = Doc_info.directory() + "/trace.dat";
      std::ofstream trace_file;
      // If the file has not been created, add the header
      trace_file.open(trace_filename, std::ios::app);
      if (!Is_trace_file_created)
      {
        add_trace_header(trace_file);
        // Next line
        trace_file << "\n";
        Is_trace_file_created = true;
      }
      add_trace(trace_file);
      // End the line and print to file
      trace_file << std::endl;

      // Close the trace file
      trace_file.close();
    }

    virtual void add_trace_header(std::ofstream& trace_file)
    {
      trace_file << "doc_number ";
      trace_file << "n_element ";
      trace_file << "dofs ";
      trace_file << "min_error ";
      trace_file << "max_error ";
      trace_file << "min_element_size ";
      trace_file << "max_element_size ";
    }

    virtual void add_trace(std::ofstream& trace_file)
    {
      // Compute info required
      double max_error = 0.0;
      double min_error = 0.0;
      double max_element_size = 0.0;
      double min_element_size = 0.0;
      unsigned n_element = 0;
      if (this->nsub_mesh() > 0)
      {
        // Set the minimums to be a large number
        min_error = 1e6;
        min_element_size = 1e6;

        // Number of elements
        // Assume submesh 0 has the bulk mesh.
        n_element = this->mesh_pt(0)->nelement();

        // Compute and set the elemental error
        Vector<double> elemental_error(n_element);
        Z2ErrorEstimator error_estimator;
        Mesh* fluid_mesh_pt = dynamic_cast<Mesh*>(mesh_pt(0));
        error_estimator.get_element_errors(fluid_mesh_pt, elemental_error);

        for (unsigned e = 0; e < n_element; e++)
        {
          max_error = std::max(max_error, elemental_error[e]);
          min_error = std::min(min_error, elemental_error[e]);
          double element_size = fluid_mesh_pt->finite_element_pt(e)->size();
          min_element_size = std::min(min_element_size, element_size);
          max_element_size = std::max(max_element_size, element_size);
        }
      }

      // Output the info
      trace_file << Doc_info.number() << " ";
      trace_file << n_element << " ";
      trace_file << this->ndof() << " ";
      trace_file << min_error << " ";
      trace_file << max_error << " ";
      trace_file << min_element_size << " ";
      trace_file << max_element_size << " ";
    }
  };
} // namespace oomph
#endif
