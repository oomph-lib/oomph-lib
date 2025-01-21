#ifndef BASE_PROBLEM_HEADER
#define BASE_PROBLEM_HEADER

#include <string>
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

  public:
    BaseProblem(std::string parameter_filename = "parameters.dat")
      : Parameters(create_parameters_from_file(parameter_filename))
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

    virtual void doc_solution() = 0;
  };
} // namespace oomph
#endif
