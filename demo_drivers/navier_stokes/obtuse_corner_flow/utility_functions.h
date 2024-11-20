#ifndef UTILITY_FUNCTIONS_HEADER
#define UTILITY_FUNCTIONS_HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "generic.h"

namespace oomph
{
  template<typename TA, typename TB>
  bool compare_matrices(TA const& A,
                        TB const& B,
                        const double& abs_err_tolerance = 1e-6)
  {
    bool matrices_are_equal = true;

    const unsigned N = A.nrow();
    const unsigned M = A.ncol();

    unsigned i = 0;
    unsigned j = 0;
    while (i < N)
    {
      j = 0;
      while (j < M)
      {
        if (A(i, j) != B(i, j))
        {
          if (std::abs(A(i, j) - B(i, j)) > abs_err_tolerance)
          {
            matrices_are_equal = false;

            printf("i: %4u, j: %4u", i, j);
            printf(", A: %10.5g", A(i, j));
            printf(", B: %10.5g", B(i, j));
            printf(", abs err: %10.5g\n", std::abs(A(i, j) - B(i, j)));
          }
        }
        j++;
      }
      i++;
    }

    return matrices_are_equal;
  }

  template<typename TA>
  void output_matrices(std::ofstream& output_stream, TA const& A)
  {
    output_stream.precision(16);
    for (unsigned i = 0; i < A.ncol(); i++)
    {
      for (unsigned j = 0; j < A.nrow(); j++)
      {
        output_stream << A(i, j) << " ";
      }
      output_stream << std::endl;
    }
  }

  template<typename TA>
  void output_matrices(TA const& A, const std::string& filename)
  {
    std::ofstream output_stream(filename);
    output_stream.precision(16);
    for (unsigned i = 0; i < A.ncol(); i++)
    {
      for (unsigned j = 0; j < A.nrow(); j++)
      {
        output_stream << A(i, j) << " ";
      }
      output_stream << std::endl;
    }
    output_stream.close();
  }

  template<typename PROBLEM_PT>
  void save_dofs_types(PROBLEM_PT const& problem_pt,
                       const std::string& filename)
  {
    std::ofstream output_stream(filename);
    problem_pt->describe_dofs(output_stream);
    output_stream.close();
  }

  CRDoubleMatrix* load_crdoublematrix(const std::string& filename,
                                      const LinearAlgebraDistribution* dist_pt,
                                      const unsigned& ncol)
  {
    Vector<double> value;
    Vector<int> column_index;
    Vector<int> row_start;

    std::ifstream file_stream(filename, std::ios::in);
    std::string line;
    int row_counter = 0;
    row_start.push_back(0);
    while (std::getline(file_stream, line))
    {
      std::istringstream ss(line);
      std::string token;

      std::getline(ss, token, ' ');
      while (row_counter < std::stoi(token))
      {
        row_start.push_back(value.size());
        row_counter++;
      }

      std::getline(ss, token, ' ');
      column_index.push_back(std::stoi(token));

      std::getline(ss, token);
      value.push_back(std::stod(token));
    }
    file_stream.close();

    row_start.push_back(value.size());

    return new CRDoubleMatrix(dist_pt, ncol, value, column_index, row_start);
  }
} // namespace oomph

#endif
