// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

#include "fpdiff.h"
#include "gzip_reader.h"
#include "oomph_definitions.h"

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <vector>

#define CHUNKSIZE 8192

namespace oomph
{
  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  enum TestStatus
  {
    Passed = 0,
    Failed = 1
  };

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  enum class Type
  {
    Number,
    String
  };

  /************************************************************************************
   * @brief
   *
   * @param value
   * @param ending
   * @return true
   * @return false
   ************************************************************************************/
  inline bool ends_with(std::string const& value, std::string const& ending)
  {
    if (ending.size() > value.size())
    {
      return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
  } // End of ends_with

  /************************************************************************************
   * @brief Add number copies of the symbol to the input string.
   *
   * @param string
   * @param symbol
   * @param number
   * @return std::string
   ************************************************************************************/
  inline std::string modify_string_inplace(std::string& text,
                                           const std::string& symbol,
                                           const unsigned& number)
  {
    for (unsigned i = 0; i < number; i++)
    {
      text += symbol;
    }
    text += " ";
    return text;
  } // End of modify_string_inplace

  /************************************************************************************
   * @brief
   *
   * @param str
   ************************************************************************************/
  inline std::vector<std::string> split_string(const std::string& text)
  {
    std::stringstream text_stream{text};
    std::string word;
    std::vector<std::string> words;
    while (std::getline(text_stream, word, ' '))
    {
      words.push_back(word);
    }
    return words;
  } // End of split_string

  /************************************************************************************
   * @brief
   *
   * @param text
   * @return std::string
   ************************************************************************************/
  inline std::string lower(const std::string& text)
  {
    std::string text_copy = text;
    std::transform(text_copy.begin(),
                   text_copy.end(),
                   text_copy.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return text_copy;
  }

  /************************************************************************************
   * @brief Distinguish between a number and a string:
   *
   * Returns integer 1 if the argument is a number,
   *                 2 if the argument is a string.
   *
   * @param text
   * @return true
   * @return false
   ************************************************************************************/
  Type get_type(const std::string& text)
  {
    std::regex is_number{
      "(^[+-]?[0-9]*[.]?[0-9]+$)|(^[+-]?[0-9]+[.]?[0-9]*$)|(^[+-]?[0-9]?"
      "[.]?[0-9]+[EeDd][+-][0-9]+$)"};
    if (std::regex_match(text, is_number))
    {
      return Type::Number;
    }
    return Type::String;
  } // End of get_type

  /************************************************************************************
   * @brief
   *
   * @param filename
   * @return std::vector<std::string>
   ************************************************************************************/
  std::vector<std::string> gzip_load(const std::string& filename)
  {
    if (!ends_with(filename, ".gz"))
    {
      throw OomphLibError("File:\n\t" + filename + "\nis not a gzip file!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    return GZipReader(filename).read_all();
  } // End of gzip_load


  /************************************************************************************
   * @brief
   *
   * @param filename
   * @return std::vector<std::string>
   ************************************************************************************/
  std::vector<std::string> load_file(const std::string& filename)
  {
    std::string line;
    std::vector<std::string> file_data;

    if (ends_with(filename, ".gz"))
    {
      file_data = std::move(gzip_load(filename));
    }
    else
    {
      std::ifstream file(filename);
      while (std::getline(file, line))
      {
        file_data.push_back(line);
      }
      file.close();
    }
    return file_data;
  } // End of load_file

  /************************************************************************************
   * @brief
   *
   * @param filename1
   * @param filename2
   * @param relative_error
   * @param small
   * @param outstream
   * @return int
   ************************************************************************************/
  int fpdiff(const std::string& filename1,
             const std::string& filename2,
             std::ostream& outstream,
             const double& relative_error,
             const double& small)
  {
    // Storage for worst case error sizes
    double max_rel_diff = 0.0;
    double max_wrong_entry = 0.0;

    // Open the files (if run as a script then open Faileds are handled higher
    // up by catching the error, otherwise it is the parent program's job to
    // handle the error and so we shouldn't do anything weird here).
    auto file1 = load_file(filename1);
    auto file2 = load_file(filename2);

    // Find the number of lines in each file
    auto n1 = file1.size();
    auto n2 = file2.size();
    auto min_lines = std::min(n1, n2);

    if (n1 < n2)
    {
      std::vector<std::string> temp(std::move(file1));
      file1 = std::move(file2);
      file2 = std::move(temp);
    }

    // Counter for the number of errors
    unsigned long long nerr = 0;

    // Counter for the number of lines
    unsigned long long count = -1;

    // Counter for the number of lines with errors
    unsigned long long nline_error = 0;

    // Loop over the lines in file1 (the file with the most lines!)
    for (const auto& line1 : file1)
    {
      // Increase the counter
      count++;

      // If we've run over the end of the file2, issue a warning and end the
      // loop
      if (count >= min_lines)
      {
        outstream << "\nWarning: files have different numbers of lines\n"
                  << "\nResults are for first " << count
                  << " lines of both files\n"
                  << std::endl;
        nerr++;
        break;
      }

      // Read the next line from file2
      auto line2 = file2[count];

      // If the lines are the same, we're done
      if (line1 == line2) continue;

      // We need to do more work

      // Split each line into its separate fields
      auto fields1 = split_string(line1);
      auto fields2 = split_string(line2);

      // Find the number of fields in each line
      auto n_field1 = fields1.size();
      auto n_field2 = fields2.size();

      if (n_field1 != n_field2)
      {
        outstream << "\n =====> line " << (count + 1)
                  << ": different number of fields\n"
                  << n_field1 << " fields: " << line1 << "\n"
                  << n_field2 << " fields: " << line2 << std::endl;
        nerr++;
        continue;
      }

      // Flag to indicate whether there has been a problem in the field
      bool encountered_problem = false;

      // Strings that will hold the output data
      std::string outputline1 = "";
      std::string outputline2 = "";
      std::string outputline3 = "";

      // Loop over the fields
      for (unsigned i = 0; i < n_field1; i++)
      {
        // Start by loading the fields into the outputlines (plus whitespace)
        outputline1 += (fields1[i] + " ");
        outputline3 += (fields2[i] + " ");

        // Find the lengths of the fields
        auto length1 = fields1[i].length();
        auto length2 = fields2[i].length();

        // Pad the shortest field so the lengths are the same
        auto field_length = length2;
        if (length1 < length2)
        {
          for (unsigned j = 0; j < (length2 - length1); j++)
          {
            outputline1 += " ";
          }
        }
        else
        {
          field_length = length1;
          for (unsigned j = 0; j < (length1 - length2); j++)
          {
            outputline3 += " ";
          }
        } // if (length1 < length2)

        // If the fields are identical, we are fine
        if (fields1[i] == fields2[i])
        {
          // Put spaces into the error line
          modify_string_inplace(outputline2, " ", field_length);
          continue;
        }

        // Find the type (numeric or string) of each field
        Type type1 = get_type(fields1[i]);
        Type type2 = get_type(fields2[i]);

        // If the data-types aren't the same issue an error
        if ((type1 != type2) || (type1 == Type::String))
        {
          std::string symbol = (type1 != type2 ? "*" : "%");
          encountered_problem = true;
          nerr++;

          // Put the appropriate symbol into the error line
          modify_string_inplace(outputline2, symbol, field_length);
          continue;
        }

        // Convert strings to floating point number
        auto string_to_double = [](const std::string& text) {
          return std::stod(
            std::regex_replace(lower(text), std::regex("d"), "e"));
        };
        double x1 = string_to_double(fields1[i]);
        double x2 = string_to_double(fields2[i]);

        // If both numbers are very small, that's fine
        if ((std::fabs(x1) <= small) && (std::fabs(x2) <= small))
        {
          // Put spaces into the error line
          modify_string_inplace(outputline2, " ", field_length);
        }
        else
        {
          // Find the relative difference based on the largest number
          // Note that this "minimises" the relative error (in some sense)
          // but means that I don't have to separately trap the cases
          // when x1, x2 are zero
          double diff = 100.0 * (std::fabs(x1 - x2) /
                                 std::max(std::fabs(x1), std::fabs(x2)));

          // If the relative error is smaller than the tolerance, that's fine
          if (diff <= relative_error)
          {
            // Put spaces into the error line
            modify_string_inplace(outputline2, " ", field_length);
          }
          // Otherwise issue an error
          else
          {
            encountered_problem = true;
            nerr++;

            // Put the appropriate symbols into the error line
            modify_string_inplace(outputline2, "-", field_length);

            // Record any changes in the worst case values
            max_rel_diff = std::max(max_rel_diff, diff);
            // max_wrong_entry =
            //   std::max({max_wrong_entry, std::fabs(x1), std::fabs(x2)});
            if (std::fabs(x1) > max_wrong_entry)
            {
              max_wrong_entry = std::fabs(x1);
            }
            else if (std::fabs(x2) > max_wrong_entry)
            {
              max_wrong_entry = std::fabs(x2);
            }
          }
        } // if ((std::fabs(x1) <= small) && (std::fabs(x2) <= small))
      } // for (unsigned i = 0; i < n_field1; i++)

      // If there has been any sort of error, print it
      if (encountered_problem)
      {
        nline_error++;
        outstream << "\n =====> line " << (count + 1) << "\n"
                  << outputline1 << "\n"
                  << outputline2 << "\n"
                  << outputline3 << std::endl;
      }
    } // for (const auto& line1 : file1)

    std::string border(80, '*');
    std::cout << border << "\nfpdiff() compared files:"
              << "\n  * " << filename1 << "\n  * " << filename2 << std::endl;

    if (nerr > 0)
    {
      outstream << "\n In files " << filename1 << " " << filename2
                << "\n number of lines processed: " << count
                << "\n number of lines containing errors: " << nline_error
                << "\n number of errors: " << nerr << " "
                << "\n largest relative error: " << max_rel_diff << " "
                << "\n largest abs value of an entry which caused an error: "
                << max_wrong_entry << " "
                << "\n========================================================"
                << "\n    Parameters used:"
                << "\n        threshold for numerical zero : " << small
                << "\n        maximum rel. difference [percent] : "
                << relative_error << "\n    Legend: "
                << "\n        *******  means differences in data type "
                   "(string vs number)"
                << "\n        -------  means real data exceeded the relative "
                   "difference maximum"
                << "\n        %%%%%%%  means that two strings are different"
                << "\n========================================================"
                << "\n\n   [FAILED]" << std::endl;
      std::cout << "Test failed!\n" << border << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      outstream << "\n\n In files " << filename1 << " " << filename2
                << "\n   [OK] for fpdiff.py parameters: - max. rel. error = "
                << relative_error << " %"
                << "\n                                  - numerical zero  = "
                << small << std::endl;
      std::cout << "Test passed!\n" << border << std::endl;
      return EXIT_SUCCESS;
    }
  } // End of fpdiff

  /************************************************************************************
   * @brief
   *
   * @param filename1:
   * @param filename2:
   * @param log_file:
   * @param relative_error:
   * @param small:
   * @return int:
   ************************************************************************************/
  int fpdiff(const std::string& filename1,
             const std::string& filename2,
             const std::string& log_file,
             const double& relative_error,
             const double& small)
  {
    // File to output fpdiff output to
    std::ofstream log_stream{log_file};

    // Run test
    int test_status = fpdiff(filename1, filename2, log_stream);

    // Close file
    log_stream.close();

    // Return test result
    return test_status;
  } // End of fpdiff
} // namespace oomph