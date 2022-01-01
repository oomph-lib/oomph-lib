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
#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <algorithm>
#include <limits.h>
#include <cstring>

#ifdef OOMPH_HAS_UNISTDH
#include <unistd.h> // for getpid()
#endif

#include "oomph_utilities.h"
#include "Vector.h"
#include "matrices.h"

namespace oomph
{
  //======start_of_ANSIEscapeCode_namespace=============================
  /// Contains an enumeration of the ANSI escape codes used for
  /// colouring text (when piped to the command line). Adapted from
  /// the guide on:
  /// https://stackoverflow.com/questions/2616906/how-do-i-output-
  /// coloured-text-to-a-linux-terminal?utm_medium=organic&utm_source=
  /// google_rich_qa&utm_campaign=google_rich_qa
  /// Here, \033 is the ESC character, ASCII 27. It is followed by [,
  /// then zero or more numbers separated by ;, and finally the letter
  /// m. The numbers describe the colour and format to switch to from
  /// that point onwards.
  //====================================================================
  namespace ANSIEscapeCode
  {
    /// Function to change text effect. NOTE: This assumes the user
    /// knows what they're doing/assigning; no error checking done here...
    void set_text_effect(std::string text_effect)
    {
      // Assign the string
      Text_effect = text_effect;
    } // End of set_text_effect

    /// Variable to decide on effects
    std::string Text_effect = "1";

    /// The code for each type of colour
    std::string Black = "\033[" + Text_effect + ";30m";
    std::string Red = "\033[" + Text_effect + ";31m";
    std::string Green = "\033[" + Text_effect + ";32m";
    std::string Yellow = "\033[" + Text_effect + ";33m";
    std::string Blue = "\033[" + Text_effect + ";34m";
    std::string Magenta = "\033[" + Text_effect + ";35m";
    std::string Cyan = "\033[" + Text_effect + ";36m";
    std::string Reset = "\033[0m";
  } // namespace ANSIEscapeCode


  //=====================================================================
  /// Namespace for debugging helpers. Currently only contains a
  /// function to prett-ify file name and line numbers (in red) to use
  /// when debugging. Makes it easy to identify where a std::cout
  /// statement was called.
  //=====================================================================
  namespace DebugHelpers
  {
    /// Return the concaternation of the initials of the input
    /// file name and line number. The make_new_line flag indicates
    /// whether the string starts with a "\n", i.e. a new line
    std::string debug_string(const std::string& filename,
                             const int& line,
                             const bool& make_new_line)
    {
      // Make a string
      std::string debug_string;

      // Temporary storage for the filename which can be edited
      std::string file = filename;

      // The delimeter
      std::string delimiter = "_";

      // Posiiton of the delimiter
      size_t pos = 0;

      // The string before the delimiter
      std::string token;

      // Get the filename (gets rid of e.g. ./ before filename)
      while ((pos = file.find("/")) != std::string::npos)
      {
        // Get the string before the delimiter
        token = file.substr(0, pos);

        // Erase the string before the delimeter
        file = file.substr(pos + delimiter.length());
      }

      // While we can find delimeters
      while ((pos = file.find(delimiter)) != std::string::npos)
      {
        // Get the string before the delimiter
        token = file.substr(0, pos);

        // Output the string
        debug_string += toupper(token.at(0));

        // Erase the delimeter
        file = file.substr(pos + delimiter.length());
      }

      // Output the string
      debug_string += toupper(file.at(0));

      // String stream
      std::ostringstream debug_stream;

      // If they want a new line
      if (make_new_line)
      {
        // Add a newline string
        debug_stream << "\n";
      }

      // Create debug string
      debug_stream << "\033[1;31m" << debug_string << line << ":\033[0m ";

      // Return the string
      return debug_stream.str();
    } // End of create_debug_string
  } // namespace DebugHelpers


  //=====================================================================
  /// Helper namespace containing function that computes second invariant
  /// of tensor
  //=====================================================================
  namespace SecondInvariantHelper
  {
    /// Compute second invariant of tensor
    double second_invariant(const DenseMatrix<double>& tensor)
    {
      // get size of tensor
      unsigned n = tensor.nrow();
      double trace_of_tensor = 0.0;
      double trace_of_tensor_squared = 0.0;
      for (unsigned i = 0; i < n; i++)
      {
        // Calculate the trace of the tensor: T_ii
        trace_of_tensor += tensor(i, i);
        // Calculate the trace of the tensor squared: T_ij T_ji
        for (unsigned j = 0; j < n; j++)
        {
          trace_of_tensor_squared += tensor(i, j) * tensor(j, i);
        }
      }

      // return the second invariant: 1.0/2.0*[(T_ii)^2 - T_ij T_ji]
      return 0.5 *
             (trace_of_tensor * trace_of_tensor - trace_of_tensor_squared);
    }

  } // namespace SecondInvariantHelper


  //==============================================
  /// Namespace for error messages for broken
  /// copy constructors and assignment operators
  //==============================================
  namespace BrokenCopy
  {
    /// Issue error message and terminate execution
    void broken_assign(const std::string& class_name)
    {
      // Write the error message into a string
      std::string error_message = "Assignment operator for class\n\n";
      error_message += class_name;
      error_message += "\n\n";
      error_message += "is deliberately broken to avoid the accidental \n";
      error_message += "use of the inappropriate C++ default.\n";
      error_message += "If you really need an assignment operator\n";
      error_message += "for this class, write it yourself...\n";

      throw OomphLibError(
        error_message, "broken_assign()", OOMPH_EXCEPTION_LOCATION);
    }


    /// Issue error message and terminate execution
    void broken_copy(const std::string& class_name)
    {
      // Write the error message into a string
      std::string error_message = "Copy constructor for class\n\n";
      error_message += class_name;
      error_message += "\n\n";
      error_message += "is deliberately broken to avoid the accidental\n";
      error_message += "use of the inappropriate C++ default.\n";
      error_message +=
        "All function arguments should be passed by reference or\n";
      error_message +=
        "constant reference. If you really need a copy constructor\n";
      error_message += "for this class, write it yourself...\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  } // namespace BrokenCopy


  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////
  /// ///////////////////////////////////////////////////////////////


  //====================================================================
  /// Namespace for global (cumulative) timings
  //====================================================================
  namespace CumulativeTimings
  {
    /// (Re-)start i-th timer
    void start(const unsigned& i)
    {
      Start_time[i] = clock();
    }

    /// Halt i-th timer
    void halt(const unsigned& i)
    {
      Timing[i] += clock() - Start_time[i];
    }

    /// Report time accumulated by i-th timer
    double cumulative_time(const unsigned& i)
    {
      return double(Timing[i]) / CLOCKS_PER_SEC;
    }

    /// Reset i-th timer
    void reset(const unsigned& i)
    {
      Timing[i] = clock_t(0.0);
    }

    /// Reset all timers
    void reset()
    {
      unsigned n = Timing.size();
      for (unsigned i = 0; i < n; i++)
      {
        Timing[i] = clock_t(0.0);
      }
    }

    /// Set number of timings that can be recorded in parallel
    void set_ntimers(const unsigned& ntimers)
    {
      Timing.resize(ntimers, clock_t(0.0));
      Start_time.resize(ntimers, clock_t(0.0));
    }

    /// Cumulative timings
    Vector<clock_t> Timing;

    /// Start times of active timers
    Vector<clock_t> Start_time;

  } // namespace CumulativeTimings


  //======================================================================
  /// Set output directory (we try to open a file in it
  /// to see if the directory exists -- if it doesn't we'll
  /// issue a warning -- or, if directory_must_exist()==true,
  /// die by throwing and OomphLibError
  //======================================================================
  void DocInfo::set_directory(const std::string& directory_)
  {
    // Try to open a file in output directory
    std::ostringstream filename;
    filename << directory_ << "/.dummy_check.dat";
    std::ofstream some_file;
    some_file.open(filename.str().c_str());
    if (!some_file.is_open())
    {
      // Construct the error message
      std::string error_message = "Problem opening output file.\n";
      error_message += "I suspect you haven't created the output directory ";
      error_message += directory_;
      error_message += "\n";

      // Issue a warning if the directory does not have to exist
      if (!Directory_must_exist)
      {
        // Create an Oomph Lib warning
        OomphLibWarning(
          error_message, "set_directory()", OOMPH_EXCEPTION_LOCATION);
      }
      // Otherwise throw an error
      else
      {
        error_message += "and the Directory_must_exist flag is true.\n";
        throw OomphLibError(
          error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }
    // Write to the dummy file
    some_file << "Dummy file, opened to check if output directory "
              << std::endl;
    some_file << "exists. Can be deleted...." << std::endl;
    some_file.close();
    // Set directory
    Directory = directory_;
  }


  // =================================================================
  /// Conversion functions for easily making strings (e.g. for filenames - to
  /// avoid stack smashing problems with cstrings and long filenames).
  // =================================================================
  namespace StringConversion
  {
    /// Convert a string to lower case (outputs a copy).
    std::string to_lower(const std::string& input)
    {
      std::string output(input);
      std::string::iterator it;
      for (it = output.begin(); it != output.end(); ++it)
      {
        ::tolower(*it);
      }

      return output;
    }

    /// Convert a string to upper case (outputs a copy).
    std::string to_upper(const std::string& input)
    {
      std::string output(input);
      std::string::iterator it;
      for (it = output.begin(); it != output.end(); ++it)
      {
        ::toupper(*it);
      }
      return output;
    }

    /// Split a string, s, into a vector of strings where ever there is
    /// an instance of delimiter (i.e. is delimiter is " " will give a list of
    /// words). Note that mutliple delimiters in a row will give empty
    /// strings.
    void split_string(const std::string& s,
                      char delim,
                      Vector<std::string>& elems)
    {
      // From http://stackoverflow.com/questions/236129/splitting-a-string-in-c
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim))
      {
        elems.push_back(item);
      }
    }

    /// Split a string, s, into a vector of strings where ever there is
    /// an instance of delimiter (i.e. is delimiter is " " will give a list of
    /// words). Note that mutliple delimiters in a row will give empty
    /// strings. Return by value.
    Vector<std::string> split_string(const std::string& s, char delim)
    {
      // From http://stackoverflow.com/questions/236129/splitting-a-string-in-c
      Vector<std::string> elems;
      split_string(s, delim, elems);
      return elems;
    }

  } // namespace StringConversion


  //====================================================================
  /// Namespace for command line arguments
  //====================================================================
  namespace CommandLineArgs
  {
    /// Number of arguments + 1
    int Argc;

    /// Arguments themselves
    char** Argv;

    /// Map to indicate an input flag as having been set
    std::map<std::string, ArgInfo<bool>> Specified_command_line_flag;

    /// Map to associate an input flag with a double -- specified via pointer
    std::map<std::string, ArgInfo<double>> Specified_command_line_double_pt;

    /// Map to associate an input flag with an int -- specified via pointer
    std::map<std::string, ArgInfo<int>> Specified_command_line_int_pt;

    /// Map to associate an input flag with an unsigned -- specified via pointer
    std::map<std::string, ArgInfo<unsigned>> Specified_command_line_unsigned_pt;

    /// Map to associate an input flag with a string -- specified via pointer
    std::map<std::string, ArgInfo<std::string>>
      Specified_command_line_string_pt;

    /// Set values
    void setup(int argc, char** argv)
    {
      Argc = argc;
      Argv = argv;
    }

    /// Doc the command line arguments
    void output()
    {
      oomph_info << "You are running the program: " << CommandLineArgs::Argv[0]
                 << std::endl;
      oomph_info << "with the following command line args: " << std::endl;
      std::stringstream str;
      for (int i = 1; i < CommandLineArgs::Argc; i++)
      {
        str << CommandLineArgs::Argv[i] << " ";
      }
      oomph_info << str.str() << std::endl;
    }


    /// Specify possible argument-free command line flag
    void specify_command_line_flag(const std::string& command_line_flag,
                                   const std::string& doc)
    {
      Specified_command_line_flag[command_line_flag] =
        ArgInfo<bool>(false, 0, doc);
    }

    /// Specify possible command line flag that specifies a double,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   double* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_double_pt[command_line_flag] =
        ArgInfo<double>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies an int,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   int* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_int_pt[command_line_flag] =
        ArgInfo<int>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies an unsigned,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   unsigned* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_unsigned_pt[command_line_flag] =
        ArgInfo<unsigned>(false, arg_pt, doc);
    }

    /// Specify possible command line flag that specifies a string,
    /// accessed via pointer
    void specify_command_line_flag(const std::string& command_line_flag,
                                   std::string* arg_pt,
                                   const std::string& doc)
    {
      Specified_command_line_string_pt[command_line_flag] =
        ArgInfo<std::string>(false, arg_pt, doc);
    }


    /// Check if command line flag has been set (value will have been
    /// assigned directly).
    bool command_line_flag_has_been_set(const std::string& flag)
    {
      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        if (it->first == flag)
        {
          return it->second.is_set;
        }
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        if (it->first == flag)
        {
          return (it->second).is_set;
        }
      }

      return false;
    }

    /// Document the values of all flags (specified or not).
    void doc_all_flags(std::ostream& outstream)
    {
      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        outstream << it->first << " " << it->second.is_set << std::endl;
      }
      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        outstream << it->first << " " << *(it->second.arg_pt) << std::endl;
      }
      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        // Quote blank strings, otherwise trying to parse the output in any way
        // will go wrong.
        std::string arg_string = *(it->second.arg_pt);
        if (arg_string == "")
        {
          arg_string = "\"\"";
        }

        outstream << it->first << " " << arg_string << std::endl;
      }
    }

    /// Document specified command line flags
    void doc_specified_flags()
    {
      oomph_info << std::endl;
      oomph_info << "Specified (and recognised) command line flags:\n";
      oomph_info << "----------------------------------------------\n";

      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        if (it->second.is_set)
        {
          oomph_info << it->first << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          oomph_info << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          oomph_info << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          oomph_info << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        if (it->second.is_set)
        {
          oomph_info << it->first << " " << *(it->second.arg_pt) << std::endl;
        }
      }

      oomph_info << std::endl;
    }


    /// Document available command line flags
    void doc_available_flags()
    {
      oomph_info << std::endl;
      oomph_info << "Available command line flags:\n";
      oomph_info << "-----------------------------\n";

      for (std::map<std::string, ArgInfo<bool>>::iterator it =
             Specified_command_line_flag.begin();
           it != Specified_command_line_flag.end();
           it++)
      {
        oomph_info << it->first << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<double>>::iterator it =
             Specified_command_line_double_pt.begin();
           it != Specified_command_line_double_pt.end();
           it++)
      {
        oomph_info << it->first << " <double> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<int>>::iterator it =
             Specified_command_line_int_pt.begin();
           it != Specified_command_line_int_pt.end();
           it++)
      {
        oomph_info << it->first << " <int> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
             Specified_command_line_unsigned_pt.begin();
           it != Specified_command_line_unsigned_pt.end();
           it++)
      {
        oomph_info << it->first << " <unsigned> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      for (std::map<std::string, ArgInfo<std::string>>::iterator it =
             Specified_command_line_string_pt.begin();
           it != Specified_command_line_string_pt.end();
           it++)
      {
        oomph_info << it->first << " <string> " << std::endl
                   << it->second.doc << std::endl
                   << std::endl;
      }

      oomph_info << std::endl;
    }


    /// Helper function to check if command line index is legal
    void check_arg_index(const int& argc, const int& arg_index)
    {
      if (arg_index >= argc)
      {
        output();
        doc_available_flags();
        std::stringstream error_stream;
        error_stream
          << "Tried to read more command line arguments than\n"
          << "specified. This tends to happen if a required argument\n"
          << "to a command line flag was omitted, e.g. by running \n\n"
          << "     ./a.out -some_double \n\n rather than\n\n"
          << "     ./a.out -some_double 1.23 \n\n"
          << "To aid the debugging I've output the available\n"
          << "command line arguments above.\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
    }


    /// Parse command line, check for recognised flags and assign
    /// associated values
    void parse_and_assign(int argc,
                          char* argv[],
                          const bool& throw_on_unrecognised_args)
    {
      // Keep looping over command line arguments
      int arg_index = 1;
      while (arg_index < argc)
      {
        bool found_match = false;

        if ((strcmp(argv[arg_index], "-help") == 0) ||
            (strcmp(argv[arg_index], "--help") == 0))
        {
          oomph_info
            << "NOTE: You entered --help or -help on the command line\n";
          oomph_info << "so I'm going to tell you about this code's\n";
          oomph_info << "available command line flags and then return.\n";
          doc_available_flags();

#ifdef OOMPH_HAS_MPI
          int flag;
          MPI_Initialized(&flag);
          if (flag != 0) MPI_Helpers::finalize();
#endif
          oomph_info << "Shutting down...\n";
          exit(0);
        }


        // Check if the flag has been previously specified as a simple argument
        // free command line argument
        for (std::map<std::string, ArgInfo<bool>>::iterator it =
               Specified_command_line_flag.begin();
             it != Specified_command_line_flag.end();
             it++)
        {
          if (it->first == argv[arg_index])
          {
            Specified_command_line_flag[argv[arg_index]].is_set = true;
            found_match = true;
            break;
          }
        }

        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) a double
          for (std::map<std::string, ArgInfo<double>>::iterator it =
                 Specified_command_line_double_pt.begin();
               it != Specified_command_line_double_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the double. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = atof(argv[arg_index]);
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) an int
          for (std::map<std::string, ArgInfo<int>>::iterator it =
                 Specified_command_line_int_pt.begin();
               it != Specified_command_line_int_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the int. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = atoi(argv[arg_index]);
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) an
          // unsigned
          for (std::map<std::string, ArgInfo<unsigned>>::iterator it =
                 Specified_command_line_unsigned_pt.begin();
               it != Specified_command_line_unsigned_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the unsigned. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = unsigned(atoi(argv[arg_index]));
              found_match = true;
              break;
            }
          }
        }


        if (!found_match)
        {
          // Check if the flag has been previously specified as a
          // command line argument that specifies (and is followed by) a string
          for (std::map<std::string, ArgInfo<std::string>>::iterator it =
                 Specified_command_line_string_pt.begin();
               it != Specified_command_line_string_pt.end();
               it++)
          {
            if (it->first == argv[arg_index])
            {
              // Next command line argument specifies the string. Set it via
              // the pointer
              arg_index++;
              check_arg_index(argc, arg_index);
              it->second.is_set = true;
              *(it->second.arg_pt) = argv[arg_index];
              found_match = true;
              break;
            }
          }
        }


        // Oh dear, we still haven't found the argument in the list.
        // Maybe it was specified wrongly -- issue warning.
        if (!found_match)
        {
          // Construct the error message
          std::string error_message = "Command line argument\n\n";
          error_message += argv[arg_index];
          error_message += "\n\nwas not recognised. This may be legal\n";
          error_message += "but seems sufficiently suspicious to flag up.\n";
          error_message += "If it should have been recognised, make sure you\n";
          error_message += "used the right number of \"-\" signs...\n";

          if (throw_on_unrecognised_args)
          {
            error_message += "Throwing an error because you requested it with";
            error_message += " throw_on_unrecognised_args option.";
            throw OomphLibError(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
          else
          {
            // Create an Oomph Lib warning
            OomphLibWarning(
              error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
          }
        }


        arg_index++;
      }
    }


    /// Parse previously specified command line, check for
    /// recognised flags and assign associated values
    void parse_and_assign(const bool& throw_on_unrecognised_args)
    {
      parse_and_assign(Argc, Argv, throw_on_unrecognised_args);
    }

  } // namespace CommandLineArgs

#ifdef OOMPH_HAS_MPI
  //========================================================================
  /// Single (global) instantiation of the mpi output modifier
  //========================================================================
  MPIOutputModifier oomph_mpi_output;

  //========================================================================
  /// Precede the output by the processor ID but output everything
  //========================================================================
  bool MPIOutputModifier::operator()(std::ostream& stream)
  {
    int my_rank = Communicator_pt->my_rank();

    if (!Output_from_single_processor)
    {
      stream << "Processor " << my_rank << ":   ";
      // Continue processing
      return true;
    }
    else
    {
      if (unsigned(my_rank) == Output_rank)
      {
        stream << "Processor " << my_rank << ":   ";
        // Continue processing
        return true;
      }
      else
      {
        return false;
      }
    }
  }

#endif

  //=============================================================================
  /// Initialize mpi. If optional boolean flag is set to false, we use
  /// MPI_COMM_WORLD itself as oomph-lib's communicator. Defaults to true.
  //=============================================================================
  void MPI_Helpers::init(int argc,
                         char** argv,
                         const bool& make_duplicate_of_mpi_comm_world)
  {
#ifdef OOMPH_HAS_MPI
    // call mpi int
    MPI_Init(&argc, &argv);


    // By default, create the oomph-lib communicator using MPI_Comm_dup so that
    // the communicator has the same group of processes but a new context
    MPI_Comm oomph_comm_world = MPI_COMM_WORLD;
    if (make_duplicate_of_mpi_comm_world)
    {
      MPI_Comm_dup(MPI_COMM_WORLD, &oomph_comm_world);
    }

    if (MPI_COMM_WORLD != oomph_comm_world)
    {
      oomph_info << "Oomph-lib communicator is a duplicate of MPI_COMM_WORLD\n";
    }
    else
    {
      oomph_info << "Oomph-lib communicator is MPI_COMM_WORLD\n";
    }

    // create the oomph-lib communicator
    // note: oomph_comm_world is deleted when the destructor of
    // Communicator_pt is called
    Communicator_pt = new OomphCommunicator(oomph_comm_world, true);

    // Change MPI error handler so that error will return
    // rather than aborting
    MPI_Comm_set_errhandler(oomph_comm_world, MPI_ERRORS_RETURN);

    // Use MPI output modifier: Each processor precedes its output
    // by its rank
    oomph_mpi_output.communicator_pt() = Communicator_pt;
    oomph_info.output_modifier_pt() = &oomph_mpi_output;
#else
    // create a serial communicator
    Communicator_pt = new OomphCommunicator;
#endif
    MPI_has_been_initialised = true;
  }

  //=============================================================================
  /// finalize mpi
  //=============================================================================
  void MPI_Helpers::finalize()
  {
    // delete the communicator
    delete Communicator_pt;

    // and call MPI_Finalize
#ifdef OOMPH_HAS_MPI
    MPI_Finalize();
#endif
  }

  //=============================================================================
  /// access to the global oomph-lib communicator
  //=============================================================================
  OomphCommunicator* MPI_Helpers::communicator_pt()
  {
#ifdef MPI
#ifdef PARANOID
    if (!MPI_has_been_initialised)
    {
      std::ostringstream error_message_stream;
      error_message_stream
        << "MPI has not been initialised.\n Call MPI_Helpers::init(...)";
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return Communicator_pt;

#else // ifndef MPI

    // If MPI is not enabled then we are unlikely to have called init, so to
    // be nice to users create a new serial comm_pt if none exits.

    if (Communicator_pt == 0)
    {
      Communicator_pt = new OomphCommunicator;
    }

    // #ifdef PARANOID
    //  if(!Communicator_pt->serial_communicator())
    //   {
    //    std::string error_msg =
    //     "MPI_Helpers has somehow ended up with a non-serial\n"
    //     + "communicator pointer even though MPI is disabled!";
    //    throw OomphLibError(error_msg.str(),
    //                        OOMPH_CURRENT_FUNCTION,
    //                        OOMPH_EXCEPTION_LOCATION);
    //   }
    // #endif

    return Communicator_pt;

#endif // end ifdef MPI
  }

  bool MPI_Helpers::MPI_has_been_initialised = false;
  OomphCommunicator* MPI_Helpers::Communicator_pt = 0;


  //====================================================================
  /// Namespace for flagging up obsolete parts of the code
  //====================================================================
  namespace ObsoleteCode
  {
    /// Flag up obsolete parts of the code
    bool FlagObsoleteCode = true;

    /// Output warning message
    void obsolete()
    {
      if (FlagObsoleteCode)
      {
        std::string junk;
        oomph_info << "\n\n--------------------------------------------\n";
        oomph_info << "You are using obsolete code " << std::endl;
        oomph_info << "--------------------------------------------\n\n";
        oomph_info << "Enter: \"s\" to suppress further messages" << std::endl;
        oomph_info << "       \"k\" to crash the code to allow a trace back in "
                      "the debugger"
                   << std::endl;

        oomph_info << "       any other key to continue\n \n";
        oomph_info << "                    [Note: Insert \n \n ";
        oomph_info << "                            "
                      "ObsoleteCode::FlagObsoleteCode=false;\n \n";
        oomph_info << "                     into your code to suppress these "
                      "messages \n";
        oomph_info << "                     altogether.] \n";

        std::cin >> junk;
        if (junk == "s")
        {
          FlagObsoleteCode = false;
        }
        if (junk == "k")
        {
          throw OomphLibError(
            "Killed", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
      }
    }


    /// Ouput a warning message with a string argument
    void obsolete(const std::string& message)
    {
      if (FlagObsoleteCode)
      {
        oomph_info << "\n\n------------------------------------" << std::endl;
        oomph_info << message << std::endl;
        oomph_info << "----------------------------------------" << std::endl;

        obsolete();
      }
    }

  } // namespace ObsoleteCode


  //====================================================================
  /// Namespace for tecplot stuff
  //====================================================================
  namespace TecplotNames
  {
    /// Tecplot colours
    Vector<std::string> colour;


    /// Setup namespace
    void setup()
    {
      colour.resize(5);
      colour[0] = "RED";
      colour[1] = "GREEN";
      colour[2] = "BLUE";
      colour[3] = "CYAN";
      colour[4] = "BLACK";
    }


  } // namespace TecplotNames


#ifdef LEAK_CHECK

  //====================================================================
  /// Namespace for leak check: Keep a running count of all instantiated
  /// objects -- add your own if you want to...
  //====================================================================
  namespace LeakCheckNames
  {
    long QuadTree_build;
    long OcTree_build;
    long QuadTreeForest_build;
    long OcTreeForest_build;
    long RefineableQElement<2> _build;
    long RefineableQElement<3> _build;
    long MacroElement_build;
    long HangInfo_build;
    long Node_build;
    long GeomReference_build;
    long AlgebraicNode_build;

    void reset()
    {
      QuadTree_build = 0;
      OcTree_build = 0;
      QuadTreeForest_build = 0;
      OcTreeForest_build = 0;
      RefineableQElement<2> _build = 0;
      RefineableQElement<3> _build = 0;
      MacroElement_build = 0;
      HangInfo_build = 0;
      Node_build = 0;
      GeomReference_build = 0;
      AlgebraicNode_build = 0;
    }

    void doc()
    {
      oomph_info << "\n Leak check: # of builds - # of deletes for the "
                    "following objects: \n\n";
      oomph_info << "LeakCheckNames::QuadTree_build "
                 << LeakCheckNames::QuadTree_build << std::endl;
      oomph_info << "LeakCheckNames::QuadTreeForest_build "
                 << LeakCheckNames::QuadTreeForest_build << std::endl;
      oomph_info << "LeakCheckNames::OcTree_build "
                 << LeakCheckNames::OcTree_build << std::endl;
      oomph_info << "LeakCheckNames::OcTreeForest_build "
                 << LeakCheckNames::OcTreeForest_build << std::endl;
      oomph_info << "LeakCheckNames::RefineableQElement<2>_build "
                 << LeakCheckNames::RefineableQElement<2> _build << std::endl;
      oomph_info << "LeakCheckNames::RefineableQElement<3>_build "
                 << LeakCheckNames::RefineableQElement<3> _build << std::endl;
      oomph_info << "LeakCheckNames::MacroElement_build "
                 << LeakCheckNames::MacroElement_build << std::endl;
      oomph_info << "LeakCheckNames::HangInfo_build "
                 << LeakCheckNames::HangInfo_build << std::endl;
      oomph_info << "LeakCheckNames::Node_build " << LeakCheckNames::Node_build
                 << std::endl;
      oomph_info << "LeakCheckNames::GeomReference_build "
                 << LeakCheckNames::GeomReference_build << std::endl;
      oomph_info << "LeakCheckNames::AlgebraicNode_build "
                 << LeakCheckNames::AlgebraicNode_build << std::endl;
      oomph_info << std::endl;
    }


  } // namespace LeakCheckNames


#endif


  /// ////////////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////////////


  //====================================================================
  /// Namespace for pause() command
  //====================================================================
  namespace PauseFlags
  {
    /// Flag to enable pausing code -- pause the code by default
    bool PauseFlag = true;

  } // namespace PauseFlags

  //======================================================================
  /// Pause and display message
  //======================================================================
  void pause(std::string message)
  {
    std::string junk;
    if (PauseFlags::PauseFlag)
    {
      oomph_info << message << "\n hit any key to continue [hit \"S\" "
                 << "to suppress further interruptions]\n";
      std::cin >> junk;
      if (junk == "S")
      {
        PauseFlags::PauseFlag = false;
      }
    }
    else
    {
      oomph_info << "\n[Suppressed pause message] \n";
    }
  }


  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////////

  //=============================================================================
  /// Helper for recording execution time.
  //=============================================================================
  namespace TimingHelpers
  {
    /// returns the time in seconds after some point in past
    double timer()
    {
#ifdef OOMPH_HAS_MPI
      if (MPI_Helpers::mpi_has_been_initialised())
      {
        return MPI_Wtime();
      }
      else
#endif
      {
        time_t t = clock();
        return double(t) / double(CLOCKS_PER_SEC);
      }
    }

    /// Returns a nicely formatted string from an input time in seconds;
    /// the format depends on the size of time, e.g.:
    /// 86510 will be printed as 1d 1m:50
    ///  3710 will be printed as 1h:01:50
    ///   700 will be printed as 11m:40
    ///    59 will be printed as 59s
    std::string convert_secs_to_formatted_string(const double& time_in_sec)
    {
      std::ostringstream ss;

      unsigned sec_within_day = unsigned(time_in_sec) % (3600 * 24);

      unsigned days = unsigned(time_in_sec) / (3600 * 24);
      unsigned hours = sec_within_day / 3600;
      unsigned minutes = (sec_within_day % 3600) / 60;
      unsigned seconds = (sec_within_day % 3600) % 60;

      if (days > 0)
      {
        ss << days << "d ";
      }

      if (hours > 0)
      {
        ss << hours << "h:";
        ss << std::setw(2) << std::setfill('0');
        ss << minutes << ":";
        ss << seconds;
      }
      else if (minutes > 0)
      {
        ss << minutes << "m:";
        ss << std::setw(2) << std::setfill('0');
        ss << seconds;
      }
      else
      {
        double seconds_double = seconds + (time_in_sec - double(seconds));
        ss << std::setprecision(1) << std::fixed << seconds_double << "s";
      }

      return ss.str();
    }
  } // end of namespace TimingHelpers


  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////


  //===============================================================
  /// Namespace with helper functions to assess total memory usage
  /// on the fly using system() -- details are very machine specific! This just
  /// provides the overall machinery with default settings for
  /// our own (linux machines). Uses the system command
  /// to spawn a command that computes the total memory usage
  /// on the machine where this is called.  [Disclaimer: works
  /// on my machine(s) -- no guarantees for any other platform;
  /// linux or not. MH]
  //===============================================================
  namespace MemoryUsage
  {
    /// Boolean to suppress synchronisation of doc memory
    /// usage on processors (via mpi barriers). True (i.e. sync is
    /// suppressed) by default because not all processors may
    /// reach the relevant doc memory usage statements
    /// causing the code to hang).
    bool Suppress_mpi_synchronisation = true;

    /// String containing system command that obtains memory usage
    /// of all processes.
    /// Default assignment for linux. [Disclaimer: works on my machine(s) --
    /// no guarantees for any other platform; linux or not. MH]
    std::string My_memory_usage_system_string = "ps aux";

    /// Bool allowing quick bypassing of ALL operations related
    /// to memory usage monitoring -- this allows the code to remain
    /// "instrumented" without incurring the heavy penalties associated
    /// with the system calls and i/o. Default setting: false.
    bool Bypass_all_memory_usage_monitoring = false;

    /// String containing name of file in which we document
    /// my memory usage -- you may want to change this to allow different
    /// processors to write to separate files (especially in mpi
    /// context). Note that file is appended to
    /// so it ought to be emptied (either manually or by calling
    /// helper function empty_memory_usage_file()
    std::string My_memory_usage_filename = "my_memory_usage.dat";

    /// Function to empty file that records my memory usage in
    /// file whose name is specified by My_memory_usage_filename.
    void empty_my_memory_usage_file()
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Open without appending and write header
      std::ofstream the_file;
      the_file.open(My_memory_usage_filename.c_str());
      the_file << "# My memory usage: \n";
      the_file.close();
    }


#ifdef OOMPH_HAS_UNISTDH

    /// Doc my memory usage, prepended by string (which allows
    /// identification from where the function is called, say) that records
    /// memory usage in file whose name is specified by
    /// My_memory_usage_filename. Data is appended to that file; wipe it with
    /// empty_my_memory_usage_file(). Note: This requires getpid() which is
    /// defined in unistd.h so if you don't have that we won't build that
    /// function!
    void doc_my_memory_usage(const std::string& prefix_string)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Write prefix
      std::ofstream the_file;
      the_file.open(My_memory_usage_filename.c_str(), std::ios::app);
      the_file << prefix_string << " ";
      the_file.close();

      // Sync all processors if in parallel (unless suppressed)
#ifdef OOMPH_HAS_MPI
      if ((MPI_Helpers::mpi_has_been_initialised()) &&
          (!Suppress_mpi_synchronisation))
      {
        MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
      }
#endif

      // Process memory usage command and write to file
      std::stringstream tmp;
      tmp << My_memory_usage_system_string << " | awk '{if ($2==" << getpid()
          << "){print $4}}' >> " << My_memory_usage_filename;
      int success = system(tmp.str().c_str());

      // Dummy command to shut up compiler warnings
      success += 1;
    }

#endif

    /// String containing system command that obtains total memory usage.
    /// Default assignment for linux. [Disclaimer: works on my machine(s) --
    /// no guarantees for any other platform; linux or not. MH]
    std::string Total_memory_usage_system_string =
      "ps aux | awk 'BEGIN{sum=0}{sum+=$4}END{print sum}'";

    /// String containing name of file in which we document total
    /// memory usage -- you may want to change this to allow different
    /// processors to write to separate files (especially in mpi
    /// context). Note that file is appended to
    /// so it ought to be emptied (either manually or by calling
    /// helper function empty_memory_usage_file()
    std::string Total_memory_usage_filename = "memory_usage.dat";

    /// Function to empty file that records total memory usage in
    /// file whose name is specified by Total_memory_usage_filename.
    void empty_total_memory_usage_file()
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Open without appending and write header
      std::ofstream the_file;
      the_file.open(Total_memory_usage_filename.c_str());
      the_file << "# Total memory usage: \n";
      the_file.close();
    }

    /// Doc total memory usage, prepended by string (which allows
    /// identification from where the function is called, say) that records
    /// memory usage in file whose name is specified by
    /// Total_memory_usage_filename. Data is appended to that file; wipe it with
    /// empty_memory_usage_file().
    void doc_total_memory_usage(const std::string& prefix_string)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Write prefix
      std::ofstream the_file;
      the_file.open(Total_memory_usage_filename.c_str(), std::ios::app);
      the_file << prefix_string << " ";
      the_file.close();

      // Sync all processors if in parallel
#ifdef OOMPH_HAS_MPI
      if ((MPI_Helpers::mpi_has_been_initialised()) &&
          (!Suppress_mpi_synchronisation))
      {
        MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
      }
#endif

      // Process memory usage command and write to file
      std::stringstream tmp;
      tmp << Total_memory_usage_system_string << " >> "
          << Total_memory_usage_filename;
      int success = system(tmp.str().c_str());

      // Dummy command to shut up compiler warnings
      success += 1;
    }


    /// Function to empty file that records total and local memory usage
    /// in appropriate files
    void empty_memory_usage_files()
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      empty_my_memory_usage_file();
      empty_total_memory_usage_file();
      empty_top_file();
    }

    /// Doc total and local memory usage, prepended by string (which
    /// allows identification from where the function is called, say). NOTE:
    /// Local memory usage only works if we have unistd.h header
    void doc_memory_usage(const std::string& prefix_string)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

#ifdef OOMPH_HAS_UNISTDH
      doc_my_memory_usage(prefix_string);
#endif

      doc_total_memory_usage(prefix_string);
    }


    /// String containing system command that runs "top" (or equivalent)
    /// "indefinitely" and writes to file specified in Top_output_filename.
    /// Default assignment for linux. [Disclaimer: works on my machine(s) --
    /// no guarantees for any other platform; linux or not. MH]
    std::string Top_system_string = "while true; do top -b -n 2 ; done  ";

    ///  String containing name of file in which we document "continuous"
    /// output from "top" (or equivalent)-- you may want to change this to
    /// allow different processors to write to separate files (especially in mpi
    /// context). Note that file is appended to
    /// so it ought to be emptied (either manually or by calling
    /// helper function empty_top_file()
    std::string Top_output_filename = "top_output.dat";

    /// Function to empty file that records continuous output from top in
    /// file whose name is specified by Top_output_filename
    void empty_top_file()
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Open without appending and write header
      std::ofstream the_file;
      the_file.open(Top_output_filename.c_str());
      the_file << "# Continuous output from top obtained with: \n";
      the_file << "# " << Top_system_string << "\n";
      the_file.close();
    }


    /// Start running top continuously and output (append) into
    /// file specified by Top_output_filename. Wipe that file  with
    /// empty_top_file() if you wish. Note that this is again quite Linux
    /// specific and unlikely to work on other operating systems. Insert
    /// optional comment into output file before starting.
    void run_continous_top(const std::string& comment)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Sync all processors if in parallel
      std::string modifier = "";

#ifdef OOMPH_HAS_MPI
      if (MPI_Helpers::mpi_has_been_initialised())
      {
        if (!Suppress_mpi_synchronisation)
        {
          MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
        }
        std::stringstream tmp;
        tmp << "_proc" << MPI_Helpers::communicator_pt()->my_rank();
        modifier = tmp.str();
      }
#endif

      // Process memory usage command and write to file
      std::stringstream tmp;

      // String stream seems unhappy about direct insertions of these
      std::string backslash = "\\";
      std::string dollar = "$";

      // Create files that spawn and kill continuous top
      tmp << "echo \"#/bin/bash\" > run_continuous_top" << modifier << ".bash; "
          << "echo \" echo " << backslash << "\" kill " << backslash << dollar
          << backslash << dollar << " " << backslash
          << "\" > kill_continuous_top" << modifier
          << ".bash; chmod a+x kill_continuous_top" << modifier << ".bash; "
          << Top_system_string << " \" >> run_continuous_top" << modifier
          << ".bash; chmod a+x run_continuous_top" << modifier << ".bash ";
      int success = system(tmp.str().c_str());

      // Add comment to annotate?
      if (comment != "")
      {
        insert_comment_to_continous_top(comment);
      }

      // Start spawning
      std::stringstream tmp2;
      tmp2 << "./run_continuous_top" << modifier << ".bash  >> "
           << Top_output_filename << " & ";
      success = system(tmp2.str().c_str());

      // Dummy command to shut up compiler warnings
      success += 1;
    }


    /// Stop running top continuously. Note that this is
    /// again quite Linux specific and unlikely to work on other operating
    /// systems. Insert optional comment into output file before stopping.
    void stop_continous_top(const std::string& comment)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      // Sync all processors if in parallel
      std::string modifier = "";

#ifdef OOMPH_HAS_MPI
      if (MPI_Helpers::mpi_has_been_initialised())
      {
        if (!Suppress_mpi_synchronisation)
        {
          MPI_Barrier(MPI_Helpers::communicator_pt()->mpi_comm());
        }
        std::stringstream tmp;
        tmp << "_proc" << MPI_Helpers::communicator_pt()->my_rank();
        modifier = tmp.str();
      }
#endif

      // Add comment to annotate?
      if (comment != "")
      {
        insert_comment_to_continous_top(comment);
      }

      // Kill
      std::stringstream tmp2;
      tmp2 << "./kill_continuous_top" << modifier << ".bash  >> "
           << Top_output_filename << " & ";
      int success = system(tmp2.str().c_str());

      // Dummy command to shut up compiler warnings
      success += 1;
    }


    /// Insert comment into running continuous top output
    void insert_comment_to_continous_top(const std::string& comment)
    {
      // bail out straight away?
      if (Bypass_all_memory_usage_monitoring) return;

      std::stringstream tmp;
      tmp << " echo \"OOMPH-LIB EVENT: " << comment << "\"  >> "
          << Top_output_filename;
      int success = system(tmp.str().c_str());

      // Dummy command to shut up compiler warnings
      success += 1;
    }


  } // end of namespace MemoryUsage


} // namespace oomph
