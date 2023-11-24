// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Non-inline and static member functions for the Oomph-lib
// exception handlers

#ifdef OOMPH_HAS_STACKTRACE
#include "stacktrace.h"
#endif
#include "oomph_definitions.h"

namespace oomph
{
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  /// =====================================================================
  /// Namespace to control level of comprehensive timings
  //======================================================================
  namespace Global_timings
  {
    /// Global boolean to switch on comprehensive timing -- can
    /// probably be declared const false when development on hector
    /// is complete
    bool Doc_comprehensive_timings = false;
  }; // namespace Global_timings

  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////
  /// //////////////////////////////////////////////////////////////////

  //=======================================================================
  /// Helper namespace for set_terminate function -- used to spawn
  /// messages from uncaught errors
  /// =======================================================================
  namespace TerminateHelper
  {
    /// Setup terminate helper
    void setup()
    {
      if (Exception_stringstream_pt != 0) delete Exception_stringstream_pt;
      Exception_stringstream_pt = new std::stringstream;
      std::set_terminate(spawn_errors_from_uncaught_errors);
    }

    /// Flush string stream of error messages (call when error has been
    /// caught)
    void suppress_exception_error_messages()
    {
      delete Exception_stringstream_pt;
      Exception_stringstream_pt = new std::stringstream;
    }

    /// Function to spawn messages from uncaught errors
    void spawn_errors_from_uncaught_errors()
    {
      (*Error_message_stream_pt) << (*Exception_stringstream_pt).str();
    }

    /// Clean up function that deletes anything dynamically allocated
    /// in this namespace
    void clean_up_memory()
    {
      // If it's a null pointer
      if (Exception_stringstream_pt != 0)
      {
        // Delete it
        delete Exception_stringstream_pt;

        // Make it a null pointer
        Exception_stringstream_pt = 0;
      }
    } // End of clean_up_memory

    /// Stream to output error messages
    std::ostream* Error_message_stream_pt = &std::cerr;

    /// String stream that records the error message
    std::stringstream* Exception_stringstream_pt = 0;
  } // namespace TerminateHelper

  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  /// =====================================================================
  /// A class for handling oomph-lib run-time exceptions quietly.
  //======================================================================
  OomphLibQuietException::OomphLibQuietException() : std::runtime_error("") {}


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //========================================================================
  /// The OomphLibException destructor actually spawns the error message
  /// created in the constructor (unless suppresed)
  //==========================================================================
  OomphLibException::~OomphLibException() throw()
  {
    if (!Suppress_error_message)
    {
      (*Exception_stream_pt) << (*Exception_stringstream_pt).str();
    }
    delete Exception_stringstream_pt;
    Exception_stringstream_pt = 0;
  }

  //========================================================================
  /// The OomphLibException constructor takes the error description,
  /// function name,  a location string provided by the
  /// OOMPH_EXCEPTION_LOCATION and an exception type "WARNING" or "ERROR"
  /// and combines them into a standard error message that is written to the
  /// exception stream. The output_width of the message can also be specified.
  /// Optionally provide a traceback of the function calls.
  //==========================================================================
  OomphLibException::OomphLibException(const std::string& error_description,
                                       const std::string& function_name,
                                       const char* location,
                                       const std::string& exception_type,
                                       std::ostream& exception_stream,
                                       const unsigned& output_width,
                                       bool list_trace_back)
    : std::runtime_error("OomphException")
  {
    // By default we shout
    Suppress_error_message = false;

    // Store exception stream
    Exception_stream_pt = &exception_stream;

    // Create storage for error message
    Exception_stringstream_pt = new std::stringstream;

    // Build an exception header string from the information passed
    // Start with a couple of new lines to space things out
    std::string exception_header = "\n\n";

    // Now add a dividing line
    for (unsigned i = 0; i < output_width; i++)
    {
      exception_header += "=";
    }
    exception_header += "\n";

    // Write the type of exception
    exception_header += "Oomph-lib ";
    exception_header += exception_type;

    // Add the function in which it occurs
    exception_header += "\n\n at ";
    exception_header += location;
    exception_header += "\n\n in ";
    exception_header += function_name;

    // Finish with two new lines
    exception_header += "\n\n";

    // and a closing line
    for (unsigned i = 0; i < (unsigned)(0.8 * output_width); i++)
    {
      exception_header += "-";
    }

    // Output the error header to the stream
    (*Exception_stringstream_pt) << exception_header << std::endl;

    // Report the error
    (*Exception_stringstream_pt) << std::endl << error_description << std::endl;

#ifdef OOMPH_HAS_STACKTRACE
    // Print the stacktrace
    if (list_trace_back)
    {
      print_stacktrace((*Exception_stringstream_pt));
    }
#endif

    // Finish off with another set of double lines
    for (unsigned i = 0; i < output_width; i++)
    {
      (*Exception_stringstream_pt) << "=";
    }
    (*Exception_stringstream_pt) << std::endl << std::endl;

    // Copy message to stream in terminate helper in case the message
    // doesn't get caught and/or doesn/t make it to the destructor
    (*TerminateHelper::Exception_stringstream_pt)
      << (*Exception_stringstream_pt).str();
  }

  //========================================================================
  /// Default output stream for OomphLibErorrs (cerr)
  //========================================================================
  std::ostream* OomphLibError::Stream_pt = &std::cerr;

  //=======================================================================
  /// Default output width for OomphLibErrors (70)
  //=======================================================================
  unsigned OomphLibError::Output_width = 70;

  //=======================================================================
  /// Default output stream for OomphLibWarnings(cerr)
  //=======================================================================
  std::ostream* OomphLibWarning::Stream_pt = &std::cerr;

  //=======================================================================
  /// Default output width for OomphLibWarnings (70)
  //=======================================================================
  unsigned OomphLibWarning::Output_width = 70;


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Namespace containing an output stream that can be used for
  /// debugging. Use at your own risk -- global data is evil!
  //=======================================================================
  namespace Global_output_stream
  {
    /// Output stream
    std::ofstream* Outfile = 0;

  } // namespace Global_output_stream


  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Namespace containing a number that can be used to annotate things for
  /// debugging. Use at your own risk -- global data is evil!
  //=======================================================================
  namespace Global_unsigned
  {
    /// The unsigned
    unsigned Number = 0;

  } // namespace Global_unsigned
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////


  //=======================================================================
  /// Namespace containing a vector of strings that can be used to
  /// to store global output modifiers. This is global data
  /// and you use it at your own risk!
  //=======================================================================
  namespace Global_string_for_annotation
  {
    /// Return the i-th string or "" if the relevant string hasn't
    /// been defined
    std::string string(const unsigned& i)
    {
      if (i < String.size())
      {
        return String[i];
      }
      else
      {
        return "";
      }
    }

    /// Storage for strings that may be used for global annotations.
    /// This is global data and you use it at your own risk!
    std::vector<std::string> String;

  } // namespace Global_string_for_annotation


  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////
  /// /////////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Single (global) instantiation of the Nullstream
  //========================================================================
  Nullstream oomph_nullstream;

  //========================================================================
  /// Single (global) instantiation of the OomphInfo object -- this
  /// is used throughout the library as a "replacement" for std::cout
  //========================================================================
  OomphInfo oomph_info;


  //========================================================================
  /// Single global instatiation of the default output modifier.
  //========================================================================
  OutputModifier default_output_modifier;

} // namespace oomph
