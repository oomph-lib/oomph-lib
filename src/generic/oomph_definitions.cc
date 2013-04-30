//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Non-inline and static member functions for the Oomph-lib
//exception handlers

#ifdef OOMPH_HAS_STACKTRACE
#include "stacktrace.h"
#endif
#include "oomph_definitions.h"

namespace oomph
{



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


///=====================================================================
/// Namespace to control level of comprehensive timings
//======================================================================
 namespace Global_timings
  {
   /// \short Global boolean to switch on comprehensive timing -- can 
   /// probably be declared const false when development on hector
   /// is complete
   bool Doc_comprehensive_timings=false;
  };


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


///=====================================================================
/// A class for handling oomph-lib run-time exceptions quietly.
//======================================================================
 OomphLibQuietException::OomphLibQuietException() : 
  std::runtime_error("")
 {}
 

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//========================================================================
/// The OomphLibException constructor takes the error description, 
/// function name,  a location string provided by the 
/// OOMPH_EXCEPTION_LOCATION and an exception type "WARNING" or "ERROR" 
/// and combines them into a standard error message that is written to the
/// exception stream. The output_width of the message can also be specified.
/// Optionally provide a traceback of the function calls.
//==========================================================================
OomphLibException::OomphLibException(const std::string &error_description,
                                     const std::string &function_name,
                                     const char *location,
                                     const std::string &exception_type,
                                     std::ostream &exception_stream,
                                     const unsigned &output_width, 
                                     bool list_trace_back) : 
 std::runtime_error("OomphException")
{
 //Build an exception header string from the information passed
 //Start with a couple of new lines to space things out
 std::string exception_header="\n\n";
 //Now add a dividing line
 for(unsigned i=0;i<output_width;i++) {exception_header += "=";}
 exception_header += "\n";
 
 //Write the type of exception
 exception_header += "Oomph-lib ";
 exception_header += exception_type;
 //Add the function in which it occurs
 exception_header += "\n\n at ";
 exception_header += location;
 exception_header += "\n\n in ";
 exception_header += function_name;
 //Finish with two new lines 
 exception_header +="\n\n";
 //and a closing line
 for(unsigned i=0;i<(unsigned)(0.8*output_width);i++)  
  {exception_header +=  "-";}
 
 //Output the error header to the stream
 exception_stream << exception_header << std::endl;
 //Report the error
 exception_stream << std::endl << error_description << std::endl;

#ifdef OOMPH_HAS_STACKTRACE
 // Print the stacktrace
 if (list_trace_back)
  {
   print_stacktrace(exception_stream);
  }
#endif

 //Finish off with another set of double lines
 for(unsigned i=0;i<output_width;i++) {exception_stream << "=";}
 exception_stream << std::endl << std::endl;
 
 //Flush the stream buffer
 exception_stream.flush();

}

//========================================================================
/// Default output stream for OomphLibErorrs (cerr)
//========================================================================
std::ostream *OomphLibError::Stream_pt = &std::cerr;

//=======================================================================
/// Default output width for OomphLibErrors (70)
//=======================================================================
unsigned OomphLibError::Output_width = 70;

//=======================================================================
/// Default output stream for OomphLibWarnings(cerr)
//=======================================================================
std::ostream *OomphLibWarning::Stream_pt = &std::cerr;

//=======================================================================
/// Default output width for OomphLibWarnings (70)
//=======================================================================
unsigned OomphLibWarning::Output_width = 70;



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Namespace containing an output stream that can be used for
/// debugging. Use at your own risk -- global data is evil!
//=======================================================================
namespace Global_output_stream
{

 /// Output stream
 std::ofstream* Outfile=0;

}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Namespace containing a number that can be used to annotate things for
/// debugging. Use at your own risk -- global data is evil!
//=======================================================================
namespace Global_unsigned
{

 /// The unsigned
 unsigned Number=0;

}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Namespace containing a vector of strings that can be used to
/// to store global output modifiers. This is global data
/// and you use it at your own risk!
//=======================================================================
namespace Global_string_for_annotation
{

 /// \short Return the i-th string or "" if the relevant string hasn't
 /// been defined
 std::string string(const unsigned& i)
 {
  if (i<String.size())
   {
    return String[i];
   }
  else
   {
    return "";
   }
 }

 /// \short Storage for strings that may be used for global annotations.
 /// This is global data and you use it at your own risk!
 std::vector<std::string> String;

}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


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

}
