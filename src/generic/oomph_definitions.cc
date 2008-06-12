//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
#include "oomph_definitions.h"


namespace oomph
{

//========================================================================
/// The OomphLibException constructor takes the error description, 
/// function name,  a location string provided by the 
/// OOMPH_EXCEPTION_LOCATION and an exception type "WARNING" or "ERROR" 
/// and combines them into a standard error message that is written to the
/// exception stream. The output_width of the message can also be specified
//==========================================================================
OomphLibException::OomphLibException(const std::string &error_description,
                                     const std::string &function_name,
                                     const char *location,
                                     const std::string &exception_type,
                                     std::ostream &exception_stream,
                                     const unsigned &output_width) : 
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
 
 //Finish off with another set of double lines
 for(unsigned i=0;i<output_width;i++) {exception_stream << "=";}
 exception_stream << std::endl << std::endl;
 
 //Flush the stream buffer
 exception_stream.flush();
}

//========================================================================
/// Default output stream for OomphLibErorrs (cout)
//========================================================================
std::ostream *OomphLibError::Stream_pt = &std::cout;

//=======================================================================
/// Default output width for OomphLibErrors (70)
//=======================================================================
unsigned OomphLibError::Output_width = 70;

//=======================================================================
/// Default output stream for OomphLibWarnings(cout)
//=======================================================================
std::ostream *OomphLibWarning::Stream_pt = &std::cout;

//=======================================================================
/// Default output width for OomphLibWarnings (70)
//=======================================================================
unsigned OomphLibWarning::Output_width = 70;


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
