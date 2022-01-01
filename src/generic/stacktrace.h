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
//====================================================================
// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

// Modified by MH to take proper C++ output stream and to produce
// code for retrieval of line numbers. Also obtains name of
// executable from oomph-lib CommandLineArgs namespace (if it's been
// set up)
//====================================================================
#ifndef _STACKTRACE_H_
#define _STACKTRACE_H_

#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <cxxabi.h>

// Include oomph-utilities to be able to obtain name of executable
#include "oomph_utilities.h"

/** Print a demangled stack backtrace of the caller function */
static inline void print_stacktrace(std::ostream& exception_stream,
                                    unsigned int max_frames = 63)
{
  exception_stream << "\nStack trace:\n";
  exception_stream << "------------\n";

  // storage array for stack trace address data
  void* addrlist[max_frames + 1];

  // retrieve current stack addresses
  int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

  if (addrlen == 0)
  {
    exception_stream << "\n<empty stack trace, possibly corrupt>\n";
    return;
  }

  // resolve addresses into strings containing "filename(function+address)",
  // this array must be free()-ed
  char** symbollist = backtrace_symbols(addrlist, addrlen);

  // allocate string which will be filled with the demangled function name
  size_t funcnamesize = 256;
  char* funcname = (char*)malloc(funcnamesize);

  // Stream to collect information that allows retrieval of line numbers
  std::stringstream address_stream;

  // iterate over the returned symbol lines. skip the first, it is the
  // address of this function.
  for (int i = 1; i < addrlen; i++)
  {
    // ------------- begin  mh addition ---------------
    char *begin_absolute = 0, *end_absolute = 0;

    // find absolute address (in square brackets)
    // ./module(function+0x15c) [0x8048a6d]
    for (char* p = symbollist[i]; *p; ++p)
    {
      if (*p == '[') begin_absolute = p;
      else if (*p == ']' && begin_absolute)
      {
        end_absolute = p;
        break;
      }
    }

    if (begin_absolute)
    {
      *begin_absolute++ = '\0';
      *end_absolute = '\0';
      if (oomph::CommandLineArgs::Argc != 0)
      {
        address_stream << "addr2line -e " << oomph::CommandLineArgs::Argv[0]
                       << " " << begin_absolute << "\n";
      }
      else
      {
        address_stream << "addr2line -e [name_of_executable] " << begin_absolute
                       << "\n";
      }
    }

    // ------------- end mh addition ---------------

    char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

    // find parentheses and +address offset surrounding the mangled name:
    // ./module(function+0x15c) [0x8048a6d]
    for (char* p = symbollist[i]; *p; ++p)
    {
      if (*p == '(') begin_name = p;
      else if (*p == '+')
        begin_offset = p;
      else if (*p == ')' && begin_offset)
      {
        end_offset = p;
        break;
      }
    }

    if (begin_name && begin_offset && end_offset && begin_name < begin_offset)
    {
      *begin_name++ = '\0';
      *begin_offset++ = '\0';
      *end_offset = '\0';

      // mangled name is now in [begin_name, begin_offset) and caller
      // offset in [begin_offset, end_offset). now apply
      // __cxa_demangle():

      int status;
      char* ret =
        abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
      if (status == 0)
      {
        funcname = ret; // use possibly realloc()-ed string
        /* fprintf(out, "  %s : %s+%s\n", */
        /*         symbollist[i], funcname, begin_offset); */
        exception_stream << symbollist[i] << " : " << funcname << "+"
                         << begin_offset << std::endl;
      }
      else
      {
        // demangling failed. Output function name as a C function with
        // no arguments.
        /* fprintf(out, "  %s : %s()+%s\n", */
        /*         symbollist[i], begin_name, begin_offset); */

        exception_stream << symbollist[i] << " : " << begin_name << "+"
                         << begin_offset << std::endl;
      }
    }
    else
    {
      // couldn't parse the line? print the whole line.
      // fprintf(out, "  %s\n", symbollist[i]);
      exception_stream << symbollist[i] << std::endl;
    }
  }

  exception_stream << "\nHere are the commmands to obtain the line numbers:\n";

  exception_stream << "--------------------------------------------------\n";
  exception_stream << address_stream.str() << std::endl;


  if (oomph::CommandLineArgs::Argc == 0)
  {
    exception_stream << "\nNOTE: Replace [name_of_executable] by the actual\n"
                     << "       name of the executable. I would have inserted\n"
                     << "       this for you if you had called \n\n"
                     << "          CommandLineArgs::setup(argc,argv);\n\n"
                     << "       in your driver code...\n\n";
  }

  free(funcname);
  free(symbollist);
}

#endif // _STACKTRACE_H_
