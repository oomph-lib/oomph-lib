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
#ifndef OOMPH_FPDIFF_HEADER
#define OOMPH_FPDIFF_HEADER

#include <fstream>
#include <ostream>
#include <string>
#include <vector>

// Config header
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

namespace oomph
{
  /************************************************************************************
   * @brief
   *
   * @param filename
   * @return std::vector<std::string>
   ************************************************************************************/
  std::vector<std::string> gzip_load(const std::string& filename);

  /************************************************************************************
   * @brief
   *
   * @param filename1
   * @param filename2
   * @param outstream
   * @param relative_error: max relative difference in percent.
   * @param small: a small number -- essentially round-off error.
   * @param details_stream
   * @return true
   * @return false
   ************************************************************************************/
  int fpdiff(const std::string& filename1,
             const std::string& filename2,
             std::ostream& outstream,
             const double& relative_error = 1.0e-01,
             const double& small = 1.0e-14);

  /************************************************************************************
   * @brief
   *
   * @param filename1
   * @param filename2
   * @param outstream
   * @param relative_error: max relative difference in percent.
   * @param small: a small number -- essentially round-off error.
   * @param details_stream
   * @return true
   * @return false
   ************************************************************************************/
  int fpdiff(const std::string& filename1,
             const std::string& filename2,
             const std::string& log_file,
             const double& relative_error = 1.0e-01,
             const double& small = 1.0e-14);
} // namespace oomph

#endif