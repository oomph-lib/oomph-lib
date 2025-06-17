// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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

#ifndef OOMPH_GZIP_READER_HEADER
#define OOMPH_GZIP_READER_HEADER

#include "oomph_zlib/zlib.h"

#include <string>
#include <memory>
#include <vector>

namespace oomph
{
  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  struct Packet
  {
    std::string data;
    bool found;
  };

  /************************************************************************************
   * @brief
   *
   * @tparam _CharT
   * @tparam _Traits
   ************************************************************************************/
  class GZipReader
  {
  public:
    GZipReader();
    GZipReader(const std::string& filename, const unsigned& chunk_size = 256);
    ~GZipReader();
    void open(const std::string& filename);
    void close();
    bool is_open() const;
    Packet getline();
    std::vector<std::string> read_all();

  private:
    gzFile Gzip_file;
    unsigned Chunk_size;
  };

} // namespace oomph

#endif
