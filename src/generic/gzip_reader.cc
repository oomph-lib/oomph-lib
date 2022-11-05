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

#include "gzip_reader.h"
#include "oomph_definitions.h"
#include "oomph_zlib/zlib.h"

#include <algorithm>
#include <fstream>
#include <vector>

namespace oomph
{
  /************************************************************************************
   * @brief Construct a new GZipReader object
   *
   ************************************************************************************/
  GZipReader::GZipReader() : Gzip_file(nullptr), Chunk_size(256) {}

  /************************************************************************************
   * @brief Construct a new GZipReader object
   *
   * @param filename:
   ************************************************************************************/
  GZipReader::GZipReader(const std::string& filename,
                         const unsigned& chunk_size)
    : Gzip_file(nullptr), Chunk_size(chunk_size)
  {
    this->open(filename);
    if (!this->is_open())
    {
      throw OomphLibError("Unable to open file '" + filename + "'!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
  }

  /************************************************************************************
   * @brief Destroy the GZipReader object
   *
   ************************************************************************************/
  GZipReader::~GZipReader()
  {
    if (this->is_open())
    {
      this->close();
    }
  }

  /************************************************************************************
   * @brief
   *
   * @param filename:
   ************************************************************************************/
  void GZipReader::open(const std::string& filename)
  {
    Gzip_file = gzopen(filename.c_str(), "rb");
  }

  /************************************************************************************
   * @brief
   *
   ************************************************************************************/
  void GZipReader::close()
  {
    if (!this->is_open())
    {
      return;
    }
    gzclose(Gzip_file);
    Gzip_file = nullptr;
  }

  /************************************************************************************
   * @brief
   *
   * @return true:
   * @return false:
   ************************************************************************************/
  bool GZipReader::is_open() const
  {
    if (Gzip_file != nullptr)
    {
      return true;
    }
    return false;
  }

  /************************************************************************************
   * @brief
   *
   * TODO: Update to using std::optional<std::string> if C++17 is allowed
   *
   * @return std::string:
   ************************************************************************************/
  Packet GZipReader::getline()
  {
    // Storage for the loaded line and flag to indicate whether we were able to
    // read something
    std::string line{};
    bool found{false};

    // Early exit if the file is closed
    if (!is_open())
    {
      return Packet{line, found};
    }

    // Allocate space for a buffer on the stack
    char buffer_pt[Chunk_size];

    // Pointer to the end of the buffer
    char* buffer_end_pt = buffer_pt + Chunk_size;

    // Flag to indicate whether we've reached the end of the line
    bool reached_end_of_line{false};

    // Repeatedly read chunks until we reach the end of the line
    while (!reached_end_of_line)
    {
      // EOF condition or error
      if (!gzgets(Gzip_file, buffer_pt, Chunk_size))
      {
        this->close();
        break;
      }

      // No EOF or error
      found = true;

      // Find the extent of the loaded string
      char* line_end_pt = std::find(buffer_pt, buffer_end_pt, '\0');

      // Append the loaded string
      line += {buffer_pt, line_end_pt};

      // Check if the last loaded non-null character is a newline character
      if ((line_end_pt != buffer_pt) && (*(line_end_pt - 1) == '\n'))
      {
        reached_end_of_line = true;
      }
    }

    // Strip the newline character, if there is one
    if (found && (!line.empty()) && (line[line.length() - 1] == '\n'))
    {
      line.erase(line.length() - 1);
    }
    return Packet{line, found};
  }

  /************************************************************************************
   * @brief
   *
   * @return std::vector<std::string>:
   ************************************************************************************/
  std::vector<std::string> GZipReader::read_all()
  {
    // Early exit if the file is closed
    if (!this->is_open())
    {
      return {};
    }

    // Storage for the loaded data
    std::vector<std::string> file_contents;
    while (this->is_open())
    {
      Packet packet = this->getline();
      if (packet.found)
      {
        file_contents.emplace_back(packet.data);
      }
    }
    return file_contents;
  }

} // namespace oomph