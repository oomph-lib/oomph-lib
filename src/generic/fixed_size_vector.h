// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

#ifndef FIXED_SIZE_VECTOR_HEADER
#define FIXED_SIZE_VECTOR_HEADER

// Generic oomph-lib header
#include "generic.h"

namespace oomph
{


  //============================================================
  /// Class for fixed-size vector storing N objects of type T.
  /// No resize or other STL functionality but (therefore)
  /// much less run-time overhead.
  //============================================================
  template<class T, unsigned N>
  class FixedSizeVector
  {
  public:
    /// Constructor; no initialisiation
    FixedSizeVector() {}

    /// Constructor with initialisation to value
    FixedSizeVector(const T& value)
    {
      fill(value);
    }

    /// Fill all entries with specified value
    void fill(const T& value)
    {
      for (unsigned i = 0; i < N; i++)
      {
        Data[i] = value;
      }
    }

    /// Size of vector
    unsigned size() const
    {
      return N;
    }

    /// Read/write access to i-th entry
    T& operator[](const unsigned& i)
    {
#ifdef RANGE_CHECKING
      if (i >= N)
      {
        std::stringstream error_message;
        error_message << "Range error: trying to access entry " << i
                      << " in a fixed-size vector of size " << N << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return Data[i];
    }


    /// Read access to i-th entry
    const T& operator[](const unsigned& i) const
    {
#ifdef RANGE_CHECKING
      if (i >= N)
      {
        std::stringstream error_message;
        error_message << "Range error: trying to access entry " << i
                      << " in a fixed-size vector of size " << N << std::endl;
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return Data[i];
    }


  private:
    /// The data, stored as a raw C-style array
    T Data[N];
  };


  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  //============================================================
  /// Class for fixed-size matrix storing NROW x NCOL objects
  /// of type T.
  ///
  /// No resize or other STL functionality but therefore
  /// much less run-time overhead.
  ///
  /// Storage is row-major:
  ///
  /// (0,0), (0,1), ..., (0,NCOL-1),
  /// (1,0), (1,1), ..., (1,NCOL-1),
  /// ...
  /// (NROW-1,NCOL-1)
  //============================================================
  template<class T, unsigned NROW, unsigned NCOL>
  class FixedSizeMatrix
  {
  public:
    /// Default constructor
    FixedSizeMatrix() {}

    /// Construct and fill all entries with specified value
    FixedSizeMatrix(const T& value)
    {
      fill(value);
    }

    // hierher copy into vector.

    /// Fill all entries with specified value
    void fill(const T& value)
    {
      for (unsigned i = 0; i < NROW * NCOL; i++)
      {
        Data[i] = value;
      }
    }

    /// Number of rows
    unsigned nrow() const
    {
      return NROW;
    }

    /// Number of columns
    unsigned ncol() const
    {
      return NCOL;
    }

    /// Read/write access to (i,j)-th entry
    T& operator()(const unsigned& i, const unsigned& j)
    {
#ifdef RANGE_CHECKING
      if ((i >= NROW) || (j >= NCOL))
      {
        std::stringstream error_message;

        error_message << "Range error: trying to access entry (" << i << ","
                      << j << ") "
                      << "in a fixed-size matrix of size " << NROW << " x "
                      << NCOL << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return Data[i * NCOL + j];
    }

    /// Read access to (i,j)-th entry
    const T& operator()(const unsigned& i, const unsigned& j) const
    {
#ifdef RANGE_CHECKING
      if ((i >= NROW) || (j >= NCOL))
      {
        std::stringstream error_message;

        error_message << "Range error: trying to access entry (" << i << ","
                      << j << ") "
                      << "in a fixed-size matrix of size " << NROW << " x "
                      << NCOL << std::endl;

        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif

      return Data[i * NCOL + j];
    }

    /// Fill all entries with specified value
    void initialise(const T& value)
    {
      for (unsigned i = 0; i < NROW * NCOL; i++)
      {
        Data[i] = value;
      }
    }

  private:
    /// Matrix entries stored in row-major order
    T Data[NROW * NCOL];
  };


} // namespace oomph

#endif
