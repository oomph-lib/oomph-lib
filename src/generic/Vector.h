// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
// A header file that is used to define the oomph-lib Vector class

// Include guards to prevent multiple inclusions of the header
#ifndef OOMPH_VECTOR_HEADER
#define OOMPH_VECTOR_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// Standard library includes
#include <vector>
#include <sstream>
#include <cmath>

// Oomph-lib error handler
#include "oomph_definitions.h"


namespace oomph
{
  //===========================================================================
  /// A slight extension to the standard template vector class so that
  /// we can include "graceful" array range checks if the RANGE_CHECKING
  /// flag is set. The generalisation to general allocators is NOT handled here,
  /// mainly because we never use it, but also because the intel and gnu
  /// compilers have different names for the internal classes, which makes
  /// writing code that works for both a pain!
  //===========================================================================
  template<class _Tp>
  class Vector : public std::vector<_Tp>
  {
  public:
    /// Typedef to make the constructors look a bit cleaner
    typedef _Tp value_type;

    /// Typedef to make the constructors look a bit cleaner
    typedef value_type& reference;

    /// Typedef to make the constructors look a bit cleaner
    typedef const value_type& const_reference;

    /// Typedef to make the constructors look a bit cleaner
    typedef size_t size_type;

// Only include this code, if we are range checking
#ifdef RANGE_CHECKING
  private:
    // Function to return a reference to a vector entry,
    // including array range checking
    reference error_checked_access(size_type __n)
    {
      // If there is an out of range error, die, but issue a warning message
      if (__n >= this->size())
      {
        // Construct an error message, as a string stream
        std::ostringstream error_message;
        if (this->size() == 0)
        {
          error_message
            << "Range Error: Vector is empty but you requested entry " << __n;
        }
        else
        {
          error_message << "Range Error: " << __n << " is not in the range (0,"
                        << this->size() - 1 << ")";
        }

        // Throw an Oomph-lib error
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);

        // This is a dummy return to keep the Intel compiler happy
        return std::vector<_Tp>::operator[](__n);
      }
      else
      {
        return std::vector<_Tp>::operator[](__n);
      }
    }

    // Function to return a constant reference to a vector entry
    // including error range checking
    const_reference error_checked_access(size_type __n) const
    {
      // If there is an out of range error, die, but issue a warning message
      if (__n >= this->size())
      {
        // Construct an error message, as a string stream
        std::ostringstream error_message;
        error_message << "Range Error: " << __n << " is not in the range (0,"
                      << this->size() - 1 << ")";

        // Throw an Oomph-lib error
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);

        // This is a dummy return to keep the Intel compiler happy
        return std::vector<_Tp>::operator[](__n);
      }
      else
      {
        return std::vector<_Tp>::operator[](__n);
      }
    }

#endif

  public:
    // Standard Constuctors (some have been omitted from the stl classes)

    /// Construct an empty vector
    Vector() : std::vector<_Tp>() {}

    /// A constructor that creates a vector of size __n.
    /// Note the use of explicit for "strong" type checking
    explicit Vector(size_type __n) : std::vector<_Tp>(__n) {}

    /// A constructor that creates a vector of size __n and
    /// initialises every entry to __value
    Vector(size_type __n, const _Tp& __value) : std::vector<_Tp>(__n, __value)
    {
    }

    /// A constructor that creates a vector with entries set by the
    /// values in the input initialiser_list
    /// Example:
    ///           Vector<int> arr{0, 20, 100, 150);
    ///           Vector<int> arr = {0, 20, 100, 150);
    Vector(std::initializer_list<_Tp> init) : std::vector<_Tp>(init) {}

    /// Copy constructor
    Vector(const Vector<_Tp>& __x) : std::vector<_Tp>(__x) {}

    // No explicit destructor is required because the base class destructor
    // handles all memory issues ~Vector() {}

    /// Iterate over all values and set to the desired value
    void initialise(const _Tp& __value)
    {
      for (typename std::vector<_Tp>::iterator it = std::vector<_Tp>::begin();
           it != std::vector<_Tp>::end();
           it++)
      {
        *it = __value;
      }
    }

#ifdef RANGE_CHECKING
    /// Overload the bracket access operator to include array-range checking
    /// if the RANGE_CHECKING flag is set
    reference operator[](size_type __n)
    {
      return error_checked_access(__n);
    }

    /// Overloaded, range-checking, bracket access operator (const version)
    const_reference operator[](size_type __n) const
    {
      return error_checked_access(__n);
    }
#endif
  };

  //==================================================================
  /// A Vector of bools cannot be created because the is no
  /// compiler-independent  implementation of the bit manipulators.
  /// Making all the constructors private should lead to compile-time
  /// errors.
  //=================================================================
  template<>
  class Vector<bool> : private std::vector<bool>
  {
  public:
    /// Typedef to make the constructors look a bit cleaner
    typedef bool value_type;

    /// Typedef to make the constructors look a bit cleaner
    typedef value_type& reference;

    /// Typedef to make the constructors look a bit cleaner
    typedef const value_type& const_reference;

    /// Typedef to make the constructors look a bit cleaner
    typedef size_t size_type;


    /// Dummy constructor to avoid compiler from warning about
    /// only-private constructors
    Vector(const double& dont_call_this_constructor)
    {
      // Throw an Oomph-lib error
      throw OomphLibError("Please use vector<bool> instead of Vector<bool>",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

  private:
    // Standard Constuctors (some have been omitted from the stl classes)

    /// Construct an empty vector
    Vector() : std::vector<bool>() {}

    /// A constructor that creates a vector of size __n.
    /// Note the use of explicit for "strong" type checking
    explicit Vector(size_type __n) : std::vector<bool>(__n) {}

    /// A constructor that creates a vector of size __n and
    /// initialises every entry to __value
    Vector(size_type __n, const bool& __value) : std::vector<bool>(__n, __value)
    {
    }

    /// A constructor that creates a vector with entries set by the
    /// values in the input initialiser_list.
    /// Example:
    ///           Vector<int> arr{true, false, false, true);
    ///           Vector<int> arr = {true, false, false, true);
    Vector(std::initializer_list<bool> init) : std::vector<bool>(init) {}

    /// Copy constructor
    Vector(const Vector<bool>& __x) : std::vector<bool>(__x) {}

    /// Iterate over all values and set to the desired value
    void initialise(const bool& __value)
    {
      for (std::vector<bool>::iterator it = std::vector<bool>::begin();
           it != std::vector<bool>::end();
           it++)
      {
        *it = __value;
      }
    }
  };


  //=================================================================
  /// Namespace for helper functions for Vector<double>
  //=================================================================
  namespace VectorHelpers
  {
    /// Check the lengths if two Vectors are the same length
    inline void check_lengths_match(const Vector<double>& a,
                                    const Vector<double>& b)
    {
#ifdef PARANOID
      if (a.size() != b.size())
      {
        std::ostringstream err;
        err << "Vectors must be the same length."
            << "len(a) = " << a.size() << ", "
            << "len(b) = " << b.size() << ".";

        throw OomphLibError(
          err.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }


    /// Probably not always best/fastest because not optimised for
    /// dimension but useful...
    inline double dot(const Vector<double>& a, const Vector<double>& b)
    {
      check_lengths_match(a, b);
      double temp = 0;
      for (unsigned i = 0, ni = a.size(); i < ni; i++)
      {
        temp += a[i] * b[i];
      }
      return temp;
    }

    /// Get the magnitude of a vector.
    inline double magnitude(const Vector<double>& a)
    {
      return (std::sqrt(dot(a, a)));
    }

    /// Get the angle between two vector.
    inline double angle(const Vector<double>& a, const Vector<double>& b)
    {
      // Notice that we use one square root operation by avoiding the
      // call to magnitude(...)
      return std::acos(dot(a, b) / std::sqrt(dot(a, a) * dot(b, b)));
    }


    /// Cross product using "proper" output (move semantics means this is
    /// ok nowadays).
    inline void cross(const Vector<double>& A,
                      const Vector<double>& B,
                      Vector<double>& C)
    {
#ifdef PARANOID
      if ((A.size() != 3) || (B.size() != 3) || (C.size() != 3))
      {
        std::ostringstream err;
        err << "Cross product only defined for vectors of length 3.\n"
            << "len(a) = " << A.size() << ", "
            << "len(b) = " << B.size() << ", "
            << "len(c) = " << C.size() << ".";

        throw OomphLibError(
          err.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      C[0] = A[1] * B[2] - A[2] * B[1];
      C[1] = A[2] * B[0] - A[0] * B[2];
      C[2] = A[0] * B[1] - A[1] * B[0];
    }

    /// Cross product using "proper" output (move semantics means this is
    /// ok This calls the other cross(...) function.
    inline Vector<double> cross(const Vector<double>& A,
                                const Vector<double>& B)
    {
      Vector<double> output(3, 0.0);
      cross(A, B, output);

      return output;
    }

  } // namespace VectorHelpers


} // namespace oomph


#endif
