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
#ifndef SPARSE_VECTOR_HEADER
#define SPARSE_VECTOR_HEADER


namespace oomph
{
  //========================================================================
  /// Sparse vector -- wrapper for map (but avoids filling map during
  /// reads of zero/empty entries). Square bracket operator is read-only.
  //========================================================================
  template<class T>
  class SparseVector
  {
  public:
    /// Constructor
    SparseVector()
    {
      // Empty instance
      Empty_pt = new T(0);
    }

    /// Destructor
    ~SparseVector()
    {
      // Delete empty instance
      delete Empty_pt;
    }

    // Initialise the bin (only called by the create_bins_of_objects()
    // method)
    void initialise(const unsigned& size)
    {
      // Create a "large" bool vector indicating all entries are empty
      Has_entry.resize(size, false);
    }

    /// Wipe storage
    void clear()
    {
      Data.clear();

      // Get current size and reset all entries
      const unsigned size = Has_entry.size();
      Has_entry.resize(size, false);
    }

    /// Square bracket access (const version)
    const T& operator[](const unsigned& i) const
    {
      typedef typename std::map<unsigned, T>::const_iterator IT;
      IT it = Data.find(i);
      if (it == Data.end())
      {
        return *Empty_pt;
      }
      return (*it).second;
    }

    /// Set value of i-th entry
    void set_value(const unsigned& i, const T& value)
    {
      Data[i] = value;
      // Mark as having entry
      Has_entry[i] = true;
    }

    /// Number of nonzero entries stored in here (specifying the
    /// size of the vector (as a mathematical object)
    /// doesn't really make sense -- the thing could be infinitely
    /// big and we wouldn't know or care)
    unsigned nnz() const
    {
      return Data.size();
    }

    /// Read-only access to underlying map
    const std::map<unsigned, T>* map_pt() const
    {
      return &Data;
    }

    /// Read/write access to underlying map -- dangerous!
    std::map<unsigned, T>* map_pt()
    {
      return &Data;
    }

    /// Check if the bin has entries
    bool has_entry(const unsigned& nbin)
    {
      return Has_entry[nbin];
    }

    /// Return vector containing all values
    /// (without reference to their specific indices)
    /// for fast direct access
    void get_all_values(Vector<T>& all_values) const
    {
      all_values.clear();
      all_values.resize(nnz());
      typedef typename std::map<unsigned, T>::const_iterator IT;
      unsigned count = 0;
      for (IT it = Data.begin(); it != Data.end(); it++)
      {
        all_values[count++] = (*it).second;
      }
    }

  private:
    /// Internal storage in map.
    std::map<unsigned, T> Data;

    /// Empty instance
    const T* Empty_pt;

    /// Keep track of the filled and empty bins
    std::vector<bool> Has_entry;
  };


#endif
