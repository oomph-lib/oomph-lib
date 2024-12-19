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
#ifndef OOMPH_MAP_MATRIX_HEADER
#define OOMPH_MAP_MATRIX_HEADER


// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

#include <map>

// oomph-lib headers
#include "Vector.h"
#include "oomph_utilities.h"

namespace oomph
{
  //================================================================
  /// MapMatrixMixed is a generalised, STL-map-based, sparse(ish) matrix
  /// class with mixed indices
  ///
  /// The matrix is indexed by indices of type KEY_TYPE_ROW and KEY_TYPE_COL
  /// and has entries of type VALUE_TYPE.
  ///
  /// Careful: If a zero entry is referenced then it is created in
  /// memory. Therefore this isn't really a practical sparse matrix scheme.
  /// Do not loop over `all' possible indices as even looking
  /// at them will inflate the matrix until it occupies as much
  /// space as a full one -- use (modification of) output routine
  /// to retrieve all nonzero entries.
  ///
  /// However, this is not a serious restriction, as the main purpose
  /// of this class is to allow non-integer indices.
  ///
  /// Example of usage:
  /// \code
  ///
  ///
  ///  // Assume we have a Vector of pointers to objects:
  ///  Vector<Rubbish*> object_pt;
  ///
  ///  [...]
  ///
  ///  // Number of objects
  ///  int nentry=object_pt.size();
  ///
  ///  // Use the pointers to the objects and associated integers as indices
  ///  // in a MapMatrixMixed whose entries are of type int
  ///  MapMatrixMixed<Rubbish*,int,int> like_a_matrix;
  ///
  ///  for (int i=1;i<nentry;i++)
  ///   {
  ///    for (int j=1;j<nentry;j++)
  ///     {
  ///      int number=100*i+j;
  ///      like_a_matrix(object_pt[i],j)=number;
  ///     }
  ///   }
  ///
  ///  oomph_info << "Matrix has nnz() " << like_a_matrix.nnz() <<
  ///          " and size() " << like_a_matrix.size() << std::endl;
  ///
  ///  oomph_info << "\n\n\n Here are the nonzero entries: i, j, a(i,j)\n";
  ///
  ///  like_a_matrix.output(oomph_info);
  ///
  ///  // Can be used like a normal matrix:
  ///  like_a_matrix(object_pt[1],20)+=13;
  ///  like_a_matrix(object_pt[1],1)+=13;
  ///
  ///  oomph_info << "\n\n\n Here are the nonzero entries: i, j, a(i,j)\n";
  /// like_a_matrix.output(oomph_info);
  ///
  /// \endcode
  ///
  //================================================================
  template<class KEY_TYPE_ROW, class KEY_TYPE_COL, class VALUE_TYPE>
  class MapMatrixMixed
  {
  public:
    /// Default (empty) constructor
    MapMatrixMixed(){};

    /// Broken assignment operator
    void operator=(const MapMatrixMixed&) = delete;

    /// Typedef to keep the code more readable
    typedef std::map<KEY_TYPE_COL, VALUE_TYPE> InnerMapMixed;

    /// Typedef to keep the code more readable
    typedef typename InnerMapMixed::iterator InnerMixedIt;

    /// Typedef to keep the code more readable const version
    typedef typename InnerMapMixed::const_iterator ConstInnerMixedIt;

    /// Typedef to keep the code more readable
    typedef std::map<KEY_TYPE_ROW, std::map<KEY_TYPE_COL, VALUE_TYPE>*>
      OuterMapMixed;

    /// Typedef to keep the code more readable
    typedef typename OuterMapMixed::iterator OuterMixedIt;

    /// Typedef to keep the code more readable const version
    typedef typename OuterMapMixed::const_iterator ConstOuterMixedIt;


    /// Copy constructor
    MapMatrixMixed(
      const MapMatrixMixed<KEY_TYPE_ROW, KEY_TYPE_COL, VALUE_TYPE>& map_mat)
    {
      // Step through the row pointers
      for (ConstOuterMixedIt it = map_mat.Row_pt.begin();
           it != map_mat.Row_pt.end();
           it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over entries in the row
          for (ConstInnerMixedIt it2 = inner_map.begin();
               it2 != inner_map.end();
               it2++)
          {
            // If the entry is nonzero: Copy
            if (it2->second != 0)
            {
              // key1, key2, value
              (*this)(it->first, it2->first) = it2->second;
            }
          }
        }
      }
    }


    /// Copy a single column into its own map
    void copy_column(const KEY_TYPE_COL& j,
                     std::map<KEY_TYPE_ROW, VALUE_TYPE>& copied_map)
    {
      // Step through the row pointers
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);
          // If the desired column of the inner map is non-zero
          if (inner_map[j] != 0)
          {
            // Set the value of the copied map to be the desired column of the
            // inner map
            copied_map[it->first] = inner_map[j];
          }
        }
      }
    }


    /// Destructor
    virtual ~MapMatrixMixed()
    {
      // Step through the pointers to row maps
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // it->second is the stored object
          delete it->second;
        }
      }
    }

    /// Wipe all entries
    void clear()
    {
      // Step through the pointers to row maps
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // it->second is the stored object: a map which can be cleared!
          it->second->clear();
        }
      }
    }


    /// Return (reference to) entry.
    /// Careful: If the entry does not exist then it is created and
    /// set to zero
    VALUE_TYPE& operator()(const KEY_TYPE_ROW& i, const KEY_TYPE_COL& j)
    {
      return *entry_pt(i, j);
    }

    /// Get an element corresponding to the key (i,j)
    ///  Searches the container for an element with a key equivalent to
    /// (i,j) and returns the element if found, otherwise the default 0 value
    /// for the value type is returned. The container is not modified.
    VALUE_TYPE get(const KEY_TYPE_ROW& i, const KEY_TYPE_COL& j) const
    {
      if (Row_pt.count(i) > 0)
      {
        // Get the pointer to the row and check the j key
        InnerMapMixed* inner_map_mixed_pt = Row_pt.find(i)->second;
        if (inner_map_mixed_pt->count(j) > 0)
        {
          return inner_map_mixed_pt->find(j)->second;
        }
        else
        {
          return VALUE_TYPE(0);
        }
      }
      else
      {
        // The key does not exist, return VALUE_TYPE(0)
        return VALUE_TYPE(0);
      }
    }


    /// Dump all non-`zero' entries to file.
    /// Output is in the format
    ///   `i', `j', `entry[i][j]'
    void output(std::ostream& outfile)
    {
      // NOTE:
      //------
      // map.first = key
      // map.second = value

      // Step through the row pointers
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over entries in the row
          for (InnerMixedIt it2 = inner_map.begin(); it2 != inner_map.end();
               it2++)
          {
            // If the entry is nonzero: Doc
            if (it2->second != 0)
            {
              // Output key1, key2, value
              outfile << it->first << " " << it2->first << " " << it2->second
                      << std::endl;
            }
          }
        }
      }
    }

    /// Work out number of non-`zero' entries
    unsigned long nnz()
    {
      // Initialise counter for # of nonzero entries
      unsigned long count = 0;

      // Step through the row pointers
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over entries in the row
          for (InnerMixedIt it2 = inner_map.begin(); it2 != inner_map.end();
               it2++)
          {
            // If the entry is nonzero: Increment counter
            if (it2->second != 0)
            {
              count++;
            }
          }
        }
      }
      return count;
    }

    /// Work out number of non-`zero' entries, const version
    unsigned long nnz() const
    {
      // Initialise counter for # of nonzero entries
      unsigned long count = 0;

      // Step through the row pointers
      for (ConstOuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over entries in the row
          for (ConstInnerMixedIt it2 = inner_map.begin();
               it2 != inner_map.end();
               it2++)
          {
            // If the entry is nonzero: Increment counter
            if (it2->second != 0)
            {
              count++;
            }
          }
        }
      }
      return count;
    }


    /// Work out total number of entries
    unsigned long size()
    {
      // Initialise counter for # of nonzero entries
      unsigned long count = 0;

      // Step through the row pointers
      for (OuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over all (!) entries in the row (incl. zero ones!)
          for (InnerMixedIt it2 = inner_map.begin(); it2 != inner_map.end();
               it2++)
          {
            count++;
          }
        }
      }
      return count;
    }

    /// Work out total number of entries const version
    unsigned long size() const
    {
      // Initialise counter for # of nonzero entries
      unsigned long count = 0;

      // Step through the row pointers
      for (ConstOuterMixedIt it = Row_pt.begin(); it != Row_pt.end(); it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMapMixed inner_map = *(it->second);

          // Loop over all (!) entries in the row (incl. zero ones!)
          for (ConstInnerMixedIt it2 = inner_map.begin();
               it2 != inner_map.end();
               it2++)
          {
            count++;
          }
        }
      }
      return count;
    }


  protected:
    /// Return pointer to entry
    VALUE_TYPE* entry_pt(const KEY_TYPE_ROW& i, const KEY_TYPE_COL& j)
    {
      // There's not a single entry in this row: Entry must be zero.
      if (Row_pt[i] == 0)
      {
        // Create row and entry in row and set the value to zero
        Row_pt[i] = new std::map<KEY_TYPE_COL, VALUE_TYPE>;
        (*Row_pt[i])[j] = VALUE_TYPE(0);
        return &(*Row_pt[i])[j];
      }
      // Simply return pointer to existing entry
      else
      {
        return &(*Row_pt[i])[j];
      }
    }


    /// Here's the generalised matrix structure: A map of pointers to
    /// the maps that hold the entries in each row.
    std::map<KEY_TYPE_ROW, std::map<KEY_TYPE_COL, VALUE_TYPE>*> Row_pt;
  };


  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////
  /// ////////////////////////////////////////////////////////////////////


  //================================================================
  /// MapMatrix is a generalised, STL-map-based, sparse(-ish) matrix
  /// class.
  ///
  /// The matrix is indexed by indices of type KEY_TYPE
  /// and has entries of type VALUE_TYPE. It is a specialisation of the
  /// class MapMatrixMixed. Please implement future functions in that class.
  ///
  /// Careful: If a zero entry is referenced then it is created in
  /// memory. Therefore this isn't really a practical sparse matrix scheme.
  /// Do not loop over `all' possible indices as even looking
  /// at them will inflate the matrix until it occupies as much
  /// space as a full one -- use (modification of) output routine
  /// to retrieve all nonzero entries.
  ///
  /// However, this is not a serious restriction, as the main purpose
  /// of this class is to allow non-integer indices.
  ///
  /// Example of usage:
  /// \code
  ///
  ///
  ///  // Assume we have a Vector of pointers to objects:
  ///  Vector<Rubbish*> object_pt;
  ///
  ///  [...]
  ///
  ///  // Number of objects
  ///  int nentry=object_pt.size();
  ///
  ///  // Use the pointers to the objects as indices
  ///  // in a MapMatrix whose entries are of type int
  ///  MapMatrix<Rubbish*,int> like_a_matrix;
  ///
  ///  for (int i=1;i<nentry;i++)
  ///   {
  ///    for (int j=1;j<nentry;j++)
  ///     {
  ///      int number=100*i+j;
  ///      like_a_matrix(object_pt[i],object_pt[j])=number;
  ///     }
  ///   }
  ///
  ///  oomph_info << "Matrix has nnz() " << like_a_matrix.nnz() <<
  ///          " and size() " << like_a_matrix.size() << std::endl;
  ///
  ///  oomph_info << "\n\n\n Here are the nonzero entries: i, j, a(i,j)\n";
  ///
  ///  like_a_matrix.output(oomph_info);
  ///
  ///  // Can be used like a normal matrix:
  ///
  ///  //Add to existing entry
  ///  like_a_matrix(object_pt[1],object_pt[2])+=13;
  ///
  ///  // Add to non-existing entry
  ///  Rubbish* temp_pt=new Rubbish(20);
  ///  like_a_matrix(object_pt[1],temp_pt)+=13;
  ///
  ///  oomph_info << "\n\n\n Here are the nonzero entries: i, j, a(i,j)\n";
  /// like_a_matrix.output(oomph_info);
  ///
  /// \endcode
  ///
  //================================================================
  template<class KEY_TYPE, class VALUE_TYPE>
  class MapMatrix : public MapMatrixMixed<KEY_TYPE, KEY_TYPE, VALUE_TYPE>
  {
  public:
    /// Default (empty) constructor
    MapMatrix(){};

    /// Typedef to keep the code more readable
    typedef std::map<KEY_TYPE, VALUE_TYPE> InnerMap;

    /// Typedef to keep the code more readable
    typedef typename InnerMap::iterator InnerIt;

    /// Typedef to keep the code more readable
    typedef typename InnerMap::const_iterator ConstInnerIt;

    /// Typedef to keep the code more readable
    typedef std::map<KEY_TYPE, std::map<KEY_TYPE, VALUE_TYPE>*> OuterMap;

    /// Typedef to keep the code more readable
    typedef typename OuterMap::iterator OuterIt;

    /// Typedef to keep the code more readable
    typedef typename OuterMap::const_iterator ConstOuterIt;

    /// Copy constructor
    MapMatrix(const MapMatrix<KEY_TYPE, VALUE_TYPE>& map_mat)
    {
      // Step through the row pointers
      for (ConstOuterIt it = map_mat.Row_pt.begin(); it != map_mat.Row_pt.end();
           it++)
      {
        // Is the row pointer nonzero, i.e. are there any entries in this row?
        if (it->second != 0)
        {
          // Identify the map that holds the entries in this row:
          InnerMap inner_map = *(it->second);

          // Loop over entries in the row
          for (ConstInnerIt it2 = inner_map.begin(); it2 != inner_map.end();
               it2++)
          {
            // If the entry is nonzero: Copy
            if (it2->second != 0)
            {
              (*this)(it->first, it2->first) = it2->second;
            }
          }
        }
      }
    }

    /// Broken assignment operator
    void operator=(const MapMatrix&) = delete;
  };

} // namespace oomph

#endif
