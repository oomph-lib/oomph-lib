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
#ifndef OOMPH_SUM_OF_MATRICES_H
#define OOMPH_SUM_OF_MATRICES_H


#include "mesh.h"
#include "matrices.h"
#include "Vector.h"
#include "oomph_utilities.h"

#include <map>


namespace oomph
{
  using namespace StringConversion;

  // =================================================================
  /// Class to store bi-directional lookup between added matrix row/col
  /// numbers to main matrix (SumOfMatrix) row/col numbers.
  // =================================================================
  class AddedMainNumberingLookup
  {
  public:
    /// Default constructor
    AddedMainNumberingLookup() {}

    /// Real constructor: construct lookup from node numbers in mesh and
    /// global equation numbers. Useful for the case when the main matrix is
    /// a Jacobian and the added matrix is a contribution only on a certain
    /// mesh.
    AddedMainNumberingLookup(const Mesh* mesh_pt, const unsigned& dof_index)
    {
      this->build(mesh_pt, dof_index);
    }

    /// Construct lookup schemes from int array (HLib's format for this
    /// data).
    AddedMainNumberingLookup(const int* lookup_array, const unsigned& length)
    {
      this->build(lookup_array, length);
    }

    /// Destructor
    ~AddedMainNumberingLookup() {}

    /// Given a main matrix row/col number get the equivalent row/col
    /// in the added matrix. Throw an error if not found.
    unsigned main_to_added(const int& main) const
    {
      int result = unsafe_main_to_added(main);
#ifdef PARANOID
      // If it's -1 then we failed to find it:
      if (result == -1)
      {
        std::string err = "Main matrix row/col number " + to_string(main) +
                          " not found in lookup.";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }

      if (result < 0)
      {
        std::string err = "Something crazy went wrong here.";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
#endif

      return unsigned(result);
    }

    /// Given a main matrix row/col number get the equivalent row/col
    /// in the added matrix. Return -1 if not found.
    int unsafe_main_to_added(const int& main) const
    {
      // Find the entry
      std::map<unsigned, unsigned>::const_iterator it =
        Main_to_added_mapping.find(unsigned(main));

      // Check the entry existed, it not then return -1.
      if (it == main_to_added_mapping_pt()->end())
      {
        return -1;
      }
      else
      {
        return it->second;
      }
    }

    /// Given a row/col number in the added matrix return the equivalent
    /// row/col number in the main matrix.
    unsigned added_to_main(const unsigned& added) const
    {
      return Added_to_main_mapping[added];
    }

    /// Construct the lookup schemes given a mesh and the degree of freedom
    /// to lookup main equation numbers for.
    void build(const Mesh* mesh_pt, const unsigned& dof_index)
    {
      construct_added_to_main_mapping(mesh_pt, dof_index);
      construct_reverse_mapping();
    }

    /// Construct lookup schemes from int array (HLib's format for this
    /// data).
    void build(const int* lookup_array, const unsigned& length)
    {
#ifdef PARANOID
      // Check for negative entries just in case (since it's an integer
      // array).
      for (unsigned j = 0; j < length; j++)
      {
        if (lookup_array[j] < 0)
        {
          std::string err = "negative entry in lookup array!";
          throw OomphLibError(
            err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
        }
      }
#endif

      // Copy array into mapping and generate the inverse lookup
      Added_to_main_mapping.assign(lookup_array, lookup_array + length);
      construct_reverse_mapping();
    }

    /// Construct lookup using node vector
    void build(const Vector<const Node*>& bem_lookup, const unsigned& dof_index)
    {
      const unsigned ni = bem_lookup.size();
      Added_to_main_mapping.assign(ni, -1);

      for (unsigned i = 0; i < ni; i++)
      {
        Added_to_main_mapping[i] = bem_lookup[i]->eqn_number(dof_index);
      }

      construct_reverse_mapping();
    }

    /// Construct an identity map (mostly for testing).
    void build_identity_map(const unsigned& n)
    {
      Added_to_main_mapping.assign(n, 0);
      for (unsigned j = 0; j < n; j++)
      {
        Added_to_main_mapping[j] = j;
      }
      construct_reverse_mapping();
    }


    // Access functions
    // ============================================================

    /// Const access function for mapping.
    const Vector<unsigned>* added_to_main_mapping_pt() const
    {
      return &Added_to_main_mapping;
    }

    /// Const access function for mapping.
    const std::map<unsigned, unsigned>* main_to_added_mapping_pt() const
    {
      return &Main_to_added_mapping;
    }

  private:
    /// Set up the lookup from added matrix row/col to main matrix.
    void construct_added_to_main_mapping(const Mesh* mesh_pt,
                                         const unsigned& dof_index)
    {
      // Basically just copy from the node data.
      Added_to_main_mapping.resize(mesh_pt->nnode());
      for (unsigned nd = 0, nnode = mesh_pt->nnode(); nd < nnode; nd++)
      {
        Added_to_main_mapping[nd] = mesh_pt->node_pt(nd)->eqn_number(dof_index);
      }
    }

    /// Set up the main to added mapping using the added to main mapping.
    void construct_reverse_mapping()
    {
#ifdef PARANOID
      if (Added_to_main_mapping.size() == 0)
      {
        std::ostringstream error_msg;
        error_msg << "Must set up Added_to_main_mapping first.";
        throw OomphLibError(
          error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Clear old data
      Main_to_added_mapping.clear();

      // Copy from Added_to_main_mapping with order reversed.
      for (unsigned j = 0; j < Added_to_main_mapping.size(); j++)
      {
        Main_to_added_mapping.insert(
          std::make_pair(Added_to_main_mapping[j], j));
      }
    }

    /// Mapping from added matrix row/col numbers to main matrix row/col
    /// numbers.
    Vector<unsigned> Added_to_main_mapping;

    /// Mapping from main matrix row/col numbers to added matrix row/col
    /// numbers. Note that we cannot use a vector here because the main
    /// matrix rows/cols mapped onto are probably not contiguous. Access
    /// times are O(log N) so if you need to iterate over all elements then
    /// use the pointer access functions and use stl iterators properly.
    std::map<unsigned, unsigned> Main_to_added_mapping;

    /// Inaccessible copy constructor
    AddedMainNumberingLookup(const AddedMainNumberingLookup& dummy) = delete;

    /// Inaccessible assignment operator
    void operator=(const AddedMainNumberingLookup& dummy) = delete;
  };


  //======================================================================
  /// Class for a matrix of the form M = S + G + H + ... where S is the main
  /// matrix and G,H etc. are matrices of size S or smaller.  This may be useful
  /// if, for example, G,H etc. are subblocks of M that must be stored in a
  /// different format to S.

  /// Maps mut be provided which gives a map from the rows/cols of the main
  /// matrix to the rows/cols of each of the added matrices.
  //======================================================================
  class SumOfMatrices : public DoubleMatrixBase,
                        public Matrix<double, SumOfMatrices>
  {
  private:
    /// Pointer to the matrix which we are adding the others to
    DoubleMatrixBase* Main_matrix_pt;

    /// List of pointers to the matrices that are added to the main matrix
    Vector<DoubleMatrixBase*> Added_matrix_pt;

    /// List of maps between row numbers of the main matrix and the
    /// added matrices.
    Vector<const AddedMainNumberingLookup*> Row_map_pt;

    /// List of maps between col numbers of the main matrix and the
    /// added matrices.
    Vector<const AddedMainNumberingLookup*> Col_map_pt;

    /// Should we delete the sub matrices when destructor is called?
    Vector<unsigned> Should_delete_added_matrix;

    /// Should we delete the main matrix when destructor is called?
    /// Default is no.
    bool Should_delete_main_matrix;

  public:
    /// Default constructor
    SumOfMatrices()
      : Main_matrix_pt(0),
        Added_matrix_pt(0),
        Row_map_pt(0),
        Col_map_pt(0),
        Should_delete_added_matrix(0),
        Should_delete_main_matrix(0)
    {
    }

    /// Constructor taking a pointer to the main matrix as input.
    SumOfMatrices(DoubleMatrixBase* main_matrix_pt)
      : Main_matrix_pt(main_matrix_pt),
        Added_matrix_pt(0),
        Row_map_pt(0),
        Col_map_pt(0),
        Should_delete_added_matrix(0),
        Should_delete_main_matrix(0)
    {
    }

    /// Broken copy constructor
    SumOfMatrices(const SumOfMatrices& matrix) = delete;

    /// Broken assignment operator
    void operator=(const SumOfMatrices&) = delete;

    /// Destructor: delete matrices as instructed by
    /// Should_delete_added_matrix vector and Should_delete_main_matrix.
    ~SumOfMatrices()
    {
      for (unsigned i_matrix = 0; i_matrix < Added_matrix_pt.size(); i_matrix++)
      {
        if (Should_delete_added_matrix[i_matrix] == 1)
        {
          delete Added_matrix_pt[i_matrix];
        }
      }

      if (Should_delete_main_matrix)
      {
        delete Main_matrix_pt;
      }
    }

    /// Access to the main matrix
    const DoubleMatrixBase* main_matrix_pt() const
    {
      return Main_matrix_pt;
    }
    DoubleMatrixBase*& main_matrix_pt()
    {
      return Main_matrix_pt;
    }

    /// Set the main matrix to be deleted by the destructor of the
    /// SumOfMatrices (default is to not delete it).
    void set_delete_main_matrix()
    {
      Should_delete_main_matrix = true;
    }


    /// Output the "bottom right" entry regardless of it being
    /// zero or not (this allows automatic detection of matrix size in
    /// e.g. matlab, python).
    void output_bottom_right_zero_helper(std::ostream& outfile) const
    {
      int last_row = this->nrow() - 1;
      int last_col = this->ncol() - 1;

      double last_value = operator()(last_row, last_col);

      if (last_value == 0.0)
      {
        outfile << last_row << " " << last_col << " " << 0.0 << std::endl;
      }
    }


    /// Output the matrix in sparse format. Note that this is going to be
    /// slow because we have to check every entry of every matrix for non-zeros.
    void sparse_indexed_output_helper(std::ostream& outfile) const
    {
      for (unsigned long i = 0; i < nrow(); i++)
      {
        for (unsigned long j = 0; j < ncol(); j++)
        {
          double entry = operator()(i, j);
          // Output if non-zero entry
          if (entry != 0.0)
          {
            outfile << i << " " << j << " " << entry << std::endl;
          }
        }
      }
    }


    /// Get a list of row/col indices and total entry for non-zeros in
    /// the matrix. e.g. for use as input to other matrix classes. Warning this
    /// is SLOW! for sparse matrices.
    void get_as_indices(Vector<int>& row,
                        Vector<int>& col,
                        Vector<double>& values)
    {
      row.clear();
      col.clear();
      values.clear();

      for (int i = 0; i < int(nrow()); i++)
      {
        for (int j = 0; j < int(ncol()); j++)
        {
          double entry = operator()(i, j);
          // Output if non-zero entry
          if (entry != 0.0)
          {
            row.push_back(i);
            col.push_back(j);
            values.push_back(entry);
          }
        }
      }
    }

    /// Add a new matrix to the sum by giving a matrix pointer and a
    /// mapping from the main matrix numbering to the added matrix's numbering.
    void add_matrix(DoubleMatrixBase* added_matrix_pt_in,
                    const AddedMainNumberingLookup* main_to_added_rows_pt,
                    const AddedMainNumberingLookup* main_to_added_cols_pt,
                    bool should_delete_matrix = false)
    {
#ifdef PARANOID
      // Check that row mapping has correct size
      if (main_to_added_rows_pt->main_to_added_mapping_pt()->size() >
          added_matrix_pt_in->nrow())
      {
        throw OomphLibError(
          "Row mapping size should be less than or equal to nrow (less than if "
          "it is a sparse matrix and there are some empty rows).",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Check that col mapping has correct size
      if (main_to_added_cols_pt->main_to_added_mapping_pt()->size() >
          added_matrix_pt_in->ncol())
      {
        throw OomphLibError(
          "Col mapping size should be less than or equal to ncol (less than if "
          "it is a sparse matrix and there are some empty cols).",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
#ifdef RANGE_CHECKING
      // Check that all entries in the row mapping are "in range" for the
      // main matrix.
      const Vector<unsigned>* rowmap_pt =
        main_to_added_rows_pt->added_to_main_mapping_pt();
      unsigned max_row =
        *std::max_element(rowmap_pt->begin(), rowmap_pt->end());
      if (max_row > main_matrix_pt()->nrow())
      {
        std::string err =
          "Trying to add a matrix with a mapping which specifices";
        err += " a max row of " + to_string(max_row) + " but the main matrix ";
        err += "only has " + to_string(main_matrix_pt()->nrow()) + " rows!";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }

      // Check that all entries in the row mapping are "in range" for the
      // main matrix.
      const Vector<unsigned>* colmap_pt =
        main_to_added_cols_pt->added_to_main_mapping_pt();
      unsigned max_col =
        *std::max_element(colmap_pt->begin(), colmap_pt->end());
      if (max_col > main_matrix_pt()->ncol())
      {
        std::string err =
          "Trying to add a matrix with a mapping which specifices";
        err += " a max col of " + to_string(max_col) + " but the main matrix ";
        err += "only has " + to_string(main_matrix_pt()->ncol()) + " cols!";
        throw OomphLibError(
          err, OOMPH_EXCEPTION_LOCATION, OOMPH_CURRENT_FUNCTION);
      }
#endif

      Added_matrix_pt.push_back(added_matrix_pt_in);
      Row_map_pt.push_back(main_to_added_rows_pt);
      Col_map_pt.push_back(main_to_added_cols_pt);
      Should_delete_added_matrix.push_back(unsigned(should_delete_matrix));
    }

    /// Access function for ith added matrix (main matrix not included in
    /// numbering).
    inline DoubleMatrixBase* added_matrix_pt(const unsigned& i) const
    {
      return Added_matrix_pt[i];
    }

    /// Access to the maps
    const AddedMainNumberingLookup* row_map_pt(const unsigned& i) const
    {
      return Row_map_pt[i];
    }
    const AddedMainNumberingLookup* col_map_pt(const unsigned& i) const
    {
      return Col_map_pt[i];
    }

    /// Return the number of rows of the main matrix.
    inline unsigned long nrow() const
    {
#ifdef PARANOID
      if (Main_matrix_pt == 0)
      {
        throw OomphLibError("Main_matrix_pt not set",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Main_matrix_pt->nrow();
    }

    /// Return the number of columns of the main matrix.
    inline unsigned long ncol() const
    {
#ifdef PARANOID
      if (Main_matrix_pt == 0)
      {
        throw OomphLibError("Main_matrix_pt not set",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Main_matrix_pt->ncol();
    }

    /// Return the number of added matrices in the sum
    inline unsigned n_added_matrix() const
    {
      return Added_matrix_pt.size();
    }

    /// Multiply: just call multiply on each of the matrices and add up
    /// the results (with appropriate bookeeping of the relative positions).
    void multiply(const DoubleVector& x, DoubleVector& soln) const;

    /// Broken operator() because it does not make sense to return
    /// anything by reference.
    double& entry(const unsigned long& i, const unsigned long& j) const
    {
      throw OomphLibError(
        "Broken write to entry: it does not make sense to write to a sum, you "
        "must write to one of the component matrices.",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Access function to get the total value of entries in position
    /// (i,j). Warning: this way of getting entries is far too slow to use
    /// inside of loops.
    double operator()(const unsigned long& i, const unsigned long& j) const
    {
      // Get contributions from all matrices
      double sum = main_matrix_pt()->operator()(i, j);
      for (unsigned i_matrix = 0; i_matrix < n_added_matrix(); i_matrix++)
      {
        int li = Row_map_pt[i_matrix]->unsafe_main_to_added(i);
        int lj = Col_map_pt[i_matrix]->unsafe_main_to_added(j);

        // If the added matrix contributes to this entry then add it
        if ((li != -1) && (lj != -1))
        {
          sum += added_matrix_pt(i_matrix)->operator()(li, lj);
        }
      }

      return sum;
    }

    /// Dummy overload of a pure virtual function. I'm not sure how best
    /// to implement this and I don't think I need it.
    virtual void multiply_transpose(const DoubleVector& x,
                                    DoubleVector& soln) const
    {
      std::ostringstream error_msg;
      error_msg << "Function not yet implemented.";
      throw OomphLibError(
        error_msg.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      // Possible implementations (not really thought through):
      // - just call multiply transpose on submatrices?
      // - do something funny with switching row and column maps?
    }
  };


} // namespace oomph
#endif
