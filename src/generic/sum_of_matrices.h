#ifndef OOMPH_SUM_OF_MATRICES_H
#define OOMPH_SUM_OF_MATRICES_H


#include "mesh.h"
#include "matrices.h"
#include "Vector.h"

#include <map>

namespace oomph
{

 // =================================================================
 /// Class to store node number to global and global to node lookups (could
 /// use mesh for node number to global but worried it might change order).
 // =================================================================
 class NodeGlobalNumbersLookup
 {
 public:
  /// Default constructor
  NodeGlobalNumbersLookup() {}

  /// Real constructor
  NodeGlobalNumbersLookup(const Mesh* mesh_pt, const unsigned& dof_index)
   {
    this->build(mesh_pt, dof_index);
   }

  /// Destructor
  ~NodeGlobalNumbersLookup() {}

  /// \short Given a global equation number in the lookup get the
  /// associated node number.
  int global_to_node(const int& global) const
   {
#ifdef PARANOID
    if(global < 0)
     {
      std::ostringstream error_msg;
      error_msg << "Pinned equation numbers do not relate to nodes.";
      throw OomphLibError(error_msg.str(),
                          "NodeGlobalNumbersLookup::global_to_node",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    // Find the entry
    std::map<unsigned, unsigned>::const_iterator
     it = Global_to_node_mapping.find(unsigned(global));

    // Check the entry existed, it not then return -1.
    if(it == global_to_node_mapping_pt()->end())
     {
      return -1;
     }
    else
     {
      return it->second;
     }
   }

  /// Convert node number (in the mesh) to global equation number.
  unsigned node_to_global(const unsigned& node) const
   {return Node_to_global_mapping[node];}

  /// Construct the lookup schemes given a mesh and the degree of freedom
  /// to lookup.
  void build(const Mesh* mesh_pt, const unsigned& dof_index)
   {
    construct_node_to_global_mapping(mesh_pt, dof_index);
    construct_reverse_mapping();
   }

  /// Construct an identity map (mostly for testing).
  void build_identity_map(const unsigned& n)
   {
    Node_to_global_mapping.assign(n, 0);
    for(unsigned j=0; j<n; j++)
     {
      Node_to_global_mapping[j] = j;
     }
    construct_reverse_mapping();
   }


  // Access functions
  // ============================================================

  /// \short Const access function for Node_to_global_mapping.
  const Vector<long>* node_to_global_mapping_pt() const
   {return &Node_to_global_mapping;}

  /// \short Const access function for Global_to_node_mapping_pt.
  const std::map<unsigned, unsigned>* global_to_node_mapping_pt() const
   {return &Global_to_node_mapping;}

 private:

  /// Set up the lookup from node number to global equation number.
  void construct_node_to_global_mapping(const Mesh* mesh_pt,
                                        const unsigned& dof_index)
   {
    // Basically just copy from the node data.
    Node_to_global_mapping.resize(mesh_pt->nnode());
    for(unsigned nd=0, nnode=mesh_pt->nnode(); nd<nnode; nd++)
     {
      Node_to_global_mapping[nd] =
       mesh_pt->node_pt(nd)->eqn_number(dof_index);
     }
   }

  /// Set up the global to node mapping using the node to global mapping.
  void construct_reverse_mapping()
   {
#ifdef PARANOID
    if(Node_to_global_mapping.size() == 0)
     {
      std::ostringstream error_msg;
      error_msg << "Must set up Node_to_global_mapping first.";
      throw OomphLibError(error_msg.str(),
                          "NodeGlobalNumbersLookup::construct_reverse_mapping",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    // Clear old data
    Global_to_node_mapping.clear();

    // Copy from Node_to_global_mapping with order reversed.
    for(unsigned j=0; j<Node_to_global_mapping.size(); j++)
     {
      std::pair<unsigned, unsigned> a = std::make_pair(Node_to_global_mapping[j], j);
      Global_to_node_mapping.insert(a);
     }
   }

  /// Mapping from node numbers to global equation numbers.
  Vector<long> Node_to_global_mapping;

  /// Mapping from global equation numbers to node numbers. Note that we
  /// cannot use a vector here because the global equation numbers are
  /// probably not contiguous. Access times are O(log N) so if you need to
  /// iterate over all elements then use the pointer access functions and
  /// use stl iterators properly.
  std::map<unsigned, unsigned> Global_to_node_mapping;

  /// Inaccessible copy constructor
  NodeGlobalNumbersLookup(const NodeGlobalNumbersLookup& dummy)
   {BrokenCopy::broken_copy("NodeGlobalNumbersLookup");}

  /// Inaccessible assignment operator
  void operator=(const NodeGlobalNumbersLookup& dummy)
   {BrokenCopy::broken_assign("NodeGlobalNumbersLookup");}

 };


//======================================================================
/// Class for a matrix of the form M = S + G + H + ... where S is the main
/// matrix and G,H etc. are matrices of size S or smaller.  This may be useful
/// if, for example, G,H etc. are subblocks of M that must be stored in a
/// different format to S.

/// Maps mut be provided which gives a map from the rows/cols of the "main"
/// matrix to the rows/cols of each of the added matrices.
//======================================================================
 class SumOfMatrices : public DoubleMatrixBase,
                       public Matrix<double,SumOfMatrices>
 {

 private:


  typedef std::map<long unsigned, long unsigned> IndexMap;

  /// Pointer to the matrix which we are adding the others to
  DoubleMatrixBase* Main_matrix_pt;

  /// List of pointers to the matrices that are added to the main matrix
  Vector<DoubleMatrixBase*> Added_matrix_pt;

  /// \short List of maps from the row numbers of the main matrix to the added matrix
  /// row numbers.
  Vector<const NodeGlobalNumbersLookup* > Main_to_individual_rows_pt;

  /// \short List of maps from the col numbers of the main matrix to the added matrix
  /// col numbers.
  Vector<const NodeGlobalNumbersLookup* > Main_to_individual_cols_pt;

  /// Should we delete the sub matrices when destructor is called?
  Vector<unsigned> Should_delete_added_matrix;

  /// \short Should we delete the main matrix when destructor is called? Default is
  /// no.
  bool Should_delete_main_matrix;

 public:

  /// Default constructor
  SumOfMatrices()
   : Main_matrix_pt(0), Added_matrix_pt(0),
     Main_to_individual_rows_pt(0), Main_to_individual_cols_pt(0),
     Should_delete_added_matrix(0), Should_delete_main_matrix(0) {}

  /// Constructor taking a pointer to the main matrix as input.
  SumOfMatrices(DoubleMatrixBase* main_matrix_pt)
   : Main_matrix_pt(main_matrix_pt), Added_matrix_pt(0),
     Main_to_individual_rows_pt(0), Main_to_individual_cols_pt(0),
     Should_delete_added_matrix(0), Should_delete_main_matrix(0) {}

  /// Broken copy constructor
  SumOfMatrices(const SumOfMatrices& matrix)
   {BrokenCopy::broken_copy("SumOfMatrices");}

  /// Broken assignment operator
  void operator=(const SumOfMatrices&)
   {BrokenCopy::broken_assign("SumOfMatrices");}

  /// Destructor: delete matrices as instructed by Should_delete_added_matrix vector
  ~SumOfMatrices()
   {
    for(unsigned i_matrix=0; i_matrix<Added_matrix_pt.size(); i_matrix++)
     {
      if(Should_delete_added_matrix[i_matrix] == 1)
       delete Added_matrix_pt[i_matrix];
     }

    if(Should_delete_main_matrix)
     delete Main_matrix_pt;
   }

  /// Access to the main matrix
  const DoubleMatrixBase* main_matrix_pt() const {return Main_matrix_pt;}
  DoubleMatrixBase*& main_matrix_pt() {return Main_matrix_pt;}

  /// \short Set the main matrix to be deleted by the destructor of the SumOfMatrices
  /// (default is to not delete it).
  void set_delete_main_matrix()
   {Should_delete_main_matrix = true;}

  /// \short Output the matrix in sparse format. Note that this is going to be
  /// slow because we have to check every entry of every matrix for non-zeros.
  void sparse_indexed_output(std::ostream &outfile,
                             const bool& force_output_final_entry=false) const
   {
    for (unsigned long i=0; i<nrow(); i++)
     {
      for (unsigned long j=0; j<ncol(); j++)
       {
        double entry = operator()(i,j);
        // Output if non-zero entry
        if(entry != 0.0)
         {
          outfile << i << " " << j << " " << entry
                  << std::endl;
         }
       }
     }

    // If there is no output for the last entry and we requested it then
    // output the zero.
    if((force_output_final_entry) &&
       (operator()(this->nrow()-1, this->ncol()-1) == 0.0))
     {
      outfile << this->nrow()-1 << " " << this->ncol()-1 << " 0"
              << std::endl;
     }
   }

  /// \short Output the matrix in sparse format to a file. Note that this
  /// is going to be slow because we have to check every entry of every
  /// matrix for non-zeros. If the boolean flag is true then output the
  /// "bottom right" entry regardless of it being zero or
  /// not (helps reading into matlab etc.).
  void sparse_indexed_output(const std::string &outfile,
                             const bool& force_output_final_entry=false) const
   {
    // Open file
    std::ofstream some_file;
    some_file.open(outfile.c_str());
    sparse_indexed_output(some_file, force_output_final_entry);
    some_file.close();
   }


  /// \short Get a list of row/col indices and total entry for non-zeros in the
  /// matrix. e.g. for use as input to other matrix classes. Warning this is
  /// SLOW! for sparse matrices.
  void get_as_indices(Vector<int>& row, Vector<int>& col,
                      Vector<double>& values)
   {
    row.clear(); col.clear(); values.clear();

    for (int i=0; i<int(nrow()); i++)
     {
      for (int j=0; j<int(ncol()); j++)
       {
        double entry = operator()(i,j);
        // Output if non-zero entry
        if(entry != 0.0)
         {
          row.push_back(i);
          col.push_back(j);
          values.push_back(entry);
         }
       }
     }
   }

  /// \short Add a new matrix to the sum by giving a matrix pointer and a mapping from
  /// the main matrix numbering to the new matrix numbering.
  void add_matrix(DoubleMatrixBase* added_matrix_pt_in,
                  const NodeGlobalNumbersLookup* main_to_individual_rows_pt,
                  const NodeGlobalNumbersLookup* main_to_individual_cols_pt,
                  bool should_delete_matrix=false)
   {
#ifdef RANGE_CHECKING
    if (main_to_individual_rows_pt->global_to_node_mapping_pt()->size()
        > added_matrix_pt_in->nrow())
     {
      throw OomphLibError("Row mapping size should be less than or equal to nrow (less than if it is a sparse matrix and there are some empty rows).",
                          "SumOfMatrices::add_matrix",
                          OOMPH_EXCEPTION_LOCATION);
     }

    if (main_to_individual_cols_pt->global_to_node_mapping_pt()->size()
        > added_matrix_pt_in->ncol())
     {
      throw OomphLibError("Col mapping size should be less than or equal to ncol (less than if it is a sparse matrix and there are some empty cols).",
                          "SumOfMatrices::add_matrix",
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    Added_matrix_pt.push_back(added_matrix_pt_in);
    Main_to_individual_rows_pt.push_back(main_to_individual_rows_pt);
    Main_to_individual_cols_pt.push_back(main_to_individual_cols_pt);
    Should_delete_added_matrix.push_back(unsigned(should_delete_matrix));
   }

  /// Access functions
  inline DoubleMatrixBase* added_matrix_pt(const unsigned& i) const
   {return Added_matrix_pt[i];}

  /// Return the number of rows of the total matrix (equal to that of the first
  /// matrix).
  inline unsigned long nrow() const
   {
#ifdef PARANOID
    if(Main_matrix_pt==0)
     {
      OomphLibError("Main_matrix_pt not set","SumOfMatrices::nrow()",
                    OOMPH_EXCEPTION_LOCATION);
     }
#endif
    return Main_matrix_pt->nrow();
   }

  /// \short Return the number of columns of the total matrix (equal to that of the
  /// first matrix).
  inline unsigned long ncol() const
   {
#ifdef PARANOID
    if(Main_matrix_pt==0)
     {
      OomphLibError("Main_matrix_pt not set","SumOfMatrices::nrow()",
                    OOMPH_EXCEPTION_LOCATION);
     }
#endif
    return Main_matrix_pt->ncol();
   }

  /// Return the number of matrices in the sum
  inline unsigned n_added_matrix() const {return Added_matrix_pt.size();}

  /// \short Multiply: just call multiply on each of the matrices and add up the
  /// results (with appropriate bookeeping of the relative positions).
  void multiply(const DoubleVector &x, DoubleVector &soln) const;

  /// \short Broken operator() because it does not make sense to return
  /// anything by reference.
  double& entry(const unsigned long& i, const unsigned long& j) const
   {
    throw OomphLibError("Broken write to entry: it does not make sense to write to a sum, you must write to one of the component matrices.",
                        "non-constant SumOfMatrices::operator()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  /// Access function to get the total value of entries in position
  /// (i,j). Warning: this way of getting entries is far too slow to use
  /// inside of loops.
  double operator()(const unsigned long &i,
                    const unsigned long &j) const
   {
    // Get contributions from all matrices
    double sum = main_matrix_pt()->operator()(i,j);
    for(unsigned i_matrix=0; i_matrix<n_added_matrix(); i_matrix++)
     {
      int li = Main_to_individual_rows_pt[i_matrix]->global_to_node(i);
      int lj = Main_to_individual_cols_pt[i_matrix]->global_to_node(j);

      // If the global numbers are in the map then add the entry
      if(( li != -1) && (lj != -1))
       {
        sum += added_matrix_pt(i_matrix)->operator()(li,lj);
       }
     }

    return sum;
   }

  /// \short Dummy overload of a pure virtual function. I'm not sure how best to
  /// implement this and I don't think I need it.
  virtual void multiply_transpose(const DoubleVector &x,
                                  DoubleVector &soln)const
   {
    std::ostringstream error_msg;
    error_msg << "Function not yet implemented.";
    throw OomphLibError(error_msg.str(),
                        "SumOfMatrices::multiply_transpose",
                        OOMPH_EXCEPTION_LOCATION);

    // Possible implementations (not really thought through):
    // - just call multiply transpose on submatrices?
    // - do something funny with switching row and column maps?
   }

 };


}
#endif
