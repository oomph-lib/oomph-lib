

#include "matrices.h"
#include "double_vector.h"
#include "sum_of_matrices.h"


namespace oomph
{

 // =================================================================
 /// \short Matrix-vector multiplication for a sumofmatrices class. Just
 /// delegate each multiplication to the appropriate class then add up the
 /// results.
 // =================================================================
 void SumOfMatrices::multiply(const DoubleVector &x, DoubleVector &soln) const
 {
  // We assume that appropriate checks and initialisation on x and soln are
  // carried out within the individual matrix multiplys.

  // Multiply for the main matrix
  Main_matrix_pt->multiply(x, soln);

  // Now add contribution for the added matrices
  for(unsigned i_matrix=0; i_matrix<Added_matrix_pt.size(); i_matrix++)
   {
    // If possible copy the matrix distribution, otherwise it isn't
    // distributed so make a serial LinearAlgebraDistribution object.
    LinearAlgebraDistribution col_dist, row_dist;
    OomphCommunicator serial_comm; // Serial communcator (does nothing)
    col_dist.build(&serial_comm, added_matrix_pt(i_matrix)->ncol(), false);
    row_dist.build(&serial_comm, added_matrix_pt(i_matrix)->nrow(), false);

    // Create temporary output DoubleVector
    DoubleVector temp_soln(row_dist);

    // Create a const iterator for the map (faster than .find() or []
    // access, const means can't change the map via the iterator).
    std::map<unsigned, unsigned>::const_iterator it;
    
    // Pull out the appropriate values into a temp vector
    //??ds not parallel
    DoubleVector temp_x(col_dist);
    for(it = Col_map_pt[i_matrix]->main_to_added_mapping_pt()->begin();
        it != Col_map_pt[i_matrix]->main_to_added_mapping_pt()->end();
        it++)
     {
        temp_x[it->second] = x[it->first];
     }
    
    // Perform the multiplication
    Added_matrix_pt[i_matrix]->multiply(temp_x, temp_soln);

    // Add result to solution vector
    //??ds not parallel
    for(it = Row_map_pt[i_matrix]->main_to_added_mapping_pt()->begin();
        it != Row_map_pt[i_matrix]->main_to_added_mapping_pt()->end();
        it++)
     {
        soln[it->first] += temp_soln[it->second];
     }
   }

 }

}
