

#include "matrices.h"
#include "double_vector.h"
#include "sum_of_matrices.h"

using namespace oomph;

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

  // Multiply for the first matrix
  Main_matrix_pt->multiply(x,soln);

  // Now add contribution for other matrices
  for(unsigned i_matrix=0; i_matrix< Added_matrix_pt.size(); i_matrix++)
   {
    // Try to cast to a distributed object to get a distribution pointer
    DistributableLinearAlgebraObject* dist_obj_pt
     = dynamic_cast < DistributableLinearAlgebraObject* >
     (added_matrix_pt(i_matrix));

    // If possible copy the matrix distribution, otherwise it isn't
    // distributed so make a serial LinearAlgebraDistribution object.
    LinearAlgebraDistribution dist;
    LinearAlgebraDistribution* dist_pt = &dist;
    OomphCommunicator* serial_comm_pt = new OomphCommunicator; // Serial communcator (does nothing)
    if(dist_obj_pt == 0)
     {
      dist_pt->build(serial_comm_pt,added_matrix_pt(i_matrix)->nrow(),false);
     }
    else
     {
      dist_pt = dist_obj_pt->distribution_pt();
     }

    // Create temporary output DoubleVector
    DoubleVector temp_soln(dist_pt);

    // Create a const iterator for the map (faster than .find() or []
    // access, const means can't change the map via the iterator)
    std::map<unsigned, unsigned>::const_iterator it;
    
    // Pull out the appropriate values into a temp vector
    //??ds not parallel
    DoubleVector temp_x(dist_pt);
    for(it = Main_to_individual_cols_pt[i_matrix]->global_to_node_mapping_pt()->begin();
        it != Main_to_individual_cols_pt[i_matrix]->global_to_node_mapping_pt()->end();
        it++)
     {
      temp_x[it->second] = x[it->first];
     }
    
    // Perform the multiplication
    Added_matrix_pt[i_matrix]->multiply(temp_x,temp_soln);

    // Add result to solution vector
    //??ds not parallel
    for(it = Main_to_individual_rows_pt[i_matrix]->global_to_node_mapping_pt()->begin();
        it != Main_to_individual_rows_pt[i_matrix]->global_to_node_mapping_pt()->end();
        it++)
     {
      soln[it->first] += temp_soln[it->second];
     }
   }

 }

}
