// This file is part of the MercuryDPM project (https://www.mercurydpm.org).
// Copyright (c), The MercuryDPM Developers Team. All rights reserved.
// License: BSD 3-Clause License; see the LICENSE file in the root directory.

#include "interpolate_from_integral_points.h"

namespace oomph
{
  OomphCommunicator* InterpolateFromIntegralPointsBase::Communicator_pt =
    new OomphCommunicator();

  void InterpolateFromIntegralPointsBase::compute_ipt_to_node_mapping()
  {
    // Clear any stored mapping
    delete_all_ipt_to_node_mapping();

    // We can delete our inverse mapping data now because we created it and own
    // it
    Can_delete_inverse_mapping_stored = true;

    // Calculate the inverse mapping.

    // Quite an inefficient way of doing things but I think this is the best I
    // can do since none of the linear solvers directly compute the inverse of
    // the matrix. Essentially we perform N solves, each time pulling out a
    // singlee column of the inverse of the matrix. This is quite inefficient
    // but for most cases we will only do it once then assign the result to all
    // other elements so it isn't so bad.

    // Get the number of nodes and the number of integral points
    const unsigned n_node = this->nnode();
    const unsigned n_ipt = this->integral_pt()->nweight();
    if (n_node != n_ipt)
    {
      throw OomphLibError(
        "InterpolateFromIntegralPoints must use an integral scheme with an "
        "equal number of nodes and integral points",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    // Check if we need to allocate storage for the integral point to node
    // mapping
    if (Inverse_mapping_pt == 0)
    {
      Inverse_mapping_pt = new DenseMatrix<double>(n_node);
    }

    // Fill in the mapping from the nodes to the integral points.
    DenseDoubleMatrix node_to_ipt_mapping(n_node);
    for (unsigned ipt = 0; ipt < n_ipt; ipt++)
    {
      // Fill in the shape functions mapping from the nodes to the integral
      // points
      Shape psi(n_node);
      this->shape_at_knot(ipt, psi);

      // Add them to the matrix
      for (unsigned n = 0; n < n_node; n++)
      {
        node_to_ipt_mapping(ipt, n) = psi(n);
      }
    }

    // Create the linear algebra distribution - do not distribute it
    LinearAlgebraDistribution* lin_alg_dist_pt =
      new LinearAlgebraDistribution(Communicator_pt, n_node, false);

    // Create the linear solver and pass to the matrix
    DenseLU* solver_pt = new DenseLU;
    node_to_ipt_mapping.linear_solver_pt() = solver_pt;

    // Create the rhs and solution vectors
    DoubleVector rhs(lin_alg_dist_pt);
    DoubleVector result(lin_alg_dist_pt);

    for (unsigned n = 0; n < n_node; n++)
    {
      // Set the value in the current row to 1
      rhs[n] = 1.0;

      // Solve for the column of the matrix inverse
      node_to_ipt_mapping.solve(rhs, result);

      // Fill in the column of the inverse mapping
      for (unsigned ipt = 0; ipt < n_ipt; ipt++)
      {
        (*Inverse_mapping_pt)(n, ipt) = result[ipt];
      }

      // Set the value in the current row back to zero
      rhs[n] = 0.0;
    }

    delete solver_pt;
    delete lin_alg_dist_pt;
  }

  void InterpolateFromIntegralPointsBase::
    set_compute_ipt_to_node_mapping_from_element(
      InterpolateFromIntegralPointsBase* const& element_pt)
  {
    if (Inverse_mapping_pt != element_pt->inverse_mapping_pt())
    {
      // Clear any stored mapping
      delete_all_ipt_to_node_mapping();

      Can_delete_inverse_mapping_stored = false;

      Inverse_mapping_pt = element_pt->inverse_mapping_pt();
    }
  }

  void InterpolateFromIntegralPointsBase::delete_all_ipt_to_node_mapping()
  {
    if ((Can_delete_inverse_mapping_stored) && (Inverse_mapping_pt))
    {
      delete Inverse_mapping_pt;
    }
    Inverse_mapping_pt = 0;
  }

} // namespace oomph
