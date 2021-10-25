// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Non-inline functions for constitutive laws and strain-energy functions

#include "constitutive_laws.h"
#include "../generic/elements.h"

namespace oomph
{
  //===============================================================
  /// This function is used to check whether a matrix is square
  //===============================================================
  bool ConstitutiveLaw::is_matrix_square(const DenseMatrix<double>& M)
  {
    // If the number rows and columns is not equal, the matrix is not square
    if (M.nrow() != M.ncol())
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  //========================================================================
  /// This function is used to check whether matrices are of equal dimension
  //========================================================================
  bool ConstitutiveLaw::are_matrices_of_equal_dimensions(
    const DenseMatrix<double>& M1, const DenseMatrix<double>& M2)
  {
    // If the numbers of rows and columns are not the same, then the
    // matrices are not of equal dimension
    if ((M1.nrow() != M2.nrow()) || (M1.ncol() != M2.ncol()))
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  //========================================================================
  /// This function is used to provide simple error (bounce) checks on the
  /// input to any calculate_second_piola_kirchhoff_stress
  //=======================================================================
  void ConstitutiveLaw::error_checking_in_input(const DenseMatrix<double>& g,
                                                const DenseMatrix<double>& G,
                                                DenseMatrix<double>& sigma)
  {
    // Test whether the undeformed metric tensor is square
    if (!is_matrix_square(g))
    {
      throw OomphLibError("Undeformed metric tensor not square",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // If the deformed metric tensor does not have the same dimension as
    // the undeformed tensor, complain
    if (!are_matrices_of_equal_dimensions(g, G))
    {
      std::string error_message = "Deformed metric tensor does  \n";
      error_message +=
        "not have same dimensions as the undeformed metric tensor\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }

    // If the stress tensor does not have the same dimensions as the others
    // complain.
    if (!are_matrices_of_equal_dimensions(g, sigma))
    {
      std::string error_message =
        "Strain tensor passed to calculate_green_strain() does  \n";
      error_message +=
        "not have same dimensions as the undeformed metric tensor\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }


  //===========================================================================
  /// The function to calculate the contravariant tensor from a covariant one
  //===========================================================================
  double ConstitutiveLaw::calculate_contravariant(
    const DenseMatrix<double>& Gdown, DenseMatrix<double>& Gup)
  {
    // Initial error checking
#ifdef PARANOID
    // Test that the matrices are of the same dimension
    if (!are_matrices_of_equal_dimensions(Gdown, Gup))
    {
      throw OomphLibError("Matrices passed to calculate_contravariant() are "
                          "not of equal dimension",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the matrix
    unsigned dim = Gdown.ncol();

    // If it's not square, I don't know what to do (yet)
#ifdef PARANOID
    if (!is_matrix_square(Gdown))
    {
      std::string error_message =
        "Matrix passed to calculate_contravariant() is not square\n";
      error_message += "non-square matrix inversion not implemented yet!\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Define the determinant of the matrix
    double det = 0.0;

    // Now the inversion depends upon the dimension of the matrix
    switch (dim)
    {
        // Zero dimensions
      case 0:
        throw OomphLibError("Zero dimensional matrix",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;

        // One dimension
      case 1:
        // The determinant is just the value of the only entry
        det = Gdown(0, 0);
        // The inverse is just the inverse of the value
        Gup(0, 0) = 1.0 / Gdown(0, 0);
        break;


        // Two dimensions
      case 2:
        // Calculate the determinant
        det = Gdown(0, 0) * Gdown(1, 1) - Gdown(0, 1) * Gdown(1, 0);
        // Calculate entries of the contravariant tensor (inverse)
        Gup(0, 0) = Gdown(1, 1) / det;
        Gup(0, 1) = -Gdown(0, 1) / det;
        Gup(1, 0) = -Gdown(1, 0) / det;
        Gup(1, 1) = Gdown(0, 0) / det;
        break;

        /// Three dimensions
      case 3:
        // Calculate the determinant of the matrix
        det = Gdown(0, 0) * Gdown(1, 1) * Gdown(2, 2) +
              Gdown(0, 1) * Gdown(1, 2) * Gdown(2, 0) +
              Gdown(0, 2) * Gdown(1, 0) * Gdown(2, 1) -
              Gdown(0, 0) * Gdown(1, 2) * Gdown(2, 1) -
              Gdown(0, 1) * Gdown(1, 0) * Gdown(2, 2) -
              Gdown(0, 2) * Gdown(1, 1) * Gdown(2, 0);

        // Calculate entries of the inverse matrix
        Gup(0, 0) =
          (Gdown(1, 1) * Gdown(2, 2) - Gdown(1, 2) * Gdown(2, 1)) / det;
        Gup(0, 1) =
          -(Gdown(0, 1) * Gdown(2, 2) - Gdown(0, 2) * Gdown(2, 1)) / det;
        Gup(0, 2) =
          (Gdown(0, 1) * Gdown(1, 2) - Gdown(0, 2) * Gdown(1, 1)) / det;
        Gup(1, 0) =
          -(Gdown(1, 0) * Gdown(2, 2) - Gdown(1, 2) * Gdown(2, 0)) / det;
        Gup(1, 1) =
          (Gdown(0, 0) * Gdown(2, 2) - Gdown(0, 2) * Gdown(2, 0)) / det;
        Gup(1, 2) =
          -(Gdown(0, 0) * Gdown(1, 2) - Gdown(0, 2) * Gdown(1, 0)) / det;
        Gup(2, 0) =
          (Gdown(1, 0) * Gdown(2, 1) - Gdown(1, 1) * Gdown(2, 0)) / det;
        Gup(2, 1) =
          -(Gdown(0, 0) * Gdown(2, 1) - Gdown(0, 1) * Gdown(2, 0)) / det;
        Gup(2, 2) =
          (Gdown(0, 0) * Gdown(1, 1) - Gdown(0, 1) * Gdown(1, 0)) / det;
        break;

      default:
        throw OomphLibError("Dimension of matrix must be 0, 1, 2 or 3\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;
    }

    // Return the determinant of the matrix
    return (det);
  }


  //===========================================================================
  /// The function to calculate the derivatives of the contravariant tensor
  /// and determinant of covariant tensor with respect to the components of
  /// the covariant tensor
  //===========================================================================
  void ConstitutiveLaw::calculate_d_contravariant_dG(
    const DenseMatrix<double>& Gdown,
    RankFourTensor<double>& d_Gup_dG,
    DenseMatrix<double>& d_detG_dG)
  {
    // Find the dimension of the matrix
    const unsigned dim = Gdown.ncol();

    // If it's not square, I don't know what to do (yet)
#ifdef PARANOID
    if (!is_matrix_square(Gdown))
    {
      std::string error_message =
        "Matrix passed to calculate_contravariant() is not square\n";
      error_message += "non-square matrix inversion not implemented yet!\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Define the determinant of the matrix
    double det = 0.0;

    // Now the inversion depends upon the dimension of the matrix
    switch (dim)
    {
        // Zero dimensions
      case 0:
        throw OomphLibError("Zero dimensional matrix",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;

        // One dimension
      case 1:
        // There is only one entry, so derivatives are easy
        d_detG_dG(0, 0) = 1.0;
        d_Gup_dG(0, 0, 0, 0) = -1.0 / (Gdown(0, 0) * Gdown(0, 0));
        break;


        // Two dimensions
      case 2:
        // Calculate the determinant
        det = Gdown(0, 0) * Gdown(1, 1) - Gdown(0, 1) * Gdown(1, 0);

        // Calculate the derivatives of the determinant
        d_detG_dG(0, 0) = Gdown(1, 1);
        // Need to use symmetry here
        d_detG_dG(0, 1) = d_detG_dG(1, 0) = -2.0 * Gdown(0, 1);
        d_detG_dG(1, 1) = Gdown(0, 0);

        // Calculate the "upper triangular" derivatives of the contravariant
        // tensor
        {
          const double det2 = det * det;
          d_Gup_dG(0, 0, 0, 0) = -Gdown(1, 1) * d_detG_dG(0, 0) / det2;
          d_Gup_dG(0, 0, 0, 1) = -Gdown(1, 1) * d_detG_dG(0, 1) / det2;
          d_Gup_dG(0, 0, 1, 1) =
            1.0 / det - Gdown(1, 1) * d_detG_dG(1, 1) / det2;

          d_Gup_dG(0, 1, 0, 0) = Gdown(0, 1) * d_detG_dG(0, 0) / det2;
          d_Gup_dG(0, 1, 0, 1) =
            -1.0 / det + Gdown(0, 1) * d_detG_dG(0, 1) / det2;
          d_Gup_dG(0, 1, 1, 1) = Gdown(0, 1) * d_detG_dG(1, 1) / det2;

          d_Gup_dG(1, 1, 0, 0) =
            1.0 / det - Gdown(0, 0) * d_detG_dG(0, 0) / det2;
          d_Gup_dG(1, 1, 0, 1) = -Gdown(0, 0) * d_detG_dG(0, 1) / det2;
          d_Gup_dG(1, 1, 1, 1) = -Gdown(0, 0) * d_detG_dG(1, 1) / det2;
        }

        // Calculate entries of the contravariant tensor (inverse)
        // Gup(0,0) = Gdown(1,1)/det;
        // Gup(0,1) = -Gdown(0,1)/det;
        // Gup(1,0) = -Gdown(1,0)/det;
        // Gup(1,1) = Gdown(0,0)/det;
        break;

        /// Three dimensions
      case 3:
        // This is not yet implemented
        throw OomphLibError(
          "Analytic derivatives of metric tensors not yet implemented in 3D\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);

        // Calculate the determinant of the matrix
        det = Gdown(0, 0) * Gdown(1, 1) * Gdown(2, 2) +
              Gdown(0, 1) * Gdown(1, 2) * Gdown(2, 0) +
              Gdown(0, 2) * Gdown(1, 0) * Gdown(2, 1) -
              Gdown(0, 0) * Gdown(1, 2) * Gdown(2, 1) -
              Gdown(0, 1) * Gdown(1, 0) * Gdown(2, 2) -
              Gdown(0, 2) * Gdown(1, 1) * Gdown(2, 0);

        // Calculate entries of the inverse matrix
        // Gup(0,0) =  (Gdown(1,1)*Gdown(2,2) - Gdown(1,2)*Gdown(2,1))/det;
        // Gup(0,1) = -(Gdown(0,1)*Gdown(2,2) - Gdown(0,2)*Gdown(2,1))/det;
        // Gup(0,2) =  (Gdown(0,1)*Gdown(1,2) - Gdown(0,2)*Gdown(1,1))/det;
        // Gup(1,0) = -(Gdown(1,0)*Gdown(2,2) - Gdown(1,2)*Gdown(2,0))/det;
        // Gup(1,1) =  (Gdown(0,0)*Gdown(2,2) - Gdown(0,2)*Gdown(2,0))/det;
        // Gup(1,2) = -(Gdown(0,0)*Gdown(1,2) - Gdown(0,2)*Gdown(1,0))/det;
        // Gup(2,0) =  (Gdown(1,0)*Gdown(2,1) - Gdown(1,1)*Gdown(2,0))/det;
        // Gup(2,1) = -(Gdown(0,0)*Gdown(2,1) - Gdown(0,1)*Gdown(2,0))/det;
        // Gup(2,2) =  (Gdown(0,0)*Gdown(1,1) - Gdown(0,1)*Gdown(1,0))/det;
        break;

      default:
        throw OomphLibError("Dimension of matrix must be 0, 1, 2 or 3\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
        break;
    }
  }


  //=========================================================================
  /// Calculate the derivatives of the contravariant
  /// 2nd Piola Kirchhoff stress tensor with respect to the deformed metric
  /// tensor. Arguments are the
  /// covariant undeformed and deformed metric tensor and the
  /// matrix in which to return the derivatives of the stress tensor
  /// The default implementation uses finite differences, but can be
  /// overloaded for constitutive laws in which an analytic formulation
  /// is possible.
  //==========================================================================
  void ConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    const DenseMatrix<double>& sigma,
    RankFourTensor<double>& d_sigma_dG,
    const bool& symmetrize_tensor)
  {
    // Initial error checking
#ifdef PARANOID
    // Test that the matrices are of the same dimension
    if (!are_matrices_of_equal_dimensions(g, G))
    {
      throw OomphLibError("Matrices passed are not of equal dimension",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the matrix (assuming that it's square)
    const unsigned dim = G.ncol();

    // Find the dimension
    // FD step
    const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Advanced metric tensor
    DenseMatrix<double> G_pls(dim, dim);
    DenseMatrix<double> sigma_pls(dim, dim);

    // Copy across the original value
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        G_pls(i, j) = G(i, j);
      }
    }

    // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
    // NOTE: We exploit the symmetry of the stress and metric tensors
    //       by incrementing G(i,j) and G(j,i) simultaenously and
    //       only fill in the "upper" triangles without copying things
    //       across the lower triangle. This is taken into account
    //       in the solid mechanics codes.
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        G_pls(i, j) += eps_fd;
        G_pls(j, i) = G_pls(i, j);

        // Get advanced stress
        this->calculate_second_piola_kirchhoff_stress(g, G_pls, sigma_pls);

        for (unsigned ii = 0; ii < dim; ii++)
        {
          for (unsigned jj = ii; jj < dim; jj++)
          {
            d_sigma_dG(ii, jj, i, j) =
              (sigma_pls(ii, jj) - sigma(ii, jj)) / eps_fd;
          }
        }

        // Reset
        G_pls(i, j) = G(i, j);
        G_pls(j, i) = G(j, i);
      }
    }

    // If we are symmetrising the tensor, do so
    if (symmetrize_tensor)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < i; j++)
        {
          for (unsigned ii = 0; ii < dim; ii++)
          {
            for (unsigned jj = 0; jj < ii; jj++)
            {
              d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
            }
          }
        }
      }
    }
  }


  //=========================================================================
  /// Calculate the derivatives of the contravariant
  /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
  /// with respect to the deformed metric tensor.
  /// Also return the derivatives of the determinant of the
  /// deformed metric tensor with respect to the deformed metric tensor.
  /// This form is appropriate
  /// for truly-incompressible materials.
  /// The default implementation uses finite differences for the
  /// derivatives that depend on the constitutive law, but not
  /// for the derivatives of the determinant, which are generic.
  //========================================================================
  void ConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    const DenseMatrix<double>& sigma,
    const double& detG,
    const double& interpolated_solid_p,
    RankFourTensor<double>& d_sigma_dG,
    DenseMatrix<double>& d_detG_dG,
    const bool& symmetrize_tensor)
  {
    // Initial error checking
#ifdef PARANOID
    // Test that the matrices are of the same dimension
    if (!are_matrices_of_equal_dimensions(g, G))
    {
      throw OomphLibError("Matrices passed are not of equal dimension",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the matrix (assuming that it's square)
    const unsigned dim = G.ncol();

    // FD step
    const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Advanced metric tensor etc.
    DenseMatrix<double> G_pls(dim, dim);
    DenseMatrix<double> sigma_dev_pls(dim, dim);
    DenseMatrix<double> Gup_pls(dim, dim);
    double detG_pls;

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        G_pls(i, j) = G(i, j);
      }
    }

    // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
    // NOTE: We exploit the symmetry of the stress and metric tensors
    //       by incrementing G(i,j) and G(j,i) simultaenously and
    //       only fill in the "upper" triangles without copying things
    //       across the lower triangle. This is taken into account
    //       in the remaining code further below.
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        G_pls(i, j) += eps_fd;
        G_pls(j, i) = G_pls(i, j);

        // Get advanced stress
        this->calculate_second_piola_kirchhoff_stress(
          g, G_pls, sigma_dev_pls, Gup_pls, detG_pls);


        // Derivative of determinant of deformed metric tensor
        d_detG_dG(i, j) = (detG_pls - detG) / eps_fd;

        // Derivatives of deviatoric stress and "upper" deformed metric
        // tensor
        for (unsigned ii = 0; ii < dim; ii++)
        {
          for (unsigned jj = ii; jj < dim; jj++)
          {
            d_sigma_dG(ii, jj, i, j) =
              (sigma_dev_pls(ii, jj) - interpolated_solid_p * Gup_pls(ii, jj) -
               sigma(ii, jj)) /
              eps_fd;
          }
        }

        // Reset
        G_pls(i, j) = G(i, j);
        G_pls(j, i) = G(j, i);
      }
    }

    // If we are symmetrising the tensor, do so
    if (symmetrize_tensor)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < i; j++)
        {
          d_detG_dG(i, j) = d_detG_dG(j, i);

          for (unsigned ii = 0; ii < dim; ii++)
          {
            for (unsigned jj = 0; jj < ii; jj++)
            {
              d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
            }
          }
        }
      }
    }
  }


  //========================================================================
  /// Calculate the derivatives of the contravariant
  /// 2nd Piola Kirchoff stress tensor with respect to the deformed metric
  /// tensor. Also return the derivatives of the generalised dilatation,
  /// \f$ d, \f$ with respect to the deformed metric tensor.
  /// This form is appropriate
  /// for near-incompressible materials.
  /// The default implementation uses finite differences.
  //=======================================================================
  void ConstitutiveLaw::calculate_d_second_piola_kirchhoff_stress_dG(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    const DenseMatrix<double>& sigma,
    const double& gen_dil,
    const double& inv_kappa,
    const double& interpolated_solid_p,
    RankFourTensor<double>& d_sigma_dG,
    DenseMatrix<double>& d_gen_dil_dG,
    const bool& symmetrize_tensor)
  {
    // Initial error checking
#ifdef PARANOID
    // Test that the matrices are of the same dimension
    if (!are_matrices_of_equal_dimensions(g, G))
    {
      throw OomphLibError("Matrices passed are not of equal dimension",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Find the dimension of the matrix (assuming that it's square)
    const unsigned dim = G.ncol();

    // FD step
    const double eps_fd = GeneralisedElement::Default_fd_jacobian_step;

    // Advanced metric tensor etc
    DenseMatrix<double> G_pls(dim, dim);
    DenseMatrix<double> sigma_dev_pls(dim, dim);
    DenseMatrix<double> Gup_pls(dim, dim);
    double gen_dil_pls;
    double inv_kappa_pls;

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        G_pls(i, j) = G(i, j);
      }
    }

    // Do FD -- only w.r.t. to upper indices, exploiting symmetry.
    // NOTE: We exploit the symmetry of the stress and metric tensors
    //       by incrementing G(i,j) and G(j,i) simultaenously and
    //       only fill in the "upper" triangles without copying things
    //       across the lower triangle. This is taken into account
    //       in the remaining code further below.
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        G_pls(i, j) += eps_fd;
        G_pls(j, i) = G_pls(i, j);

        // Get advanced stress
        this->calculate_second_piola_kirchhoff_stress(
          g, G_pls, sigma_dev_pls, Gup_pls, gen_dil_pls, inv_kappa_pls);

        // Derivative of generalised dilatation
        d_gen_dil_dG(i, j) = (gen_dil_pls - gen_dil) / eps_fd;

        // Derivatives of deviatoric stress and "upper" deformed metric
        // tensor
        for (unsigned ii = 0; ii < dim; ii++)
        {
          for (unsigned jj = ii; jj < dim; jj++)
          {
            d_sigma_dG(ii, jj, i, j) =
              (sigma_dev_pls(ii, jj) - interpolated_solid_p * Gup_pls(ii, jj) -
               sigma(ii, jj)) /
              eps_fd;
          }
        }

        // Reset
        G_pls(i, j) = G(i, j);
        G_pls(j, i) = G(j, i);
      }
    }

    // If we are symmetrising the tensor, do so
    if (symmetrize_tensor)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < i; j++)
        {
          d_gen_dil_dG(i, j) = d_gen_dil_dG(j, i);

          for (unsigned ii = 0; ii < dim; ii++)
          {
            for (unsigned jj = 0; jj < ii; jj++)
            {
              d_sigma_dG(ii, jj, i, j) = d_sigma_dG(jj, ii, j, i);
            }
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// Calculate the contravariant 2nd Piola Kirchhoff
  /// stress tensor. Arguments are the
  /// covariant undeformed (stress-free) and deformed metric
  /// tensors, g and G, and the matrix in which to return the stress tensor.
  //=====================================================================
  void GeneralisedHookean::calculate_second_piola_kirchhoff_stress(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    DenseMatrix<double>& sigma)
  {
    // Error checking
#ifdef PARANOID
    error_checking_in_input(g, G, sigma);
#endif

    // Find the dimension of the problem
    unsigned dim = G.nrow();

    // Calculate the contravariant Deformed metric tensor
    DenseMatrix<double> Gup(dim);
    // We don't need the Jacobian so cast the function to void
    (void)calculate_contravariant(G, Gup);

    // Premultiply some constants
    double C1 = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));
    double C2 = 2.0 * (*Nu_pt) / (1.0 - 2.0 * (*Nu_pt));

    // Strain tensor
    DenseMatrix<double> strain(dim, dim);

    // Upper triangle
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        strain(i, j) = 0.5 * (G(i, j) - g(i, j));
      }
    }

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < i; j++)
      {
        strain(i, j) = strain(j, i);
      }
    }

    // Compute upper triangle of stress
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        // Initialise this component of sigma
        sigma(i, j) = 0.0;
        for (unsigned k = 0; k < dim; k++)
        {
          for (unsigned l = 0; l < dim; l++)
          {
            sigma(i, j) += C1 *
                           (Gup(i, k) * Gup(j, l) + Gup(i, l) * Gup(j, k) +
                            C2 * Gup(i, j) * Gup(k, l)) *
                           strain(k, l);
          }
        }
      }
    }

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < i; j++)
      {
        sigma(i, j) = sigma(j, i);
      }
    }
  }

  //===========================================================================
  /// Calculate the deviatoric part of the contravariant
  /// 2nd Piola Kirchoff stress tensor. Also return the contravariant
  /// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
  /// the inverse of the bulk modulus \f$ \kappa\f$. This form is appropriate
  /// for near-incompressible materials for which
  /// \f$ \sigma^{ij} = -p G^{ij} + \overline{ \sigma^{ij}}  \f$
  /// where the "pressure" \f$ p \f$ is determined from
  /// \f$ p / \kappa - d =0 \f$.
  //===========================================================================
  void GeneralisedHookean::calculate_second_piola_kirchhoff_stress(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    DenseMatrix<double>& sigma_dev,
    DenseMatrix<double>& Gup,
    double& gen_dil,
    double& inv_kappa)
  {
    // Find the dimension of the problem
    unsigned dim = G.nrow();

    // Assign memory for the determinant of the deformed metric tensor
    double detG = 0.0;

    // Compute deviatoric stress by calling the incompressible
    // version of this function
    calculate_second_piola_kirchhoff_stress(g, G, sigma_dev, Gup, detG);

    // Calculate the inverse of the "bulk" modulus
    inv_kappa =
      (1.0 - 2.0 * (*Nu_pt)) * (1.0 + (*Nu_pt)) / ((*E_pt) * (*Nu_pt));

    // Finally compute the generalised dilatation (i.e. the term that
    // must be zero if \kappa \to \infty
    gen_dil = 0.0;
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        gen_dil += Gup(i, j) * 0.5 * (G(i, j) - g(i, j));
      }
    }
  }

  //======================================================================
  /// Calculate the deviatoric part
  /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
  /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
  /// Also return the contravariant deformed metric
  /// tensor and the determinant of the deformed metric tensor.
  /// This form is appropriate
  /// for truly-incompressible materials for which
  /// \f$ \sigma^{ij} = - p G^{ij} +\overline{ \sigma^{ij}}  \f$
  /// where the "pressure" \f$ p \f$ is determined by
  /// \f$ \det G_{ij} - \det g_{ij} = 0 \f$.
  //======================================================================
  void GeneralisedHookean::calculate_second_piola_kirchhoff_stress(
    const DenseMatrix<double>& g,
    const DenseMatrix<double>& G,
    DenseMatrix<double>& sigma_dev,
    DenseMatrix<double>& Gup,
    double& detG)
  {
    // Error checking
#ifdef PARANOID
    error_checking_in_input(g, G, sigma_dev);
#endif

    // Find the dimension of the problem
    unsigned dim = G.nrow();

    // Calculate the contravariant Deformed metric tensor
    detG = calculate_contravariant(G, Gup);

    // Premultiply the appropriate physical constant
    double C1 = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));

    // Strain tensor
    DenseMatrix<double> strain(dim, dim);

    // Upper triangle
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        strain(i, j) = 0.5 * (G(i, j) - g(i, j));
      }
    }

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < i; j++)
      {
        strain(i, j) = strain(j, i);
      }
    }

    // Compute upper triangle of stress
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = i; j < dim; j++)
      {
        // Initialise this component of sigma
        sigma_dev(i, j) = 0.0;
        for (unsigned k = 0; k < dim; k++)
        {
          for (unsigned l = 0; l < dim; l++)
          {
            sigma_dev(i, j) += C1 *
                               (Gup(i, k) * Gup(j, l) + Gup(i, l) * Gup(j, k)) *
                               strain(k, l);
          }
        }
      }
    }

    // Copy across
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < i; j++)
      {
        sigma_dev(i, j) = sigma_dev(j, i);
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Calculate the contravariant 2nd Piola Kirchhoff
  /// stress tensor. Arguments are the
  /// covariant undeformed and deformed metric tensor and the
  /// matrix in which to return the stress tensor.
  /// Uses correct 3D invariants for 2D (plane strain) problems.
  //=======================================================================
  void IsotropicStrainEnergyFunctionConstitutiveLaw::
    calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                            const DenseMatrix<double>& G,
                                            DenseMatrix<double>& sigma)
  {
// Error checking
#ifdef PARANOID
    error_checking_in_input(g, G, sigma);
#endif

    // Find the dimension of the problem
    unsigned dim = g.nrow();

#ifdef PARANOID
    if (dim == 1)
    {
      std::string function_name =
        "IsotropicStrainEnergyFunctionConstitutiveLaw::";
      function_name += "calculate_second_piola_kirchhoff_stress()";

      throw OomphLibError("Check constitutive equations carefully when dim=1",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Calculate the contravariant undeformed and deformed metric tensors
    // and get the determinants of the metric tensors
    DenseMatrix<double> gup(dim), Gup(dim);
    double detg = calculate_contravariant(g, gup);
    double detG = calculate_contravariant(G, Gup);

    // Calculate the strain invariants
    Vector<double> I(3, 0.0);
    // The third strain invaraint is the volumetric change
    I[2] = detG / detg;
    // The first and second are a bit more complex --- see G&Z
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        I[0] += gup(i, j) * G(i, j);
        I[1] += g(i, j) * Gup(i, j);
      }
    }

    // If 2D we assume plane strain: In this case the 3D tensors have
    // a 1 on the diagonal and zeroes in the off-diagonals of their
    // third rows and columns. Only effect: Increase the first two
    // invariants by one; rest of the computation can just be performed
    // over the 2d set of coordinates.
    if (dim == 2)
    {
      I[0] += 1.0;
      I[1] += 1.0;
    }

    // Second strain invariant is multiplied by the third.
    I[1] *= I[2];

    // Calculate the derivatives of the strain energy function wrt the
    // strain invariants
    Vector<double> dWdI(3, 0.0);
    Strain_energy_function_pt->derivatives(I, dWdI);


    // Only bother to compute the tensor B^{ij} (Green & Zerna notation)
    // if the derivative wrt the second strain invariant is non-zero
    DenseMatrix<double> Bup(dim, dim, 0.0);
    if (std::fabs(dWdI[1]) > 0.0)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          Bup(i, j) = I[0] * gup(i, j);
          for (unsigned r = 0; r < dim; r++)
          {
            for (unsigned s = 0; s < dim; s++)
            {
              Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
            }
          }
        }
      }
    }

    // Now set the values of the functions phi, psi and p (Green & Zerna
    // notation) Note that the Green & Zerna stress \tau^{ij} is
    // s^{ij}/sqrt(I[2]), where s^{ij} is the desired second Piola-Kirchhoff
    // stress tensor so we multiply their constants by sqrt(I[2])
    double phi = 2.0 * dWdI[0];
    double psi = 2.0 * dWdI[1];
    double p = 2.0 * dWdI[2] * I[2];

    // Put it all together to get the stress
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        sigma(i, j) = phi * gup(i, j) + psi * Bup(i, j) + p * Gup(i, j);
      }
    }
  }

  //===========================================================================
  /// Calculate the deviatoric part
  /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
  /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
  /// Also return the contravariant deformed metric
  /// tensor and the determinant of the deformed metric tensor.
  /// Uses correct 3D invariants for 2D (plane strain) problems.
  /// This is the version for the pure incompressible formulation.
  //============================================================================
  void IsotropicStrainEnergyFunctionConstitutiveLaw::
    calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                            const DenseMatrix<double>& G,
                                            DenseMatrix<double>& sigma_dev,
                                            DenseMatrix<double>& Gup,
                                            double& detG)
  {
// Error checking
#ifdef PARANOID
    error_checking_in_input(g, G, sigma_dev);
#endif

    // Find the dimension of the problem
    unsigned dim = g.nrow();


#ifdef PARANOID
    if (dim == 1)
    {
      std::string function_name =
        "IsotropicStrainEnergyFunctionConstitutiveLaw::";
      function_name += "calculate_second_piola_kirchhoff_stress()";

      throw OomphLibError("Check constitutive equations carefully when dim=1",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Calculate the contravariant undeformed and deformed metric tensors
    DenseMatrix<double> gup(dim);
    // Don't need this determinant
    (void)calculate_contravariant(g, gup);
    // These are passed back
    detG = calculate_contravariant(G, Gup);

    // Calculate the strain invariants
    Vector<double> I(3, 0.0);
    // The third strain invaraint must be one (incompressibility)
    I[2] = 1.0;
    // The first and second are a bit more complex
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        I[0] += gup(i, j) * G(i, j);
        I[1] += g(i, j) * Gup(i, j);
      }
    }

    // If 2D we assume plane strain: In this case the 3D tensors have
    // a 1 on the diagonal and zeroes in the off-diagonals of their
    // third rows and columns. Only effect: Increase the first two
    // invariants by one; rest of the computation can just be performed
    // over the 2d set of coordinates.
    if (dim == 2)
    {
      I[0] += 1.0;
      I[1] += 1.0;
    }

    // Calculate the derivatives of the strain energy function wrt the
    // strain invariants
    Vector<double> dWdI(3, 0.0);
    Strain_energy_function_pt->derivatives(I, dWdI);

    // Only bother to compute the tensor B^{ij} (Green & Zerna notation)
    // if the derivative wrt the second strain invariant is non-zero
    DenseMatrix<double> Bup(dim, dim, 0.0);
    if (std::fabs(dWdI[1]) > 0.0)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          Bup(i, j) = I[0] * gup(i, j);
          for (unsigned r = 0; r < dim; r++)
          {
            for (unsigned s = 0; s < dim; s++)
            {
              Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
            }
          }
        }
      }
    }

    // Now set the values of the functions phi and psi (Green & Zerna notation)
    double phi = 2.0 * dWdI[0];
    double psi = 2.0 * dWdI[1];
    // Calculate the trace/dim of the first two terms of the stress tensor
    // phi g^{ij} + psi B^{ij} (see Green & Zerna)
    double K;
    // In two-d, we cannot use the strain invariants directly
    // but can use symmetry of the tensors involved
    if (dim == 2)
    {
      K = 0.5 * ((I[0] - 1.0) * phi +
                 psi * (Bup(0, 0) * G(0, 0) + Bup(1, 1) * G(1, 1) +
                        2.0 * Bup(0, 1) * G(0, 1)));
    }
    // In three-d we can make use of the strain invariants, see Green & Zerna
    else
    {
      K = (I[0] * phi + 2.0 * I[1] * psi) / 3.0;
    }

    // Put it all together to get the stress, subtracting the trace of the
    // first two terms to ensure that the stress is deviatoric, which means
    // that the computed pressure is the mechanical pressure
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        sigma_dev(i, j) = phi * gup(i, j) + psi * Bup(i, j) - K * Gup(i, j);
      }
    }
  }

  //===========================================================================
  /// Calculate the deviatoric part of the contravariant
  /// 2nd Piola Kirchoff stress tensor. Also return the contravariant
  /// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
  /// the inverse of the bulk modulus \f$ \kappa\f$.
  /// Uses correct 3D invariants for 2D (plane strain) problems.
  /// This is the version for the near-incompressible formulation.
  //===========================================================================
  void IsotropicStrainEnergyFunctionConstitutiveLaw::
    calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                            const DenseMatrix<double>& G,
                                            DenseMatrix<double>& sigma_dev,
                                            DenseMatrix<double>& Gup,
                                            double& gen_dil,
                                            double& inv_kappa)
  {
// Error checking
#ifdef PARANOID
    error_checking_in_input(g, G, sigma_dev);
#endif

    // Find the dimension of the problem
    unsigned dim = g.nrow();

#ifdef PARANOID
    if (dim == 1)
    {
      std::string function_name =
        "IsotropicStrainEnergyFunctionConstitutiveLaw::";
      function_name += "calculate_second_piola_kirchhoff_stress()";

      throw OomphLibError("Check constitutive equations carefully when dim=1",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Calculate the contravariant undeformed and deformed metric tensors
    // and get the determinants of the metric tensors
    DenseMatrix<double> gup(dim);
    double detg = calculate_contravariant(g, gup);
    double detG = calculate_contravariant(G, Gup);

    // Calculate the strain invariants
    Vector<double> I(3, 0.0);
    // The third strain invaraint is the volumetric change
    I[2] = detG / detg;
    // The first and second are a bit more complex --- see G&Z
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        I[0] += gup(i, j) * G(i, j);
        I[1] += g(i, j) * Gup(i, j);
      }
    }

    // If 2D we assume plane strain: In this case the 3D tensors have
    // a 1 on the diagonal and zeroes in the off-diagonals of their
    // third rows and columns. Only effect: Increase the first two
    // invariants by one; rest of the computation can just be performed
    // over the 2d set of coordinates.
    if (dim == 2)
    {
      I[0] += 1.0;
      I[1] += 1.0;
    }

    // Second strain invariant is multiplied by the third.
    I[1] *= I[2];

    // Calculate the derivatives of the strain energy function wrt the
    // strain invariants
    Vector<double> dWdI(3, 0.0);
    Strain_energy_function_pt->derivatives(I, dWdI);

    // Only bother to calculate the tensor B^{ij} (Green & Zerna notation)
    // if the derivative wrt the second strain invariant is non-zero
    DenseMatrix<double> Bup(dim, dim, 0.0);
    if (std::fabs(dWdI[1]) > 0.0)
    {
      for (unsigned i = 0; i < dim; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          Bup(i, j) = I[0] * gup(i, j);
          for (unsigned r = 0; r < dim; r++)
          {
            for (unsigned s = 0; s < dim; s++)
            {
              Bup(i, j) -= gup(i, r) * gup(j, s) * G(r, s);
            }
          }
        }
      }
    }

    // Now set the values of the functions phi and psi (Green & Zerna notation)
    // but multiplied by sqrt(I[2]) to recover the second Piola-Kirchhoff stress
    double phi = 2.0 * dWdI[0];
    double psi = 2.0 * dWdI[1];

    // Calculate the trace/dim of the first two terms of the stress tensor
    // phi g^{ij} + psi B^{ij} (see Green & Zerna)
    double K;
    // In two-d, we cannot use the strain invariants directly,
    // but we can use symmetry of the tensors involved
    if (dim == 2)
    {
      K = 0.5 * ((I[0] - 1.0) * phi +
                 psi * (Bup(0, 0) * G(0, 0) + Bup(1, 1) * G(1, 1) +
                        2.0 * Bup(0, 1) * G(0, 1)));
    }
    // In three-d we can make use of the strain invariants
    else
    {
      K = (I[0] * phi + 2.0 * I[1] * psi) / 3.0;
    }

    // Choose inverse kappa to be one...
    inv_kappa = 1.0;

    //...then the generalised dilation is the same as p  in Green & Zerna's
    // notation, but multiplied by sqrt(I[2]) with the addition of the
    // terms that are subtracted to make the other part of the stress
    // deviatoric
    gen_dil = 2.0 * dWdI[2] * I[2] + K;

    // Calculate the deviatoric part of the stress by subtracting
    // the computed trace/dim
    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned j = 0; j < dim; j++)
      {
        sigma_dev(i, j) = phi * gup(i, j) + psi * Bup(i, j) - K * Gup(i, j);
      }
    }
  }

} // namespace oomph
