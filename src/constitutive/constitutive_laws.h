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
// Header file for ConstitutiveLaw objects that will be used in all
// elasticity-type elements

#ifndef OOMPH_CONSTITUTIVE_LAWS_HEADER
#define OOMPH_CONSTITUTIVE_LAWS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// OOMPH-LIB includes
#include "generic/oomph_utilities.h"
#include "generic/matrices.h"

namespace oomph
{
  //=====================================================================
  /// \short Base class for strain energy functions to be used in solid
  /// mechanics computations.
  //====================================================================
  class StrainEnergyFunction
  {
  public:
    /// Constructor takes no arguments
    StrainEnergyFunction() {}

    /// Empty virtual destructor
    virtual ~StrainEnergyFunction() {}


    /// Return the strain energy in terms of the strain tensor
    virtual double W(const DenseMatrix<double>& gamma)
    {
      std::string error_message =
        "The strain-energy function as a function of the strain-tensor,\n";
      error_message +=
        "gamma, is not implemented for this strain energy function.\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      return 0.0;
    }


    /// Return the strain energy in terms of the strain invariants
    virtual double W(const Vector<double>& I)
    {
      std::string error_message =
        "The strain-energy function as a function of the strain\n ";
      error_message +=
        "invariants, I1, I2, I3, is not implemented for this strain\n ";
      error_message += "energy function\n";

      throw OomphLibError(
        error_message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      return 0.0;
    }


    /// \short Return the derivatives of the strain energy function with
    /// respect to the components of the strain tensor (default is to use
    /// finite differences).
    virtual void derivative(const DenseMatrix<double>& gamma,
                            DenseMatrix<double>& dWdgamma)
    {
      throw OomphLibError(
        "Sorry, the FD setup of dW/dgamma hasn't been implemented yet",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    /// \short Return the derivatives of the strain energy function with
    /// respect to the strain invariants. Default version is to use finite
    /// differences
    virtual void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      // Calculate the derivatives of the strain-energy-function wrt the strain
      // invariants
      double FD_Jstep = 1.0e-8; // Usual comments about global stuff
      double energy = W(I);

      // Loop over the strain invariants
      for (unsigned i = 0; i < 3; i++)
      {
        // Store old value
        double I_prev = I[i];
        // Increase ith strain invariant
        I[i] += FD_Jstep;
        // Get the new value of the strain energy
        double energy_new = W(I);
        // Calculate the value of the derivative
        dWdI[i] = (energy_new - energy) / FD_Jstep;
        // Reset value of ith strain invariant
        I[i] = I_prev;
      }
    }

    /// \short Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode.
    virtual bool requires_incompressibility_constraint() = 0;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// \short MooneyRivlin strain-energy function.
  /// with constitutive parameters C1 and C2:
  /// \f[
  /// W = C_1 (I_0 - 3) + C_2 (I_1 - 3)
  /// \f]
  /// where incompressibility (\f$ I_2 \equiv 1\f$) is assumed.
  //====================================================================
  class MooneyRivlin : public StrainEnergyFunction
  {
  public:
    /// Constructor takes the pointer to the value of the constants
    MooneyRivlin(double* c1_pt, double* c2_pt)
      : StrainEnergyFunction(), C1_pt(c1_pt), C2_pt(c2_pt)
    {
    }

    /// Empty Virtual destructor
    virtual ~MooneyRivlin() {}

    /// Return the strain energy in terms of strain tensor
    double W(const DenseMatrix<double>& gamma)
    {
      return StrainEnergyFunction::W(gamma);
    }

    /// Return the strain energy in terms of the strain invariants
    double W(const Vector<double>& I)
    {
      return (*C1_pt) * (I[0] - 3.0) + (*C2_pt) * (I[1] - 3.0);
    }


    /// \short Return the derivatives of the strain energy function with
    /// respect to the strain invariants
    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      dWdI[0] = (*C1_pt);
      dWdI[1] = (*C2_pt);
      dWdI[2] = 0.0;
    }

    /// \short Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. True
    bool requires_incompressibility_constraint()
    {
      return true;
    }


  private:
    /// Pointer to first Mooney Rivlin constant
    double* C1_pt;

    /// Pointer to second Mooney Rivlin constant
    double* C2_pt;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// \short Generalisation of Mooney Rivlin constitutive law to compressible
  /// media as suggested on p. 553 of Fung, Y.C. & Tong, P. "Classical and
  /// Computational Solid Mechanics" World Scientific (2001).
  /// Input parameters are Young's modulus E, Poisson ratio nu and
  /// the Mooney-Rivlin constant C1. In the small-deformation-limit
  /// the behaviour becomes equivalent to that of linear elasticity
  /// with the same E and nu.
  ///
  /// Note that there's a factor of 2 difference between C1 and the Mooney
  /// Rivlin C1!
  //====================================================================
  class GeneralisedMooneyRivlin : public StrainEnergyFunction
  {
  public:
    /// \short Constructor takes the pointers to the constitutive parameters:
    /// Poisson's ratio, the Mooney-Rivlin parameter. Young's modulus is set
    /// to 1, implying that it has been used to scale the stresses
    GeneralisedMooneyRivlin(double* nu_pt, double* c1_pt)
      : StrainEnergyFunction(),
        Nu_pt(nu_pt),
        C1_pt(c1_pt),
        E_pt(new double(1.0)),
        Must_delete_e(true)
    {
    }

    /// \short Constructor takes the pointers to the constitutive parameters:
    /// Poisson's ratio, the Mooney-Rivlin parameter and Young's modulus
    GeneralisedMooneyRivlin(double* nu_pt, double* c1_pt, double* e_pt)
      : StrainEnergyFunction(),
        Nu_pt(nu_pt),
        C1_pt(c1_pt),
        E_pt(e_pt),
        Must_delete_e(false)
    {
    }


    /// Virtual destructor
    virtual ~GeneralisedMooneyRivlin()
    {
      if (Must_delete_e) delete E_pt;
    }

    /// Return the strain energy in terms of strain tensor
    double W(const DenseMatrix<double>& gamma)
    {
      return StrainEnergyFunction::W(gamma);
    }


    /// Return the strain energy in terms of the strain invariants
    double W(const Vector<double>& I)
    {
      double G = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));
      return 0.5 * ((*C1_pt) * (I[0] - 3.0) + (G - (*C1_pt)) * (I[1] - 3.0) +
                    ((*C1_pt) - 2.0 * G) * (I[2] - 1.0) +
                    (1.0 - (*Nu_pt)) * G * (I[2] - 1.0) * (I[2] - 1.0) /
                      (2.0 * (1.0 - 2.0 * (*Nu_pt))));
    }


    /// \short Return the derivatives of the strain energy function with
    /// respect to the strain invariants
    void derivatives(Vector<double>& I, Vector<double>& dWdI)
    {
      double G = (*E_pt) / (2.0 * (1.0 + (*Nu_pt)));
      dWdI[0] = 0.5 * (*C1_pt);
      dWdI[1] = 0.5 * (G - (*C1_pt));
      dWdI[2] = 0.5 * ((*C1_pt) - 2.0 * G +
                       2.0 * (1.0 - (*Nu_pt)) * G * (I[2] - 1.0) /
                         (2.0 * (1.0 - 2.0 * (*Nu_pt))));
    }


    /// \short Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. False.
    bool requires_incompressibility_constraint()
    {
      return false;
    }

  private:
    /// Poisson's ratio
    double* Nu_pt;

    /// Mooney-Rivlin parameter
    double* C1_pt;

    /// Young's modulus
    double* E_pt;

    /// \short Boolean flag to indicate if storage for elastic modulus
    /// must be deleted in destructor
    bool Must_delete_e;
  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //===========================================================================
  /// A class for constitutive laws for elements that solve
  /// the equations of solid mechanics based upon the principle of virtual
  /// displacements. In that formulation, the information required from a
  /// constitutive law is the (2nd Piola-Kirchhoff) stress tensor
  /// \f$ \sigma^{ij} \f$ as a function of the (Green) strain
  /// \f$ \gamma^{ij} \f$:
  /// \f[
  /// \sigma^{ij} = \sigma^{ij}(\gamma_{ij}).
  /// \f]
  /// The Green strain is defined as
  /// \f[
  /// \gamma_{ij} = \frac{1}{2} (G_{ij} - g_{ij}), \ \ \ \ \ \ \ \ \ \ \ (1)
  /// \f]
  /// where \f$G_{ij} \f$ and \f$ g_{ij}\f$ are the metric tensors
  /// in the deformed and undeformed (stress-free) configurations, respectively.
  /// A specific ConstitutiveLaw needs to be implement the pure
  /// virtual function
  /// \code
  /// ConstitutiveLaw::calculate_second_piola_kirchhoff_stress(...)
  /// \endcode
  /// Equation (1) shows that the strain may be calculated from the
  /// undeformed and deformed metric tensors.  Frequently, these tensors are
  /// also required in the constitutive law itself.
  /// To avoid unnecessary re-computation of these quantities, we
  /// pass the deformed and undeformed metric tensor to
  /// \c calculate_second_piola_kirchhoff_stress(...)
  /// rather than the strain tensor itself.
  ///
  /// The functional form of the constitutive equation is different
  /// for compressible/incompressible/near-incompressible behaviour
  /// and we provide interfaces that are appropriate for all of these cases.
  /// -# \b Compressible \b Behaviour: \n If the material is compressible,
  ///    the stress can be computed from the deformed and undeformed
  ///    metric tensors,
  ///    \f[
  ///   \sigma^{ij} = \sigma^{ij}(\gamma_{ij}) =
  ///   \sigma^{ij}\bigg( \frac{1}{2} (G_{ij} - g_{ij})\bigg),
  ///   \f]
  ///   using the interface
  ///   \code
  ///   // 2nd Piola Kirchhoff stress tensor
  ///   DenseMatrix<double> sigma(DIM,DIM);
  ///
  ///   // Metric tensor in the undeformed (stress-free) configuration
  ///   DenseMatrix<double> g(DIM,DIM);
  ///
  ///   // Metric tensor in the deformed  configuration
  ///   DenseMatrix<double> G(DIM,DIM);
  ///
  ///   // Compute stress from the two metric tensors:
  ///   calculate_second_piola_kirchhoff_stress(g,G,sigma);
  ///   \endcode
  ///   \n \n  \n
  /// -# \b Incompressible \b Behaviour: \n If the material is incompressible,
  ///    its deformation is constrained by the condition that
  ///    \f[
  ///    \det G_{ij} - \det g_{ij}= 0
  ///    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (2)
  ///    \f]
  ///    which ensures that the volume of infinitesimal material
  ///    elements remains constant during the deformation. This
  ///    condition is typically enforced by a Lagrange multiplier which
  ///    plays the role of a pressure. In such cases, the
  ///    stress tensor has form
  ///    \f[
  ///    \sigma^{ij} = -p G^{ij} +
  ///    \overline{\sigma}^{ij}\big(\gamma_{kl}\big),
  ///    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)
  ///    \f]
  ///    where only the deviatoric part of the stress tensor,
  ///    \f$ \overline{\sigma}^{ij}, \f$ depends directly on the
  ///    strain. The pressure \f$ p \f$ needs to be determined
  ///    independently from (2).
  ///    Given the deformed and undeformed metric tensors,
  ///    the computation of the stress tensor \f$ \sigma^{ij} \f$
  ///    for an incompressible
  ///    material therefore requires the computation of the following
  ///    quantities:
  ///    - The deviatoric stress \f$ \overline{\sigma}^{ij} \f$
  ///    - The contravariant deformed metric tensor  \f$ G^{ij} \f$
  ///    - The determinant of the deformed
  ///      metric tensors, \f$ \det G_{ij}, \f$ which
  ///      is required in equation (2) whose solution determines the pressure.
  ///    .
  ///   \n
  ///   These quantities can be obtained from the following interface \n
  ///   \code
  ///   // Deviatoric part of the 2nd Piola Kirchhoff stress tensor
  ///   DenseMatrix<double> sigma_dev(DIM,DIM);
  ///
  ///   // Metric tensor in the undeformed (stress-free) configuration
  ///   DenseMatrix<double> g(DIM,DIM);
  ///
  ///   // Metric tensor in the deformed  configuration
  ///   DenseMatrix<double> G(DIM,DIM);
  ///
  ///   // Determinant of the deformed metric tensor
  ///   double Gdet;
  ///
  ///   // Contravariant deformed metric tensor
  ///   DenseMatrix<double> G_contra(DIM,DIM);
  ///
  ///   // Compute stress from the two metric tensors:
  ///   calculate_second_piola_kirchhoff_stress(g,G,sigma_dev,G_contra,Gdet);
  ///   \endcode
  ///   \n \n \n
  /// -# \b Nearly \b Incompressible \b Behaviour: \n If the material is nearly
  ///    incompressible, it is advantageous to split the stress into
  ///    its deviatoric and hydrostatic parts by writing the
  ///    constitutive law in the form
  ///    \f[
  ///    \sigma^{ij} = -p G^{ij} +
  ///    \overline{\sigma}^{ij}\big(\gamma_{kl}\big),
  ///    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (3)
  ///    \f]
  ///    where the deviatoric part of the stress tensor,
  ///    \f$ \overline{\sigma}^{ij}, \f$ depends on the
  ///    strain. This form of the constitutive
  ///    law is identical to that of the incompressible
  ///    case and it involves a pressure \f$ p \f$ which needs to be
  ///    determined from an additional equation. In the
  ///    incompressible case, this equation was given by the incompressibility
  ///    constraint (2). Here, we need to augment the constitutive law (3) by
  ///    a separate equation for the pressure. Generally this takes the
  ///    form
  ///    \f[
  ///    p = - \kappa \ d
  ///    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (4)
  ///    \f]
  ///    where \f$ \kappa \f$ is the "bulk modulus", a material property
  ///    that needs to be specified by the constitutive law.
  ///    \f$ d \f$ is the (generalised) dilatation, i.e. the relative change
  ///    in the volume of an infinitesimal material element (or some
  ///    suitable generalised quantitiy that is related to it). As the
  ///    material approaches incompressibility, \f$ \kappa \to \infty\f$, so
  ///    that infinitely large pressures would be required to achieve any change
  ///    in volume. To facilitate the implementation of (4) as the equation for
  ///    the pressure, we re-write it in the form \f[ p \ \frac{1}{\kappa} +
  ///    d\big(g_{ij},G_{ij}\big) = 0 \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (5) \f]
  ///    which only involves quantities that remain finite
  ///    as we approach true incompressibility.
  ///    \n
  ///    Given the deformed and undeformed metric tensors,
  ///    the computation of the stress tensor \f$ \sigma^{ij} \f$
  ///    for a nearly incompressible
  ///    material therefore requires the computation of the following
  ///    quantities:
  ///    - The deviatoric stress \f$ \overline{\sigma}^{ij} \f$
  ///    - The contravariant deformed metric tensor  \f$ G^{ij} \f$
  ///    - The generalised dilatation \f$ d \f$
  ///    - The inverse of the bulk modulus \f$ \kappa \f$
  ///    .
  ///   \n
  ///   These quantities can be obtained from the following interface \n
  ///   \code
  ///   // Deviatoric part of the 2nd Piola Kirchhoff stress tensor
  ///   DenseMatrix<double> sigma_dev(DIM,DIM);
  ///
  ///   // Metric tensor in the undeformed (stress-free) configuration
  ///   DenseMatrix<double> g(DIM,DIM);
  ///
  ///   // Metric tensor in the deformed  configuration
  ///   DenseMatrix<double> G(DIM,DIM);
  ///
  ///   // Contravariant deformed metric tensor
  ///   DenseMatrix<double> G_contra(DIM,DIM);
  ///
  ///   // Inverse of the bulk modulus
  ///   double inv_kappa;
  ///
  ///   // Generalised dilatation
  ///   double gen_dil;
  ///
  ///   // Compute stress from the two metric tensors:
  ///   calculate_second_piola_kirchhoff_stress(g,G,sigma_dev,G_contra,inv_kappa,gen_dil);
  ///   \endcode
  //==========================================================================
  class ConstitutiveLaw
  {
  protected:
    /// \short Test whether a matrix is square
    bool is_matrix_square(const DenseMatrix<double>& M);

    /// \short Test whether two matrices are of equal dimensions
    bool are_matrices_of_equal_dimensions(const DenseMatrix<double>& M1,
                                          const DenseMatrix<double>& M2);

    /// \short Check for errors in the input,
    /// i.e. check that the dimensions of the arrays are all consistent
    void error_checking_in_input(const DenseMatrix<double>& g,
                                 const DenseMatrix<double>& G,
                                 DenseMatrix<double>& sigma);

    /// \short Calculate a contravariant tensor from a covariant tensor,
    /// and return the determinant of the covariant tensor.
    double calculate_contravariant(const DenseMatrix<double>& Gcov,
                                   DenseMatrix<double>& Gcontra);

    /// \short Calculate the derivatives of the contravariant tensor
    /// and the derivatives of the determinant of the covariant tensor
    /// with respect to the components of the covariant tensor
    void calculate_d_contravariant_dG(const DenseMatrix<double>& Gcov,
                                      RankFourTensor<double>& dGcontra_dG,
                                      DenseMatrix<double>& d_detG_dG);


  public:
    /// Empty constructor
    ConstitutiveLaw() {}


    /// Empty virtual destructor
    virtual ~ConstitutiveLaw() {}


    /// \short Calculate the contravariant 2nd Piola Kirchhoff
    /// stress tensor. Arguments are the
    /// covariant undeformed and deformed metric tensor and the
    /// matrix in which to return the stress tensor
    virtual void calculate_second_piola_kirchhoff_stress(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      DenseMatrix<double>& sigma) = 0;

    /// \short Calculate the derivatives of the contravariant
    /// 2nd Piola Kirchhoff stress tensor with respect to the deformed metric
    /// tensor. Arguments are the
    /// covariant undeformed and deformed metric tensor, the current value of
    /// the stress tensor and the
    /// rank four tensor in which to return the derivatives of the stress tensor
    /// The default implementation uses finite differences, but can be
    /// overloaded for constitutive laws in which an analytic formulation
    /// is possible.
    /// If the boolean flag symmetrize_tensor is false, only the
    /// "upper  triangular" entries of the tensor will be filled in. This is
    /// a useful efficiency when using the derivatives in Jacobian calculations.
    virtual void calculate_d_second_piola_kirchhoff_stress_dG(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      RankFourTensor<double>& d_sigma_dG,
      const bool& symmetrize_tensor = true);


    /// \short Calculate the deviatoric part
    /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// Also return the contravariant deformed metric
    /// tensor and the determinant of the deformed metric tensor.
    /// This form is appropriate
    /// for truly-incompressible materials for which
    /// \f$ \sigma^{ij} = - p G^{ij} +\overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined by
    /// \f$ \det G_{ij} - \det g_{ij} = 0 \f$.
    virtual void calculate_second_piola_kirchhoff_stress(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      DenseMatrix<double>& sigma_dev,
      DenseMatrix<double>& G_contra,
      double& Gdet)
    {
      throw OomphLibError(
        "Incompressible formulation not implemented for this constitutive law",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Calculate the derivatives of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// with respect to the deformed metric tensor.
    /// Also return the derivatives of the determinant of the
    /// deformed metric tensor with respect to the deformed metric tensor.
    /// This form is appropriate
    /// for truly-incompressible materials.
    /// The default implementation uses finite differences for the
    /// derivatives that depend on the constitutive law, but not
    /// for the derivatives of the determinant, which are generic.
    //// If the boolean flag symmetrize_tensor is false, only the
    /// "upper  triangular" entries of the tensor will be filled in. This is
    /// a useful efficiency when using the derivatives in Jacobian calculations.
    virtual void calculate_d_second_piola_kirchhoff_stress_dG(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      const double& detG,
      const double& interpolated_solid_p,
      RankFourTensor<double>& d_sigma_dG,
      DenseMatrix<double>& d_detG_dG,
      const bool& symmetrize_tensor = true);


    /// \short Calculate the deviatoric part of the contravariant
    /// 2nd Piola Kirchoff stress tensor. Also return the contravariant
    /// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
    /// the inverse of the bulk modulus \f$ \kappa\f$. This form is appropriate
    /// for near-incompressible materials for which
    /// \f$ \sigma^{ij} = -p G^{ij} + \overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined from
    /// \f$ p / \kappa - d =0 \f$.
    virtual void calculate_second_piola_kirchhoff_stress(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      DenseMatrix<double>& sigma_dev,
      DenseMatrix<double>& Gcontra,
      double& gen_dil,
      double& inv_kappa)
    {
      throw OomphLibError(
        "Near-incompressible formulation not implemented for constitutive law",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// \short Calculate the derivatives of the contravariant
    /// 2nd Piola Kirchoff stress tensor with respect to the deformed metric
    /// tensor. Also return the derivatives of the generalised dilatation,
    /// \f$ d, \f$ with respect to the deformed metric tensor.
    /// This form is appropriate
    /// for near-incompressible materials.
    /// The default implementation uses finite differences.
    /// If the boolean flag symmetrize_tensor is false, only the
    /// "upper  triangular" entries of the tensor will be filled in. This is
    /// a useful efficiency when using the derivatives in Jacobian calculations.
    virtual void calculate_d_second_piola_kirchhoff_stress_dG(
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      const double& gen_dil,
      const double& inv_kappa,
      const double& interpolated_solid_p,
      RankFourTensor<double>& d_sigma_dG,
      DenseMatrix<double>& d_gen_dil_dG,
      const bool& symmetrize_tensor = true);


    /// \short Pure virtual function in which the user must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode.
    virtual bool requires_incompressibility_constraint() = 0;
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //========================================================================
  /// Class for a "non-rational" extension of classical linear elasticity
  /// to large displacements:
  /// \f[ \sigma^{ij} = E^{ijkl} \gamma_{kl} \f]
  /// where
  /// \f[ E^{ijkl} = \frac{E}{(1+\nu)} \left( \frac{\nu}{(1-2\nu)} G^{ij} G^{kl}
  ///                + \frac{1}{2} \left(
  ///                G^{ik} G^{jl} + G^{il} G^{jk} \right) \right) \f]
  /// For small strains \f$ (| G_{ij} - g_{ij} | \ll 1)\f$ this approaches
  /// the version appropriate for linear elasticity, obtained
  /// by replacing \f$ G^{ij}\f$ with \f$ g^{ij}\f$.
  ///
  /// We provide three versions of \c calculate_second_piola_kirchhoff_stress():
  /// -# If \f$ \nu \ne 1/2 \f$ (and not close to \f$ 1/2 \f$), the
  ///    constitutive law can be used directly in the above form, using
  ///    the deformed and undeformed metric tensors as input.
  /// -# If the material is incompressible (\f$ \nu = 1/2 \f$),
  ///    the first term in the above expression for  \f$ E^{ijkl} \f$
  ///    is singular. We re-write the constitutive equation for this
  ///    case as
  ///    \f[ \sigma^{ij} = -p G^{ij}
  ///                + \frac{E}{3} \left(
  ///                G^{ik} G^{jl} + G^{il} G^{jk} \right) \gamma_{kl} \f]
  ///    where the pressure \f$ p \f$ needs to be determined independently
  ///    via the incompressibility constraint.
  ///    In this case, the stress returned by
  ///    \c calculate_second_piola_kirchhoff_stress()
  ///    contains only the deviatoric part of the 2nd Piola Kirchhoff stress,
  ///    \f[ \overline{\sigma}^{ij} =
  ///                \frac{E}{3} \left(
  ///                G^{ik} G^{jl} + G^{il} G^{jk} \right) \gamma_{kl}. \f]
  ///    The function also returns the contravariant metric tensor
  ///    \f$ G^{ij}\f$ (since it is needed to form the complete stress
  ///    tensor), and the determinant of the deformed covariant metric
  ///    tensor \f$ {\tt detG} = \det G_{ij} \f$ (since it is needed
  ///    in the equation that enforces the incompressibility).
  /// -# If  \f$ \nu \approx 1/2 \f$, the original form of the
  ///     constitutive equation could be used, but the resulting
  ///     equations tend to be ill-conditioned since they contain
  ///     the product of the large "bulk modulus"
  ///     \f[ \kappa = \frac{E\nu}{(1+\nu)(1-2\nu)} \f]
  ///     and the small "generalised dilatation"
  ///     \f[ d = \frac{1}{2} G^{ij} (G_{ij}-g_{ij}). \f]
  ///     [\f$ d \f$ represents the actual dilatation in the small
  ///     strain limit; for large deformations it doesn't have
  ///     any sensible interpretation (or does it?). It is simply
  ///     the term that needs to go to zero as \f$ \kappa \to \infty\f$.]
  ///     In this case, the stress returned by
  ///     \c calculate_second_piola_kirchhoff_stress()
  ///     contains only the deviatoric part of the 2nd Piola Kirchhoff stress,
  ///     \f[ \overline{\sigma}^{ij} =
  ///                  \frac{E}{3} \left(
  ///                 G^{ik} G^{jl} + G^{il} G^{jk} \right) \gamma_{kl}. \f]
  ///     The function also returns the contravariant metric tensor
  ///     \f$ G^{ij}\f$ (since it is needed to form the complete stress
  ///     tensor), the inverse of the bulk modulus, and the generalised
  ///     dilatation (since they are needed in the equation
  ///     that determines the pressure).
  ///
  //=========================================================================
  class GeneralisedHookean : public ConstitutiveLaw
  {
  public:
    /// The constructor takes the pointers to values of material parameters:
    /// Poisson's ratio and Young's modulus.
    GeneralisedHookean(double* nu_pt, double* e_pt)
      : ConstitutiveLaw(), Nu_pt(nu_pt), E_pt(e_pt), Must_delete_e(false)
    {
    }

    /// The constructor takes the pointers to value of
    /// Poisson's ratio . Young's modulus is set to E=1.0,
    /// implying that all stresses have been non-dimensionalised
    /// on on it.
    GeneralisedHookean(double* nu_pt)
      : ConstitutiveLaw(),
        Nu_pt(nu_pt),
        E_pt(new double(1.0)),
        Must_delete_e(true)
    {
    }


    /// Virtual destructor
    virtual ~GeneralisedHookean()
    {
      if (Must_delete_e) delete E_pt;
    }

    /// \short Calculate the contravariant 2nd Piola Kirchhoff
    /// stress tensor. Arguments are the
    /// covariant undeformed and deformed metric tensor and the
    /// matrix in which to return the stress tensor
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma);


    /// \short Calculate the deviatoric part
    /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// Also return the contravariant deformed metric
    /// tensor and the determinant of the deformed metric tensor.
    /// This form is appropriate
    /// for truly-incompressible materials for which
    /// \f$ \sigma^{ij} = - p G^{ij} +\overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined by
    /// \f$ \det G_{ij} - \det g_{ij} = 0 \f$.
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma_dev,
                                                 DenseMatrix<double>& G_contra,
                                                 double& Gdet);


    /// \short Calculate the deviatoric part of the contravariant
    /// 2nd Piola Kirchoff stress tensor. Also return the contravariant
    /// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
    /// the inverse of the bulk modulus \f$ \kappa\f$. This form is appropriate
    /// for near-incompressible materials for which
    /// \f$ \sigma^{ij} = -p G^{ij} + \overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined from
    /// \f$ p / \kappa - d =0 \f$.
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma_dev,
                                                 DenseMatrix<double>& Gcontra,
                                                 double& gen_dil,
                                                 double& inv_kappa);


    /// \short Pure virtual function in which the writer must declare if the
    /// constitutive equation requires an incompressible formulation
    /// in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. False.
    bool requires_incompressibility_constraint()
    {
      return false;
    }

  private:
    /// Poisson ratio
    double* Nu_pt;

    /// Young's modulus
    double* E_pt;

    /// \short Boolean flag to indicate if storage for elastic modulus
    /// must be deleted in destructor
    bool Must_delete_e;
  };


  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////


  //=====================================================================
  /// A class for constitutive laws derived from strain-energy functions.
  /// Theory is in Green and Zerna.
  //=====================================================================
  class IsotropicStrainEnergyFunctionConstitutiveLaw : public ConstitutiveLaw
  {
  private:
    /// Pointer to the strain energy function
    StrainEnergyFunction* Strain_energy_function_pt;

  public:
    /// Constructor takes a pointer to the strain energy function
    IsotropicStrainEnergyFunctionConstitutiveLaw(
      StrainEnergyFunction* const& strain_energy_function_pt)
      : ConstitutiveLaw(), Strain_energy_function_pt(strain_energy_function_pt)
    {
    }

    /// \short Calculate the contravariant 2nd Piola Kirchhoff
    /// stress tensor. Arguments are the
    /// covariant undeformed and deformed metric tensor and the
    /// matrix in which to return the stress tensor.
    /// Uses correct 3D invariants for 2D (plane strain) problems.
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma);


    /// \short Calculate the deviatoric part
    /// \f$ \overline{ \sigma^{ij}}\f$  of the contravariant
    /// 2nd Piola Kirchhoff stress tensor \f$ \sigma^{ij}\f$.
    /// Also return the contravariant deformed metric
    /// tensor and the determinant of the deformed metric tensor.
    /// This form is appropriate
    /// for truly-incompressible materials for which
    /// \f$ \sigma^{ij} = - p G^{ij} +\overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined by
    /// \f$ \det G_{ij} - \det g_{ij} = 0 \f$.
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma_dev,
                                                 DenseMatrix<double>& G_contra,
                                                 double& Gdet);


    /// \short Calculate the deviatoric part of the contravariant
    /// 2nd Piola Kirchoff stress tensor. Also return the contravariant
    /// deformed metric tensor, the generalised dilatation, \f$ d, \f$ and
    /// the inverse of the bulk modulus \f$ \kappa\f$. This form is appropriate
    /// for near-incompressible materials for which
    /// \f$ \sigma^{ij} = -p G^{ij} + \overline{ \sigma^{ij}}  \f$
    /// where the "pressure" \f$ p \f$ is determined from
    /// \f$ p / \kappa - d =0 \f$.
    void calculate_second_piola_kirchhoff_stress(const DenseMatrix<double>& g,
                                                 const DenseMatrix<double>& G,
                                                 DenseMatrix<double>& sigma_dev,
                                                 DenseMatrix<double>& Gcontra,
                                                 double& gen_dil,
                                                 double& inv_kappa);


    /// \short State if the constitutive equation requires an incompressible
    /// formulation in which the volume constraint is enforced explicitly.
    /// Used as a sanity check in PARANOID mode. This is determined
    /// by interrogating the associated strain energy function.
    bool requires_incompressibility_constraint()
    {
      return Strain_energy_function_pt->requires_incompressibility_constraint();
    }
  };

} // namespace oomph

#endif
