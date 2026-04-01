#ifndef OOMPH_PLASTICITY_CONSTITUTIVE_HEADER
#define OOMPH_PLASTICITY_CONSTITUTIVE_HEADER

#include "constitutive/constitutive_laws.h"
#include "generic/matrix_helpers.h"
#include "generic/elements.h"

/// Within this file, Hashiguchi (2019) refers to
///  Hashiguchi, K. Multiplicative Hyperelastic-Based Plasticity for Finite
///  Elastoplastic Deformation/Sliding: A Comprehensive Review. Arch Computat
///  Methods Eng 26, 597–637 (2019). https://doi.org/10.1007/s11831-018-9256-5
namespace oomph
{

  // ===========================================================================
  /// \short A purely virtual class for isotropic hardening laws to be used /
  /// stored by the class PlasticConstitutiveLaw.
  // ===========================================================================
  class IsotropicHardeningLaw
  {
  public:
    // =========================================================================
    /// \short Constructor: Initializes the isotropic hardening law.
    ///
    /// \param isotropic_hardening_factor_pt Pointer to the parameter
    /// that describes the relation between the derivatives of
    /// isotropic hardening and the plastic multiplier (lambda).
    /// For mathematical details, see Eq. (152) of Hashiguchi (2019).
    // =========================================================================
    IsotropicHardeningLaw(double* isotropic_hardening_factor_pt)
      : Isotropic_hardening_factor_pt(isotropic_hardening_factor_pt)
    {
    }

    // =========================================================================
    /// \short A virtual destructor to prevent potential memory leaks.
    // =========================================================================
    virtual ~IsotropicHardeningLaw() = default;

    // =========================================================================
    /// \short Computes the (isotropic) yield function.
    /// This is a convenience wrapper that calls the purely virtual
    /// version of yield_function, disabling the derivative computation.
    /// This function is included to enable a simpler function call
    /// when derivatives are not needed.
    ///
    /// \param h The isotropic hardening value.
    /// \returns The isotropic yield surface stress.
    // =========================================================================
    virtual double yield_function(const double& h)
    {
      double a;
      bool compute = false;
      return yield_function(h, a, compute);
    }

    // =========================================================================
    /// \short Computes the (isotropic) yield function and its derivative.
    /// This purely virtual function must be overridden by derived classes
    /// to define the specific yield surface stress and, if requested,
    /// its derivative with respect to the hardening parameter.
    ///
    /// \param[in] h The isotropic hardening value.
    /// \param[out] dfdh The derivative of the yield surface stress w.r.t. h.
    /// \param[in] compute_derivative A flag to indicate whether dfdh
    /// should be computed.
    /// \returns The isotropic yield surface stress.
    // =========================================================================
    virtual double yield_function(const double& h,
                                  double& dfdh,
                                  const bool& compute_derivative) = 0;

    // =========================================================================
    /// \returns The pointer to the isotropic_hardening_factor.
    // =========================================================================
    double*& isotropic_hardening_factor_pt()
    {
      return Isotropic_hardening_factor_pt;
    }

    // =========================================================================
    /// \returns The isotropic_hardening_factor.
    // =========================================================================
    double isotropic_hardening_factor()
    {
      return *Isotropic_hardening_factor_pt;
    }

  private:
    /// A pointer to the isotropic hardening factor. It describes the
    /// relation between the derivatives of isotropic hardening and the plastic
    /// multiplier (lambda). For mathematical details, see Eq. (152) of
    /// Hashiguchi (2019).
    double* Isotropic_hardening_factor_pt;
  };

  // ===========================================================================
  /// \short An exponential isotropic hardening class, where the hardening is
  /// described by Eq. (151) of Hashiguchi (2019):
  /// \f[ f(H) = f_0 \left( 1 + h_1 \left[ 1 - \exp (-h_2 H) \right] \right) \f]
  // ===========================================================================
  class ExponentialIsotropicHardeningLaw : public IsotropicHardeningLaw
  {
  public:
    // =========================================================================
    /// \short Constructor: Initializes the exponential isotropic hardening law.
    ///
    /// \param isotropic_hardening_factor_pt \see IsotropicHardeningLaw
    /// \param f0_pt The initial isotropic yield stress.
    /// \param h1_pt A constant describing the evolution of the isotropic stress
    /// with the isotropic hardening value.
    /// \param h2_pt A constant describing the evolution of the isotropic stress
    /// with the isotropic hardening value.
    // =========================================================================
    ExponentialIsotropicHardeningLaw(double* isotropic_hardening_factor_pt,
                                     double* f0_pt,
                                     double* h1_pt,
                                     double* h2_pt)
      : IsotropicHardeningLaw(isotropic_hardening_factor_pt),
        F0_pt(f0_pt),
        H1_pt(h1_pt),
        H2_pt(h2_pt)
    {
    }

    // =========================================================================
    /// \short Implements the exponential hardening law.
    /// For a description of the equation \see ExponentialIsotropicHardeningLaw.
    /// For a description of the arguments \see
    /// IsotropicHardeningLaw::yield_function
    // =========================================================================
    double yield_function(const double& h,
                          double& dfdh,
                          const bool& compute_derivative) override
    {
      const double exp_val = std::exp(-(*H2_pt) * h);
      if (compute_derivative)
      {
        dfdh = (*F0_pt) * (*H1_pt) * (*H2_pt) * exp_val;
      }

      return (*F0_pt) * (1 + (*H1_pt) * (1 - exp_val));
    }

    // =========================================================================
    /// \short Access function to the initial isotropic yield stress pointer.
    // =========================================================================
    double* f0_pt() const
    {
      return F0_pt;
    }

    // =========================================================================
    /// \short Access function to the h1 parameter pointer.
    // =========================================================================
    double* h1_pt() const
    {
      return H1_pt;
    }

    // =========================================================================
    /// \short Access function to the h2 parameter pointer.
    // =========================================================================
    double* h2_pt() const
    {
      return H2_pt;
    }

  private:
    /// A pointer to the initial isotropic yield stress.
    double* F0_pt;

    /// Pointer to a value describing the evolution of the isotropic stress
    /// with the isotropic hardening value.
    double* H1_pt;

    /// Pointer to a value describing the evolution of the isotropic stress
    /// with the isotropic hardening value.
    double* H2_pt;
  };

  // ===========================================================================
  /// \short A purely virtual class for the yield criterion.
  // ===========================================================================
  class YieldCriterion
  {
  public:
    // =========================================================================
    /// \short A virtual destructor to prevent potential memory leaks.
    // =========================================================================
    virtual ~YieldCriterion() = default;

    // =========================================================================
    /// \short A purely virtual function to compute the scalar yield function
    /// value and its derivative w.r.t. M (if requested).
    ///
    /// \param[in] M The Mandel stress.
    /// \param[out] dsigmadM The derivative of the yield function value.
    /// \param[in] compute_derivative Whether to compute dsigmadM.
    /// \returns The yield function value.
    // =========================================================================
    virtual double surface_function(const DenseMatrix<double>& M,
                                    DenseMatrix<double>& dsigmadM,
                                    const bool& compute_derivative) = 0;

    // =========================================================================
    /// \short A purely virtual function to compute the second derivative of the
    /// scalar yield function value w.r.t. M.
    ///
    /// \param[in] f The yield function value.
    /// \param[in] dsigmadM The derivative of the yield function value.
    /// \param[out] ddsigmadMdM The second derivative of the yield function
    /// value.
    // =========================================================================
    virtual void surface_function_second_derivative(
      const double& f,
      const DenseMatrix<double>& dsigmadM,
      RankFourTensor<double>& ddsigmadMdM) = 0;
  };


  // ===========================================================================
  /// \short A class for the von Mises yield criterion (see Eq. 150 of
  /// Hashiguchi (2019)):
  /// \f[ f(M) = \sqrt{\frac{3.0}{2.0}} ||M^\prime|| \f]
  // ===========================================================================
  class VonMisesYieldCriterion : public YieldCriterion
  {
  public:
    // =========================================================================
    /// \short Empty constructor: This class does not have any internal
    /// variables.
    // =========================================================================
    VonMisesYieldCriterion() : YieldCriterion() {}

    // =========================================================================
    /// \short Computes the von Mises yield function. For a description of the
    /// arguments \see YieldCriterion::surface_function.
    // =========================================================================
    double surface_function(const DenseMatrix<double>& M,
                            DenseMatrix<double>& dsigmadM,
                            const bool& compute_derivative) override
    {
      // Calculate the trace of M
      double trace = 0.0;
      for (unsigned i = 0; i < M.nrow(); i++)
      {
        trace += M(i, i);
      }
      double mean_stress = trace / 3.0;

      // Get the magnitude of the deviatoric part
      double dev_mag_sq = 0.0;
      for (unsigned i = 0; i < M.nrow(); i++)
      {
        for (unsigned j = 0; j < M.ncol(); j++)
        {
          double val = M(i, j);
          if (i == j) val -= mean_stress; // Subtract mean from diagonal
          dev_mag_sq += val * val;
        }
      }

      // The function value
      double f = std::sqrt(3.0 / 2.0 * dev_mag_sq);

      if (compute_derivative)
      {
        dsigmadM.resize(M.nrow(), M.ncol(), 0.0);
        if (f > 1e-15)
        {
          // The derivative reads as
          // df/dM_ij = 3./2. S_ij / f, where
          // S_ij = M_ij - 1./3. M_kk \delta_ij
          // Note: S_kk = 0
          double scale_factor = 3.0 / 2.0 / f;
          for (unsigned i = 0; i < M.nrow(); i++)
          {
            for (unsigned j = 0; j < M.ncol(); j++)
            {
              double val = M(i, j);
              if (i == j) val -= mean_stress;
              dsigmadM(i, j) = val * scale_factor;
            }
          }
        }
        else
        {
          dsigmadM.initialise(0.0);
        }
      }

      // Return von Mises equivalent Stress
      return f;
    }

    // =========================================================================
    /// \short Computes the second derivative of the von Mises yield function.
    /// For a description of the arguments \see
    /// YieldCriterion::surface_function_second_derivative.
    // =========================================================================
    void surface_function_second_derivative(
      const double& f,
      const DenseMatrix<double>& dsigmadM,
      RankFourTensor<double>& ddsigmadMdM) override
    {
      const unsigned int nrow = dsigmadM.nrow();
      const unsigned int ncol = dsigmadM.ncol();
      ddsigmadMdM.resize(nrow, ncol, nrow, ncol);

      // If f is too small, e.g., 0, the derivative would explode, so we skip
      // it.
      if (f < 1e-15)
      {
        ddsigmadMdM.initialise(0.0);
        return;
      }

      // The hessian is
      // d^2f/(dM_ij dM_kl) = 1 / f *
      //           [ 3 / 2 (\delta_ik \delta_jl - 1 / 3 \delta_ij \delta_kl)
      //             - df/dM_kl * df/dM_ij]
      double invf = 1.0 / f;

      // We will split this in three parts to avoid the if statements
      // The first part - df/dM_kl * df/dM_ij / f
      for (unsigned i = 0; i < nrow; i++)
      {
        for (unsigned j = 0; j < ncol; j++)
        {
          double dSigmavmdM_ij_invf = -dsigmadM(i, j) * invf;
          for (unsigned k = 0; k < nrow; k++)
          {
            for (unsigned l = 0; l < ncol; l++)
            {
              ddsigmadMdM(i, j, k, l) = dsigmadM(k, l) * dSigmavmdM_ij_invf;
            }
          }
        }
      }

      // Now the first delta part
      // 3 / 2 (\delta_ik \delta_jl) / f
      double delta_ik_delta_jl_factor = 1.5 * invf;
      for (unsigned i = 0; i < nrow; i++)
      {
        for (unsigned j = 0; j < ncol; j++)
        {
          ddsigmadMdM(i, j, i, j) += delta_ik_delta_jl_factor;
        }
      }

      // And the second delta part
      // - \delta_ij \delta_kl / 2 / f
      double delta_ij_delta_kl_factor = invf / 2;
      const unsigned int min_dim = std::min(nrow, ncol);
      for (unsigned i = 0; i < min_dim; i++)
      {
        for (unsigned j = 0; j < min_dim; j++)
        {
          ddsigmadMdM(i, i, j, j) -= delta_ij_delta_kl_factor;
        }
      }
    }
  };

  // ===========================================================================
  /// \short A class which describes the behaviour of the normal-yield ratio
  /// \f$ R \f$ defined in Eq. (65) of Hashiguchi (2019).
  ///
  /// This class implements the whole evolution as described in sections 7.1
  /// (Eqs. (65) to (69)) and 8.6 (Eqs. (111) and (112)) of Hashiguchi (2019).
  // ===========================================================================
  class NormalYieldRatioLaw
  {
  public:
    // =========================================================================
    /// \short Constructor: Initializes the normal-yield ratio law.
    ///
    /// \param[in] normal_yield_ratio_elastic_pt The elastic/initial value
    /// of the normal-yield ratio.
    /// \param[in] u_pt An evolution constant. This is \f$ u \f$ from Eq. (68)
    /// or \f$ \bar{u} \f$ from Eq. (111).
    /// \param[in] elastic_core_u_pt Another evolution constant. See Eq. (111).
    /// \param[in] regularization_constant_pt A regularization constant, which
    /// is used to prevent division by zero.
    // =========================================================================
    NormalYieldRatioLaw(double* normal_yield_ratio_elastic_pt,
                        double* u_pt,
                        double* elastic_core_u_pt,
                        double* regularization_constant_pt)
      : Normal_yield_ratio_elastic_pt(normal_yield_ratio_elastic_pt),
        U_pt(u_pt),
        Elastic_core_u_pt(elastic_core_u_pt),
        Regularization_constant_pt(regularization_constant_pt)
    {
    }

    // =========================================================================
    /// \short A convenience constructor.
    ///
    /// For a description of the parameters see the called constructor.
    // =========================================================================
    NormalYieldRatioLaw(double* normal_yield_ratio_elastic_pt,
                        double* u_pt,
                        double* elastic_core_u_pt)
      : NormalYieldRatioLaw(normal_yield_ratio_elastic_pt,
                            u_pt,
                            elastic_core_u_pt,
                            &FiniteElement::Tolerance_for_singular_jacobian)
    {
    }

    // =========================================================================
    /// \short A convenience function to compute r (Eq. (68)) in the plastic
    /// case without computing any derivative.
    ///
    /// \param[in] u The value of the evolution constant.
    /// \param[in] delta_lambda The increase of the plastic multiplier since a
    /// reference step.
    /// \param[in] r_prev The value of r at the reference step.
    /// \returns The new value for r.
    // =========================================================================
    double compute_r_plastic(const double& u,
                             double delta_lambda,
                             const double& r_prev)
    {
      double derivative = 0;
      const bool compute_derivative = 0;
      return compute_r_plastic(
        u, delta_lambda, r_prev, derivative, derivative, compute_derivative);
    }

    // =========================================================================
    /// \short A function to compute r (Eq. (68)) in the plastic case
    ///
    /// \param[in] u The value of the evolution constant.
    /// \param[in] delta_lambda The increase of the plastic multiplier since a
    /// reference step.
    /// \param[in] r_prev The value of r at the reference step.
    /// \param[out] drdlambda The derivative of r w.r.t. lambda.
    /// \param[out] drdu The derivative of r w.r.t. u.
    /// \param[in] compute_derivative Whether the derivatives should be
    /// computed.
    /// \returns The new value for r.
    // =========================================================================
    double compute_r_plastic(const double& u,
                             double delta_lambda,
                             const double& R_prev,
                             double& drdlambda,
                             double& drdu,
                             const bool& compute_derivative)
    {
      // Initialize derivatives
      if (compute_derivative)
      {
        drdlambda = 0.0;
        drdu = 0.0;
      }

      // Retreive the minimum value of R
      const double Re = (*Normal_yield_ratio_elastic_pt);

      // Early exit if Re is 1.0 (the maximum value allowed for R)
      if (std::abs(1.0 - Re) < *Regularization_constant_pt) return 1.0;

      // Precompute some constants
      const double OneMinusRe = 1.0 - Re;
      const double preFactor = (2.0 * OneMinusRe) / MathematicalConstants::Pi;
      const double preFactorArg =
        MathematicalConstants::Pi / (2.0 * OneMinusRe);


      // Compute the argument of cos^{-1}
      // We do not need a max around R_prev - Re, since R_prev >= Re by
      // definition
      double inner_arg = std::cos(preFactorArg * (R_prev - Re)) *
                         std::exp(-u * preFactorArg * (delta_lambda));

      // Hard clamp strictly for the safety of std::acos to prevent NaNs
      // from standard floating-point overshoot.
      if (inner_arg > 1.0) inner_arg = 1.0;
      if (inner_arg < -1.0) inner_arg = -1.0;

      // Compute R - Eq. (68)
      double R = preFactor * std::acos(inner_arg) + Re;

      // Calculate dR/dLambda and dR/du for Jacobian analytically
      if (compute_derivative)
      {
        // d/dx(acos(f(x))) = -1 / sqrt(1 - f(x)^2) * df/dx
        double denom_sq = 1.0 - inner_arg * inner_arg;
        double inv_sqrt =
          -1.0 / std::sqrt(denom_sq + (*Regularization_constant_pt));

        // df/dLambda = f * (-u * preFactorArg)
        // dRdLambda = preFactor * inv_sqrt * df/dLambda
        // preFactor * preFactorArg = 1
        // Hence: dR/dLambda = - u * f * inv_sqrt
        drdlambda = -u * inner_arg * inv_sqrt;

        // df/du = f * (-preFactorArg * delta_lambda)
        // dR/du = preFactor * inv_sqrt * dfdu
        //       = - f * inv_sqrt * delta_lambda
        drdu = -inner_arg * inv_sqrt * delta_lambda;
      }

      return R;
    }

    // =========================================================================
    /// \short A function to compute the evolution variable u using Eqs.
    /// (111) and (112).
    ///
    /// \param[in] rc The value of rc (Eq. (105)).
    /// \param[in] barbar_N The value of bar_bar_N (Eq. (107)).
    /// \param[in] hat_bar_Nc The value of hat_bar_Nc (Eq. (110)).
    /// The following input arguments will only be used if compute_derivative is
    /// true.
    /// \param[in] dbarbar_N_dbarbar_M The derivative of bar_bar_N w.r.t.
    /// bar_bar_M
    /// \param[in] dbarbar_M_dbar_Mk The derivative of bar_bar_M w.r.t. bar_Mk.
    /// \param[in] dbarbar_M_dbar_Mc The derivative of bar_bar_M w.r.t. bar_Mc.
    /// \param[in] dc_sigma_dbar_Mk The derivative of c_sigma w.r.t. bar_Mk.
    /// \param[in] dc_sigma_dbar_Mc The derivative of c_sigma w.r.t. bar_Mc.
    /// \param[in] dbarbar_M_dr The derivative of bar_bar_M w.r.t. R.
    /// \param[in] dhatbar_Nc_dhatbar_Mc The derivative of bar_bar_M w.r.t.
    /// \param[in] drc_dhatbar_Mc The derivative of r_c w.r.t. hat_bar_Mc
    /// \param[in] drc_dh The derivative of rc w.r.t. h.
    /// The following output arguments will only be resized and initialised, if
    /// compute_derivative is true.
    /// \param[out] du_dbar_M The derivative of u w.r.t. bar_M.
    /// \param[out] du_dbar_Mk The derivative of u w.r.t. bar_Mk.
    /// \param[out] du_dbar_Mc The derivative of u w.r.t. bar_Mc.
    /// \param[out] du_dh The derivative of u w.r.t. h.
    /// \param[out] du_dr The derivative of u w.r.t. r.
    /// \param[in] compute_derivative Whether the derivatives should be
    /// computed.
    /// \returns The new value for u.
    // =========================================================================
    double compute_u_with_elastic_core(
      const double& rc,
      const DenseMatrix<double>& barbar_N,
      const DenseMatrix<double>& hat_bar_Nc,
      const RankFourTensor<double>& dbarbar_N_dbarbar_M,
      const double& dbarbar_M_dbar_Mk,
      const double& dbarbar_M_dbar_Mc,
      const DenseMatrix<double>& dbarbar_M_dr,
      const RankFourTensor<double>& dhatbar_Nc_dhatbar_Mc,
      const DenseMatrix<double>& drc_dhatbar_Mc,
      const double& drc_dh,
      DenseMatrix<double>& du_dbar_M,
      DenseMatrix<double>& du_dbar_Mk,
      DenseMatrix<double>& du_dbar_Mc,
      double& du_dh,
      double& du_dr,
      bool compute_derivative

    )
    {
      DenseMatrix<double> dc_dbar_M, dc_dbar_Mk, dc_dbar_Mc;
      double dcdr = 0;

      double c_sigma = compute_c_sigma(barbar_N,
                                       hat_bar_Nc,
                                       dbarbar_N_dbarbar_M,
                                       dbarbar_M_dbar_Mk,
                                       dbarbar_M_dbar_Mc,
                                       dbarbar_M_dr,
                                       dhatbar_Nc_dhatbar_Mc,
                                       dc_dbar_M,
                                       dc_dbar_Mk,
                                       dc_dbar_Mc,
                                       dcdr,
                                       compute_derivative);

      return compute_u_with_elastic_core_from_c(rc,
                                                c_sigma,
                                                drc_dhatbar_Mc,
                                                drc_dh,
                                                dc_dbar_M,
                                                dc_dbar_Mk,
                                                dc_dbar_Mc,
                                                dcdr,
                                                du_dbar_M,
                                                du_dbar_Mk,
                                                du_dbar_Mc,
                                                du_dh,
                                                du_dr,
                                                compute_derivative);
    }

    // =========================================================================
    /// \returns The elastic/initial value of the normal-yield ratio.
    // =========================================================================
    double get_re()
    {
      return *Normal_yield_ratio_elastic_pt;
    }

    // =========================================================================
    /// \returns The evolution parameter u.
    // =========================================================================
    double get_u()
    {
      return *U_pt;
    }

    // =========================================================================
    /// \short Access to the elastic/initial value of the normal-yield ratio Re,
    /// as in Eq. (68).
    // =========================================================================
    double*& normal_yield_ratio_elastic_pt()
    {
      return Normal_yield_ratio_elastic_pt;
    }

    // =========================================================================
    /// \short Access to the evolution parameter u parameter. (bar{u} in Eq.
    /// (111); u in Eq. (68), if there was no elastic core)
    // =========================================================================
    double*& u_pt()
    {
      return U_pt;
    }

    // =========================================================================
    /// \short Access to the elastic core u parameter. (u_c in Eq. (111))
    // =========================================================================
    double*& elastic_core_u_pt()
    {
      return Elastic_core_u_pt;
    }

    // =========================================================================
    /// \short Access to the regularization constant.
    // =========================================================================
    double*& regularization_constant_pt()
    {
      return Regularization_constant_pt;
    }

  protected:
    // =========================================================================
    /// \short A function to compute the evolution variable u from c using Eq.
    /// (111):
    /// \f[ u = \bar{u} \mathrm{exp}
    ///           \left(u_\mathrm{c} \mathcal{R}_\mathrm{c0} C_\sigma\right) \f]
    ///
    /// \param[in] rc The value of \f$ R_\mathrm{c} \f$ (Eq. 105).
    /// \param[in] c_sigma The value of \f$ C_\sigma \f$ (Eq. 112).
    /// \param[in] r_prev The value of r at the reference step.
    /// The following input arguments will only be used if compute_derivative is
    /// true.
    /// \param[in] drc_dhatbar_Mc The derivative of r_c w.r.t. hat_bar_Mc
    /// \param[in] drc_dh The derivative of rc w.r.t. h.
    /// \param[in] dc_sigma_dbar_M The derivative of c_sigma w.r.t. bar_M.
    /// \param[in] dc_sigma_dbar_Mk The derivative of c_sigma w.r.t. bar_Mk.
    /// \param[in] dc_sigma_dbar_Mc The derivative of c_sigma w.r.t. bar_Mc.
    /// \param[in] dc_sigma_dr The derivative of c_sigma w.r.t. r.
    /// The following output arguments will only be resized and initialised, if
    /// compute_derivative is true.
    /// \param[out] du_dbar_M The derivative of u w.r.t. bar_M.
    /// \param[out] du_dbar_Mk The derivative of u w.r.t. bar_Mk.
    /// \param[out] du_dbar_Mc The derivative of u w.r.t. bar_Mc.
    /// \param[out] du_dh The derivative of u w.r.t. h.
    /// \param[out] du_dr The derivative of u w.r.t. r.
    /// \param[in] compute_derivative Whether the derivatives should be
    /// computed.
    /// \returns The new value for u.
    // =========================================================================
    double compute_u_with_elastic_core_from_c(
      const double& rc,
      const double& c_sigma,
      const DenseMatrix<double>& drc_dhatbar_Mc,
      const double& drc_dh,
      const DenseMatrix<double>& dc_sigma_dbar_M,
      const DenseMatrix<double>& dc_sigma_dbar_Mk,
      const DenseMatrix<double>& dc_sigma_dbar_Mc,
      const double& dc_sigma_dr,
      DenseMatrix<double>& du_dbar_M,
      DenseMatrix<double>& du_dbar_Mk,
      DenseMatrix<double>& du_dbar_Mc,
      double& du_dh,
      double& du_dr,
      bool compute_derivative)
    {
      // Compute Eq. 111
      const double uc = (*Elastic_core_u_pt);
      const double u_out = (*U_pt) * std::exp(uc * rc * c_sigma);

      // Compute derivative of Eq. 111 with simple chain rules
      if (compute_derivative)
      {
        const double prefactor_dc_sigma = uc * u_out * rc;
        const double prefactor_drc = uc * u_out * c_sigma;

        const unsigned DIM = drc_dhatbar_Mc.ncol();
        du_dbar_M.resize(DIM, DIM);
        du_dbar_Mk.resize(DIM, DIM);
        du_dbar_Mc.resize(DIM, DIM);

        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            du_dbar_M(i, j) = prefactor_dc_sigma * dc_sigma_dbar_M(i, j);
            du_dbar_Mk(i, j) = prefactor_dc_sigma * dc_sigma_dbar_Mk(i, j) -
                               prefactor_drc * drc_dhatbar_Mc(i, j);
            du_dbar_Mc(i, j) = prefactor_dc_sigma * dc_sigma_dbar_Mc(i, j) +
                               prefactor_drc * drc_dhatbar_Mc(i, j);
          }
        }

        du_dh = prefactor_drc * drc_dh;
        du_dr = prefactor_dc_sigma * dc_sigma_dr;
      }

      return u_out;
    }

    // =========================================================================
    /// \short A function to compute the c_sigma using Eq. (112):
    /// \f[ C_\sigma = \bar{\bar{N}} \colon \hat{\bar{N}}_\mathrm{c} \f]
    ///
    /// \param[in] bar_bar_N The value of bar_bar_N (Eq. 107).
    /// \param[in] hat_bar_Nc The value of hat_bar_Nc (Eq. 110).
    /// The following input arguments will only be used if compute_derivative is
    /// true.
    /// \param[in] dbarbar_N_dbarbar_M The derivative of bar_bar_N w.r.t.
    /// bar_bar_M
    /// \param[in] dbarbar_M_dbar_Mk The derivative of bar_bar_M w.r.t. bar_Mk.
    /// \param[in] dbarbar_M_dbar_Mc The derivative of bar_bar_M w.r.t. bar_Mc.
    /// \param[in] dbarbar_M_dr The derivative of bar_bar_M w.r.t. R.
    /// \param[in] dhatbar_Nc_dhatbar_Mc The derivative of bar_bar_M w.r.t.
    /// hat_bar_Mc.
    /// The following output arguments will only be resized and initialised, if
    /// compute_derivative is true.
    /// \param[out] dc_dbar_M The derivative of c_sigma w.r.t. bar_M.
    /// \param[out] dc_dbar_Mk The derivative of c_sigma w.r.t. bar_Mk.
    /// \param[out] dc_dbar_Mc The derivative of c_sigma w.r.t. bar_Mc.
    /// \param[out] dcdr The derivative of c_sigma w.r.t. h.
    /// \param[in] compute_derivative Whether the derivatives should be
    /// computed.
    /// \returns The new value of c_sigma.
    // =========================================================================
    double compute_c_sigma(const DenseMatrix<double>& bar_bar_N,
                           const DenseMatrix<double>& hat_bar_Nc,
                           const RankFourTensor<double>& dbarbar_N_dbarbar_M,
                           const double& dbarbar_M_dbar_Mk,
                           const double& dbarbar_M_dbar_Mc,
                           const DenseMatrix<double>& dbarbar_M_dr,
                           const RankFourTensor<double>& dhatbar_Nc_dhatbar_Mc,
                           DenseMatrix<double>& dc_dbar_M,
                           DenseMatrix<double>& dc_dbar_Mk,
                           DenseMatrix<double>& dc_dbar_Mc,
                           double& dcdr,
                           bool compute_derivative)
    {
      // The reduction is c_sigma = bar_bar_N(i, j) hat_bar_N_c(i, j). The
      // derivative can be obtained by applying product and chain rules.
      if (compute_derivative)
      {
        const unsigned DIM = bar_bar_N.ncol();
        dc_dbar_M.resize(DIM, DIM);
        dc_dbar_Mk.resize(DIM, DIM);
        dc_dbar_Mc.resize(DIM, DIM);
        dcdr = 0.0;

        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            double term1_sum = 0.0;
            double term2_sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              for (unsigned int b = 0; b < DIM; b++)
              {
                term1_sum += dbarbar_N_dbarbar_M(a, b, i, j) * hat_bar_Nc(a, b);
                term2_sum +=
                  bar_bar_N(a, b) * dhatbar_Nc_dhatbar_Mc(a, b, i, j);
              }
            }

            dcdr += term1_sum * dbarbar_M_dr(i, j);

            dc_dbar_M(i, j) = term1_sum;
            dc_dbar_Mk(i, j) = term1_sum * dbarbar_M_dbar_Mk - term2_sum;
            dc_dbar_Mc(i, j) = term1_sum * dbarbar_M_dbar_Mc + term2_sum;
          }
        }
      }

      return MatrixHelpers::reduce(bar_bar_N, hat_bar_Nc);
    }

  private:
    /// Pointer to the elastic/initial value of the normal-yield ratio Re, as in
    /// Eq. (68).
    double* Normal_yield_ratio_elastic_pt;

    /// Pointer the evolution parameter u parameter. (bar{u} in Eq. (111); u in
    /// Eq. (68), if there was no elastic core).
    double* U_pt;

    /// Pointer to the elastic core u parameter. (u_c in Eq. (111)).
    double* Elastic_core_u_pt;

    /// Pointer to a regularization constant regularizing the division by 0.
    double* Regularization_constant_pt;
  };

  // ===========================================================================
  /// \short A container class to store the pointers to laws for the individual
  /// plasticity parameters as well as to some base parameters. An instance of
  /// this class is required by \ref PlasticEquationsBase.
  // ===========================================================================
  class PlasticConstitutiveLaw
  {
  public:
    /// The isotropic hardening law, \see IsotropicHardeningLaw.
    IsotropicHardeningLaw* isotropic_hardening_law_pt;

    /// The yield criterion, \see YieldCriterion.
    YieldCriterion* yield_criterion_pt;

    /// The kinematic hardening law describing the relationship between the
    /// kinematic hardening deformation gradient tensor and the associated
    /// stress. See, e.g., Eq. (33).
    ConstitutiveLaw* kinematic_hardening_law_pt;

    /// The elastic core law describing the relationship between the elastic
    /// core deformation gradient tensor and the associated stress. See, e.g.,
    /// Eqs. (96) and (97).
    ConstitutiveLaw* elastic_core_law_pt;

    /// The normal yield ratio law, \see NormalYieldRatioLaw.
    NormalYieldRatioLaw* normal_yield_ratio_law_pt;

    /// A pointer to eta_p, see. Eq. (114).
    double* eta_p_pt;

    /// A pointer to eta_pk, see. Eq. (114).
    double* kinematic_hardening_eta_pt = nullptr;

    /// A pointer to b_k, see. Eq. (114).
    double* kinematic_hardening_b_pt = nullptr;

    /// A pointer to eta_pc, see. Eq. (114).
    double* elastic_core_eta_pt = nullptr;

    /// A pointer to X, see. Eq. (114).
    double* elastic_core_x_pt = nullptr;
  };
} // namespace oomph
#endif