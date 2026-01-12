#ifndef OOMPH_PLASTICITY_CONSTITUTIVE_HEADER
#define OOMPH_PLASTICITY_CONSTITUTIVE_HEADER

#include "constitutive_laws.h"
#include "generic/matrix_helpers.h"

namespace oomph
{
  class PlasticConstitutiveLaw
  {
  public:
    /*!
     * Computes the von mises stress and its derivative wrt. M
     */
    double compute_yield_surface_function(const DenseMatrix<double>& M,
                                          DenseMatrix<double>& dSigmavmdM,
                                          bool computeDerivative)
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

      if (computeDerivative)
      {
        dSigmavmdM.resize(M.nrow(), M.ncol(), 0.0);
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
              dSigmavmdM(i, j) = val * scale_factor;
            }
          }
        }
        else
        {
          dSigmavmdM.initialise(0.0);
        }
      }

      // Return Von Mises Equivalent Stress
      return f;
    }

    /*!
     * \brief computes the second derivative of the con mises yield surface
     * function
     * \param[in] f: The yield surface function value
     * \param[in] dSigmavmdM: Its first derivative
     * \param[out] ddSigmavmdMdM: Its second derivative
     */
    void compute_ddyield_surface_functiondMdM(
      const double& f,
      const DenseMatrix<double>& dSigmavmdM,
      RankFourTensor<double>& ddSigmavmdMdM)
    {
      const unsigned int nrow = dSigmavmdM.nrow();
      const unsigned int ncol = dSigmavmdM.ncol();
      ddSigmavmdMdM.resize(nrow, ncol, nrow, ncol);

      if (f < 1e-15)
      {
        ddSigmavmdMdM.initialise(0.0);
        return;
      }

      // The hessian is
      // d^2f/(dM_ij dM_kl) = 1 / f *
      //           [ 3 / 2 (\delta_ik \delta_jl - 1 / 3 \delta_ij \delta_kl)
      //             - df/dM_kl * df/dM_ij]
      double invf = 1 / f;

      // We will split this in three parts to avoid the if statements
      // The first part - df/dM_kl * df/dM_ij / f
      for (unsigned i = 0; i < nrow; i++)
      {
        for (unsigned j = 0; j < ncol; j++)
        {
          double dSigmavmdM_ij_invf = -dSigmavmdM(i, j) * invf;
          for (unsigned k = 0; k < nrow; k++)
          {
            for (unsigned l = 0; l < ncol; l++)
            {
              ddSigmavmdMdM(i, j, k, l) = dSigmavmdM(k, l) * dSigmavmdM_ij_invf;
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
          ddSigmavmdMdM(i, j, i, j) += delta_ik_delta_jl_factor;
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
          ddSigmavmdMdM(i, i, j, j) -= delta_ij_delta_kl_factor;
        }
      }
    }

    double isotropic_hardening_yield_function(const double& H)
    {
      double a;
      return isotropic_hardening_yield_function(H, a, false);
    }

    double isotropic_hardening_yield_function(const double& H,
                                              double& dfdH,
                                              bool computeDerivative)
    {
      const double exp_val = std::exp(-isotropic_hardening_h2 * H);
      if (computeDerivative)
      {
        dfdH = isotropic_hardening_f0 * isotropic_hardening_h1 *
               isotropic_hardening_h2 * exp_val;
      }

      return isotropic_hardening_f0 *
             (1 + isotropic_hardening_h1 * (1 - exp_val));
    }

    void mandel_like_kinematic_hardening_variable(
      const DenseMatrix<double>& Fpks,
      const DenseMatrix<double>& intermediateMetric,
      DenseMatrix<double>& Mbark,
      RankFourTensor<double>& dbar_Mk_dFpks,
      bool computeDerivative)
    {
      const unsigned int nrow = Fpks.nrow();
      const unsigned int ncol = Fpks.ncol();
      Mbark.resize(nrow, nrow);

      for (unsigned int i = 0; i < nrow; i++)
      {
        for (unsigned int j = 0; j < nrow; j++)
        {
          double FpksFpksT_ij = 0;
          for (unsigned int m = 0; m < ncol; m++)
          {
            FpksFpksT_ij += Fpks(i, m) * Fpks(j, m);
          }
          Mbark(i, j) = kinematic_hardening_stress_c *
                        (FpksFpksT_ij - intermediateMetric(i, j));
        }
      }

      if (!computeDerivative) return;

      dbar_Mk_dFpks.resize(nrow, nrow, nrow, ncol, 0.0);
      dbar_Mk_dFpks.initialise(0.0);

      // dM_ik / dF_kl = c * (F_jl \delta_ik + F_il \delta_jk)
      for (unsigned int i = 0; i < nrow; i++)
      {
        for (unsigned int j = 0; j < nrow; j++)
        {
          for (unsigned int l = 0; l < ncol; l++)
          {
            // Contribution i == k
            dbar_Mk_dFpks(i, j, i, l) +=
              kinematic_hardening_stress_c * Fpks(j, l);

            // Contribution j == k
            dbar_Mk_dFpks(i, j, j, l) +=
              kinematic_hardening_stress_c * Fpks(i, l);
          }
        }
      }
    }

    void mandel_like_elastic_core_variable(
      const DenseMatrix<double>& Fpcs,
      const DenseMatrix<double>& intermediateMetric,
      DenseMatrix<double>& Mbarc,
      RankFourTensor<double>& dbar_Mk_dFpcs,
      bool computeDerivative)
    {
      const unsigned int nrow = Fpcs.nrow();
      const unsigned int ncol = Fpcs.ncol();
      Mbarc.resize(nrow, nrow);

      for (unsigned int i = 0; i < nrow; i++)
      {
        for (unsigned int j = 0; j < nrow; j++)
        {
          double FpcsFpcsT_ij = 0;
          for (unsigned int m = 0; m < ncol; m++)
          {
            FpcsFpcsT_ij += Fpcs(i, m) * Fpcs(j, m);
          }
          Mbarc(i, j) =
            elastic_core_stress_c * (FpcsFpcsT_ij - intermediateMetric(i, j));
        }
      }

      if (!computeDerivative) return;

      dbar_Mk_dFpcs.resize(nrow, nrow, nrow, ncol, 0.0);
      dbar_Mk_dFpcs.initialise(0.0);

      // dM_ik / dF_kl = c * (F_jl \delta_ik + F_il \delta_jk)
      for (unsigned int i = 0; i < nrow; i++)
      {
        for (unsigned int j = 0; j < nrow; j++)
        {
          for (unsigned int l = 0; l < ncol; l++)
          {
            // Contribution i == k
            dbar_Mk_dFpcs(i, j, i, l) += elastic_core_stress_c * Fpcs(j, l);

            // Contribution j == k
            dbar_Mk_dFpcs(i, j, j, l) += elastic_core_stress_c * Fpcs(i, l);
          }
        }
      }
    }

    /*!
     * \brief computes the yield ratio based on the equation f(Mbarbar') = R
     * F(H)
     */
    double compute_R_elastic(const DenseMatrix<double>& bar_M,
                             const DenseMatrix<double>& bar_Mk,
                             const DenseMatrix<double>& bar_Mc,
                             const double& H)
    {
      // Compute tilde_bar_M and hat_bar_Mc
      DenseMatrix<double> tilde_bar_M(bar_M.nrow(), bar_M.ncol(), 0.0),
        hat_bar_Mc(bar_M.nrow(), bar_M.ncol(), 0.0);

      // Variables to store the trace (sum of diagonals)
      double trace_tilde = 0.0;
      double trace_hat = 0.0;

      for (unsigned int i = 0; i < bar_M.nrow(); i++)
      {
        for (unsigned int j = 0; j < bar_M.ncol(); j++)
        {
          tilde_bar_M(i, j) = bar_M(i, j) - bar_Mc(i, j);
          hat_bar_Mc(i, j) = bar_Mc(i, j) - bar_Mk(i, j);

          // Accumulate trace
          if (i == j)
          {
            trace_tilde += tilde_bar_M(i, j);
            trace_hat += hat_bar_Mc(i, j);
          }
        }
      }

      // Divide by 3 here even for DIM=2, because stress is 3D
      double mean_tilde = trace_tilde / 3;
      double mean_hat = trace_hat / 3;

      double term_MM = 0.0; // tilde_bar_M : tilde_bar_M
      double term_MMc = 0.0; // tilde_bar_M : hat_bar_Mc
      double term_McMc = 0.0; // hat_bar_Mc  : hat_bar_Mc (needed for norm)

      for (unsigned int i = 0; i < bar_M.nrow(); i++)
      {
        for (unsigned int j = 0; j < bar_M.ncol(); j++)
        {
          if (i == j)
          {
            // Compute deviatoric part
            tilde_bar_M(i, j) -= mean_tilde;
            hat_bar_Mc(i, j) -= mean_hat;
          }

          double m_val = tilde_bar_M(i, j);
          double mc_val = hat_bar_Mc(i, j);

          // Compute the contractions
          term_MM += m_val * m_val;
          term_MMc += m_val * mc_val;
          term_McMc += mc_val * mc_val;
        }
      }

      // Finally compute R
      double Fh = isotropic_hardening_yield_function(H);

      double denominator = 2. / 3. * Fh * Fh - term_McMc;
      double enumerator =
        term_MMc + std::sqrt(term_MMc * term_MMc + denominator * term_MM);

      return enumerator / denominator;
    }


    double isotropic_hardening_f0 = 500e6;

    double isotropic_hardening_h1 = 0.8;

    double isotropic_hardening_h2 = 0.50;

    double eta_p = 0.0;

    /*!
     * \brief R^\text{e}
     */
    double normal_yield_ratio_elastic = 0.5;

    /*!
     * normal_yield_ratio_constant_u
     */
    double normal_yield_ratio_constant_u = 200;

    /*!
     * \brief f^\text{H} from the formulas
     */
    double isotropic_hardening_factor = std::sqrt(3.0 / 2.0);

    /*!
     * \brief \eta^\text{pk}
     */
    double kinematic_hardening_eta = 0.0;

    /*!
     * \brief b^\text{pk}
     */
    double kinematic_hardening_b = 2.0;

    /*!
     * \brief C^\text{k}
     */
    double kinematic_hardening_stress_c = 1.0;

    /*!
     * \brief \eta^\text{pc}
     */
    double elastic_core_eta = 0.0;

    /*!
     * \brief X
     */
    double elastic_core_x = 1.0;

    /*!
     * \brief C^\text{c}
     */
    double elastic_core_stress_c = 1.0;

    /*!
     * \brief u_c
     */
    double elasic_core_u = 1.0;
  };
} // namespace oomph
#endif