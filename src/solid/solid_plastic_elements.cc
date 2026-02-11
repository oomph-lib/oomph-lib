#include "solid/solid_plastic_elements.h"

using namespace oomph;

template<unsigned DIM>
const std::vector<std::string> PlasticEquations<DIM>::Plastic_data_names{
  "Fe", "Fpks", "Fpcs", "H", "Lambda", "R"};

template class PlasticEquations<2>;
template class QPlasticPVDElement<2, 2>;

template class PlasticEquations<3>;
template class QPlasticPVDElement<3, 2>;

/*
 * \details Computes the Caucystress using the plastic deformation gradient
 */
template<unsigned DIM>
void PlasticEquations<DIM>::get_cauchy_stress(unsigned ipt,
                                              DenseMatrix<double>& sigma)
{
  DenseMatrix<double> F(DIM);
  compute_deformation_gradient_tensor(ipt, F);

  // calculate J = det(F)
  double J = MatrixHelpers::determinant<DIM>(F);

  // calculate C = F^T * F
  DenseMatrix<double> C(DIM, DIM, 0.0);
  MatrixHelpers::transpose_multiply(F, F, C);

  // compute Cp = Fp^T * Fp
  DenseMatrix<double> invFp(DIM, DIM);
  get_inv_fp_matrix(ipt, invFp);

  DenseMatrix<double> Fp(DIM, DIM);
  MatrixHelpers::invert_matrix<DIM>(invFp, Fp);

  DenseMatrix<double> Cp(DIM, DIM, 0.0);
  MatrixHelpers::transpose_multiply(Fp, Fp, Cp);


  // compute second Piola–Kirchhoff stress tensor
  DenseMatrix<double> S_PK2(DIM, DIM, 0.0);

  // This calculates S = 2 * dPsi(C, Cp) / dC
  this->Constitutive_law_pt->calculate_second_piola_kirchhoff_stress(
    Cp, C, S_PK2);

  // compute cauchy stress
  // Standard formula: sigma = (1/J) * F * S_PK2 * F^T
  DenseMatrix<double> tmp(DIM, DIM, 0.0);
  DenseMatrix<double> Ft;
  MatrixHelpers::transpose(F, Ft);
  MatrixHelpers::multiply(S_PK2, Ft, tmp); // tmp = S * F^T

  MatrixHelpers::multiply(F, tmp, sigma); // sigma = F * tmp

  // Divide by Jacobian
  for (unsigned i = 0; i < DIM; i++)
  {
    for (unsigned j = 0; j < DIM; j++)
    {
      sigma(i, j) /= J;
    }
  }
}

/*!
 * \brief computes R as a dependent variable R(H)
 */
template<unsigned DIM>
double PlasticEquations<DIM>::compute_r_plastic(const double& u,
                                                const double& delta_lambda,
                                                const double& R_prev,
                                                double& dRdLambda,
                                                double& dRdu,
                                                bool computeDerivative)
{
  const double Re =
    this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic;

  // Initialize derivatives
  if (computeDerivative)
  {
    dRdLambda = 0.0;
    dRdu = 0.0;
  }

  if (std::abs(1.0 - Re) < 1.0e-12) return 1.0;

  const double OneMinusRe = 1.0 - Re;
  const double preFactor = (2.0 * OneMinusRe) / MathematicalConstants::Pi;
  const double preFactorArg = MathematicalConstants::Pi / (2.0 * OneMinusRe);

  // Compute the argument of cos^{-1}
  // We do not need a max around R_prev - Re, since R_prev >= Re by definition
  double inner_arg = std::cos(preFactorArg * (R_prev - Re)) *
                     std::exp(-u * preFactorArg * delta_lambda);

  // Limit the arg. This has two effects, it prevents the acos to take wierd
  // numbers and makes the derivative always not infinite.
  const double TOLERANCE = 1.0e-12;
  const double MAX_ARG = std::sqrt(1.0 - TOLERANCE);

  if (inner_arg > MAX_ARG) inner_arg = MAX_ARG;
  if (inner_arg < -MAX_ARG) inner_arg = -MAX_ARG;

  // Compute R
  double R = preFactor * std::acos(inner_arg) + Re;

  if (computeDerivative)
  {
    // Calculate dR/dLambda and dR/du for Jacobian
    double denom_sq = 1.0 - inner_arg * inner_arg;

    // d/dx(acos(f(x))) = -1 / sqrt(1 - f(x)^2) * df/dx
    double inv_sqrt = -1.0 / std::sqrt(denom_sq);

    // df/dLambda = f * (-u * preFactorArg)
    double dWdLambda = inner_arg * (-u * preFactorArg);
    dRdLambda = preFactor * inv_sqrt * dWdLambda;

    // df/du = f * (-preFactorArg * delta_lambda)
    double dWdu = inner_arg * (-preFactorArg * delta_lambda);
    dRdu = preFactor * inv_sqrt * dWdu;
  }

  return R;
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_mandel_stress_elastic(
  const DenseMatrix<double>& invFp,
  const DenseMatrix<double>& FtF,
  DenseMatrix<double>& bar_M,
  RankFourTensor<double>& dbar_M_dinvFp,
  bool computeDerivative)
{
  // Compute GF^{p, -1} and Ce
  DenseMatrix<double> FtF_invFp(DIM, DIM, 0.0), Ce(DIM, DIM, 0.0);
  MatrixHelpers::multiply(FtF, invFp, FtF_invFp);
  MatrixHelpers::transpose_multiply(invFp, FtF_invFp, Ce);

  // Now compute bar_M
  DenseMatrix<double> S(DIM); // Second Piola–Kirchhoff stress tensors
  this->Constitutive_law_pt->calculate_second_piola_kirchhoff_stress(
    this->unity, Ce, S);

  // The mandel stress
  MatrixHelpers::multiply(Ce, S, bar_M);

  if (!computeDerivative) return;

  // The derivative dSdC
  RankFourTensor<double> dSdC(DIM, DIM, DIM, DIM, 0.0);
  this->Constitutive_law_pt->calculate_d_second_piola_kirchhoff_stress_dG(
    this->unity, Ce, S, dSdC);

  // Now compude dMdC
  // dM_{ij} / dC_{kl} = delta_{ik} * S_{lj} + C_{im} * (dS_{mj} / dC_{kl})
  RankFourTensor<double> dMbardCe(DIM, DIM, DIM, DIM, 0.0);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      for (unsigned int k = 0; k < DIM; k++)
      {
        for (unsigned int l = 0; l < DIM; l++)
        {
          // Term1: delta_{ik} * S_{lj}
          double term1 = 0.0;
          if (i == k)
          {
            term1 = S(l, j);
          }

          // Term 2: C_{im} * dSdC_{mjkl}
          double term2 = 0.0;
          for (unsigned int m = 0; m < DIM; m++)
          {
            term2 += Ce(i, m) * dSdC(m, j, k, l);
          }

          // Sum contributions
          dMbardCe(i, j, k, l) = term1 + term2;
        }
      }
    }
  }

  // Finally dMdFpinv
  // dMdFpinv_{abcd} = (GF^{p, -1})_{ci} *
  //                                     (dM_{ab} / dC_{di} + dM_{ab} /
  //                                     dC_{id})
  dbar_M_dinvFp.resize(DIM, DIM, DIM, DIM, 0.0);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      for (unsigned int k = 0; k < DIM; k++)
      {
        for (unsigned int l = 0; l < DIM; l++)
        {
          dbar_M_dinvFp(i, j, k, l) = 0;
          for (unsigned int a = 0; a < DIM; a++)
          {
            dbar_M_dinvFp(i, j, k, l) +=
              FtF_invFp(k, a) * (dMbardCe(i, j, l, a) + dMbardCe(i, j, a, l));
          }
        }
      }
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_mandellike_kinematic_hardening(
  const DenseMatrix<double>& Fpks,
  DenseMatrix<double>& bar_Mk,
  RankFourTensor<double>& dbar_MkdFpks,
  bool computeDerivative)
{
  bar_Mk.resize(DIM, DIM, 0.0);
  // We now always build all plastic data. Do we need to check if we have pinned
  // the data instead?
  // if (!Plastic_data_has_been_built[Fpks_INDEX])
  // {
  //   bar_Mk.initialise(0.0);

  //   if (computeDerivative)
  //   {
  //     dbar_MkdFpks.resize(DIM, DIM, DIM, DIM, 0.0);
  //   }

  //   return;
  // }

  this->Plastic_consitutive_law_pt->mandel_like_kinematic_hardening_variable(
    Fpks, this->unity, bar_Mk, dbar_MkdFpks, computeDerivative);
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_mandellike_elastic_core(
  const DenseMatrix<double>& Fpcs,
  DenseMatrix<double>& bar_Mc,
  RankFourTensor<double>& dbar_McdFpcs,
  bool computeDerivative)
{
  bar_Mc.resize(DIM, DIM, 0.0);
  // We now always build all plastic data. Do we need to check if we have pinned
  // the data instead?
  // if (!Plastic_data_has_been_built[Fpcs_INDEX])
  // {
  //   bar_Mc.initialise(0.0);

  //   if (computeDerivative)
  //   {
  //     dbar_McdFpcs.resize(DIM, DIM, DIM, DIM, 0.0);
  //   }

  //   return;
  // }

  this->Plastic_consitutive_law_pt->mandel_like_elastic_core_variable(
    Fpcs, this->unity, bar_Mc, dbar_McdFpcs, computeDerivative);
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_mandel_stress_total(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& bar_Mk,
  const DenseMatrix<double>& bar_Mc,
  const double& R,
  DenseMatrix<double>& barbar_M,
  double& dbarbar_M_dMk,
  double& dbarbar_M_dMc,
  DenseMatrix<double>& dbarbar_M_dR,
  bool computeJacobian)
{
  barbar_M.resize(DIM, DIM);
  if (computeJacobian) dbarbar_M_dR.resize(DIM, DIM);
  // We now always build all plastic data. Do we need to check if we have pinned
  // the data instead?
  // if (Plastic_data_has_been_built[Fpcs_INDEX] &&
  //     Plastic_data_has_been_built[Fpks_INDEX])
  // {
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      barbar_M(i, j) = bar_M(i, j) - R * bar_Mk(i, j) - (1 - R) * bar_Mc(i, j);

      if (computeJacobian)
      {
        dbarbar_M_dR(i, j) = -bar_Mk(i, j) + bar_Mc(i, j);
      }
    }
  }

  if (computeJacobian)
  {
    dbarbar_M_dMk = -R;
    dbarbar_M_dMc = R - 1;
  }

  return;
  // }
  // else if (Plastic_data_has_been_built[Fpks_INDEX])
  // {
  //   for (unsigned int i = 0; i < DIM; i++)
  //   {
  //     for (unsigned int j = 0; j < DIM; j++)
  //     {
  //       barbar_M(i, j) = bar_M(i, j) - bar_Mk(i, j);
  //     }
  //   }

  //   if (computeJacobian)
  //   {
  //     dbarbar_M_dMk = -1.0;
  //     dbarbar_M_dMc = 0.0;
  //   }
  // }
  // else
  // {
  //   for (unsigned int i = 0; i < DIM; i++)
  //   {
  //     for (unsigned int j = 0; j < DIM; j++)
  //     {
  //       barbar_M(i, j) = bar_M(i, j);
  //     }
  //   }

  //   if (computeJacobian)
  //   {
  //     dbarbar_M_dMk = 0.0;
  //     dbarbar_M_dMc = 0.0;
  //   }
  // }

  // dbarbar_M_dR.initialise(0.0);
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_barbar_N(
  const DenseMatrix<double>& barbar_M,
  const double& f,
  const DenseMatrix<double>& dfdM,
  DenseMatrix<double>& barbar_N,
  RankFourTensor<double>& dbarbar_N_dbarbar_M,
  bool computeDerivative)
{
  barbar_N.resize(DIM, DIM, 0.0);

  MatrixHelpers::symmetrize(dfdM, barbar_N);
  double nMag = MatrixHelpers::magnitude(barbar_N);

  // Safety check for zero length
  if (nMag < 1.0e-15)
  {
    barbar_N.initialise(0.0);
    if (computeDerivative)
    {
      dbarbar_N_dbarbar_M.resize(DIM, DIM, DIM, DIM, 0.0);
    }
    return;
  }

  // Compute barbar_N
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      barbar_N(i, j) /= nMag;
    }
  }

  if (!computeDerivative) return;
  // The derivative is
  // dN_{ij}/dM_{kl} = [\delta_{is}\delta_{jt} - N_{ij}N_{st}] *
  // [ ddf/(dM_{st}dM_kl) + ddf/(dM_{ts}dM_kl) ] / (2 nMag)

  // We split the computation in two loop blocks - this saves if statements and
  // loop iterations

  // First, compute N_{st} * [ ddf/(dM_{st}dM_kl) + ddf/(dM_{ts}dM_kl) ]
  RankFourTensor<double> ddfdMdM;
  this->Plastic_consitutive_law_pt->compute_ddyield_surface_functiondMdM(
    f, dfdM, ddfdMdM);
  DenseMatrix<double> N_ddfdMdM(DIM, DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      N_ddfdMdM(i, j) = 0;
      for (unsigned int s = 0; s < DIM; s++)
      {
        for (unsigned int t = 0; t < DIM; t++)
        {
          N_ddfdMdM(i, j) +=
            barbar_N(s, t) * (ddfdMdM(s, t, i, j) + ddfdMdM(t, s, i, j));
        }
      }
    }
  }

  // Now compute the whole thing. We can now immediately contract the delta term
  // i.e. we do not need to iterate through s and t again
  dbarbar_N_dbarbar_M.resize(DIM, DIM, DIM, DIM, 0.0);
  double prefactor = 1 / (2 * nMag);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      for (unsigned int k = 0; k < DIM; k++)
      {
        for (unsigned int l = 0; l < DIM; l++)
        {
          double deltaFactor = ddfdMdM(i, j, k, l) + ddfdMdM(j, i, k, l);

          dbarbar_N_dbarbar_M(i, j, k, l) =
            prefactor * (deltaFactor - barbar_N(i, j) * N_ddfdMdM(k, l));
        }
      }
    }
  }
}

/*
 * \details
 * computes \bar{L}^\text{p} / \dot{\bar{\lambda}} =
 *                                            barbar_N + eta [\bar{M}, barbar_N]
 * and its derivatives wrt. to barbarM and \bar_{M}
 */
template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_bar_Lp(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& barbar_N,
  const RankFourTensor<double>& dbarbar_N_dbarbar_M,
  DenseMatrix<double>& bar_Lp,
  RankFourTensor<double>& dbar_Lp_dbarbar_M,
  RankFourTensor<double>& dbar_Lp_dbar_M,
  bool computeDerivative)
{
  const double eta = this->Plastic_consitutive_law_pt->eta_p;

  // Lp = barbar_N (if eta = 0)
  bar_Lp.resize(DIM, DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      bar_Lp(i, j) = barbar_N(i, j);
    }
  }

  // Add the commutator part, if neccesary
  if (eta != 0.0)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum = 0.0;

        // Compute (M*N - N*M)_ij
        for (unsigned int k = 0; k < DIM; k++)
        {
          sum += bar_M(i, k) * barbar_N(k, j) - barbar_N(i, k) * bar_M(k, j);
        }

        bar_Lp(i, j) += eta * sum;
      }
    }
  }

  // The derivative.
  if (computeDerivative)
  {
    // We do not need to initialize the value, it is done in the next for loop
    dbar_Lp_dbarbar_M.resize(DIM, DIM, DIM, DIM);
    dbar_Lp_dbar_M.resize(DIM, DIM, DIM, DIM);

    // First part: dN/dM
    // Use one loop only
    const unsigned long size4D = DIM * DIM * DIM * DIM;
    for (unsigned long idx = 0; idx < size4D; idx++)
    {
      double val = dbarbar_N_dbarbar_M.raw_direct_access(idx);
      dbar_Lp_dbarbar_M.raw_direct_access(idx) = val;
      dbar_Lp_dbar_M.raw_direct_access(idx) = val;
    }

    // The derivative of the second part [bar_M, barbar_N]
    // [bar_M, barbar_N]_{ij} = bar_M_{ia} barbar_N{aj} - barbar_N{ia} bar_M{aj}
    // d[bar_M, barbar_N]_{ij}/dbar_M{kl} =
    //      \delta_{ik}barbar_N{lj} + bar_M_{ia} dbarbar_N{aj}/dbar_M{kl}
    //    - dbarbar_N{ia}/dbar_M{kl} bar_M{aj} - barbar_N{ik}\delta{jl}
    // and d[bar_M, barbar_N]_{ij}/dbarbar_M{kl} =
    //      bar_M_{ia} dbarbar_N{aj}/dbarbar_M{kl}
    //      - dbarbar_N{ia}/dbarbar_M{kl} bar_M{aj}
    //
    // Note that dbarbar_N/dbar_M = dbarbar_N/dbarbar_M
    if (eta != 0.0)
    {
      // First we consider only the parts with dbar_M/dbar_M
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          // \delta_{ik}barbar_N{lj}
          for (unsigned int l = 0; l < DIM; l++)
          {
            dbar_Lp_dbar_M(i, j, i, l) += eta * barbar_N(l, j);
          }

          // - barbar_N{ik}\delta{jl}
          for (unsigned int k = 0; k < DIM; k++)
          {
            dbar_Lp_dbar_M(i, j, k, j) -= eta * barbar_N(i, k);
          }
        }
      }

      // Now the remainder
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          for (unsigned int k = 0; k < DIM; k++)
          {
            for (unsigned int l = 0; l < DIM; l++)
            {
              double sum = 0.0;
              for (unsigned int a = 0; a < DIM; a++)
              {
                // bar_M_{ia} dbarbar_N{aj}/dbar_M{kl}
                //                            dbarbar_N{ia}/dbar_M{kl} bar_M{aj}
                sum += bar_M(i, a) * dbarbar_N_dbarbar_M(a, j, k, l) -
                       dbarbar_N_dbarbar_M(i, a, k, l) * bar_M(a, j);
              }
              dbar_Lp_dbarbar_M(i, j, k, l) += eta * sum;
              dbar_Lp_dbar_M(i, j, k, l) += eta * sum;
            }
          }
        }
      }
    }
  }
}

/*
 * \details
 * computes \bar{L}^\text{pkd} / \dot{\bar{\lambda}} =
 *                       bar_Mk + etapk [\bar{M}, \bar{M}^\text{k}] / b^\text{k}
 * and its derivatives wrt. to \bar{M}^\text{k} and \bar{M}
 */
template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_bar_Lpkd(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& bar_Mk,
  DenseMatrix<double>& bar_Lpkd,
  RankFourTensor<double>& dbar_Lpkd_dbar_M,
  RankFourTensor<double>& dbar_Lpkd_dbar_Mk,
  bool computeDerivative)
{
  bar_Lpkd.resize(DIM, DIM);
  const double inv_bk =
    1.0 / this->Plastic_consitutive_law_pt->kinematic_hardening_b;
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      bar_Lpkd(i, j) = inv_bk * bar_Mk(i, j);
    }
  }

  const double eta = this->Plastic_consitutive_law_pt->kinematic_hardening_eta;
  const double invbk_Eta = eta * inv_bk;
  if (eta != 0.0)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum = 0.0;
        for (unsigned int k = 0; k < DIM; k++)
        {
          sum += bar_M(i, k) * bar_Mk(k, j) - bar_Mk(i, k) * bar_M(k, j);
        }
        bar_Lpkd(i, j) += invbk_Eta * sum;
      }
    }
  }

  if (computeDerivative)
  {
    dbar_Lpkd_dbar_M.resize(DIM, DIM, DIM, DIM, 0.0);
    dbar_Lpkd_dbar_Mk.resize(DIM, DIM, DIM, DIM, 0.0);
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        dbar_Lpkd_dbar_Mk(i, j, i, j) = inv_bk;
      }
    }

    // The derivative of the second part [bar_M, bar_Mk]
    // [bar_M, bar_Mk]_{ij} = bar_M_{ia} bar_Mk{aj} - bar_Mk{ia} bar_M{aj}
    // d[bar_M, bar_Mk]_{ij}/dbar_M{kl} =
    //      \delta_{ik}bar_Mk{lj} - bar_Mk{ik}\delta{jl}
    // and d[bar_M, bar_Mk]_{ij}/dbar_Mk{kl} =
    //      bar_M_{ik}\delta{jl} - \delta{ik}bar_M{lj}
    if (eta != 0.0)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          // The contribution of delta jl
          for (unsigned int k = 0; k < DIM; k++)
          {
            dbar_Lpkd_dbar_M(i, j, k, j) -= invbk_Eta * bar_Mk(i, k);
            dbar_Lpkd_dbar_Mk(i, j, k, j) += invbk_Eta * bar_M(i, k);
          }
          // The contribution of delta ik
          for (unsigned int l = 0; l < DIM; l++)
          {
            dbar_Lpkd_dbar_M(i, j, i, l) += invbk_Eta * bar_Mk(l, j);
            dbar_Lpkd_dbar_Mk(i, j, i, l) -= invbk_Eta * bar_M(l, j);
          }
        }
      }
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_hat_bar_Nc(
  const double& f_Mc,
  const DenseMatrix<double>& df_Mc_dMc,
  DenseMatrix<double>& hat_bar_Nc,
  RankFourTensor<double>& dhat_bar_Nc_dMc,
  bool computeDerivative)
{
  hat_bar_Nc.resize(DIM, DIM);

  double nMag = MatrixHelpers::magnitude(df_Mc_dMc);
  // Safety check for zero length
  if (nMag < 1.0e-15)
  {
    hat_bar_Nc.initialise(0.0);
    if (computeDerivative)
    {
      dhat_bar_Nc_dMc.resize(DIM, DIM, DIM, DIM, 0.0);
    }
    return;
  }


  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      hat_bar_Nc(i, j) = df_Mc_dMc(i, j) / nMag;
    }
  }

  if (!computeDerivative) return;

  // The derivative is (A = df_Mc_dMc)
  // dN_ij/dA_mn = (\delta_im \delta_jn - N_ij N_mn) / nMag
  // For dN_ij_dhat_bar_Mc_{mn} one must just use the chain rule:
  // dN_ij/dhat_bar_Mc_{kl} = (dN_ij/dA_kl + N_ij N_mn dA_mn / dM_kl) / nMag
  RankFourTensor<double> dfdMcdMc;
  this->Plastic_consitutive_law_pt->compute_ddyield_surface_functiondMdM(
    f_Mc, df_Mc_dMc, dfdMcdMc);

  DenseMatrix<double> Nc_ddfdMcdMc(DIM, DIM);
  for (unsigned int k = 0; k < DIM; k++)
  {
    for (unsigned int l = 0; l < DIM; l++)
    {
      double sum = 0.0;
      for (unsigned int m = 0; m < DIM; m++)
      {
        for (unsigned int n = 0; n < DIM; n++)
        {
          sum += hat_bar_Nc(m, n) * dfdMcdMc(m, n, k, l);
        }
      }

      Nc_ddfdMcdMc(k, l) = sum;
    }
  }

  dhat_bar_Nc_dMc.resize(DIM, DIM, DIM, DIM, 0.0);
  double nMag_inv = 1 / nMag;
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      for (unsigned int k = 0; k < DIM; k++)
      {
        for (unsigned int l = 0; l < DIM; l++)
        {
          dhat_bar_Nc_dMc(i, j, k, l) =
            nMag_inv *
            (dfdMcdMc(i, j, k, l) - hat_bar_Nc(i, j) * Nc_ddfdMcdMc(k, l));
        }
      }
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::compute_bar_Lpcd(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& hat_bar_Nc,
  const double& Rc,
  const RankFourTensor<double>& dhat_bar_Nc_dhat_bar_Mc,
  const DenseMatrix<double>& dRc_dMc,
  const double& dRc_dH,
  DenseMatrix<double>& bar_Lpcd,
  RankFourTensor<double>& dbar_Lpcd_dbar_M,
  RankFourTensor<double>& dbar_Lpcd_dhat_bar_Mc,
  DenseMatrix<double>& dbar_Lpcd_dh,
  bool computeDerivative)
{
  bar_Lpcd.resize(DIM, DIM);

  const double invX = 1 / this->Plastic_consitutive_law_pt->elastic_core_x;
  const double Rc_by_X = Rc * invX;
  const double eta = this->Plastic_consitutive_law_pt->elastic_core_eta;
  const double eta_prefactor = Rc_by_X * eta;

  // Really, we compute Lpcd / Rc; Rc is muliplited to it after the derivative
  // computation
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      double sum = hat_bar_Nc(i, j);
      if (eta != 0.0)
      {
        for (unsigned int a = 0; a < DIM; a++)
        {
          sum += eta * (bar_M(i, a) * hat_bar_Nc(a, j) -
                        hat_bar_Nc(i, a) * bar_M(a, j));
        }
      }

      bar_Lpcd(i, j) = invX * sum;
    }
  }

  if (computeDerivative)
  {
    // dLpcs_ij/dbar_M_kl =
    //          Rc/X * eta (\delta_ik \hat\bar{Nc}_lj - \hat\bar{Nc}_ik
    //          \delta_jl)
    dbar_Lpcd_dbar_M.resize(DIM, DIM, DIM, DIM, 0.0);
    dbar_Lpcd_dhat_bar_Mc.resize(DIM, DIM, DIM, DIM);

    // dLpcd_dh = (\hat\bar{N} - eta [\bar{M}, \hat\bar{N}]) / X * dRc/dh
    dbar_Lpcd_dh.resize(DIM, DIM);
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        dbar_Lpcd_dh(i, j) = bar_Lpcd(i, j) * dRc_dH;

        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            // The partial derivative
            double dNdMc = dhat_bar_Nc_dhat_bar_Mc(i, j, k, l);

            if (eta != 0.0)
            {
              // d/dMc_kl [M, Nc]_ij = M_ia dNcdMc_ajkl - M_aj dNcdMc_iakl
              for (unsigned int a = 0; a < DIM; a++)
              {
                dNdMc +=
                  eta * (bar_M(i, a) * dhat_bar_Nc_dhat_bar_Mc(a, j, k, l) -
                         bar_M(a, j) * dhat_bar_Nc_dhat_bar_Mc(i, a, k, l));
              }
            }

            // the total derivative
            dbar_Lpcd_dhat_bar_Mc(i, j, k, l) =
              Rc_by_X * dNdMc + bar_Lpcd(i, j) * dRc_dMc(k, l);
          }
        }
      }
    }


    if (eta != 0.0)
    {
      // Derivative wrt. \bar M
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          for (unsigned int k = 0; k < DIM; k++)
          {
            dbar_Lpcd_dbar_M(i, j, k, j) -= eta_prefactor * hat_bar_Nc(i, k);
          }

          for (unsigned int l = 0; l < DIM; l++)
          {
            dbar_Lpcd_dbar_M(i, j, i, l) += eta_prefactor * hat_bar_Nc(l, j);
          }
        }
      }
    }
  }

  // finally multiply with Rc
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      bar_Lpcd(i, j) = Rc * bar_Lpcd(i, j);
    }
  }
}

template<unsigned DIM>
double oomph::PlasticEquations<DIM>::compute_c_sigma(
  const DenseMatrix<double>& bar_bar_N,
  const DenseMatrix<double>& hat_bar_Nc,
  const RankFourTensor<double>& dbarbar_N_dbarbar_M,
  const double& dbarbar_M_dbar_Mk,
  const double& dbarbar_M_dbar_Mc,
  const DenseMatrix<double>& dbarbarM_dR,
  const RankFourTensor<double>& dhatbar_Nc_dhatbar_Mc,
  DenseMatrix<double>& dC_dbar_M,
  DenseMatrix<double>& dC_dbar_Mk,
  DenseMatrix<double>& dC_dbar_Mc,
  double& dCdR,
  bool computeDerivative)
{
  if (computeDerivative)
  {
    dC_dbar_M.resize(DIM, DIM);
    dC_dbar_Mk.resize(DIM, DIM);
    dC_dbar_Mc.resize(DIM, DIM);
    dCdR = 0.0;

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
            term2_sum += bar_bar_N(a, b) * dhatbar_Nc_dhatbar_Mc(a, b, i, j);
          }
        }

        dCdR += term1_sum * dbarbarM_dR(i, j);

        dC_dbar_M(i, j) = term1_sum;
        dC_dbar_Mk(i, j) = term1_sum * dbarbar_M_dbar_Mk - term2_sum;
        dC_dbar_Mc(i, j) = term1_sum * dbarbar_M_dbar_Mc + term2_sum;
      }
    }
  }

  return MatrixHelpers::reduce(bar_bar_N, hat_bar_Nc);
}

template<unsigned DIM>
double oomph::PlasticEquations<DIM>::compute_u(
  const double& u_in,
  const double& Rc,
  const double& c_sigma,
  const DenseMatrix<double>& dRc_dhatbar_Mc,
  const double& dRc_dh,
  const DenseMatrix<double>& dc_sigma_dbar_M,
  const DenseMatrix<double>& dc_sigma_dbar_Mk,
  const DenseMatrix<double>& dc_sigma_dbar_Mc,
  const double& dc_sigma_dR,
  DenseMatrix<double>& du_dbar_M,
  DenseMatrix<double>& du_dbar_Mk,
  DenseMatrix<double>& du_dbar_Mc,
  double& du_dh,
  double& du_dR,
  bool computeDerivative)
{
  const double uc = this->Plastic_consitutive_law_pt->elasic_core_u;
  const double u_out = u_in * std::exp(uc * Rc * c_sigma);

  if (computeDerivative)
  {
    const double prefactor_dc_sigma = uc * u_out * Rc;
    const double prefactor_drc = uc * u_out * c_sigma;

    du_dbar_M.resize(DIM, DIM);
    du_dbar_Mk.resize(DIM, DIM);
    du_dbar_Mc.resize(DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        du_dbar_M(i, j) = prefactor_dc_sigma * dc_sigma_dbar_M(i, j);
        du_dbar_Mk(i, j) = prefactor_dc_sigma * dc_sigma_dbar_Mk(i, j) -
                           prefactor_drc * dRc_dhatbar_Mc(i, j);
        du_dbar_Mc(i, j) = prefactor_dc_sigma * dc_sigma_dbar_Mc(i, j) +
                           prefactor_drc * dRc_dhatbar_Mc(i, j);
      }
    }

    du_dh = prefactor_drc * dRc_dh;
    du_dR = prefactor_dc_sigma * dc_sigma_dR;
  }

  return u_out;
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::initialise_solve(const unsigned ipt)
{
  for (unsigned data_type = 0; data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
       data_type++)
  {
    Data* data_pt = Plastic_data_pt[ipt][data_type];

    // If there is not enough history, just initiallise all data with their
    // initial condition
    if (data_pt->time_stepper_pt()->ntstorage() < 2)
    {
      set_intial_condition(ipt);
      return;
    }

    // Otherwise, use the value from the last timestep
    const unsigned nval = data_pt->nvalue();
    for (unsigned i = 0; i < nval; i++)
    {
      data_pt->set_value(i, data_pt->value(1, i));
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquations<DIM>::set_intial_condition(const unsigned int ipt)
{
  for (unsigned data_type = 0; data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
       data_type++)
  {
    Data* data_pt = Plastic_data_pt[ipt][data_type];

    // Set the deformation gradient tensors to unity
    if (data_type == invFp_INDEX || data_type == Fpks_INDEX ||
        data_type == Fpcs_INDEX)
    {
      const unsigned nval = data_pt->nvalue();
      for (unsigned i = 0; i < nval; i++)
      {
        // Check if diagonal element (Assuming row-major storage for DIMxDIM)
        if (i % (DIM + 1) == 0) data_pt->set_value(i, 1.0);
        else
          data_pt->set_value(i, 0.0);
      }
    }

    // Set R to Re
    else if (data_type == R_INDEX)
    {
      double Re = 0.0;
      if (this->Plastic_consitutive_law_pt)
      {
        Re = this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic;
      }
      data_pt->set_value(0, Re);
    }

    // Set Lambda and H to 0
    else
    {
      const unsigned nval = data_pt->nvalue();
      for (unsigned i = 0; i < nval; i++)
      {
        data_pt->set_value(i, 0.0);
      }
    }
  }
}

/*!
 * \details this function checks two conditions. If both are true, there is
 * plastic deformation The conditions are:
 *        1. f(M) > F(H)
 *        2. varbarN : LL : sym(Ce barL) > 0
 */
template<unsigned DIM>
bool PlasticEquations<DIM>::is_there_plastic_deformation(const unsigned int ipt)
{
  // Yet another version. This time Eq. 193 from
  // https://doi.org/10.1007/s11831-018-9256-5

  // Get the plastic quantities
  DenseMatrix<double> invFp(DIM, DIM, 0.0), Fpks(DIM, DIM, 0.0),
    Fpcs(DIM, DIM, 0.0);

  get_inv_fp_matrix(ipt, invFp);

  get_fpks_matrix(ipt, Fpks);

  get_fpcs_matrix(ipt, Fpcs);

  // Get R
  double R = get_r(ipt);

  // Retrieve deformation gradient tensor
  DenseMatrix<double> F(DIM), C_Total(DIM, DIM, 0.0);
  compute_deformation_gradient_tensor(0, ipt, F);
  MatrixHelpers::transpose_multiply(F, F, C_Total);

  // Now compute the elastic Mandel stress
  DenseMatrix<double> bar_M(DIM, DIM, 0.0);
  compute_mandel_stress_elastic(invFp, C_Total, bar_M);

  // Compute bar_Mk and bar_Mc
  DenseMatrix<double> bar_Mk(DIM, DIM, 0.0), bar_Mc(DIM, DIM, 0.0);
  compute_mandellike_kinematic_hardening(Fpks, bar_Mk);
  compute_mandellike_elastic_core(Fpcs, bar_Mc);

  DenseMatrix<double> barbar_M(DIM, DIM, 0.0);
  compute_mandel_stress_total(bar_M, bar_Mk, bar_Mc, R, barbar_M);

  // Now compute the previous value - the only thing that has changed if F
  // Get previous F
  DenseMatrix<double> F_prev(DIM), C_Total_prev(DIM, DIM, 0.0);
  compute_deformation_gradient_tensor(1, ipt, F_prev);
  MatrixHelpers::transpose_multiply(F_prev, F_prev, C_Total_prev);

  // Previous bar_M
  DenseMatrix<double> bar_M_prev(DIM, DIM, 0.0);
  compute_mandel_stress_elastic(invFp, C_Total_prev, bar_M_prev);

  // Compute Mbarbar_prev
  DenseMatrix<double> barbar_M_prev(DIM);
  compute_mandel_stress_total(bar_M_prev, bar_Mk, bar_Mc, R, barbar_M_prev);

  // Now we can calculate deltaMTrial = bar_M - bar_M_prev
  DenseMatrix<double> deltaMTrial(DIM, DIM, 0.0);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      deltaMTrial(i, j) = bar_M(i, j) - bar_M_prev(i, j);
    }
  }

  // In case only Fpks has been active
  // We now always build all plastic data. Do we need to check if we have pinned
  // the data instead?
  // if (Plastic_data_has_been_built[Fpks_INDEX] &&
  //     !Plastic_data_has_been_built[Fpcs_INDEX])
  // {
  //   bar_Mc.resize(DIM, DIM);
  //   for (unsigned int i = 0; i < DIM; i++)
  //   {
  //     for (unsigned int j = 0; j < DIM; j++)
  //     {
  //       bar_Mc(i, j) = bar_Mk(i, j);
  //     }
  //   }
  // }

  // Compute Nbarbar_prev
  DenseMatrix<double> Nbarbar_prev_nsym(DIM);
  double yieldSurfaceStress_prev =
    this->Plastic_consitutive_law_pt->compute_yield_surface_function(
      barbar_M_prev, Nbarbar_prev_nsym, true);
  MatrixHelpers::normalise(Nbarbar_prev_nsym);

  // Compute Mbarbar_e = barbar_M - Re Mbark_prev - (1 - Re) Mbarc_prev
  double Re = 1.0;

  Re = this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic;

  DenseMatrix<double> Mbarbar_e(DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      Mbarbar_e(i, j) =
        bar_M(i, j) - Re * bar_Mk(i, j) - (1 - Re) * bar_Mc(i, j);
    }
  }

  double h_var = get_lambda(ipt) *
                 this->Plastic_consitutive_law_pt->isotropic_hardening_factor;
  bool plasticDeformation = false;
  if (MatrixHelpers::reduce(Nbarbar_prev_nsym, deltaMTrial) >= 0)
  {
    // This is the case of forward loading

    // Check yield surface
    plasticDeformation =
      this->Plastic_consitutive_law_pt->compute_yield_surface_function(
        Mbarbar_e, GeneralisedElement::Dummy_matrix, false) -
        Re *
          this->Plastic_consitutive_law_pt->isotropic_hardening_yield_function(
            h_var) >
      0;
  }
  else
  {
    // Inverse loading

    // Compute barbar_N
    DenseMatrix<double> Nbarbar_nsym(DIM);
    double yieldSurfaceStress =
      this->Plastic_consitutive_law_pt->compute_yield_surface_function(
        barbar_M, Nbarbar_nsym, true);
    MatrixHelpers::normalise(Nbarbar_nsym);

    if (MatrixHelpers::reduce(Nbarbar_nsym, deltaMTrial) <= 0)
    {
      plasticDeformation = false;
    }
    else
    {
      plasticDeformation =
        this->Plastic_consitutive_law_pt->compute_yield_surface_function(
          Mbarbar_e, GeneralisedElement::Dummy_matrix, false) -
          Re * this->Plastic_consitutive_law_pt
                 ->isotropic_hardening_yield_function(h_var) >
        0;
    }
  }

  if (!plasticDeformation)
  {
    // Compute R from the yield surface condition
    double R =
      this->Plastic_consitutive_law_pt->compute_yield_surface_function(
        barbar_M, GeneralisedElement::Dummy_matrix, false) /
      this->Plastic_consitutive_law_pt->isotropic_hardening_yield_function(
        get_lambda(ipt) *
        this->Plastic_consitutive_law_pt->isotropic_hardening_factor);

    set_r(ipt,
          std::max(
            R, this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic));
  }

  return plasticDeformation;
}

/*!
 * \details This function computes the residuals of the plastic data. This is
 * done in two ways.
 *
 * 1. If the timestepper is unsteady, the residual is the difference of the time
 * derivatives of the plastic variables, e.g.,
 *    r_{Fp} = \dot{F}^{\text{p} -1}_\text{timestepper}
 *           - \dot{F}^{\text{p} -1}_\text{analytical}.
 *
 * 2. If the timestepper is steady, the residual is computed based on a finite
 * increment of the plastic variables compared to the previous values. If the
 * time stepper does not save any previous value, the default value is assumed
 * as the previous values. The increments are computed assuming backwards Euler
 * integration, e.g.,
 * \Delta F^{\text{p} -1} = \dot{F}^{\text{p} -1} / \dot{\lambda} \Delta\lambda.
 */
template<unsigned DIM>
void PlasticEquations<DIM>::fill_in_generic_residual_and_jacobian_plastic(
  DoubleVector& residuals,
  DenseMatrix<double>& jacobian,
  const unsigned& ipt,
  const DenseMatrix<double>& C,
  const unsigned& flag)
{
  // Is the time-stepper steady?
  Plastic_data_pt[ipt][invFp_INDEX]->time_stepper_pt()->is_steady();

  // Retreive the time-step weights
  double inv_Fp_time_step_weight = 1.0;
  double dot_or_delta_lambda_time_step_weight = 1.0;
  double fpks_time_step_weight = 1.0;
  double fpcs_time_step_weight = 1.0;
  double r_time_step_weight = 1.0;

  if (!Plastic_data_pt[ipt][invFp_INDEX]->time_stepper_pt()->is_steady())
  {
    inv_Fp_time_step_weight =
      Plastic_data_pt[ipt][invFp_INDEX]->time_stepper_pt()->weight(1, 0);
  }
  if (!Plastic_data_pt[ipt][Lambda_INDEX]->time_stepper_pt()->is_steady())
  {
    dot_or_delta_lambda_time_step_weight =
      Plastic_data_pt[ipt][Lambda_INDEX]->time_stepper_pt()->weight(1, 0);
  }
  if (!Plastic_data_pt[ipt][Fpks_INDEX]->time_stepper_pt()->is_steady())
  {
    fpks_time_step_weight =
      Plastic_data_pt[ipt][Fpks_INDEX]->time_stepper_pt()->weight(1, 0);
  }
  if (!Plastic_data_pt[ipt][Fpcs_INDEX]->time_stepper_pt()->is_steady())
  {
    fpcs_time_step_weight =
      Plastic_data_pt[ipt][Fpcs_INDEX]->time_stepper_pt()->weight(1, 0);
  }
  if (!Plastic_data_pt[ipt][R_INDEX]->time_stepper_pt()->is_steady())
  {
    r_time_step_weight =
      Plastic_data_pt[ipt][R_INDEX]->time_stepper_pt()->weight(1, 0);
  }

  // Retreive dot_or_delta_lambda
  const double dot_or_delta_lambda = this->get_dot_or_delta_lambda(ipt);

  // Compute R
  double R = get_r(ipt);


  // Will be used to help with the computations
  DenseMatrix<double> tmp(DIM);

  // Get invFp
  DenseMatrix<double> invFp(DIM, DIM, 0.0);
  get_inv_fp_matrix(ipt, invFp);

  // Now compute the elastic Mandel stress
  DenseMatrix<double> bar_M(DIM, DIM, 0.0);
  RankFourTensor<double> dbar_M_dinv_Fp;
  compute_mandel_stress_elastic(invFp, C, bar_M, dbar_M_dinv_Fp, flag);

  // Get Fpks and Fpkc
  DenseMatrix<double> Fpks, Fpcs;

  get_fpks_matrix(ipt, Fpks);

  get_fpcs_matrix(ipt, Fpcs);

  // Compute bar_Mk and bar_Mc
  DenseMatrix<double> bar_Mk(DIM, DIM, 0.0), bar_Mc(DIM, DIM, 0.0);
  RankFourTensor<double> dbar_Mk_dFpks, dbar_Mc_dFpcs;
  compute_mandellike_kinematic_hardening(Fpks, bar_Mk, dbar_Mk_dFpks, flag);
  compute_mandellike_elastic_core(Fpcs, bar_Mc, dbar_Mc_dFpcs, flag);

  DenseMatrix<double> barbar_M(DIM);
  double dbarbar_M_dMk, dbarbar_M_dMc;
  DenseMatrix<double> dbarbar_M_dR;
  compute_mandel_stress_total(bar_M,
                              bar_Mk,
                              bar_Mc,
                              R,
                              barbar_M,
                              dbarbar_M_dMk,
                              dbarbar_M_dMc,
                              dbarbar_M_dR,
                              flag);

  // barbar_N
  DenseMatrix<double> barbar_N(DIM);
  DenseMatrix<double> dfdM(DIM);
  double yieldSurfaceStress =
    this->Plastic_consitutive_law_pt->compute_yield_surface_function(
      barbar_M, dfdM, true);

  RankFourTensor<double> dbarbarN_dbarbar_M;
  compute_barbar_N(
    barbar_M, yieldSurfaceStress, dfdM, barbar_N, dbarbarN_dbarbar_M, flag);

  // L_p / dot_lambda
  DenseMatrix<double> bar_Lp(DIM);
  RankFourTensor<double> dbar_Lp_dbarbar_M;
  RankFourTensor<double> dbar_Lp_dbar_M;
  compute_bar_Lp(bar_M,
                 barbar_N,
                 dbarbarN_dbarbar_M,
                 bar_Lp,
                 dbar_Lp_dbarbar_M,
                 dbar_Lp_dbar_M,
                 flag);

  // Apply the chain rule, if neccesary
  RankFourTensor<double> dbar_Lp_dinv_Fp;
  RankFourTensor<double> dbar_Lp_dFpks;
  RankFourTensor<double> dbar_Lp_dFpcs;
  DenseMatrix<double> dbar_Lp_dR;

  if (flag)
  {
    dbar_Lp_dinv_Fp.resize(DIM, DIM, DIM, DIM);
    dbar_Lp_dFpks.resize(DIM, DIM, DIM, DIM);
    dbar_Lp_dFpcs.resize(DIM, DIM, DIM, DIM);
    dbar_Lp_dR.resize(DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            double sum_invFp = 0.0;
            double sum_Fpks = 0.0;
            double sum_Fpcs = 0.0;
            for (unsigned int p = 0; p < DIM; p++)
            {
              for (unsigned int q = 0; q < DIM; q++)
              {
                sum_invFp +=
                  dbar_Lp_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                const double dbar_Lp_dbarbar_M_ijpq =
                  dbar_Lp_dbarbar_M(i, j, p, q);

                sum_Fpks += dbar_Lp_dbarbar_M_ijpq * dbar_Mk_dFpks(p, q, k, l);

                sum_Fpcs += dbar_Lp_dbarbar_M_ijpq * dbar_Mc_dFpcs(p, q, k, l);
              }
            }
            dbar_Lp_dinv_Fp(i, j, k, l) = sum_invFp;

            dbar_Lp_dFpks(i, j, k, l) = sum_Fpks * dbarbar_M_dMk;

            dbar_Lp_dFpcs(i, j, k, l) = sum_Fpcs * dbarbar_M_dMc;
          }
        }
      }
    }

    // now we do the R derivative: dbar_Lp_dbar_M : dbarbar_M_dR
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum_invFp_ij = 0;
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            sum_invFp_ij += dbar_Lp_dbarbar_M(i, j, k, l) * dbarbar_M_dR(k, l);
          }
        }

        dbar_Lp_dR(i, j) = sum_invFp_ij;
      }
    }
  }

  // Finally minus_dotinvFp
  DenseMatrix<double> minus_dotinvFp;
  MatrixHelpers::multiply(invFp, bar_Lp, minus_dotinvFp);

  DenseMatrix<double> dot_or_delta_invFp_fromTimeStepper(DIM);
  get_dot_or_delta_inv_fp_matrix(ipt, dot_or_delta_invFp_fromTimeStepper);

  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      residuals[this->plastic_inv_fp_eqn_number(ipt, i, j)] =
        dot_or_delta_invFp_fromTimeStepper(i, j) +
        dot_or_delta_lambda * minus_dotinvFp(i, j);
    }
  }

  if (flag)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        const unsigned int row_eq = this->plastic_inv_fp_eqn_number(ipt, i, j);

        // Contribution from the time stepper
        jacobian(row_eq, row_eq) += inv_Fp_time_step_weight;

        // Contribution from
        // d(inv_Fp bar_Lp)_ij/dinvFp_kl = \delta_{ik} bar_Lp_{lj}
        //                                     + inv_Fp_ia, dL_{aj}/dinv_FP_{kl}
        for (unsigned int l = 0; l < DIM; l++)
        {
          // First the delta term:
          jacobian(row_eq, this->plastic_inv_fp_eqn_number(ipt, i, l)) +=
            dot_or_delta_lambda * bar_Lp(l, j);

          // Now the remainder
          for (unsigned int k = 0; k < DIM; k++)
          {
            double sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              sum += invFp(i, a) * dbar_Lp_dinv_Fp(a, j, k, l);
            }
            // Derivative wrt. invFp
            const unsigned int col_eq_inv_Fp =
              this->plastic_inv_fp_eqn_number(ipt, k, l);

            jacobian(row_eq, col_eq_inv_Fp) += dot_or_delta_lambda * sum;
          }
        }

        // Contribution from
        // d(inv_Fp bar_Lp)_ij/dinvFp_kl = inv_Fp_ia dL_{aj}/dFPks_{kl}
        for (unsigned int l = 0; l < DIM; l++)
        {
          // Now the remainder
          for (unsigned int k = 0; k < DIM; k++)
          {
            double sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              sum += invFp(i, a) * dbar_Lp_dFpks(a, j, k, l);
            }
            // Derivative wrt. invFp
            const unsigned int col_eq_fpks =
              this->plastic_fpks_eqn_number(ipt, k, l);

            jacobian(row_eq, col_eq_fpks) += dot_or_delta_lambda * sum;
          }
        }

        // Contribution from
        // d(inv_Fp bar_Lp)_ij/dFpcs_kl = inv_Fp_ia, dL_{aj}/dFPcs_{kl}
        for (unsigned int l = 0; l < DIM; l++)
        {
          // Now the remainder
          for (unsigned int k = 0; k < DIM; k++)
          {
            double sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              sum += invFp(i, a) * dbar_Lp_dFpcs(a, j, k, l);
            }
            // Derivative wrt. invFp
            const unsigned int col_eq_fpcs =
              this->plastic_fpcs_eqn_number(ipt, k, l);

            jacobian(row_eq, col_eq_fpcs) += dot_or_delta_lambda * sum;
          }
        }


        // Derivative wrt. lambda
        const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
        jacobian(row_eq, lambda_col) +=
          minus_dotinvFp(i, j) * dot_or_delta_lambda_time_step_weight;

        // Derivative wrt. R
        const unsigned int R_col = this->plastic_r_eqn_number(ipt);
        double sum_R = 0.0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          // invFp(i,a) * dbar_Lp_dR(a,j)
          sum_R += invFp(i, a) * dbar_Lp_dR(a, j);
        }
        jacobian(row_eq, R_col) += dot_or_delta_lambda * sum_R;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // The residual for lambda and R                                            //
  //////////////////////////////////////////////////////////////////////////////
  double h_var = this->get_lambda(ipt) *
                 this->Plastic_consitutive_law_pt->isotropic_hardening_factor;

  const unsigned int row_lamda = this->plastic_lambda_eqn_number(ipt);

  double disotropic_f_dH = 0.0;

  double isotropic_yield_stress =
    this->Plastic_consitutive_law_pt->isotropic_hardening_yield_function(
      h_var, disotropic_f_dH, flag);
  residuals[row_lamda] = yieldSurfaceStress - R * isotropic_yield_stress;

  if (flag)
  {
    // Contribution from lambda
    jacobian(row_lamda, row_lamda) -=
      disotropic_f_dH *
      this->Plastic_consitutive_law_pt->isotropic_hardening_factor * R;


    // Contribution from Fpks:
    // df/dM_ab * dM_ab/dFpks_ij
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum_invFp_ij = 0;
        double sum_Fpks_ij = 0;
        double sum_Fpcs_ij = 0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          for (unsigned int b = 0; b < DIM; b++)
          {
            const double dfdM_ab = dfdM(a, b);
            sum_invFp_ij += dfdM_ab * dbar_M_dinv_Fp(a, b, i, j);

            sum_Fpks_ij += dfdM_ab * dbar_Mk_dFpks(a, b, i, j);

            sum_Fpcs_ij += dfdM_ab * dbar_Mc_dFpcs(a, b, i, j);
          }
        }

        const unsigned int invFp_ij_col =
          this->plastic_inv_fp_eqn_number(ipt, i, j);
        jacobian(row_lamda, invFp_ij_col) += sum_invFp_ij;

        const unsigned int Fpks_ij_col =
          this->plastic_fpks_eqn_number(ipt, i, j);
        jacobian(row_lamda, Fpks_ij_col) += sum_Fpks_ij * dbarbar_M_dMk;

        const unsigned int Fpcs_ij_col =
          this->plastic_fpcs_eqn_number(ipt, i, j);
        jacobian(row_lamda, Fpcs_ij_col) += sum_Fpcs_ij * dbarbar_M_dMc;
      }
    }

    const unsigned int r_col = this->plastic_r_eqn_number(ipt);
    jacobian(row_lamda, r_col) -= isotropic_yield_stress;

    // yieldSurfaceStress = f(barbarM) also depends on R
    double sum_df_dR = 0.0;
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        sum_df_dR += dfdM(i, j) * dbarbar_M_dR(i, j);
      }
    }
    jacobian(row_lamda, r_col) += sum_df_dR;
  }


  //////////////////////////////////////////////////////////////////////////////
  // The residual for Fpks: \dot{Fpks} = (\bar Lp - \bar Lpkd) Fpks           //
  //////////////////////////////////////////////////////////////////////////////

  DenseMatrix<double> bar_Lpkd(DIM, DIM);
  RankFourTensor<double> dbar_Lpkd_dbar_M, dbar_Lpkd_dbar_Mk;
  compute_bar_Lpkd(
    bar_M, bar_Mk, bar_Lpkd, dbar_Lpkd_dbar_M, dbar_Lpkd_dbar_Mk, flag);

  DenseMatrix<double> dotFpks(DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      double sum_ij = 0.0;
      for (unsigned int a = 0; a < DIM; a++)
      {
        sum_ij += (bar_Lp(i, a) - bar_Lpkd(i, a)) * Fpks(a, j);
      }
      dotFpks(i, j) = sum_ij;
    }
  }

  DenseMatrix<double> dot_or_delta_Fpks_fromTimeStepper;
  get_dot_or_delta_fpks_matrix(ipt, dot_or_delta_Fpks_fromTimeStepper);

  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      residuals[this->plastic_fpks_eqn_number(ipt, i, j)] =
        dot_or_delta_Fpks_fromTimeStepper(i, j) -
        dot_or_delta_lambda * dotFpks(i, j);
    }
  }

  if (flag)
  {
    // Precompute some things
    RankFourTensor<double> dbar_Lpkd_dinvFp(DIM, DIM, DIM, DIM),
      dbar_Lpkd_dFpks(DIM, DIM, DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            double sum_invFp = 0.0;
            double sum_Fpks = 0.0;
            for (unsigned int p = 0; p < DIM; p++)
            {
              for (unsigned int q = 0; q < DIM; q++)
              {
                sum_invFp +=
                  dbar_Lpkd_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                sum_Fpks +=
                  dbar_Lpkd_dbar_Mk(i, j, p, q) * dbar_Mk_dFpks(p, q, k, l);
              }
            }
            dbar_Lpkd_dinvFp(i, j, k, l) = sum_invFp;
            dbar_Lpkd_dFpks(i, j, k, l) = sum_Fpks * dbarbar_M_dMk;
          }
        }
      }
    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        const unsigned int row_eq_fpks =
          this->plastic_fpks_eqn_number(ipt, i, j);

        // Contribution from the time stepper
        jacobian(row_eq_fpks, row_eq_fpks) += fpks_time_step_weight;

        // Derivative wrt. Fpks
        for (unsigned int k = 0; k < DIM; k++)
        {
          // The delta term \delta_{jl} (Lp - Lpkd)_{ik}
          jacobian(row_eq_fpks, this->plastic_fpks_eqn_number(ipt, k, j)) -=
            dot_or_delta_lambda * (bar_Lp(i, k) - bar_Lpkd(i, k));

          for (unsigned int l = 0; l < DIM; l++)
          {
            const unsigned int col_eq_fpks =
              this->plastic_fpks_eqn_number(ipt, k, l);

            const unsigned int col_eq_invFp =
              this->plastic_inv_fp_eqn_number(ipt, k, l);

            unsigned int col_eq_fpcs = 0;
            col_eq_fpcs = this->plastic_fpcs_eqn_number(ipt, k, l);

            // Remainder FPks term:
            // Fpks_{aj} dbar_Lp_ia / dFpks_kl
            // - Fpks_{aj} dbar_Lpks_ia / dFpks_kl
            //
            // invFp term:
            // Fpks_{aj} dbar_Lp_ia / dinvFp_kl
            // - Fpks_{aj} dbar_Lpks_ia / dinvFp_kl
            //
            // Fpcs term:
            // Fpks_{aj} dbar_Lp_{ia} / dFpcs_{kl}
            double fpks_sum = 0.0;
            double invFp_sum = 0.0;
            double fpcs_sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              fpks_sum += Fpks(a, j) * (dbar_Lp_dFpks(i, a, k, l) -
                                        dbar_Lpkd_dFpks(i, a, k, l));

              invFp_sum += Fpks(a, j) * (dbar_Lp_dinv_Fp(i, a, k, l) -
                                         dbar_Lpkd_dinvFp(i, a, k, l));

              fpcs_sum += Fpks(a, j) * (dbar_Lp_dFpcs(i, a, k, l));
            }
            jacobian(row_eq_fpks, col_eq_fpks) -=
              dot_or_delta_lambda * fpks_sum;
            jacobian(row_eq_fpks, col_eq_invFp) -=
              dot_or_delta_lambda * invFp_sum;
            jacobian(row_eq_fpks, col_eq_fpcs) -=
              dot_or_delta_lambda * fpcs_sum;
          }
        }


        // Derivative wrt. lambda
        const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
        jacobian(row_eq_fpks, lambda_col) -=
          dotFpks(i, j) * dot_or_delta_lambda_time_step_weight;

        // Derivative wrt. R
        const unsigned int R_col = this->plastic_r_eqn_number(ipt);
        double sum_R = 0.0;
        for (unsigned int l = 0; l < DIM; l++)
        {
          sum_R += Fpks(l, j) * dbar_Lp_dR(i, l);
        }
        jacobian(row_eq_fpks, R_col) -= dot_or_delta_lambda * sum_R;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // The residual for Fpcs: \dot{Fpcs} = (\bar Lp - \bar Lpcd) Fpcs           //
  //////////////////////////////////////////////////////////////////////////////
  double u = this->Plastic_consitutive_law_pt->normal_yield_ratio_constant_u;
  DenseMatrix<double> du_dbar_M, du_dbar_Mk, du_dbar_Mc;
  double du_dh = 0, du_dR = 0;

  if (flag)
  {
    du_dbar_M.resize(DIM, DIM, 0.0);
    du_dbar_Mk.resize(DIM, DIM, 0.0);
    du_dbar_Mc.resize(DIM, DIM, 0.0);
  }

  // First, we need to compute Mbarhat_c
  DenseMatrix<double> hat_bar_Mc(DIM, DIM, 0.0);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      hat_bar_Mc(i, j) = bar_Mc(i, j) - bar_Mk(i, j);
    }
  }

  // Yield function of hat_bar_Mc
  DenseMatrix<double> df_Mc_dMc(DIM, DIM, 0.0);
  double f_hat_bar_Mc =
    Plastic_consitutive_law_pt->compute_yield_surface_function(
      hat_bar_Mc, df_Mc_dMc, true);

  // hat_bar_Mc
  DenseMatrix<double> hat_bar_Nc;
  RankFourTensor<double> dhat_bar_Nc_dMc;
  compute_hat_bar_Nc(
    f_hat_bar_Mc, df_Mc_dMc, hat_bar_Nc, dhat_bar_Nc_dMc, flag);

  // Now compute Rc
  double Rc = f_hat_bar_Mc / isotropic_yield_stress;
  double dRcdH;
  DenseMatrix<double> dRc_dMc;

  if (flag)
  {
    dRcdH = Rc / isotropic_yield_stress * disotropic_f_dH;
    dRc_dMc.resize(DIM, DIM);
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        dRc_dMc(i, j) = df_Mc_dMc(i, j) / isotropic_yield_stress;
      }
    }
  }

  // Compute bar_Lpcs
  DenseMatrix<double> bar_Lpcd;
  RankFourTensor<double> dbar_Lpcd_dbar_M, dbar_Lpcd_dhat_bar_Mc;
  DenseMatrix<double> dbar_Lpcd_dh;
  compute_bar_Lpcd(bar_M,
                   hat_bar_Nc,
                   Rc,
                   dhat_bar_Nc_dMc,
                   dRc_dMc,
                   dRcdH,
                   bar_Lpcd,
                   dbar_Lpcd_dbar_M,
                   dbar_Lpcd_dhat_bar_Mc,
                   dbar_Lpcd_dh,
                   flag);

  // Compute dot_Fpcs = (bar_Lp - bar_Lpcd) Fpcs
  DenseMatrix<double> dotFpcs(DIM, DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      double sum = 0.0;
      for (unsigned int a = 0; a < DIM; a++)
      {
        sum += (bar_Lp(i, a) - bar_Lpcd(i, a)) * Fpcs(a, j);
      }
      dotFpcs(i, j) = sum;
    }
  }

  DenseMatrix<double> dot_or_delta_Fpcs_fromTimeStepper;
  get_dot_or_delta_fpcs_matrix(ipt, dot_or_delta_Fpcs_fromTimeStepper);


  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      residuals[this->plastic_fpcs_eqn_number(ipt, i, j)] =
        dot_or_delta_Fpcs_fromTimeStepper(i, j) -
        dot_or_delta_lambda * dotFpcs(i, j);
    }
  }

  if (flag)
  {
    // Precompute some things
    RankFourTensor<double> dbar_Lpcd_dinvFp(DIM, DIM, DIM, DIM),
      dbar_Lpcd_dFpks(DIM, DIM, DIM, DIM), dbar_Lpcd_dFpcs(DIM, DIM, DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            double sum_invFp = 0.0;
            double sum_Fpks = 0.0;
            double sum_Fpcs = 0.0;
            for (unsigned int p = 0; p < DIM; p++)
            {
              for (unsigned int q = 0; q < DIM; q++)
              {
                sum_invFp +=
                  dbar_Lpcd_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                sum_Fpks -=
                  dbar_Lpcd_dhat_bar_Mc(i, j, p, q) * dbar_Mk_dFpks(p, q, k, l);

                sum_Fpcs +=
                  dbar_Lpcd_dhat_bar_Mc(i, j, p, q) * dbar_Mc_dFpcs(p, q, k, l);
              }
            }
            dbar_Lpcd_dinvFp(i, j, k, l) = sum_invFp;
            dbar_Lpcd_dFpks(i, j, k, l) = sum_Fpks;
            dbar_Lpcd_dFpcs(i, j, k, l) = sum_Fpcs;
          }
        }
      }
    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        const unsigned int row_eq_fpcs =
          this->plastic_fpcs_eqn_number(ipt, i, j);

        // Contribution from the time stepper
        jacobian(row_eq_fpcs, row_eq_fpcs) += fpcs_time_step_weight;

        // Derivative wrt. invFp, Fpks, Fpcs
        for (unsigned int k = 0; k < DIM; k++)
        {
          // The delta term \delta_{jl} (Lp - Lpcd)_{ik}
          jacobian(row_eq_fpcs, this->plastic_fpcs_eqn_number(ipt, k, j)) -=
            dot_or_delta_lambda * (bar_Lp(i, k) - bar_Lpcd(i, k));

          for (unsigned int l = 0; l < DIM; l++)
          {
            const unsigned int col_eq_fpks =
              this->plastic_fpks_eqn_number(ipt, k, l);

            const unsigned int col_eq_invFp =
              this->plastic_inv_fp_eqn_number(ipt, k, l);

            const unsigned int col_eq_fpcs =
              this->plastic_fpcs_eqn_number(ipt, k, l);

            const unsigned int r_col = this->plastic_r_eqn_number(ipt);

            // Remainder FPcs term:
            // Fpcs_{aj} dbar_Lp_ia / dFpcs_kl
            // - Fpcs_{aj} dbar_Lpcs_ia / dFpcs_kl
            //
            // FPks term:
            // Fpcs_{aj} dbar_Lp_ia / dFpks_kl
            // - Fpcs_{aj} dbar_Lpcs_ia / dFpks_kl
            //
            // invFp term:
            // Fpcs_{aj} dbar_Lp_ia / dinvFp_kl
            // - Fpcs_{aj} dbar_Lpcs_ia / dinvFp_kl

            double fpcs_sum = 0.0;
            double fpks_sum = 0.0;
            double invFp_sum = 0.0;
            double sum_R = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              fpcs_sum += Fpcs(a, j) * (dbar_Lp_dFpcs(i, a, k, l) -
                                        dbar_Lpcd_dFpcs(i, a, k, l));

              fpks_sum += Fpcs(a, j) * (dbar_Lp_dFpks(i, a, k, l) -
                                        dbar_Lpcd_dFpks(i, a, k, l));

              invFp_sum += Fpcs(a, j) * (dbar_Lp_dinv_Fp(i, a, k, l) -
                                         dbar_Lpcd_dinvFp(i, a, k, l));

              sum_R += Fpcs(l, j) * dbar_Lp_dR(i, l);
            }

            jacobian(row_eq_fpcs, col_eq_fpcs) -=
              dot_or_delta_lambda * fpcs_sum;
            jacobian(row_eq_fpcs, col_eq_fpks) -=
              dot_or_delta_lambda * fpks_sum;
            jacobian(row_eq_fpcs, col_eq_invFp) -=
              dot_or_delta_lambda * invFp_sum;
            jacobian(row_eq_fpcs, r_col) -= dot_or_delta_lambda * sum_R;
          }
        }

        // Now the lambda contribution:
        // - dotFpcs(i, j) - dot_lambda dotFpcs_dh * dh_dlambda
        const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
        jacobian(row_eq_fpcs, lambda_col) -=
          dotFpcs(i, j) * dot_or_delta_lambda_time_step_weight;

        // The second part: (dot_Fpcs/dh)_{i, j} = - dLpcd_ia/h Fpcs_aj;
        double sum_lambda = 0.0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          sum_lambda += dbar_Lpcd_dh(i, a) * Fpcs(a, j);
        }
        jacobian(row_eq_fpcs, lambda_col) +=
          dot_or_delta_lambda * sum_lambda *
          this->Plastic_consitutive_law_pt->isotropic_hardening_factor;
      }
    }
  }

  // The remainder is relevant for R below.
  DenseMatrix<double> dC_dbar_M, dC_dbar_Mk, dC_dbar_Mc;
  double dCdR;

  double c_sigma = compute_c_sigma(barbar_N,
                                   hat_bar_Nc,
                                   dbarbarN_dbarbar_M,
                                   dbarbar_M_dMk,
                                   dbarbar_M_dMc,
                                   dbarbar_M_dR,
                                   dhat_bar_Nc_dMc,
                                   dC_dbar_M,
                                   dC_dbar_Mk,
                                   dC_dbar_Mc,
                                   dCdR,
                                   flag);

  u = compute_u(u,
                Rc,
                c_sigma,
                dRc_dMc,
                dRcdH,
                dC_dbar_M,
                dC_dbar_Mk,
                dC_dbar_Mc,
                dCdR,
                du_dbar_M,
                du_dbar_Mk,
                du_dbar_Mc,
                du_dh,
                du_dR,
                flag);

  //////////////////////////////////////////////////////////////////////////////
  // The residual for R //
  //////////////////////////////////////////////////////////////////////////////
  const double delta_lambda = get_delta_lambda(ipt);
  double R_prev = get_r(1, ipt);
  double dRdLambda;
  double dRdu;
  double computed_R =
    compute_r_plastic(u, delta_lambda, R_prev, dRdLambda, dRdu, flag);

  const unsigned int row_r = this->plastic_r_eqn_number(ipt);
  residuals[row_r] = get_r(ipt) - computed_R;

  if (flag)
  {
    jacobian(row_r, row_r) += 1 - dRdu * du_dR;

    // Contribution from lambda
    const unsigned lambda_col = this->plastic_lambda_eqn_number(ipt);
    jacobian(row_r, lambda_col) -=
      dRdLambda +
      dRdu * du_dh *
        this->Plastic_consitutive_law_pt->isotropic_hardening_factor;

    // Now the remaining derivatives of u:
    // If Fpcs has not been built, this would all be 0.
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        const unsigned int col_bar_M_ij =
          this->plastic_inv_fp_eqn_number(ipt, i, j);

        const unsigned int col_bar_Mk_ij =
          this->plastic_fpks_eqn_number(ipt, i, j);

        const unsigned int col_bar_Mc_ij =
          this->plastic_fpcs_eqn_number(ipt, i, j);

        double contrib_M = 0.0;
        double contrib_Mk = 0.0;
        double contrib_Mc = 0.0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          for (unsigned int b = 0; b < DIM; b++)
          {
            contrib_M += du_dbar_M(a, b) * dbar_M_dinv_Fp(a, b, i, j);
            contrib_Mk += du_dbar_Mk(a, b) * dbar_Mk_dFpks(a, b, i, j);
            contrib_Mc += du_dbar_Mc(a, b) * dbar_Mc_dFpcs(a, b, i, j);
          }
        }

        jacobian(row_r, col_bar_M_ij) -= dRdu * contrib_M;
        jacobian(row_r, col_bar_Mk_ij) -= dRdu * contrib_Mk;
        jacobian(row_r, col_bar_Mc_ij) -= dRdu * contrib_Mc;
      }
    }
  }
}

template<unsigned DIM>
void PlasticEquations<DIM>::compute_deformation_gradient_tensor(
  unsigned int t, unsigned int intpt, DenseMatrix<double>& F) const
{
  // Find out how many nodes there are
  const unsigned n_node = this->nnode();

  // Find out how many position types of dof there are
  const unsigned n_position_type = this->nnodal_position_type();

  // Set up memory for the shape functions
  Shape psi(n_node, n_position_type);
  DShape dpsidxi(n_node, n_position_type, DIM);

  // Call the derivatives of the shape functions
  double J = this->dshape_lagrangian_at_knot(intpt, psi, dpsidxi);

  // Initialise to zero
  F.initialise(0.0);

  // Calculate displacements and derivatives and lagrangian coordinates
  for (unsigned l = 0; l < n_node; l++)
  {
    // Check if the node has enough history for the requested time t.
    // If not, we will take the lagrangian positions to compute F.
    // \todo Really, we could also just use unity in that case.
    bool has_enough_history =
      (t < this->node_pt(l)->position_time_stepper_pt()->ntstorage());

    // Loop over positional dofs
    for (unsigned k = 0; k < n_position_type; k++)
    {
      // Loop over displacement components (deformed position)
      for (unsigned i = 0; i < DIM; i++)
      {
        double coordinate_val;

        if (has_enough_history)
        {
          coordinate_val = this->nodal_position_gen(t, l, k, i);
        }
        else
        {
          coordinate_val = this->lagrangian_position_gen(l, k, i);
        }

        // Loop over derivative directions to build F_ij = d(x_i)/d(Xi_j)
        for (unsigned j = 0; j < DIM; j++)
        {
          F(i, j) += coordinate_val * dpsidxi(l, k, j);
        }
      }
    }
  }
}