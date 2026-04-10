#include "solid/solid_plastic_elements.h"

using namespace oomph;

template<unsigned DIM>
const std::vector<std::string> PlasticEquationsBase<DIM>::Plastic_data_names{
  "Fe", "invBpks", "invBpcs", "H", "Lambda", "R"};

namespace oomph
{
  template class PlasticEquationsBase<2>;
  template class QPlasticPVDElement<2, 2>;

  template class PlasticEquationsBase<3>;
  template class QPlasticPVDElement<3, 2>;
} // namespace oomph

template<unsigned DIM>
void PlasticEquationsBase<DIM>::get_cauchy_stress(
  const unsigned& ipt, DenseMatrix<double>& sigma) const
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
  MatrixHelpers::multiply_transpose(S_PK2, F, tmp); // tmp = S * F^T

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

////////////////////////////////////////////////////////////////////////////////
// Mathematical functions involved in computing the plastic residuals or
// similar
////////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void PlasticEquationsBase<DIM>::
  compute_total_right_cauchy_green_deformation_tensor(
    const unsigned int& t,
    const unsigned int& ipt,
    DenseMatrix<double>& G) const
{
  Vector<double> s(DIM);

  // Assign the values of s
  for (unsigned i = 0; i < DIM; ++i)
  {
    s[i] = this->integral_pt()->knot(ipt, i);
  }

  // FIRST WE CALCULATE G (C)

  // Find out how many nodes there are
  unsigned n_node = this->nnode();

  // Find out how many positional dofs there are
  unsigned n_position_type = this->nnodal_position_type();

  // Set up memory for the shape functions
  Shape psi(n_node, n_position_type);
  DShape dpsidxi(n_node, n_position_type, DIM);

  // Call the derivatives of the shape functions (ignore Jacobian)
  (void)this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

  // Storage for Lagrangian coordinates (initialised to zero)
  Vector<double> interpolated_xi(DIM, 0.0);

  // Calculate interpolated values of the derivative of global position
  // wrt lagrangian coordinates
  DenseMatrix<double> interpolated_G(DIM);

  // Initialise to zero
  for (unsigned i = 0; i < DIM; i++)
  {
    for (unsigned j = 0; j < DIM; j++)
    {
      interpolated_G(i, j) = 0.0;
    }
  }

  // Calculate displacements and derivatives
  for (unsigned l = 0; l < n_node; l++)
  {
    // Check if the node has enough history for the requested time t.
    // If not, we will take the lagrangian positions to compute F.
    bool has_enough_history =
      (t < this->node_pt(l)->position_time_stepper_pt()->ntstorage());

    // Loop over positional dofs
    for (unsigned k = 0; k < n_position_type; k++)
    {
      // Loop over displacement components (deformed position)
      for (unsigned i = 0; i < DIM; i++)
      {
        // Calculate the lagrangian coordinates and the accelerations
        interpolated_xi[i] +=
          this->lagrangian_position_gen(l, k, i) * psi(l, k);

        // Loop over derivative directions
        for (unsigned j = 0; j < DIM; j++)
        {
          if (has_enough_history)
          {
            interpolated_G(j, i) +=
              this->nodal_position_gen(t, l, k, i) * dpsidxi(l, k, j);
          }
          else
          {
            interpolated_G(j, i) +=
              this->lagrangian_position_gen(l, k, i) * dpsidxi(l, k, j);
          }
        }
      }
    }
  }

  // Assign values of G
  for (unsigned i = 0; i < DIM; i++)
  {
    // Do upper half of matrix
    for (unsigned j = i; j < DIM; j++)
    {
      // Initialise G(i,j) to zero
      G(i, j) = 0.0;
      // Now calculate the dot product
      for (unsigned k = 0; k < DIM; k++)
      {
        G(i, j) += interpolated_G(i, k) * interpolated_G(j, k);
      }
    }
    // Matrix is symmetric so just copy lower half
    for (unsigned j = 0; j < i; j++)
    {
      G(i, j) = G(j, i);
    }
  }

  // Push forward G using the isotropic growth term
  double gamma;
  this->get_isotropic_growth(ipt, s, interpolated_xi, gamma);
  double diag_entry = pow(gamma, 2.0 / double(DIM));
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      G(i, j) = G(i, j) / diag_entry;
    }
  }
}

template<unsigned DIM>
void PlasticEquationsBase<DIM>::compute_deformation_gradient_tensor(
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

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_mandel_stress_elastic(
  const DenseMatrix<double>& invFp,
  const DenseMatrix<double>& C,
  DenseMatrix<double>& bar_M,
  RankFourTensor<double>& dbar_M_dinvFp,
  const bool& compute_derivative)
{
  // Compute GF^{p, -1} and Ce
  DenseMatrix<double> FtF_invFp(DIM, DIM, 0.0), Ce(DIM, DIM, 0.0);
  MatrixHelpers::multiply(C, invFp, FtF_invFp);
  MatrixHelpers::transpose_multiply(invFp, FtF_invFp, Ce);

  // Now compute bar_M
  DenseMatrix<double> S(DIM); // Second Piola–Kirchhoff stress tensors
  this->Constitutive_law_pt->calculate_second_piola_kirchhoff_stress(
    this->unity, Ce, S);

  // The mandel stress
  MatrixHelpers::multiply(Ce, S, bar_M);

  if (!compute_derivative) return;

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
void oomph::PlasticEquationsBase<DIM>::compute_mandellike_kinematic_hardening(
  const DenseMatrix<double>& invBpks,
  DenseMatrix<double>& bar_Mk,
  RankFourTensor<double>& dbar_Mk_dinvBpks,
  const bool& compute_derivative)
{
  bar_Mk.resize(DIM, DIM, 0.0);

  this->Plastic_constitutive_law_pt->kinematic_hardening_law_pt
    ->calculate_second_piola_kirchhoff_stress(invBpks, this->unity, bar_Mk);

  if (compute_derivative)
  {
    this->Plastic_constitutive_law_pt->kinematic_hardening_law_pt
      ->calculate_d_second_piola_kirchhoff_stress_dg(
        invBpks, this->unity, bar_Mk, dbar_Mk_dinvBpks);
  }
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_mandellike_elastic_core(
  const DenseMatrix<double>& invBpcs,
  DenseMatrix<double>& bar_Mc,
  RankFourTensor<double>& dbar_Mc_dinvBpcs,
  const bool& compute_derivative)
{
  bar_Mc.resize(DIM, DIM, 0.0);

  this->Plastic_constitutive_law_pt->elastic_core_law_pt
    ->calculate_second_piola_kirchhoff_stress(invBpcs, this->unity, bar_Mc);

  if (compute_derivative)
  {
    this->Plastic_constitutive_law_pt->elastic_core_law_pt
      ->calculate_d_second_piola_kirchhoff_stress_dg(
        invBpcs, this->unity, bar_Mc, dbar_Mc_dinvBpcs);
  }
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_mandel_stress_total(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& bar_Mk,
  const DenseMatrix<double>& bar_Mc,
  const double& r,
  DenseMatrix<double>& barbar_M,
  double& dbarbar_M_dMk,
  double& dbarbar_M_dMc,
  DenseMatrix<double>& dbarbar_M_dr,
  bool compute_jacobian)
{
  // Check the plastic model
  const bool solve_core =
    Plastic_model_type >= PlasticModel::ExtendedSubloadingSurface;
  const bool solve_r =
    Plastic_model_type >= PlasticModel::SubloadingSurface;

  barbar_M.resize(DIM, DIM);
  if (compute_jacobian) dbarbar_M_dr.resize(DIM, DIM);

  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      if (solve_core)
      {
        barbar_M(i, j) =
          bar_M(i, j) - r * bar_Mk(i, j) - (1 - r) * bar_Mc(i, j);
      }
      else
      {
        barbar_M(i, j) = bar_M(i, j) - bar_Mk(i, j);
      }

      if (compute_jacobian && solve_r)
      {
        if (solve_core)
        {
          dbarbar_M_dr(i, j) = -bar_Mk(i, j) + bar_Mc(i, j);
        }
        else
        {
          dbarbar_M_dr(i, j) = 0.0;
        }
      }
    }
  }

  if (compute_jacobian)
  {
    if (solve_core)
    {
      dbarbar_M_dMk = -r;
      dbarbar_M_dMc = r - 1.0;
    }
    else
    {
      dbarbar_M_dMk = -1.0;
      dbarbar_M_dMc = 0.0;
    }
  }

  return;
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_barbar_N(
  const DenseMatrix<double>& barbar_M,
  const double& f,
  const DenseMatrix<double>& dfdM,
  DenseMatrix<double>& barbar_N,
  RankFourTensor<double>& dbarbar_N_dbarbar_M,
  const bool& compute_derivative)
{
  barbar_N.resize(DIM, DIM, 0.0);

  MatrixHelpers::symmetrize(dfdM, barbar_N);
  double nMag = MatrixHelpers::magnitude(barbar_N);

  // Safety check for zero length
  if (nMag < 1.0e-15)
  {
    barbar_N.initialise(0.0);
    if (compute_derivative)
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

  if (!compute_derivative) return;
  // The derivative is
  // dN_{ij}/dM_{kl} = [\delta_{is}\delta_{jt} - N_{ij}N_{st}] *
  // [ ddf/(dM_{st}dM_kl) + ddf/(dM_{ts}dM_kl) ] / (2 nMag)

  // We split the computation in two loop blocks - this saves if statements and
  // loop iterations

  // First, compute N_{st} * [ ddf/(dM_{st}dM_kl) + ddf/(dM_{ts}dM_kl) ]
  RankFourTensor<double> ddfdMdM;
  this->Plastic_constitutive_law_pt->yield_criterion_pt
    ->surface_function_second_derivative(f, dfdM, ddfdMdM);
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
void oomph::PlasticEquationsBase<DIM>::compute_bar_Lp(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& barbar_N,
  const RankFourTensor<double>& dbarbar_N_dbarbar_M,
  DenseMatrix<double>& bar_Lp,
  RankFourTensor<double>& dbar_Lp_dbarbar_M,
  RankFourTensor<double>& dbar_Lp_dbar_M,
  const bool& compute_derivative)
{
  const double eta = (*this->Plastic_constitutive_law_pt->eta_p_pt);

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
  if (compute_derivative)
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

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_bar_Lpkd(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& bar_Mk,
  DenseMatrix<double>& bar_Lpkd,
  RankFourTensor<double>& dbar_Lpkd_dbar_M,
  RankFourTensor<double>& dbar_Lpkd_dbar_Mk,
  const bool& compute_derivative)
{
  // \bar{L}^\text{pkd} / \dot{\bar{\lambda}} =
  //                     bar_Mk + etapk [\bar{M}, \bar{M}^\text{k}] / b^\text{k}
  bar_Lpkd.resize(DIM, DIM);
  const double inv_bk =
    1.0 / (*this->Plastic_constitutive_law_pt->kinematic_hardening_b_pt);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      bar_Lpkd(i, j) = inv_bk * bar_Mk(i, j);
    }
  }

  const double eta =
    *this->Plastic_constitutive_law_pt->kinematic_hardening_eta_pt;
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

  if (compute_derivative)
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
void oomph::PlasticEquationsBase<DIM>::compute_hat_bar_Nc(
  const double& f_Mc,
  const DenseMatrix<double>& df_Mc_dMc,
  DenseMatrix<double>& hat_bar_Nc,
  RankFourTensor<double>& dhat_bar_Nc_dMc,
  const bool& compute_derivative)
{
  hat_bar_Nc.resize(DIM, DIM);

  double n_mag = MatrixHelpers::magnitude(df_Mc_dMc);
  // Safety check for zero length
  if (n_mag < 1.0e-15)
  {
    hat_bar_Nc.initialise(0.0);
    if (compute_derivative)
    {
      dhat_bar_Nc_dMc.resize(DIM, DIM, DIM, DIM, 0.0);
    }
    return;
  }


  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      hat_bar_Nc(i, j) = df_Mc_dMc(i, j) / n_mag;
    }
  }

  if (!compute_derivative) return;

  // The derivative is (A = df_Mc_dMc)
  // dN_ij/dA_mn = (\delta_im \delta_jn - N_ij N_mn) / nMag
  // For dN_ij_dhat_bar_Mc_{mn} one must just use the chain rule:
  // dN_ij/dhat_bar_Mc_{kl} = (dN_ij/dA_kl + N_ij N_mn dA_mn / dM_kl) / nMag
  RankFourTensor<double> dfdMcdMc;
  this->Plastic_constitutive_law_pt->yield_criterion_pt
    ->surface_function_second_derivative(f_Mc, df_Mc_dMc, dfdMcdMc);

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
  double n_mag_inv = 1 / n_mag;
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      for (unsigned int k = 0; k < DIM; k++)
      {
        for (unsigned int l = 0; l < DIM; l++)
        {
          dhat_bar_Nc_dMc(i, j, k, l) =
            n_mag_inv *
            (dfdMcdMc(i, j, k, l) - hat_bar_Nc(i, j) * Nc_ddfdMcdMc(k, l));
        }
      }
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::compute_bar_Lpcd(
  const DenseMatrix<double>& bar_M,
  const DenseMatrix<double>& hat_bar_Nc,
  const double& rc,
  const RankFourTensor<double>& dhat_bar_Nc_dhat_bar_Mc,
  const DenseMatrix<double>& drc_dMc,
  const double& drc_dh,
  DenseMatrix<double>& bar_Lpcd,
  RankFourTensor<double>& dbar_Lpcd_dbar_M,
  RankFourTensor<double>& dbar_Lpcd_dhat_bar_Mc,
  DenseMatrix<double>& dbar_Lpcd_dh,
  const bool& compute_derivative)
{
  bar_Lpcd.resize(DIM, DIM);

  const double inv_x =
    1 / (*this->Plastic_constitutive_law_pt->elastic_core_x_pt);
  const double rc_by_x = rc * inv_x;
  const double eta = (*this->Plastic_constitutive_law_pt->elastic_core_eta_pt);
  const double eta_prefactor = rc_by_x * eta;

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

      bar_Lpcd(i, j) = inv_x * sum;
    }
  }

  if (compute_derivative)
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
        dbar_Lpcd_dh(i, j) = bar_Lpcd(i, j) * drc_dh;

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
              rc_by_x * dNdMc + bar_Lpcd(i, j) * drc_dMc(k, l);
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
      bar_Lpcd(i, j) = rc * bar_Lpcd(i, j);
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::initialise_solve(const unsigned ipt)
{
  for (unsigned data_type: active_plastic_data_indices)
  {
    Data* data_pt = Plastic_data_pt[ipt][data_type];

    TimeStepper* data_time_stepper = data_pt->time_stepper_pt();

    // If there is not enough history, just initiallise all data with their
    // initial condition
    if (data_pt->time_stepper_pt()->ntstorage() < 2)
    {
      set_intial_condition(ipt);
      return;
    }
    // If it is steady initialise to the previous value
    else if (data_time_stepper->is_steady())
    {
      const unsigned nval = data_pt->nvalue();
      for (unsigned i = 0; i < nval; i++)
      {
        data_pt->set_value(i, data_pt->value(1, i));
      }
    }
    // Else, initialise such that the derivative is zero
    else
    {
      const unsigned nval = data_pt->nvalue();
      for (unsigned i = 0; i < nval; i++)
      {
        data_pt->set_value(
          i, compute_internal_data_from_time_derivative(data_pt, i, 0.0));
      }

      // Make sure that R has a valid value
      enforce_boundaries_of_r(ipt);
    }
  }
}

template<unsigned DIM>
void oomph::PlasticEquationsBase<DIM>::set_intial_condition(
  const unsigned int ipt)
{
  for (unsigned data_type: active_plastic_data_indices)
  {
    Data* data_pt = Plastic_data_pt[ipt][data_type];

    // Set the deformation gradient tensors to unity
    if (data_type == invFp_INDEX || data_type == invBpks_INDEX ||
        data_type == invBpcs_INDEX)
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
      double re = 0.0;
      if (this->Plastic_constitutive_law_pt)
      {
        re = this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt
               ->get_re();
      }
      data_pt->set_value(0, re);
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

template<unsigned DIM>
bool PlasticEquationsBase<DIM>::is_there_plastic_deformation(
  const unsigned int ipt)
{
  // This implements Eq. 193 from
  // https://doi.org/10.1007/s11831-018-9256-5
  //
  // For the initial subloading surface model, bar_Mc = bar_Mk is assumed
  // For the Conventional model, R = Re = 1.0 is assumed

  // Get the plastic quantities
  DenseMatrix<double> invFp(DIM, DIM, 0.0), invBpks(DIM, DIM, 0.0),
    invBpcs(DIM, DIM, 0.0);

  get_inv_fp_matrix(ipt, invFp);
  get_invBpks_matrix(ipt, invBpks);

  if (Plastic_model_type >= PlasticModel::ExtendedSubloadingSurface)
    get_invBpcs_matrix(ipt, invBpcs);

  // Get R
  double r = 1.0;
  if (Plastic_model_type >= PlasticModel::SubloadingSurface)
    r = get_r(ipt);

  DenseMatrix<double> C_total(DIM, DIM, 0.0);
  compute_total_right_cauchy_green_deformation_tensor(0, ipt, C_total);

  // Now compute the elastic Mandel stress
  DenseMatrix<double> bar_M(DIM, DIM, 0.0);
  compute_mandel_stress_elastic(invFp, C_total, bar_M);

  // Compute bar_Mk and bar_Mc
  DenseMatrix<double> bar_Mk(DIM, DIM, 0.0), bar_Mc(DIM, DIM, 0.0);
  compute_mandellike_kinematic_hardening(invBpks, bar_Mk);

  if (Plastic_model_type >= PlasticModel::ExtendedSubloadingSurface)
  {
    compute_mandellike_elastic_core(invBpcs, bar_Mc);
  }

  // Get the total mandel stress
  DenseMatrix<double> barbar_M(DIM, DIM, 0.0);
  // Takes care internally of the plasticity model
  compute_mandel_stress_total(bar_M, bar_Mk, bar_Mc, r, barbar_M);

  // Now compute the previous value - the only thing that has changed if C
  DenseMatrix<double> C_total_prev(DIM, DIM, 0.0);
  compute_total_right_cauchy_green_deformation_tensor(1, ipt, C_total_prev);

  // Previous bar_M
  DenseMatrix<double> bar_M_prev(DIM, DIM, 0.0);
  compute_mandel_stress_elastic(invFp, C_total_prev, bar_M_prev);

  // Compute Mbarbar_prev
  DenseMatrix<double> barbar_M_prev(DIM, DIM, 0.0);
  // Takes care internally of the plasticity model
  compute_mandel_stress_total(bar_M_prev, bar_Mk, bar_Mc, r, barbar_M_prev);

  // Now we can calculate deltaMTrial = bar_M - bar_M_prev
  DenseMatrix<double> delta_M_trial(DIM, DIM, 0.0);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      delta_M_trial(i, j) = bar_M(i, j) - bar_M_prev(i, j);
    }
  }

  // Compute Nbarbar_prev
  DenseMatrix<double> barbar_N_prev_nsym(DIM, DIM, 0.0);
  (void)this->Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
    barbar_M_prev, barbar_N_prev_nsym, true);
  MatrixHelpers::normalise(barbar_N_prev_nsym);

  double re = 1.0;
  if (Plastic_model_type >= PlasticModel::SubloadingSurface)
    re = this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt->get_re();

  // Compute Mbarbar_e = barbar_M - Re Mbark_prev - (1 - Re) Mbarc_prev
  DenseMatrix<double> Mbarbar_e(DIM, DIM, 0.0);
  if (Plastic_model_type >= PlasticModel::ExtendedSubloadingSurface)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        Mbarbar_e(i, j) =
          bar_M(i, j) - re * bar_Mk(i, j) - (1.0 - re) * bar_Mc(i, j);
      }
    }
  }
  else
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        Mbarbar_e(i, j) = bar_M(i, j) - bar_Mk(i, j);
      }
    }
  }

  double h_var = get_lambda(ipt) *
                 this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
                   ->isotropic_hardening_factor();
  bool plastic_deformation = false;
  if (MatrixHelpers::reduce(barbar_N_prev_nsym, delta_M_trial) >= 0.0)
  {
    // This is the case of forward loading

    // Check yield surface
    plastic_deformation =
      this->Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
        Mbarbar_e, GeneralisedElement::Dummy_matrix, false) -
        re * this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
               ->yield_function(h_var) >
      0;
  }
  else
  {
    // Inverse loading

    // Compute barbar_N
    DenseMatrix<double> barbar_N_nsym(DIM, DIM, 0.0);
    (void)this->Plastic_constitutive_law_pt->yield_criterion_pt
      ->surface_function(barbar_M, barbar_N_nsym, true);
    MatrixHelpers::normalise(barbar_N_nsym);

    if (MatrixHelpers::reduce(barbar_N_nsym, delta_M_trial) > 0.0)
    {
      plastic_deformation =
        this->Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
          Mbarbar_e, GeneralisedElement::Dummy_matrix, false) -
          re * this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
                 ->yield_function(h_var) >
        0;
    }
    else
    {
      plastic_deformation = false;
    }
  }

  if (!plastic_deformation &&
      Plastic_model_type > PlasticModel::Conventional)
  {
    // Compute R from the yield surface condition
    r = this->Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
          barbar_M, GeneralisedElement::Dummy_matrix, false) /
        this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
          ->yield_function(
            get_lambda(ipt) *
            this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
              ->isotropic_hardening_factor());

    set_r(ipt,
          std::max(r,
                   this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt
                     ->get_re()));
  }

  return plastic_deformation;
}

template<unsigned DIM>
void PlasticEquationsBase<DIM>::fill_in_generic_residual_and_jacobian_plastic(
  DoubleVector& residuals,
  DenseMatrix<double>& jacobian,
  const unsigned& ipt,
  const DenseMatrix<double>& C,
  const unsigned& flag)
{
  // Check the plastic model
  const bool solve_core =
    Plastic_model_type >= PlasticModel::ExtendedSubloadingSurface;
  const bool solve_r =
    Plastic_model_type >= PlasticModel::SubloadingSurface;

  // Retreive the derivative weights associated with the current value of a
  // plastic variable.
  double inv_Fp_time_step_weight = 1.0;
  double dot_or_delta_lambda_time_step_weight = 1.0;
  double invBpks_time_step_weight = 1.0;
  double invBpcs_time_step_weight = 1.0;

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
  if (!Plastic_data_pt[ipt][invBpks_INDEX]->time_stepper_pt()->is_steady())
  {
    invBpks_time_step_weight =
      Plastic_data_pt[ipt][invBpks_INDEX]->time_stepper_pt()->weight(1, 0);
  }


  if (solve_core &&
      !Plastic_data_pt[ipt][invBpcs_INDEX]->time_stepper_pt()->is_steady())
  {
    invBpcs_time_step_weight =
      Plastic_data_pt[ipt][invBpcs_INDEX]->time_stepper_pt()->weight(1, 0);
  }

  // Retreive dot_or_delta_lambda
  const double dot_or_delta_lambda = this->get_dot_or_delta_lambda(ipt);

  // Compute R
  double r = 1.0;
  if (solve_r) r = get_r(ipt);


  // Will be used to help with the computations
  DenseMatrix<double> tmp(DIM);

  // Get invFp
  DenseMatrix<double> invFp(DIM, DIM, 0.0);
  get_inv_fp_matrix(ipt, invFp);

  // Now compute the elastic Mandel stress
  DenseMatrix<double> bar_M(DIM, DIM, 0.0);
  RankFourTensor<double> dbar_M_dinv_Fp;
  compute_mandel_stress_elastic(invFp, C, bar_M, dbar_M_dinv_Fp, flag);

  // Get invBpks and Fpkc
  DenseMatrix<double> invBpks, invBpcs;

  get_invBpks_matrix(ipt, invBpks);

  if (solve_core) get_invBpcs_matrix(ipt, invBpcs);

  // Compute bar_Mk and bar_Mc (kinematic hardening and elastic core mandel-like
  // quantity)
  DenseMatrix<double> bar_Mk(DIM, DIM, 0.0), bar_Mc(DIM, DIM, 0.0);
  RankFourTensor<double> dbar_Mk_dinvBpks, dbar_Mc_dinvBpcs;
  compute_mandellike_kinematic_hardening(
    invBpks, bar_Mk, dbar_Mk_dinvBpks, flag);
  if (solve_core)
    compute_mandellike_elastic_core(invBpcs, bar_Mc, dbar_Mc_dinvBpcs, flag);

  DenseMatrix<double> barbar_M(DIM);
  double dbarbar_M_dMk, dbarbar_M_dMc;
  DenseMatrix<double> dbarbar_M_dr;
  // This function internally takes care of the plasticity model
  compute_mandel_stress_total(bar_M,
                              bar_Mk,
                              bar_Mc,
                              r,
                              barbar_M,
                              dbarbar_M_dMk,
                              dbarbar_M_dMc,
                              dbarbar_M_dr,
                              flag);

  // barbar_N
  DenseMatrix<double> barbar_N(DIM);
  DenseMatrix<double> dfdM(DIM);
  double yield_surface_stress =
    this->Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
      barbar_M, dfdM, true);

  RankFourTensor<double> dbarbarN_dbarbar_M;
  compute_barbar_N(
    barbar_M, yield_surface_stress, dfdM, barbar_N, dbarbarN_dbarbar_M, flag);

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
  RankFourTensor<double> dbar_Lp_dinvBpks;
  RankFourTensor<double> dbar_Lp_dinvBpcs;
  DenseMatrix<double> dbar_Lp_dr;

  if (flag)
  {
    dbar_Lp_dinv_Fp.resize(DIM, DIM, DIM, DIM);
    dbar_Lp_dinvBpks.resize(DIM, DIM, DIM, DIM);
    if (solve_r) dbar_Lp_dr.resize(DIM, DIM);
    if (solve_core) dbar_Lp_dinvBpcs.resize(DIM, DIM, DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            double sum_invFp = 0.0;
            double sum_invBpks = 0.0;
            double sum_invBpcs = 0.0;
            for (unsigned int p = 0; p < DIM; p++)
            {
              for (unsigned int q = 0; q < DIM; q++)
              {
                sum_invFp +=
                  dbar_Lp_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                const double dbar_Lp_dbarbar_M_ijpq =
                  dbar_Lp_dbarbar_M(i, j, p, q);

                sum_invBpks +=
                  dbar_Lp_dbarbar_M_ijpq * dbar_Mk_dinvBpks(p, q, k, l);

                if (solve_core)
                  sum_invBpcs +=
                    dbar_Lp_dbarbar_M_ijpq * dbar_Mc_dinvBpcs(p, q, k, l);
              }
            }
            dbar_Lp_dinv_Fp(i, j, k, l) = sum_invFp;

            dbar_Lp_dinvBpks(i, j, k, l) = sum_invBpks * dbarbar_M_dMk;

            if (solve_core)
              dbar_Lp_dinvBpcs(i, j, k, l) = sum_invBpcs * dbarbar_M_dMc;
          }
        }
      }
    }

    // now we do the R derivative: dbar_Lp_dbar_M : dbarbar_M_dR
    if (solve_r)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          double sum_invFp_ij = 0;
          for (unsigned int k = 0; k < DIM; k++)
          {
            for (unsigned int l = 0; l < DIM; l++)
            {
              sum_invFp_ij +=
                dbar_Lp_dbarbar_M(i, j, k, l) * dbarbar_M_dr(k, l);
            }
          }
          dbar_Lp_dr(i, j) = sum_invFp_ij;
        }
      }
    }
  }

  // Finally minus_dotinvFp
  DenseMatrix<double> minus_dotinvFp;
  MatrixHelpers::multiply(invFp, bar_Lp, minus_dotinvFp);

  DenseMatrix<double> dot_or_delta_invFp_from_time_stepper(DIM);
  get_dot_or_delta_inv_fp_matrix(ipt, dot_or_delta_invFp_from_time_stepper);

  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      residuals[this->plastic_inv_fp_eqn_number(ipt, i, j)] =
        dot_or_delta_invFp_from_time_stepper(i, j) +
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
        // d(inv_Fp bar_Lp)_ij/dinvFp_kl = inv_Fp_ia dL_{aj}/dinvBpks_{kl}
        for (unsigned int l = 0; l < DIM; l++)
        {
          // Now the remainder
          for (unsigned int k = 0; k < DIM; k++)
          {
            double sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              sum += invFp(i, a) * dbar_Lp_dinvBpks(a, j, k, l);
            }
            // Derivative wrt. invFp
            const unsigned int col_eq_invBpks =
              this->plastic_invBpks_eqn_number(ipt, k, l);

            jacobian(row_eq, col_eq_invBpks) += dot_or_delta_lambda * sum;
          }
        }

        // Contribution from
        // d(inv_Fp bar_Lp)_ij/dinvBpcs_kl = inv_Fp_ia, dL_{aj}/dinvBpcs_{kl}
        if (solve_core)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            // Now the remainder
            for (unsigned int k = 0; k < DIM; k++)
            {
              double sum = 0.0;
              for (unsigned int a = 0; a < DIM; a++)
              {
                sum += invFp(i, a) * dbar_Lp_dinvBpcs(a, j, k, l);
              }
              // Derivative wrt. invFp
              const unsigned int col_eq_invBpcs =
                this->plastic_invBpcs_eqn_number(ipt, k, l);

              jacobian(row_eq, col_eq_invBpcs) += dot_or_delta_lambda * sum;
            }
          }
        }


        // Derivative wrt. lambda
        const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
        jacobian(row_eq, lambda_col) +=
          minus_dotinvFp(i, j) * dot_or_delta_lambda_time_step_weight;

        // Derivative wrt. R
        if (solve_r)
        {
          const unsigned int R_col = this->plastic_r_eqn_number(ipt);
          double sum_r = 0.0;
          for (unsigned int a = 0; a < DIM; a++)
          {
            // invFp(i,a) * dbar_Lp_dR(a,j)
            sum_r += invFp(i, a) * dbar_Lp_dr(a, j);
          }
          jacobian(row_eq, R_col) += dot_or_delta_lambda * sum_r;
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // The residual for lambda and R                                            //
  //////////////////////////////////////////////////////////////////////////////
  double h_var = this->get_lambda(ipt) *
                 this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
                   ->isotropic_hardening_factor();

  const unsigned int row_lambda = this->plastic_lambda_eqn_number(ipt);

  double disotropic_f_dh = 0.0;

  double isotropic_yield_stress =
    this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
      ->yield_function(h_var, disotropic_f_dh, flag);
  residuals[row_lambda] = yield_surface_stress - r * isotropic_yield_stress;

  if (flag)
  {
    // Contribution from lambda
    jacobian(row_lambda, row_lambda) -=
      disotropic_f_dh *
      this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
        ->isotropic_hardening_factor() *
      r;


    // Contribution from invBpks:
    // df/dM_ab * dM_ab/dinvBpks_ij
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum_invFp_ij = 0;
        double sum_invBpks_ij = 0;
        double sum_invBpcs_ij = 0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          for (unsigned int b = 0; b < DIM; b++)
          {
            const double dfdM_ab = dfdM(a, b);
            sum_invFp_ij += dfdM_ab * dbar_M_dinv_Fp(a, b, i, j);

            sum_invBpks_ij += dfdM_ab * dbar_Mk_dinvBpks(a, b, i, j);

            if (solve_core)
              sum_invBpcs_ij += dfdM_ab * dbar_Mc_dinvBpcs(a, b, i, j);
          }
        }

        const unsigned int invFp_ij_col =
          this->plastic_inv_fp_eqn_number(ipt, i, j);
        jacobian(row_lambda, invFp_ij_col) += sum_invFp_ij;

        const unsigned int invBpks_ij_col =
          this->plastic_invBpks_eqn_number(ipt, i, j);
        jacobian(row_lambda, invBpks_ij_col) += sum_invBpks_ij * dbarbar_M_dMk;

        if (solve_core)
        {
          const unsigned int invBpcs_ij_col =
            this->plastic_invBpcs_eqn_number(ipt, i, j);
          jacobian(row_lambda, invBpcs_ij_col) += sum_invBpcs_ij * dbarbar_M_dMc;
        }
      }
    }

    if (solve_r)
    {
      const unsigned int r_col = this->plastic_r_eqn_number(ipt);
      jacobian(row_lambda, r_col) -= isotropic_yield_stress;

      // yieldSurfaceStress = f(barbarM) also depends on R
      double sum_df_dr = 0.0;
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          sum_df_dr += dfdM(i, j) * dbarbar_M_dr(i, j);
        }
      }
      jacobian(row_lambda, r_col) += sum_df_dr;
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // The residual for invBpks: -\dot{invBpks} = invBpks (\bar Lp - \bar Lpkd) //
  // + (\bar Lp - \bar Lpkd)^T invBpks                                        //
  //////////////////////////////////////////////////////////////////////////////

  DenseMatrix<double> bar_Lpkd(DIM, DIM);
  RankFourTensor<double> dbar_Lpkd_dbar_M, dbar_Lpkd_dbar_Mk;
  compute_bar_Lpkd(
    bar_M, bar_Mk, bar_Lpkd, dbar_Lpkd_dbar_M, dbar_Lpkd_dbar_Mk, flag);

  DenseMatrix<double> minus_dotinvBpks(DIM);
  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      double sum_ij = 0.0;
      for (unsigned int a = 0; a < DIM; a++)
      {
        sum_ij += invBpks(i, a) * (bar_Lp(a, j) - bar_Lpkd(a, j)) +
                  (bar_Lp(a, i) - bar_Lpkd(a, i)) * invBpks(a, j);
      }
      minus_dotinvBpks(i, j) = sum_ij;
    }
  }

  DenseMatrix<double> dot_or_delta_invBpks_from_time_stepper;
  get_dot_or_delta_invBpks_matrix(ipt, dot_or_delta_invBpks_from_time_stepper);

  for (unsigned int i = 0; i < DIM; i++)
  {
    for (unsigned int j = 0; j < DIM; j++)
    {
      residuals[this->plastic_invBpks_eqn_number(ipt, i, j)] =
        dot_or_delta_invBpks_from_time_stepper(i, j) +
        dot_or_delta_lambda * minus_dotinvBpks(i, j);
    }
  }

  if (flag)
  {
    // Precompute some things
    RankFourTensor<double> dbar_Lpkd_dinvFp(DIM, DIM, DIM, DIM),
      dbar_Lpkd_dinvBpks(DIM, DIM, DIM, DIM);

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            double sum_invFp = 0.0;
            double sum_invBpks = 0.0;
            for (unsigned int p = 0; p < DIM; p++)
            {
              for (unsigned int q = 0; q < DIM; q++)
              {
                sum_invFp +=
                  dbar_Lpkd_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                sum_invBpks +=
                  dbar_Lpkd_dbar_Mk(i, j, p, q) * dbar_Mk_dinvBpks(p, q, k, l);
              }
            }
            dbar_Lpkd_dinvFp(i, j, k, l) = sum_invFp;
            dbar_Lpkd_dinvBpks(i, j, k, l) = sum_invBpks;
          }
        }
      }
    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        const unsigned int row_eq_invBpks =
          this->plastic_invBpks_eqn_number(ipt, i, j);

        // Contribution from the time stepper
        jacobian(row_eq_invBpks, row_eq_invBpks) += invBpks_time_step_weight;

        // Derivative wrt. invBpks

        // The delta terms
        //  \delta_{ik}\delta{al} (Lp - Lpkd)_{aj}
        //  \delta_{ak}\delta{jl} (Lp - Lpkd)_{ai}
        for (unsigned int a = 0; a < DIM; a++)
        {
          // \delta_{ik}\delta{al} (Lp - Lpkd)_{aj}
          jacobian(row_eq_invBpks,
                   this->plastic_invBpks_eqn_number(ipt, i, a)) +=
            dot_or_delta_lambda * (bar_Lp(a, j) - bar_Lpkd(a, j));

          // \delta_{ak}\delta{jl} (Lp - Lpkd)_{ai}
          jacobian(row_eq_invBpks,
                   this->plastic_invBpks_eqn_number(ipt, a, j)) +=
            dot_or_delta_lambda * (bar_Lp(a, i) - bar_Lpkd(a, i));
        }

        // The remainder
        for (unsigned int k = 0; k < DIM; k++)
        {
          for (unsigned int l = 0; l < DIM; l++)
          {
            const unsigned int col_eq_invBpks =
              this->plastic_invBpks_eqn_number(ipt, k, l);

            const unsigned int col_eq_invFp =
              this->plastic_inv_fp_eqn_number(ipt, k, l);

            // Remainder invBpks term:
            // invBpks_{ia} dbar_Lp_ia / dinvBpks_kl
            // - invBpks_{ia} dbar_Lpks_ia / dinvBpks_kl
            // dbar_Lp_ai / dinvBpks_kl * invBpks_{aj}
            // - dbar_Lpks_ai / dinvBpks_kl * invBpks_{aj}
            //
            // invFp term:
            // invBpks_{ia} dbar_Lp_ia / dinvFp_kl
            // - invBpks_{ia} dbar_Lpks_ia / dinvFp_kl
            // dbar_Lp_ai / dinvFp_kl invBpks_{aj}
            // - dbar_Lpks_ai / dinvFp_kl invBpks_{aj}
            //
            // invBpcs term:
            // invBpks_{ia} dbar_Lp_{ia} / dinvBpcs_{kl}
            // + invBpks_{ia} * dbar_Lpcd_{ia} / dinvBpcs_{kl}
            // dbar_Lp_{ai} / dinvBpcs_{kl} * invBpks_{aj}
            // - dbar_Lpcd_{ai} / dinvBpcs_{kl} * invBpks_{aj}
            double invBpks_sum = 0.0;
            double invFp_sum = 0.0;
            double invBpcs_sum = 0.0;
            for (unsigned int a = 0; a < DIM; a++)
            {
              invBpks_sum += invBpks(i, a) * (dbar_Lp_dinvBpks(a, j, k, l) -
                                              dbar_Lpkd_dinvBpks(a, j, k, l));
              invBpks_sum += (dbar_Lp_dinvBpks(a, i, k, l) -
                              dbar_Lpkd_dinvBpks(a, i, k, l)) *
                             invBpks(a, j);

              invFp_sum += invBpks(i, a) * (dbar_Lp_dinv_Fp(a, j, k, l) -
                                            dbar_Lpkd_dinvFp(a, j, k, l));

              invFp_sum += invBpks(a, j) * (dbar_Lp_dinv_Fp(a, i, k, l) -
                                            dbar_Lpkd_dinvFp(a, i, k, l));

              if (solve_core)
              {
                invBpcs_sum += invBpks(i, a) * (dbar_Lp_dinvBpcs(a, j, k, l));
                invBpcs_sum += invBpks(a, j) * (dbar_Lp_dinvBpcs(a, i, k, l));
              }
            }
            jacobian(row_eq_invBpks, col_eq_invBpks) +=
              dot_or_delta_lambda * invBpks_sum;
            jacobian(row_eq_invBpks, col_eq_invFp) +=
              dot_or_delta_lambda * invFp_sum;
            if (solve_core)
            {
              unsigned int col_eq_invBpcs = 0;
              col_eq_invBpcs = this->plastic_invBpcs_eqn_number(ipt, k, l);
              jacobian(row_eq_invBpks, col_eq_invBpcs) +=
                dot_or_delta_lambda * invBpcs_sum;
            }
          }
        }


        // Derivative wrt. lambda
        const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
        jacobian(row_eq_invBpks, lambda_col) +=
          minus_dotinvBpks(i, j) * dot_or_delta_lambda_time_step_weight;

        // Derivative wrt. R
        if (solve_r)
        {
          const unsigned int R_col = this->plastic_r_eqn_number(ipt);
          double sum_r = 0.0;
          for (unsigned int l = 0; l < DIM; l++)
          {
            sum_r += invBpks(i, l) * dbar_Lp_dr(l, j) +
                     dbar_Lp_dr(l, i) * invBpks(l, j);
          }
          jacobian(row_eq_invBpks, R_col) += dot_or_delta_lambda * sum_r;
        }
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // The residual for invBpcs: -\dot{invBpcs} = invBpcs (\bar Lp - \bar Lpkd) //
  // + (\bar Lp - \bar Lpkd)^T invBpcs                                        //
  //////////////////////////////////////////////////////////////////////////////
  double yield_ratio_u =
    this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt->get_u();
  DenseMatrix<double> du_dbar_M, du_dbar_Mk, du_dbar_Mc;
  double du_dh = 0, du_dr = 0;

  if (flag)
  {
    du_dbar_M.resize(DIM, DIM, 0.0);
    du_dbar_Mk.resize(DIM, DIM, 0.0);
    if (solve_core) du_dbar_Mc.resize(DIM, DIM, 0.0);
  }

  if (solve_core)
  {
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
      Plastic_constitutive_law_pt->yield_criterion_pt->surface_function(
        hat_bar_Mc, df_Mc_dMc, true);

    // hat_bar_Mc
    DenseMatrix<double> hat_bar_Nc;
    RankFourTensor<double> dhat_bar_Nc_dMc;
    compute_hat_bar_Nc(
      f_hat_bar_Mc, df_Mc_dMc, hat_bar_Nc, dhat_bar_Nc_dMc, flag);

    // Now compute Rc
    double rc = f_hat_bar_Mc / isotropic_yield_stress;
    double drcdh = 0;
    DenseMatrix<double> drc_dMc;

    if (flag)
    {
      drcdh = rc / isotropic_yield_stress * disotropic_f_dh;
      drc_dMc.resize(DIM, DIM);
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          drc_dMc(i, j) = df_Mc_dMc(i, j) / isotropic_yield_stress;
        }
      }
    }

    // Compute bar_Lpcs
    DenseMatrix<double> bar_Lpcd;
    RankFourTensor<double> dbar_Lpcd_dbar_M, dbar_Lpcd_dhat_bar_Mc;
    DenseMatrix<double> dbar_Lpcd_dh;
    compute_bar_Lpcd(bar_M,
                     hat_bar_Nc,
                     rc,
                     dhat_bar_Nc_dMc,
                     drc_dMc,
                     drcdh,
                     bar_Lpcd,
                     dbar_Lpcd_dbar_M,
                     dbar_Lpcd_dhat_bar_Mc,
                     dbar_Lpcd_dh,
                     flag);

    // Compute dot_invBpcs = (bar_Lp - bar_Lpcd) invBpcs
    DenseMatrix<double> minus_dotinvBpcs(DIM, DIM);
    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        double sum = 0.0;
        for (unsigned int a = 0; a < DIM; a++)
        {
          sum += invBpcs(i, a) * (bar_Lp(a, j) - bar_Lpcd(a, j)) +
                 (bar_Lp(a, i) - bar_Lpcd(a, i)) * invBpcs(a, j);
        }
        minus_dotinvBpcs(i, j) = sum;
      }
    }

    DenseMatrix<double> dot_or_delta_invBpcs_from_time_stepper;
    get_dot_or_delta_invBpcs_matrix(ipt, dot_or_delta_invBpcs_from_time_stepper);


    for (unsigned int i = 0; i < DIM; i++)
    {
      for (unsigned int j = 0; j < DIM; j++)
      {
        residuals[this->plastic_invBpcs_eqn_number(ipt, i, j)] =
          dot_or_delta_invBpcs_from_time_stepper(i, j) +
          dot_or_delta_lambda * minus_dotinvBpcs(i, j);
      }
    }

    if (flag)
    {
      // Precompute some things
      RankFourTensor<double> dbar_Lpcd_dinvFp(DIM, DIM, DIM, DIM),
        dbar_Lpcd_dinvBpks(DIM, DIM, DIM, DIM),
        dbar_Lpcd_dinvBpcs(DIM, DIM, DIM, DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          for (unsigned int k = 0; k < DIM; k++)
          {
            for (unsigned int l = 0; l < DIM; l++)
            {
              double sum_invFp = 0.0;
              double sum_invBpks = 0.0;
              double sum_invBpcs = 0.0;
              for (unsigned int p = 0; p < DIM; p++)
              {
                for (unsigned int q = 0; q < DIM; q++)
                {
                  sum_invFp +=
                    dbar_Lpcd_dbar_M(i, j, p, q) * dbar_M_dinv_Fp(p, q, k, l);

                  sum_invBpks -= dbar_Lpcd_dhat_bar_Mc(i, j, p, q) *
                                 dbar_Mk_dinvBpks(p, q, k, l);

                  sum_invBpcs += dbar_Lpcd_dhat_bar_Mc(i, j, p, q) *
                                 dbar_Mc_dinvBpcs(p, q, k, l);
                }
              }
              dbar_Lpcd_dinvFp(i, j, k, l) = sum_invFp;
              dbar_Lpcd_dinvBpks(i, j, k, l) = sum_invBpks;
              dbar_Lpcd_dinvBpcs(i, j, k, l) = sum_invBpcs;
            }
          }
        }
      }

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          const unsigned int row_eq_invBpcs =
            this->plastic_invBpcs_eqn_number(ipt, i, j);

          // Contribution from the time stepper
          jacobian(row_eq_invBpcs, row_eq_invBpcs) += invBpcs_time_step_weight;

          // Derivative wrt. invFp, invBpks, invBpcs
          // The delta terms
          //  \delta_{ik}\delta{al} (Lp - Lpcd)_{aj}
          //  \delta_{ak}\delta{jl} (Lp - Lpcd)_{ai}
          for (unsigned int a = 0; a < DIM; a++)
          {
            // \delta_{ik}\delta{al} (Lp - Lpcd)_{aj}
            jacobian(row_eq_invBpcs,
                     this->plastic_invBpcs_eqn_number(ipt, i, a)) +=
              dot_or_delta_lambda * (bar_Lp(a, j) - bar_Lpcd(a, j));

            // \delta_{ak}\delta{jl} (Lp - Lpkd)_{ai}
            jacobian(row_eq_invBpcs,
                     this->plastic_invBpcs_eqn_number(ipt, a, j)) +=
              dot_or_delta_lambda * (bar_Lp(a, i) - bar_Lpcd(a, i));
          }

          // The remainder
          for (unsigned int k = 0; k < DIM; k++)
          {
            for (unsigned int l = 0; l < DIM; l++)
            {
              const unsigned int col_eq_invBpks =
                this->plastic_invBpks_eqn_number(ipt, k, l);

              const unsigned int col_eq_invFp =
                this->plastic_inv_fp_eqn_number(ipt, k, l);

              const unsigned int col_eq_invBpcs =
                this->plastic_invBpcs_eqn_number(ipt, k, l);

              const unsigned int r_col = this->plastic_r_eqn_number(ipt);

              // Remainder invBpcs term:
              // invBpcs_{aj} dbar_Lp_ia / dinvBpcs_kl
              // - invBpcs_{aj} dbar_Lpcs_ia / dinvBpcs_kl
              // dbar_Lp_ai / dinvBpcs_kl * invBpcs_{aj}
              // - invBpcs__ai / dinvBpcs_kl * invBpcs_{aj}
              //
              // invBpks term:
              // invBpcs_{aj} dbar_Lp_ia / dinvBpks_kl
              // - invBpcs_{aj} dbar_Lpcs_ia / dinvBpks_kl
              //
              // invFp term:
              // invBpcs_{aj} dbar_Lp_ia / dinvFp_kl
              // - invBpcs_{aj} dbar_Lpcs_ia / dinvFp_kl

              double invBpcs_sum = 0.0;
              double invBpks_sum = 0.0;
              double invFp_sum = 0.0;
              for (unsigned int a = 0; a < DIM; a++)
              {
                invBpcs_sum += invBpcs(i, a) * (dbar_Lp_dinvBpcs(a, j, k, l) -
                                                dbar_Lpcd_dinvBpcs(a, j, k, l));
                invBpcs_sum += (dbar_Lp_dinvBpcs(a, i, k, l) -
                                dbar_Lpcd_dinvBpcs(a, i, k, l)) *
                               invBpcs(a, j);

                invBpks_sum += invBpcs(i, a) * (dbar_Lp_dinvBpks(a, j, k, l) -
                                                dbar_Lpcd_dinvBpks(a, j, k, l));
                invBpks_sum += (dbar_Lp_dinvBpks(a, i, k, l) -
                                dbar_Lpcd_dinvBpks(a, i, k, l)) *
                               invBpcs(a, j);

                invFp_sum += invBpcs(i, a) * (dbar_Lp_dinv_Fp(a, j, k, l) -
                                              dbar_Lpcd_dinvFp(a, j, k, l));
                invFp_sum +=
                  (dbar_Lp_dinv_Fp(a, i, k, l) - dbar_Lpcd_dinvFp(a, i, k, l)) *
                  invBpcs(a, j);
              }

              jacobian(row_eq_invBpcs, col_eq_invBpcs) +=
                dot_or_delta_lambda * invBpcs_sum;
              jacobian(row_eq_invBpcs, col_eq_invBpks) +=
                dot_or_delta_lambda * invBpks_sum;
              jacobian(row_eq_invBpcs, col_eq_invFp) +=
                dot_or_delta_lambda * invFp_sum;
            }
          }

          // Derivative wrt. R
          const unsigned int r_col = this->plastic_r_eqn_number(ipt);
          double sum_r = 0.0;
          for (unsigned int a = 0; a < DIM; a++)
          {
            sum_r += invBpcs(i, a) * dbar_Lp_dr(a, j) +
                     dbar_Lp_dr(a, i) * invBpcs(a, j);
          }
          jacobian(row_eq_invBpcs, r_col) += dot_or_delta_lambda * sum_r;

          // Now the lambda contribution:
          // - dotinvBpcs(i, j) - dot_lambda dotinvBpcs_dh * dh_dlambda
          const unsigned int lambda_col = this->plastic_lambda_eqn_number(ipt);
          jacobian(row_eq_invBpcs, lambda_col) +=
            minus_dotinvBpcs(i, j) * dot_or_delta_lambda_time_step_weight;

          // The second part: (dot_invBpcs/dh)_{i, j} = - dLpcd_ia/h invBpcs_aj;
          double sum_lambda = 0.0;
          for (unsigned int a = 0; a < DIM; a++)
          {
            sum_lambda += invBpcs(i, a) * dbar_Lpcd_dh(a, j) +
                          dbar_Lpcd_dh(a, i) * invBpcs(a, j);
          }
          jacobian(row_eq_invBpcs, lambda_col) -=
            dot_or_delta_lambda * sum_lambda *
            this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
              ->isotropic_hardening_factor();
        }
      }
    }

    // The remainder is relevant for R below.
    yield_ratio_u = this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt
                      ->compute_u_with_elastic_core(rc,
                                                    barbar_N,
                                                    hat_bar_Nc,
                                                    dbarbarN_dbarbar_M,
                                                    dbarbar_M_dMk,
                                                    dbarbar_M_dMc,
                                                    dbarbar_M_dr,
                                                    dhat_bar_Nc_dMc,
                                                    drc_dMc,
                                                    drcdh,
                                                    du_dbar_M,
                                                    du_dbar_Mk,
                                                    du_dbar_Mc,
                                                    du_dh,
                                                    du_dr,
                                                    flag);
  }

  //////////////////////////////////////////////////////////////////////////////
  // The residual for R //
  //////////////////////////////////////////////////////////////////////////////
  if (solve_r)
  {
    const double delta_lambda = get_delta_lambda(ipt);
    double r_prev = get_r(1, ipt);
    double dr_dlambda;
    double drdu;
    double computed_r =
      this->Plastic_constitutive_law_pt->normal_yield_ratio_law_pt
        ->compute_r_plastic(
          yield_ratio_u, delta_lambda, r_prev, dr_dlambda, drdu, flag);

    const unsigned int row_r = this->plastic_r_eqn_number(ipt);
    residuals[row_r] = get_r(ipt) - computed_r;

    if (flag)
    {
      jacobian(row_r, row_r) += 1 - drdu * du_dr;

      // Contribution from lambda
      const unsigned lambda_col = this->plastic_lambda_eqn_number(ipt);
      jacobian(row_r, lambda_col) -=
        dr_dlambda +
        drdu * du_dh *
          this->Plastic_constitutive_law_pt->isotropic_hardening_law_pt
            ->isotropic_hardening_factor();

      // Now the remaining derivatives of u:
      // If invBpcs has not been built, this would all be 0.
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          const unsigned int col_bar_M_ij =
            this->plastic_inv_fp_eqn_number(ipt, i, j);

          const unsigned int col_bar_Mk_ij =
            this->plastic_invBpks_eqn_number(ipt, i, j);

          double contrib_M = 0.0;
          double contrib_Mk = 0.0;
          double contrib_Mc = 0.0;
          for (unsigned int a = 0; a < DIM; a++)
          {
            for (unsigned int b = 0; b < DIM; b++)
            {
              contrib_M += du_dbar_M(a, b) * dbar_M_dinv_Fp(a, b, i, j);
              contrib_Mk += du_dbar_Mk(a, b) * dbar_Mk_dinvBpks(a, b, i, j);
              if (solve_core)
                contrib_Mc += du_dbar_Mc(a, b) * dbar_Mc_dinvBpcs(a, b, i, j);
            }
          }

          jacobian(row_r, col_bar_M_ij) -= drdu * contrib_M;
          jacobian(row_r, col_bar_Mk_ij) -= drdu * contrib_Mk;


          if (solve_core)
          {
            const unsigned int col_bar_Mc_ij = this->plastic_invBpcs_eqn_number(ipt, i, j);
            jacobian(row_r, col_bar_Mc_ij) -= drdu * contrib_Mc;
          }
        }
      }
    }
  }
}