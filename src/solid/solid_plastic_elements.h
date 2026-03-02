#ifndef OOMPH_PLASTICITY_SOLID_ELEMENTS_HEADER
#define OOMPH_PLASTICITY_SOLID_ELEMENTS_HEADER

#include "solid_elements.h"
#include "constitutive/plastic_constitutive_laws.h"
#include "refineable_solid_elements.h"
#include "generic/interpolate_from_integral_points.h"

// We initially implement the simplified plasticity model with no elastic core
// and no plastic dissipation

// At some point it may be worth using composition to store the plastic data
// within a sub element This will give us the ability to manage the data with
// the oomph-lib functions such as add_internal_data and
// assign_internal_eqn_numbers etc without the normal oomph-lib behaviour
// trampling over it all. This will be useful because we can then delete a good
// portion of the code below and it will allow us to handle pinning and
// unpinning our plastic data without having to keep track of the eqn numbers
// ourselves. The downside is that we will no longer benefit from the oomph-lib
// automatic time shifting.

namespace oomph
{
  template<unsigned DIM>
  class PlasticEquationsBase : public virtual InterpolateFromIntegralPointsBase,
                               public virtual PVDEquationsBase<DIM>
  {
  protected:
    // If finite difference should be used for the plastic solve;
    bool Plastic_solve_use_fd = false;

    // A pointer to the plastic constitutive law
    PlasticConstitutiveLaw* Plastic_consitutive_law_pt;

    // We use an enum to define the indices in the vectors at which the
    // different types of plastic data are stored. This simplifies things
    // considerably when searching for plastic data.
    enum Plastic_Variables_Indexes
    {
      invFp_INDEX,
      invBpks_INDEX,
      Fpcs_INDEX,
      Lambda_INDEX,
      R_INDEX,
      NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE,
      H_INDEX = NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE,
      NUMBER_OF_PLASTIC_VARIABLE_TYPES
    };

  private:
    const static std::vector<std::string> Plastic_data_names;
    // We store a vector of indices of the plastic data in internal data
    // so we can assign timesteppers etc more easily
    // We resize this every time we build a new set of plastic internal data but
    // since we only do this once when the element is constructed this shouldn't
    // be a problem - we can't do this at construction of this class because
    // the number of integral points is defined by derived classes

    // We store a separate set of plastic: internal data indices, pinned status,
    // and eqn numbers per integral point in the element.

    // [Number of ipts, number of types of plastic data]
    // Pointer to the plastic data at the given integral point and of the given
    // type.
    Vector<Vector<Data*>> Plastic_data_pt;

    // [Number of ipts, number of types of plastic data, number of variables in
    // that data type] Different to the data pinned status. That is always true
    // because we don't want to use the data in the normal way.
    Vector<Vector<std::vector<bool>>> Plastic_data_pinned_status;

    // [Number of ipts, number of types of plastic data, number of variables in
    // that data type]
    Vector<Vector<Vector<int>>> Plastic_data_eqn_number;

    // [Number of ipts, number of data to be finite differenced]
    // We store pointers to all the plastic data which are to be finite
    // differenced
    Vector<Vector<double*>> Plastic_dof_data_pt;

    // Keeps track of if data for the  plastic dof numbers has been allocated
    // We do it within this function because this allows us to reliably zero the
    // counters Num_plastic_Dofs and Num_plastic_residuals
    bool Plastic_dof_nunbers_has_been_resized = false;

    // The solver tolerance for the plastic newton solve
    double Plastic_Newton_Solver_Tolerance = 1.0e-8;

    // Pin the plastic dof of type data_type at the integral point, if index<0
    // then pin all plastic data of that type, otherwise pin the value specified
    void pin_plastic_dof(const unsigned& ipt,
                         const unsigned& data_type,
                         const int& index = -1)
    {
      if (index < 0)
      {
        const unsigned nvalue =
          Plastic_data_pinned_status[ipt][data_type].size();
        for (unsigned i = 0; i < nvalue; i++)
        {
          Plastic_data_pinned_status[ipt][data_type][i] = true;
        }
      }
      else
      {
#ifdef PARANOID
        // Check if the index makes sense
        if (index > Plastic_data_pinned_status[ipt][data_type].size() - 1)
        {
          throw OomphLibError("Plastic data index is too large",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
        Plastic_data_pinned_status[ipt][data_type][index] = true;
      }
    }

    void unpin_plastic_dof(const unsigned& ipt,
                           const unsigned& data_type,
                           const int& index = -1)
    {
      if (index < 0)
      {
        const unsigned nvalue =
          Plastic_data_pinned_status[ipt][data_type].size();
        for (unsigned i = 0; i < nvalue; i++)
        {
          Plastic_data_pinned_status[ipt][data_type][i] = false;
        }
      }
      else
      {
#ifdef PARANOID
        // Check if the index makes sense
        if (index > Plastic_data_pinned_status[ipt][data_type].size() - 1)
        {
          throw OomphLibError("Plastic data index is too large",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
        }
#endif
        Plastic_data_pinned_status[ipt][data_type][index] = false;
      }
    }

    // Change to not have an argument - in fact if we're always building all
    // plastic data then we don't need thi Collapse all building plastic data
    // into a single function
    void resize_plastic_dof_numbers()
    {
      if (Plastic_dof_nunbers_has_been_resized) return;

      const unsigned nipt = this->integral_pt()->nweight();

      Plastic_data_pt.resize(nipt);
      Plastic_data_pinned_status.resize(nipt);
      Plastic_data_eqn_number.resize(nipt);
      Plastic_dof_data_pt.resize(nipt);
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Plastic_data_pt[ipt].resize(NUMBER_OF_PLASTIC_VARIABLE_TYPES);
        Plastic_data_pinned_status[ipt].resize(
          NUMBER_OF_PLASTIC_VARIABLE_TYPES);
        Plastic_data_eqn_number[ipt].resize(NUMBER_OF_PLASTIC_VARIABLE_TYPES);
      }
      Plastic_dof_nunbers_has_been_resized = true;
    }

  protected:
    // We need a dummy double matrix for the plastic residual fill in
    DenseMatrix<double> unity;

    // Assign plastic eqn numbers for each of the integral points individually
    // Add pointers to all values which are not pinned to Plastic_dof_data_pt
    void assign_plastic_eqn_numbers()
    {
      for (unsigned ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        Plastic_dof_data_pt[ipt].clear();

        unsigned eqn_count = 0;
        for (unsigned data_type = 0;
             data_type < NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE;
             data_type++)
        {
          Data* data_pt = Plastic_data_pt[ipt][data_type];
          const unsigned nvalue = data_pt->nvalue();
          for (unsigned i = 0; i < nvalue; i++)
          {
            Plastic_data_eqn_number[ipt][data_type][i] = -1;
            if (!Plastic_data_pinned_status[ipt][data_type][i])
            {
              Plastic_data_eqn_number[ipt][data_type][i] = eqn_count++;
              Plastic_dof_data_pt[ipt].push_back(data_pt->value_pt(i));
            }
          }
        }
      }
    }

    void construct_inv_fp_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(DIM * DIM);
        Plastic_data_pt[ipt][invFp_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        for (unsigned i = 0; i < DIM * DIM; i++)
        {
          // Pin the plastic degree of freedom
          Plastic_data_pt[ipt][invFp_INDEX]->pin(i);
          // By default the plastic data is not pinned
          Plastic_data_pinned_status[ipt][invFp_INDEX].push_back(false);
          // But we have to initialise the eqn number to something so it may as
          // well be a safe value
          Plastic_data_eqn_number[ipt][invFp_INDEX].push_back(-1);
          // We initialise the plastic deformation gradient tensors to the
          // identity
          if (i % (DIM + 1) == 0)
          {
            data_pt->set_value(i, 1.0);
          }
        }
      }
    }

    void construct_invBpks_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(DIM * DIM);
        Plastic_data_pt[ipt][invBpks_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        for (unsigned i = 0; i < DIM * DIM; i++)
        {
          // Pin the plastic degree of freedom
          Plastic_data_pt[ipt][invBpks_INDEX]->pin(i);
          // By default the plastic data is not pinned
          Plastic_data_pinned_status[ipt][invBpks_INDEX].push_back(false);
          // But we have to initialise the eqn number to something so it may as
          // well be a safe value
          Plastic_data_eqn_number[ipt][invBpks_INDEX].push_back(-1);
          // We initialise the plastic deformation gradient tensors to the
          // identity
          if (i % (DIM + 1) == 0)
          {
            data_pt->set_value(i, 1.0);
          }
        }
      }
    }

    void construct_fpcs_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(DIM * DIM);
        Plastic_data_pt[ipt][Fpcs_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        for (unsigned i = 0; i < DIM * DIM; i++)
        {
          // Pin the plastic degree of freedom
          Plastic_data_pt[ipt][Fpcs_INDEX]->pin(i);
          // By default the plastic data is not pinned
          Plastic_data_pinned_status[ipt][Fpcs_INDEX].push_back(false);
          // But we have to initialise the eqn number to something so it may as
          // well be a safe value
          Plastic_data_eqn_number[ipt][Fpcs_INDEX].push_back(-1);
          // We initialise the plastic deformation gradient tensors to the
          // identity
          if (i % (DIM + 1) == 0)
          {
            data_pt->set_value(i, 1.0);
          }
        }
      }
    }

    void construct_h_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(1);
        Plastic_data_pt[ipt][H_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        // Pin the plastic degree of freedom
        Plastic_data_pt[ipt][H_INDEX]->pin(0);
        // By default the plastic data is not pinned
        Plastic_data_pinned_status[ipt][H_INDEX].push_back(false);
        // But we have to initialise the eqn number to something so it may as
        // well be a safe value
        Plastic_data_eqn_number[ipt][H_INDEX].push_back(-1);
        // What default value to set?
        data_pt->set_value(0, 0.0);
      }
    }

    void construct_r_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(1);
        Plastic_data_pt[ipt][R_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        // Pin the plastic degree of freedom
        Plastic_data_pt[ipt][R_INDEX]->pin(0);
        // By default the plastic data is not pinned
        Plastic_data_pinned_status[ipt][R_INDEX].push_back(false);
        // But we have to initialise the eqn number to something so it may as
        // well be a safe value
        Plastic_data_eqn_number[ipt][R_INDEX].push_back(-1);
        // What default value to set?
        data_pt->set_value(0, 0.0);
      }
    }

    void construct_lambda_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(1);
        Plastic_data_pt[ipt][Lambda_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        // Pin the plastic degree of freedom
        Plastic_data_pt[ipt][Lambda_INDEX]->pin(0);
        // By default the plastic data is not pinned
        Plastic_data_pinned_status[ipt][Lambda_INDEX].push_back(false);
        // But we have to initialise the eqn number to something so it may as
        // well be a safe value
        Plastic_data_eqn_number[ipt][Lambda_INDEX].push_back(-1);
        data_pt->set_value(0, 0.0);
      }
    }

    void construct_plastic_data()
    {
      resize_plastic_dof_numbers();

      construct_inv_fp_internal_data();
      construct_invBpks_internal_data();
      construct_fpcs_internal_data();
      construct_r_internal_data();
      construct_lambda_internal_data();

      // No longer really neccesary
      construct_h_internal_data();

      // We assign the equation numbers now because the user will likely want
      // all plastic data unpinned. If they pin any then they will need to call
      // assign_plastic_eqn_numbers again
      assign_plastic_eqn_numbers();
    }

    unsigned plastic_inv_fp_eqn_number(const unsigned& ipt,
                                       const unsigned& i,
                                       const unsigned& j) const
    {
      return Plastic_data_eqn_number[ipt][invFp_INDEX][i * DIM + j];
    }

    unsigned plastic_invBpks_eqn_number(const unsigned& ipt,
                                     const unsigned& i,
                                     const unsigned& j) const
    {
      return Plastic_data_eqn_number[ipt][invBpks_INDEX][i * DIM + j];
    }

    unsigned plastic_fpcs_eqn_number(const unsigned& ipt,
                                     const unsigned& i,
                                     const unsigned& j) const
    {
      return Plastic_data_eqn_number[ipt][Fpcs_INDEX][i * DIM + j];
    }

    unsigned plastic_h_eqn_number(const unsigned& ipt) const
    {
      return Plastic_data_eqn_number[ipt][H_INDEX][0];
    }

    unsigned plastic_r_eqn_number(const unsigned& ipt) const
    {
      return Plastic_data_eqn_number[ipt][R_INDEX][0];
    }

    unsigned plastic_lambda_eqn_number(const unsigned& ipt) const
    {
      return Plastic_data_eqn_number[ipt][Lambda_INDEX][0];
    }

    /*!
     * \brief Computes the deformation gradient tensor at a given integration
     * point
     * \param[in] t The time, for which F should be computed. Note, dotF is
     * computed only for t=0
     * \param[in] inpt The integration point
     * \param[out] F The resulting deformation gradient
     * \param[out] dotF The first time derivative of F
     */
    void compute_deformation_gradient_tensor(unsigned int t,
                                             unsigned int intpt,
                                             DenseMatrix<double>& F) const;

    void compute_deformation_gradient_tensor(unsigned int intpt,
                                             DenseMatrix<double>& F) const
    {
      compute_deformation_gradient_tensor(0, intpt, F);
    }

    void compute_mandel_stress_elastic(const DenseMatrix<double>& invFp,
                                       const DenseMatrix<double>& FtF,
                                       DenseMatrix<double>& Mbar,
                                       RankFourTensor<double>& dMbardinvFp,
                                       bool computeDerivative);

    void compute_mandel_stress_elastic(const DenseMatrix<double>& invFp,
                                       const DenseMatrix<double>& FtF,
                                       DenseMatrix<double>& Mbar)
    {
      return compute_mandel_stress_elastic(
        invFp,
        FtF,
        Mbar,
        PlasticEquationsBase<DIM>::Dummy_rankfourtensor,
        false);
    }

    void compute_mandellike_kinematic_hardening(
      const DenseMatrix<double>& invBpks,
      DenseMatrix<double>& bar_Mk,
      RankFourTensor<double>& dbar_MkdinvBpks,
      bool computeDerivative);

    void compute_mandellike_kinematic_hardening(const DenseMatrix<double>& invBpks,
                                                DenseMatrix<double>& bar_Mk)
    {
      return compute_mandellike_kinematic_hardening(
        invBpks, bar_Mk, PlasticEquationsBase<DIM>::Dummy_rankfourtensor, false);
    }

    void compute_mandellike_elastic_core(const DenseMatrix<double>& Fpcs,
                                         DenseMatrix<double>& bar_Mc,
                                         RankFourTensor<double>& dbar_McdFpcs,
                                         bool computeDerivative);

    void compute_mandellike_elastic_core(const DenseMatrix<double>& Fpcs,
                                         DenseMatrix<double>& bar_Mc)
    {
      return compute_mandellike_elastic_core(
        Fpcs, bar_Mc, PlasticEquationsBase<DIM>::Dummy_rankfourtensor, false);
    }

    void compute_mandel_stress_total(const DenseMatrix<double>& bar_M,
                                     const DenseMatrix<double>& bar_Mk,
                                     const DenseMatrix<double>& bar_Mc,
                                     const double& R,
                                     DenseMatrix<double>& barbar_M,
                                     double& dbarbar_M_dMk,
                                     double& dbarbar_M_dMc,
                                     DenseMatrix<double>& dbarbar_M_dR,
                                     bool computeJacobian);

    void compute_mandel_stress_total(const DenseMatrix<double>& bar_M,
                                     const DenseMatrix<double>& bar_Mk,
                                     const DenseMatrix<double>& bar_Mc,
                                     const double& R,
                                     DenseMatrix<double>& barbar_M)
    {
      double ddR = 0.0;
      compute_mandel_stress_total(bar_M,
                                  bar_Mk,
                                  bar_Mc,
                                  R,
                                  barbar_M,
                                  ddR,
                                  ddR,
                                  GeneralisedElement::Dummy_matrix,
                                  false);
    }

    void compute_barbar_N(const DenseMatrix<double>& barbar_M,
                          const double& f,
                          const DenseMatrix<double>& dfdM,
                          DenseMatrix<double>& barbar_N,
                          RankFourTensor<double>& dbarbar_N_dbarba_M,
                          bool computeDerivative);

    void compute_barbar_N(const DenseMatrix<double>& barbar_M,
                          const double& f,
                          const DenseMatrix<double>& dfdM,
                          DenseMatrix<double>& barbar_N)
    {
      compute_barbar_N(barbar_M,
                       f,
                       dfdM,
                       barbar_N,
                       GeneralisedElement::Dummy_rankfourtensor,
                       false);
    }

    /*!
     * \brief This function computes \bar{L}^\text{p}/\dot{\bar{\lambda}}
     */
    void compute_bar_Lp(const DenseMatrix<double>& bar_M,
                        const DenseMatrix<double>& barbar_N,
                        const RankFourTensor<double>& dbarbar_N_dbarbar_M,
                        DenseMatrix<double>& bar_Lp,
                        RankFourTensor<double>& dbar_Lp_dbarbar_M,
                        RankFourTensor<double>& dbar_Lp_dbar_M,
                        bool computeDerivative);

    void compute_bar_Lp(const DenseMatrix<double>& bar_M,
                        const DenseMatrix<double>& barbar_N,
                        DenseMatrix<double>& bar_Lp)
    {
      compute_bar_Lp(bar_M,
                     barbar_N,
                     GeneralisedElement::Dummy_rankfourtensor,
                     bar_Lp,
                     GeneralisedElement::Dummy_rankfourtensor,
                     GeneralisedElement::Dummy_rankfourtensor,
                     false);
    }

    /*!
     * \brief This function computes \bar{L}^\text{pkd}/\dot{\bar{\lambda}}
     */
    void compute_bar_Lpkd(const DenseMatrix<double>& bar_M,
                          const DenseMatrix<double>& bar_Mk,
                          DenseMatrix<double>& bar_Lpkd,
                          RankFourTensor<double>& dbar_Lpkd_dbar_M,
                          RankFourTensor<double>& dbar_Lpkd_dbar_Mk,
                          bool computeDerivative);

    void compute_bar_Lpkd(const DenseMatrix<double>& bar_M,
                          const DenseMatrix<double>& bar_Mk,
                          DenseMatrix<double>& bar_Lpkd)
    {
      compute_bar_Lpkd(bar_M,
                       bar_Mk,
                       bar_Lpkd,
                       GeneralisedElement::Dummy_rankfourtensor,
                       GeneralisedElement::Dummy_rankfourtensor,
                       false);
    }

    void compute_hat_bar_Nc(const double& f_Mc,
                            const DenseMatrix<double>& df_Mc_dMc,
                            DenseMatrix<double>& hat_bar_Nc,
                            RankFourTensor<double>& dhar_bar_Nc_dMc,
                            bool computederivative);

    void compute_hat_bar_Nc(const double& f_Mc,
                            const DenseMatrix<double>& df_Mc_dMc,
                            DenseMatrix<double>& hat_bar_Nc)
    {
      compute_hat_bar_Nc(f_Mc,
                         df_Mc_dMc,
                         hat_bar_Nc,
                         GeneralisedElement::Dummy_rankfourtensor,
                         false);
    }

    /*!
     * \brief This function computes \bar{L}^\text{pcd}/\dot{\bar{\lambda}}
     */
    void compute_bar_Lpcd(const DenseMatrix<double>& bar_M,
                          const DenseMatrix<double>& hat_bar_Nc,
                          const double& Rc,
                          const RankFourTensor<double>& dhat_bar_Nc_dhat_bar_Mc,
                          const DenseMatrix<double>& dRc_dMc,
                          const double& dRc_dH,
                          DenseMatrix<double>& bar_Lpcd,
                          RankFourTensor<double>& dbar_Lpcd_dbar_M,
                          RankFourTensor<double>& dbar_Lpcd_dhat_bar_Mc,
                          DenseMatrix<double>& dbar_Lpcd_dh,
                          bool computeDerivative);

    void compute_bar_Lpcd(const DenseMatrix<double>& bar_M,
                          const DenseMatrix<double>& hat_bar_Nc,
                          const double& Rc,
                          DenseMatrix<double>& bar_Lpcd)
    {
      double a;
      compute_bar_Lpcd(bar_M,
                       hat_bar_Nc,
                       Rc,
                       GeneralisedElement::Dummy_rankfourtensor,
                       GeneralisedElement::Dummy_matrix,
                       a,
                       bar_Lpcd,
                       GeneralisedElement::Dummy_rankfourtensor,
                       GeneralisedElement::Dummy_rankfourtensor,
                       GeneralisedElement::Dummy_matrix,
                       false);
    }

    double compute_c_sigma(const DenseMatrix<double>& bar_bar_N,
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
                           bool computeDerivative);

    double compute_u(const double& u_in,
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
                     bool computeDerivative);

    /*!
     * \brief Initialised the plastic parameters after the global deformation
     * has been updated
     * \details Assumes, all but the elastic parameters are equal to the ones
     * from the last time step. invFp is updated assuming all deformation was
     * purely elastic.
     */
    void initialise_solve(const unsigned int ipt);

    /*!
     * \brief Sets the values of the current plastic variables to their
     * respective initial values
     */
    void set_intial_condition(const unsigned int ipt);

    /*!
     * \brief Determins, if there is plastic deformation and the plastic solve
     * routine has to be called.
     * \details Checks if the elastic stress is contained by the yield surface
     * and if the plastic multiplier is positive.
     */
    bool is_there_plastic_deformation(const unsigned int ipt);

  public:
    Data* plastic_data_pt(const unsigned& ipt, const unsigned& data_type) const
    {
      return Plastic_data_pt[ipt][data_type];
    }

    double* plastic_dof_data_pt(const unsigned& ipt, const unsigned& ndof) const
    {
      return Plastic_dof_data_pt[ipt][ndof];
    }

    void enable_plastic_solve_by_fd()
    {
      Plastic_solve_use_fd = true;
    }
    void disable_plastic_solve_by_fd()
    {
      Plastic_solve_use_fd = false;
    }

    const bool get_plastic_solve_by_fd() const
    {
      return Plastic_solve_use_fd;
    }

    /*!
     * \brief returns the Cauchy stress for an integration point
     */
    void get_cauchy_stress(unsigned ipt, DenseMatrix<double>& sigma);

    double compute_r_plastic(const double& u,
                             const double& delta_lambda,
                             const double& R_prev)
    {
      double derivative;
      return compute_r_plastic(
        u, delta_lambda, R_prev, derivative, derivative, 0);
    }

    double compute_r_plastic(const double& u,
                             const double& delta_lambda,
                             const double& R_prev,
                             double& dRdLambda,
                             double& dRdu,
                             bool computeDerivative);

    std::string detail_plastic_dofs() const
    {
      std::stringstream str_str;
      for (unsigned ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        str_str << "Integral point " << ipt << ": ";
        for (unsigned data_type = 0;
             data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
             data_type++)
        {
          str_str << "[" << Plastic_data_names[data_type] << "]:";
          const unsigned n_value = Plastic_data_pt[ipt][data_type]->nvalue();
          if (data_type < NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE)
          {
            for (unsigned i = 0; i < n_value; i++)
            {
              str_str << " " << Plastic_data_eqn_number[ipt][data_type][i]
                      << " " << plastic_data_pt(ipt, data_type)->value(i)
                      << ", ";
            }
            str_str << ", ";
          }
          else
          {
            str_str << " not solved through newton solve";
          }
        }
        str_str << std::endl;
      }
      return str_str.str();
    }

    // Serialise all the plastic data from the integral points in order
    void serialise_all_plastic_data(Vector<double>& data,
                                    const unsigned& t = 0) const
    {
      for (unsigned ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        for (unsigned data_type = 0;
             data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
             data_type++)
        {
          const unsigned n_value = Plastic_data_pt[ipt][data_type]->nvalue();
          for (unsigned i = 0; i < n_value; i++)
          {
            data.push_back(plastic_data_pt(ipt, data_type)->value(t, i));
          }
        }
      }
    }

    // Serialise the plastic data from a specific integral point
    void serialise_plastic_data(Vector<double>& data,
                                const unsigned& ipt,
                                const unsigned& t = 0) const
    {
      for (unsigned data_type = 0; data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
           data_type++)
      {
        const unsigned n_value = Plastic_data_pt[ipt][data_type]->nvalue();
        for (unsigned i = 0; i < n_value; i++)
        {
          data.push_back(plastic_data_pt(ipt, data_type)->value(t, i));
        }
      }
    }

    // Serialise the plastic data interpolated at a given point in the element
    void interpolate_plastic_data_serialised(Vector<double>& data,
                                             Vector<double>& s,
                                             const unsigned& t = 0)
    {
      const unsigned n_node = this->nnode();
      Shape psi(n_node);
      shape(s, psi);

      const unsigned n_ipt = this->integral_pt()->nweight();
      Vector<Vector<double>> ipt_data(n_ipt);
      for (unsigned ipt = 0; ipt < n_ipt; ipt++)
      {
        serialise_plastic_data(ipt_data[ipt], ipt);
      }

      const unsigned n_data = ipt_data[0].size();
#ifdef PARANOID
      for (unsigned ipt = 1; ipt < n_ipt; ipt++)
      {
        if (ipt_data[ipt].size() != n_data)
        {
          throw OomphLibError("Size of data at ipts does not match",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION)
        }
      }
#endif

      data.resize(n_data, 0.0);
      for (unsigned ipt = 0; ipt < n_ipt; ipt++)
      {
        for (unsigned l = 0; l < n_node; l++)
        {
          const double interp_weight =
            this->integral_point_to_node_weight(ipt, l) * psi[l];
          for (unsigned i = 0; i < n_data; i++)
          {
            data[i] += interp_weight * ipt_data[ipt][i];
          }
        }
      }
    }

    // Assign the plastic data to a specific integral point
    void assign_plastic_data_serialised(const Vector<double>& data,
                                        const unsigned& ipt,
                                        const unsigned& t = 0)
    {
      unsigned data_count = 0;
      for (unsigned data_type = 0; data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
           data_type++)
      {
        if (data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES)
        {
          const unsigned n_value = Plastic_data_pt[ipt][data_type]->nvalue();
          for (unsigned i = 0; i < n_value; i++)
          {
            plastic_data_pt(ipt, data_type)
              ->set_value(t, i, data[data_count++]);
          }
        }
      }
    }

    /// Return the plastic constitutive law pointer
    PlasticConstitutiveLaw*& plastic_constitutive_law_pt()
    {
      return Plastic_consitutive_law_pt;
    }

    void assign_default_values_based_on_constitutive_law()
    {
      double Re = this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic;
      for (unsigned int ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        set_r(ipt, Re);
      }
    }

    /// Return the plastic newton solver tolerance
    double& plastic_newton_solver_tolerance()
    {
      return Plastic_Newton_Solver_Tolerance;
    }

    // Get the number of plastic data at the given ipt which is not pinned and
    // is to be solved for
    unsigned get_num_plastic_dofs(const unsigned& ipt)
    {
      return Plastic_dof_data_pt[ipt].size();
    }

    // Access the plasticity variable values
    double get_inv_fp(const unsigned& ipt,
                      const unsigned& i,
                      const unsigned& j) const
    {
      return get_inv_fp(0, ipt, i, j);
    }

    // Access the plasticity variable values
    double get_inv_fp(const unsigned& t,
                      const unsigned& ipt,
                      const unsigned& i,
                      const unsigned& j) const
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      return Plastic_data_pt[ipt][invFp_INDEX]->value(t, i * DIM + j);
    }

    /*!
     * \brief returns the elastic deformation gradient as a matrix
     */
    void get_inv_fp_matrix(const unsigned& ipt,
                           DenseMatrix<double>& invFp) const
    {
      get_inv_fp_matrix(0, ipt, invFp);
    }

    /*!
     * \brief returns the elastic deformation gradient as a matrix
     */
    void get_inv_fp_matrix(const unsigned& t,
                           const unsigned& ipt,
                           DenseMatrix<double>& invFp) const
    {
      invFp.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          invFp(i, j) =
            Plastic_data_pt[ipt][invFp_INDEX]->value(t, i * DIM + j);
        }
      }
    }

    // Access the plasticity variable values
    double get_dot_inv_fp(const unsigned& ipt,
                          const unsigned& i,
                          const unsigned& j) const
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      return this->dinternal_data_dt(Plastic_data_pt[ipt][invFp_INDEX],
                                     i * DIM + j);
    }

    /*!
     * \brief returns the time derivative of the elastic deformation gradient as
     * a matrix
     * \details Uses the data's time stepper to compute the derivative
     */
    void get_dot_inv_fp_matrix(const unsigned& ipt,
                               DenseMatrix<double>& dotinvFp)
    {
      dotinvFp.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          dotinvFp(i, j) = this->dinternal_data_dt(
            Plastic_data_pt[ipt][invFp_INDEX], i * DIM + j);
        }
      }
    }

    /*!
     * \details Returns the difference between the current and the previous
     * timestep. If the timestepper does not contain history, it will return
     * the difference between the current timestep and the initial value.
     */
    void get_delta_inv_fp_matrix(const unsigned& ipt,
                                 DenseMatrix<double>& deltainvFp)
    {
      const TimeStepper* invFp_time_stepper_pt =
        Plastic_data_pt[ipt][invFp_INDEX]->time_stepper_pt();

      deltainvFp.resize(DIM);

      if (invFp_time_stepper_pt->ntstorage() > 1)
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            deltainvFp(i, j) =
              get_inv_fp(0, ipt, i, j) - get_inv_fp(1, ipt, i, j);
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            deltainvFp(i, j) = get_inv_fp(0, ipt, i, j);

            if (i == j) deltainvFp(i, j) -= 1;
          }
        }
      }
    }

    void get_dot_or_delta_inv_fp_matrix(const unsigned& ipt,
                                        DenseMatrix<double>& dot_or_delta_invFp)
    {
      const TimeStepper* invFp_time_stepper_pt =
        Plastic_data_pt[ipt][invFp_INDEX]->time_stepper_pt();

      if (invFp_time_stepper_pt->is_steady())
      {
        return get_delta_inv_fp_matrix(ipt, dot_or_delta_invFp);
      }

      return get_dot_inv_fp_matrix(ipt, dot_or_delta_invFp);
    }

    double get_invBpks(const unsigned& ipt,
                    const unsigned& i,
                    const unsigned& j) const
    {
      return get_invBpks(0, ipt, i, j);
    }

    double get_invBpks(const unsigned int t,
                    const unsigned& ipt,
                    const unsigned& i,
                    const unsigned& j) const
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      return Plastic_data_pt[ipt][invBpks_INDEX]->value(t, i * DIM + j);
    }

    void get_invBpks_matrix(const unsigned& ipt, DenseMatrix<double>& invBpks)
    {
      return get_invBpks_matrix(0, ipt, invBpks);
    }

    void get_invBpks_matrix(const unsigned int t,
                         const unsigned int& ipt,
                         DenseMatrix<double>& invBpks)
    {
      invBpks.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          invBpks(i, j) = Plastic_data_pt[ipt][invBpks_INDEX]->value(t, i * DIM + j);
        }
      }
    }

    /*!
     * \brief returns the time derivative of the kinematic hardening deformation
     * gradient as a matrix
     * \details Uses the data's time stepper to compute the derivative
     */
    void get_dot_invBpks_matrix(const unsigned& ipt, DenseMatrix<double>& dotinvBpks)
    {
      dotinvBpks.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          dotinvBpks(i, j) = this->dinternal_data_dt(
            Plastic_data_pt[ipt][invBpks_INDEX], i * DIM + j);
        }
      }
    }

    /*!
     * \details Returns the difference between the current and the previous
     * timestep. If the timestepper does not contain history, it will return
     * the difference between the current timestep and the initial value.
     */
    void get_delta_invBpks_matrix(const unsigned& ipt,
                               DenseMatrix<double>& delta_invBpks)
    {
      const TimeStepper* invBpks_time_stepper_pt =
        Plastic_data_pt[ipt][invBpks_INDEX]->time_stepper_pt();

      delta_invBpks.resize(DIM);

      if (invBpks_time_stepper_pt->ntstorage() > 1)
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_invBpks(i, j) = get_invBpks(0, ipt, i, j) - get_invBpks(1, ipt, i, j);
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_invBpks(i, j) = get_invBpks(0, ipt, i, j);

            if (i == j) delta_invBpks(i, j) -= 1;
          }
        }
      }
    }

    void get_dot_or_delta_invBpks_matrix(const unsigned& ipt,
                                      DenseMatrix<double>& dot_or_delta_invBpks)
    {
      const TimeStepper* invBpks_time_stepper_pt =
        Plastic_data_pt[ipt][invBpks_INDEX]->time_stepper_pt();

      if (invBpks_time_stepper_pt->is_steady())
      {
        return get_delta_invBpks_matrix(ipt, dot_or_delta_invBpks);
      }

      return get_dot_invBpks_matrix(ipt, dot_or_delta_invBpks);
    }

    double get_fpcs(const unsigned& ipt,
                    const unsigned& i,
                    const unsigned& j) const
    {
      return get_fpcs(0, ipt, i, j);
    }

    double get_fpcs(const unsigned& t,
                    const unsigned& ipt,
                    const unsigned& i,
                    const unsigned& j) const
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      return Plastic_data_pt[ipt][Fpcs_INDEX]->value(t, i * DIM + j);
    }

    void get_fpcs_matrix(const unsigned int& ipt, DenseMatrix<double>& Fpcs)
    {
      get_fpcs_matrix(0, ipt, Fpcs);
    }

    void get_fpcs_matrix(const unsigned int t,
                         const unsigned int& ipt,
                         DenseMatrix<double>& Fpcs)
    {
      Fpcs.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          Fpcs(i, j) = Plastic_data_pt[ipt][Fpcs_INDEX]->value(t, i * DIM + j);
        }
      }
    }

    /*!
     * \brief returns the time derivative of the elastic core gradient as a
     * matrix
     * \details Uses the data's time stepper to compute the derivative
     */
    void get_dot_fpcs_matrix(const unsigned& ipt, DenseMatrix<double>& dotFpcs)
    {
      dotFpcs.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          dotFpcs(i, j) = this->dinternal_data_dt(
            Plastic_data_pt[ipt][Fpcs_INDEX], i * DIM + j);
        }
      }
    }

    /*!
     * \details Returns the difference between the current and the previous
     * timestep. If the timestepper does not contain history, it will return
     * the difference between the current timestep and the initial value.
     */
    void get_delta_fpcs_matrix(const unsigned& ipt,
                               DenseMatrix<double>& delta_Fpcs)
    {
      const TimeStepper* fpcs_time_stepper_pt =
        Plastic_data_pt[ipt][Fpcs_INDEX]->time_stepper_pt();

      delta_Fpcs.resize(DIM);

      if (fpcs_time_stepper_pt->ntstorage() > 1)
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_Fpcs(i, j) = get_fpcs(0, ipt, i, j) - get_fpcs(1, ipt, i, j);
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_Fpcs(i, j) = get_fpcs(0, ipt, i, j);

            if (i == j) delta_Fpcs(i, j) -= 1;
          }
        }
      }
    }

    void get_dot_or_delta_fpcs_matrix(const unsigned& ipt,
                                      DenseMatrix<double>& dot_or_delta_Fpcs)
    {
      const TimeStepper* fpcs_time_stepper_pt =
        Plastic_data_pt[ipt][Fpcs_INDEX]->time_stepper_pt();

      if (fpcs_time_stepper_pt->is_steady())
      {
        return get_delta_fpcs_matrix(ipt, dot_or_delta_Fpcs);
      }

      return get_dot_fpcs_matrix(ipt, dot_or_delta_Fpcs);
    }

    double get_h(const unsigned& ipt) const
    {
      return get_h(0, ipt);
    }

    double get_h(const unsigned t, const unsigned& ipt) const
    {
      return Plastic_data_pt[ipt][H_INDEX]->value(t, 0);
    }

    double get_dot_h(const unsigned& ipt) const
    {
      return this->dinternal_data_dt(Plastic_data_pt[ipt][H_INDEX], 0);
    }

    void set_h(const unsigned& ipt, const double& val) const
    {
      return Plastic_data_pt[ipt][H_INDEX]->set_value(0, val);
    }

    double get_r(const unsigned& ipt) const
    {
      return get_r(0, ipt);
    }

    double get_r(const unsigned t, const unsigned& ipt) const
    {
      if (Plastic_data_pt[ipt][R_INDEX]->time_stepper_pt()->ntstorage() < t + 1)
        return this->Plastic_consitutive_law_pt->normal_yield_ratio_elastic;
      return Plastic_data_pt[ipt][R_INDEX]->value(t, 0);
    }

    double get_dot_r(const unsigned& ipt) const
    {
      return this->dinternal_data_dt(Plastic_data_pt[ipt][R_INDEX], 0);
    }

    double get_lambda(const unsigned& ipt) const
    {
      return get_lambda(0, ipt);
    }

    double get_lambda(const unsigned t, const unsigned& ipt) const
    {
      return Plastic_data_pt[ipt][Lambda_INDEX]->value(t, 0);
    }

    double get_dot_lambda(const unsigned& ipt) const
    {
      return this->dinternal_data_dt(Plastic_data_pt[ipt][Lambda_INDEX], 0);
    }

    double get_dot_or_delta_lambda(const unsigned& ipt) const
    {
      const TimeStepper* lambda_time_stepper_pt =
        Plastic_data_pt[ipt][Lambda_INDEX]->time_stepper_pt();
      if (!lambda_time_stepper_pt->is_steady())
        return this->dinternal_data_dt(Plastic_data_pt[ipt][Lambda_INDEX], 0);

      return get_delta_lambda(ipt);
    }

    /*!
     * \details Returns the difference between the current and the previous
     * timestep. If the timestepper does not contain history, it will return
     * the difference between the current timestep and the initial value.
     */
    double get_delta_lambda(const unsigned& ipt) const
    {
      const TimeStepper* lambda_time_stepper_pt =
        Plastic_data_pt[ipt][Lambda_INDEX]->time_stepper_pt();
      if (lambda_time_stepper_pt->ntstorage() > 1)
        return get_lambda(0, ipt) - get_lambda(1, ipt);

      return get_lambda(ipt);
    }

    void set_lambda(const unsigned& ipt, const double& val) const
    {
      return Plastic_data_pt[ipt][Lambda_INDEX]->set_value(0, val);
    }

    // Set the plasticity variable values
    void set_inv_fp(const unsigned& ipt,
                    const unsigned& i,
                    const unsigned& j,
                    const double& val)
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      Plastic_data_pt[ipt][invFp_INDEX]->set_value(i * DIM + j, val);
    }

    /*!
     * \brief returns the elastic deformation gradient as a matrix
     */
    void set_inv_fp_matrix(const unsigned& ipt,
                           const DenseMatrix<double>& invFp)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          Plastic_data_pt[ipt][invFp_INDEX]->set_value(i * DIM + j,
                                                       invFp(i, j));
        }
      }
    }

    void set_invBpks(const unsigned& ipt,
                  const unsigned& i,
                  const unsigned& j,
                  const double& val)
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      Plastic_data_pt[ipt][invBpks_INDEX]->set_value(i * DIM + j, val);
    }
    void set_fpcs(const unsigned& ipt,
                  const unsigned& i,
                  const unsigned& j,
                  const double& val)
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      Plastic_data_pt[ipt][Fpcs_INDEX]->set_value(i * DIM + j, val);
    }
    void set_r(const unsigned& ipt, const double& val)
    {
      Plastic_data_pt[ipt][R_INDEX]->set_value(0, val);
    }

    virtual void fill_in_generic_residual_and_jacobian_plastic(
      DoubleVector& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& ipt,
      const DenseMatrix<double>& C,
      const unsigned& flag);

    void fill_in_residuals_plastic(DoubleVector& residuals,
                                   const unsigned& ipt,
                                   const DenseMatrix<double>& C)
    {
      fill_in_generic_residual_and_jacobian_plastic(
        residuals, GeneralisedElement::Dummy_matrix, ipt, C, 0);
    }

    // By default we finite difference for the jacobian
    virtual void fill_in_jacobian_plastic(DoubleVector& residuals,
                                          DenseMatrix<double>& jacobian,
                                          const unsigned& ipt,
                                          const DenseMatrix<double>& C)
    {
      if (Plastic_solve_use_fd)
      {
        fill_in_residuals_plastic(residuals, ipt, C);

        DoubleVector test_residuals(residuals);

        // Begin finite differencing wrt the plastic dofs
        for (unsigned local_unknown = 0;
             local_unknown < get_num_plastic_dofs(ipt);
             local_unknown++)
        {
          const double record_val = *Plastic_dof_data_pt[ipt][local_unknown];
          *Plastic_dof_data_pt[ipt][local_unknown] +=
            FiniteElement::Default_fd_jacobian_step;

          fill_in_residuals_plastic(test_residuals, ipt, C);
          for (unsigned local_eqn = 0; local_eqn < get_num_plastic_dofs(ipt);
               local_eqn++)
          {
            jacobian(local_eqn, local_unknown) =
              (test_residuals[local_eqn] - residuals[local_eqn]) /
              FiniteElement::Default_fd_jacobian_step;
          }
          *Plastic_dof_data_pt[ipt][local_unknown] = record_val;
        }
      }
      else
      {
        fill_in_generic_residual_and_jacobian_plastic(
          residuals, jacobian, ipt, C, 1);
      }
    }

    // Change this to solve for a specific integral point - this will make
    // finite differencing wrt C much easier
    void plastic_newton_solve(const unsigned& ipt, const DenseMatrix<double>& C)
    {
      // // std::cout << detail_plastic_dofs() << std::endl;
      // // std::cout << "Entering plastic Newton solve" << std::endl;
      // for (unsigned ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      // {
      // Plastic deformation depends on the previous solution. So if the
      // intermediate steps in newton solve were outide, it should be
      // reconsidered.
      initialise_solve(ipt);

      if (!is_there_plastic_deformation(ipt))
      {
        // oomph_info << "The deformation is purely elastic, continue"
        //            << std::endl;
        return;
      }

      // std::cout << "\tIntegral point: " << ipt << std::endl;
      // For now we just build the linear algebra distribution and the solver
      // each time we solve:

      OomphCommunicator* communicator_pt = new OomphCommunicator();

      // Create the linear algebra distribution - do not distribute it
      LinearAlgebraDistribution* lin_alg_dist_pt =
        new LinearAlgebraDistribution(
          communicator_pt, this->get_num_plastic_dofs(ipt), false);

      // Create the linear solver and pass to the matrix
      DenseLU* solver_pt = new DenseLU;

      DoubleVector residuals(lin_alg_dist_pt);
      DenseDoubleMatrix jacobian(this->get_num_plastic_dofs(ipt));

      jacobian.linear_solver_pt() = solver_pt;

      double maxres;
      unsigned int nIter = 0;
      do
      {
        // Get the residuals only
        this->fill_in_residuals_plastic(residuals, ipt, C);

        maxres = residuals.max();

        // oomph_info << "Maximum residual: " << maxres << std::endl;
        // if (Plastic_data_has_been_built[R_INDEX])
        // {
        //   oomph_info << "R = " << get_r(ipt) << std::endl;
        // }

        // need to initialise the jacobian to 0
        jacobian.initialise(0.0);
        // get jacobian
        fill_in_jacobian_plastic(residuals, jacobian, ipt, C);

        // std::cout << "Residuals are " << std::endl
        //           << "   invFp: "
        //           << this->plastic_inv_fp_eqn_number(ipt, 0, 0) << " to "
        //           << this->plastic_inv_fp_eqn_number(ipt, DIM - 1, DIM - 1)
        //           << std::endl;
        // if (Plastic_data_has_been_built[invBpks_INDEX])
        // {
        //   std::cout << "   invBpks: " << this->plastic_invBpks_eqn_number(ipt, 0,
        //   0)
        //             << " to "
        //             << this->plastic_invBpks_eqn_number(ipt, DIM - 1, DIM - 1)
        //             << std::endl;
        // }
        // if (Plastic_data_has_been_built[R_INDEX])
        // {
        //   std::cout << "   R: " << this->plastic_r_eqn_number(ipt)
        //             << std::endl;
        // }
        // if (Plastic_data_has_been_built[H_INDEX])
        // {
        //   std::cout << "   H: " << this->plastic_h_eqn_number(ipt)
        //             << std::endl;
        // }
        // if (Plastic_data_has_been_built[Lambda_INDEX])
        // {
        //   std::cout << "   YieldSurface: " <<
        //   this->plastic_lambda_eqn_number(ipt)
        //             << std::endl;
        // }

        // std::cout << "\t\tResiduals:" << std::endl;
        // for (unsigned i = 0; i < this->get_num_plastic_dofs(ipt); i++)
        // {
        //   std::cout << "\t\t\t" << residuals[i] << std::endl;
        // }
        // std::cout << "\t\tJacobian:" << std::endl;
        // for (unsigned i = 0; i < this->get_num_plastic_dofs(ipt); i++)
        // {
        //   std::cout << "\t\t\t";
        //   for (unsigned j = 0; j < this->get_num_plastic_dofs(ipt); j++)
        //   {
        //     std::cout << jacobian(i, j) << " ";
        //   }
        //   std::cout << std::endl;
        // }

        // compute delta for plastic degrees of freedom
        DoubleVector resid(residuals);
        jacobian.solve(resid, residuals);

        // update values
        double* dx_pt = residuals.values_pt();
        for (unsigned i = 0; i < this->get_num_plastic_dofs(ipt); i++)
        {
          *Plastic_dof_data_pt[ipt][i] -= dx_pt[i];
        }

        // Because h is not solved for, it will be updated here
        set_h(ipt,
              get_lambda(ipt) *
                this->Plastic_consitutive_law_pt->isotropic_hardening_factor);

        nIter++;
      } while (maxres > Plastic_Newton_Solver_Tolerance);

      delete lin_alg_dist_pt;
      delete solver_pt;
      delete communicator_pt;
      // }
    }

    void plastic_newton_solve()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
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
          // Loop over positional dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            // Loop over displacement components (deformed position)
            for (unsigned i = 0; i < DIM; i++)
            {
              // Loop over derivative directions
              for (unsigned j = 0; j < DIM; j++)
              {
                interpolated_G(j, i) +=
                  this->nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
              }
            }
          }
        }
        // Declare and calculate the deformed metric tensor
        DenseMatrix<double> G(DIM);

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

        // THEN WE SOLVE THE PLASTIC EQUATIONS
        plastic_newton_solve(ipt, G);
      }
    }

    PlasticEquationsBase()
      : InterpolateFromIntegralPointsBase(), PVDEquationsBase<DIM>()
    {
      this->unity.resize(DIM, DIM, 0.0);
      for (unsigned int i = 0; i < DIM; i++) this->unity(i, i) = 1;

      // Construct all plastic data
      construct_plastic_data();
    }

    ~PlasticEquationsBase() {}

    virtual void set_internal_data_time_stepper(
      const unsigned& i,
      TimeStepper* const& time_stepper_pt,
      const bool& preserve_existing_data) override
    {
      PVDEquationsBase<DIM>::set_internal_data_time_stepper(
        i, time_stepper_pt, preserve_existing_data);

      // We need to reassign plastic eqn numbers, after the data storage has
      // changed to update possibly old pointers
      assign_plastic_eqn_numbers();
    }

    // Assign the time-stepper for the plastic data only.
    // It should be preffered to call the function
    // set_internal_data_time_stepper
    void assign_plastic_timestepper(TimeStepper* time_stepper_pt,
                                    const bool& preserve_existing_data)
    {
      for (unsigned ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        for (unsigned data_type = 0;
             data_type < NUMBER_OF_PLASTIC_VARIABLE_TYPES;
             data_type++)
        {
          // Could swap the loops, slightly more efficient but we only do this
          // once per element
          Data* data_pt = Plastic_data_pt[ipt][data_type];
          // We want to preserve the existing data hence true
          data_pt->set_time_stepper(time_stepper_pt, preserve_existing_data);
        }
      }

      // Need to reassign eqn numbers because pointers may have changed.
      assign_plastic_eqn_numbers();
    }

    void check_initial_condition()
    {
      // Now check the initial values
      unsigned ipt = 0;
      DenseMatrix<double> invFp;
      get_inv_fp_matrix(0, ipt, invFp);
      oomph_info << "invFp(0): " << std::endl
                 << MatrixHelpers::format(invFp) << std::endl;

      get_inv_fp_matrix(1, ipt, invFp);
      oomph_info << "invFp(1): " << std::endl
                 << MatrixHelpers::format(invFp) << std::endl;

      get_inv_fp_matrix(2, ipt, invFp);
      oomph_info << "invFp(2): " << std::endl
                 << MatrixHelpers::format(invFp) << std::endl;

      DenseMatrix<double> invBpks;
      get_inv_fp_matrix(0, ipt, invBpks);
      oomph_info << "invBpks(0): " << std::endl
                 << MatrixHelpers::format(invBpks) << std::endl;

      get_inv_fp_matrix(1, ipt, invBpks);
      oomph_info << "invBpks(1): " << std::endl
                 << MatrixHelpers::format(invBpks) << std::endl;

      get_inv_fp_matrix(2, ipt, invBpks);
      oomph_info << "invBpks(2): " << std::endl
                 << MatrixHelpers::format(invBpks) << std::endl;

      oomph_info << "H(0) = " << get_h(0, ipt) << " H(1) = " << get_h(1, ipt)
                 << " H(2) = " << get_h(2, ipt) << std::endl;

      oomph_info << "lambda(0) = " << get_lambda(0, ipt)
                 << " lambda(1) = " << get_lambda(1, ipt)
                 << " lambda(2) = " << get_lambda(2, ipt) << std::endl;

      oomph_info << "R(0) = " << get_r(0, ipt) << " R(1) = " << get_r(1, ipt)
                 << " R(2) = " << get_r(2, ipt) << std::endl;
    }

    void calculate_g(const unsigned& ipt,
                     const double& diag_entry,
                     const DenseMatrix<double>& G,
                     DenseMatrix<double>& g) const override
    {
      // Compute the undeformed coordinates from the deformed ones
      // Solve the plastic equations
      // Could we just assume that the plastic data is consistent with G?
      // This will be the case after the combined (plastic and elastic) Newton
      // solve has converged so shouldn't affect outputting
      // this->plastic_newton_solve(ipt, G);

      DenseMatrix<double> invFp(DIM, DIM, 0.0);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          invFp(i, j) = get_inv_fp(ipt, i, j);
        }
      }

      // Calculate Fp
      DenseMatrix<double> Fp(DIM, DIM, 0.0);
      MatrixHelpers::invert_matrix<DIM>(invFp, Fp);

      g.resize(DIM, DIM, 0.0);
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = 0; j < DIM; j++)
        {
          for (unsigned k = 0; k < DIM; k++)
          {
            g(i, j) += Fp(k, i) * Fp(k, j);
          }
          g(i, j) *= diag_entry;
        }
      }
    }

    /*!
     * \brief compute the time derivative of an internal data item
     */
    double dinternal_data_dt(Data* data_pt, const unsigned value_idx) const
    {
      // Number of timsteps (past & present)
      const TimeStepper* data_time_stepper_pt = data_pt->time_stepper_pt();
      const unsigned n_time = data_time_stepper_pt->ntstorage();

      double dxdt = 0.0;

      // If the timestepper is not steady
      if (!data_time_stepper_pt->is_steady())
      {
        // Loop over the additional time storage and add the appropriate
        // contributions
        for (unsigned t = 0; t < n_time; t++)
        {
          dxdt +=
            data_time_stepper_pt->weight(1, t) * data_pt->value(t, value_idx);
        }
      }

      return dxdt;
    }
  };


  template<unsigned DIM>
  class PlasticEquations : public virtual PlasticEquationsBase<DIM>,
                           public virtual PVDEquations<DIM>
  {
  public:
    PlasticEquations() : PlasticEquationsBase<DIM>(), PVDEquations<DIM>() {}
    /// Return the derivatives of the 2nd Piola Kirchhoff stress tensor,
    /// as calculated from the constitutive law: Pass the interpolation point,
    /// the diagonal value of g, the metric tensors in the stress free and
    /// current configurations and the current value of the the stress tensor.
    inline virtual void get_d_stress_dG_upper(
      const unsigned& ipt,
      const double& diag_entry,
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      RankFourTensor<double>& d_sigma_dG) override
    {
      PVDEquations<DIM>::get_d_stress_dG_upper(g, G, sigma, d_sigma_dG);

      // Now compute the rest
      RankFourTensor<double> d_sigma_dg;
      this->Constitutive_law_pt->calculate_d_second_piola_kirchhoff_stress_dg(
        g, G, sigma, d_sigma_dg, false);

      //////////////////////////777
      // d_g_dG
      RankFourTensor<double> d_g_dG(DIM);

      // Record the plastic data values
      const unsigned num_plastic_dof = this->get_num_plastic_dofs(ipt);
      Vector<double> plastic_values_prior_to_fd(num_plastic_dof, 0.0);
      for (unsigned i = 0; i < num_plastic_dof; i++)
      {
        plastic_values_prior_to_fd[i] = *this->plastic_dof_data_pt(ipt, i);
      }

      // Perform finite differencing g wrt F
      DenseMatrix<double> g_new(DIM, DIM, 0.0);
      DenseMatrix<double> G_test(DIM, DIM);
      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned j = i; j < DIM; j++)
        {
          G_test(i, j) = G(i, j);
          G_test(j, i) = G(j, i);
        }
      }

      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i; j < DIM; j++)
        {
          const double saved_value = G(i, j);
          G_test(i, j) += FiniteElement::Default_fd_jacobian_step;
          if (i != j)
          {
            G_test(j, i) = G_test(i, j);
          }

          g_new.initialise(0.0);

          // Need to call plastic solve to update the internal variables to the
          // new G
          this->plastic_newton_solve(ipt, G_test);
          this->calculate_g(ipt, diag_entry, G_test, g_new);

          // We can reduce this to only the upper (or lower?) triangular
          // elements
          for (unsigned n = 0; n < DIM; n++)
          {
            for (unsigned m = 0; m < DIM; m++)
            {
              d_g_dG(n, m, i, j) = (g_new(n, m) - g(n, m)) /
                                   FiniteElement::Default_fd_jacobian_step;
            }
          }
          G_test(i, j) = saved_value;
          if (i != j)
          {
            G_test(j, i) = saved_value;
          }
        }
      }

      // Restore the values of the plastic variables
      for (unsigned i = 0; i < num_plastic_dof; i++)
      {
        *this->plastic_dof_data_pt(ipt, i) = plastic_values_prior_to_fd[i];
      }

      ////////////////
      // Put it all together
      // Loop Output Stress (i, j): Upper Triangle Only
      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i; j < DIM; j++)
        {
          // Loop Input Metric G (n, m): Upper Triangle Only
          for (unsigned n = 0; n < DIM; n++)
          {
            for (unsigned m = n; m < DIM; m++)
            {
              double sum = 0.0;

              // Summation Loop (k, l): Upper Triangle Only
              for (unsigned k = 0; k < DIM; k++)
              {
                for (unsigned l = k; l < DIM; l++)
                {
                  sum += d_sigma_dg(i, j, k, l) * d_g_dG(k, l, n, m);
                }
              }

              // Assign to the target tensor
              d_sigma_dG(i, j, n, m) += sum;
            }
          }
        }
      }
    }
  };

  template<unsigned DIM>
  class PlasticEquationswithPressure
    : public virtual PlasticEquationsBase<DIM>,
      public virtual PVDEquationsWithPressure<DIM>
  {
  public:
    PlasticEquationswithPressure()
      : PlasticEquationsBase<DIM>(), PVDEquationsWithPressure<DIM>()
    {
    }

    inline virtual void get_d_stress_dG_upper(
      const unsigned& ipt,
      const double& diag_entry,
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      const double& gen_dil,
      const double& inv_kappa,
      const double& interpolated_solid_p,
      RankFourTensor<double>& d_sigma_dG,
      DenseMatrix<double>& d_gen_dil_dG) override
    {
    }

    inline virtual void get_d_stress_dG_upper(
      const unsigned& ipt,
      const double& diag_entry,
      const DenseMatrix<double>& g,
      const DenseMatrix<double>& G,
      const DenseMatrix<double>& sigma,
      const double& detG,
      const double& interpolated_solid_p,
      RankFourTensor<double>& d_sigma_dG,
      DenseMatrix<double>& d_detG_dG) override
    {
    }
  };


  template<unsigned DIM, unsigned NNODE>
  class QPlasticPVDElement : public virtual SolidQElement<DIM, NNODE>,
                             public virtual PlasticEquations<DIM>
  {
  public:
    QPlasticPVDElement() : SolidQElement<DIM, NNODE>(), PlasticEquations<DIM>()
    {
      this->compute_ipt_to_node_mapping();
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // This is called at least twice per Newton solve but we only want one
      // So we only all when we compute residuals, NOT when we compute jacobian
      // since jacobian computation always follows a residual computation
      PlasticEquationsBase<DIM>::plastic_newton_solve();

      PVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Fill in contribution to Jacobian (either by FD or analytically,
    /// control this via evaluate_jacobian_by_fd()
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      PVDEquations<DIM>::fill_in_generic_contribution_to_residuals_pvd(
        residuals, jacobian, 1);
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      PVDEquations<DIM>::output(outfile);
    }

    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      PVDEquations<DIM>::output(outfile, n_plot);
    }

    /// C-style output function
    void output(FILE* file_pt)
    {
      PVDEquations<DIM>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PVDEquations<DIM>::output(file_pt, n_plot);
    }
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<QPlasticPVDElement<2, NNODE_1D>>
    : public virtual SolidQElement<1, NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<QPlasticPVDElement<3, NNODE_1D>>
    : public virtual SolidQElement<2, NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<2, NNODE_1D>() {}
  };


  template<unsigned DIM>
  class QPlasticPVDElementWithPressure
    : public virtual SolidQElement<DIM, 3>,
      public virtual PlasticEquationswithPressure<DIM>
  {
  public:
    QPlasticPVDElementWithPressure()
      : SolidQElement<DIM, 3>(), PlasticEquationswithPressure<DIM>()
    {
      this->compute_ipt_to_node_mapping();
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // This is called at least twice per Newton solve but we only want one
      // So we only all when we compute residuals, NOT when we compute jacobian
      // since jacobian computation always follows a residual computation
      PlasticEquationsBase<DIM>::plastic_newton_solve();

      PVDEquationsWithPressure<DIM>::
        fill_in_generic_contribution_to_residuals_pvd(
          residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Fill in contribution to Jacobian (either by FD or analytically,
    /// control this via evaluate_jacobian_by_fd()
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      PVDEquationsWithPressure<
        DIM>::fill_in_generic_contribution_to_residuals_pvd(residuals,
                                                            jacobian,
                                                            1);
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      PVDEquationsWithPressure<DIM>::output(outfile);
    }

    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      PVDEquationsWithPressure<DIM>::output(outfile, n_plot);
    }

    /// C-style output function
    void output(FILE* file_pt)
    {
      PVDEquationsWithPressure<DIM>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      PVDEquationsWithPressure<DIM>::output(file_pt, n_plot);
    }
  };


  template<>
  class FaceGeometry<QPlasticPVDElementWithPressure<2>>
    : public virtual SolidQElement<1, 3>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<1, 3>() {}
  };

  template<>
  class FaceGeometry<QPlasticPVDElementWithPressure<3>>
    : public virtual SolidQElement<2, 3>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<2, 3>() {}
  };


  // We define in here the additional functions required to build child elements
  // We need to pass any flags to the child as well as build the data at the
  // integral points, to do this we interpolate from the integral points to the
  // nodes first then interpolate to the child integral points
  template<unsigned DIM>
  class RefineablePlasticEquations : public virtual PlasticEquations<DIM>,
                                     public virtual RefineableSolidElement
  {
  public:
    RefineablePlasticEquations()
      : PlasticEquations<DIM>(), RefineableSolidElement()
    {
    }

    // We need to get the plastic data from the children at the integral points
    void rebuild_from_sons(Mesh*& mesh_pt) override
    {
      const unsigned n_ipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < n_ipt; ipt++)
      {
        // get the local coordinate of the integral point
        Vector<double> s(DIM);
        for (unsigned i = 0; i < DIM; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        // Determine what child element the integral point is within
        // and get the local cooordinate within that element
        unsigned child_num = 0;
        Vector<double> s_child(DIM, 0.0);

        if (DIM == 2)
        {
          using namespace QuadTreeNames;

          child_num += (s[0] >= 0.0);
          child_num += (s[1] >= 0.0) * 2;

          s_child[0] = (s[0] >= 0.0) ? (s[0] - 0.5) * 2.0 : (s[0] + 0.5) * 2.0;
          s_child[1] = (s[1] >= 0.0) ? (s[1] - 0.5) * 2.0 : (s[1] + 0.5) * 2.0;
        }
        else if (DIM == 3)
        {
          child_num += (s[0] >= 0.0);
          child_num += (s[1] >= 0.0) * 2;
          child_num += (s[2] >= 0.0) * 4;

          s_child[0] = (s[0] >= 0.0) ? (s[0] - 0.5) * 2.0 : (s[0] + 0.5) * 2.0;
          s_child[1] = (s[1] >= 0.0) ? (s[1] - 0.5) * 2.0 : (s[1] + 0.5) * 2.0;
          s_child[2] = (s[2] >= 0.0) ? (s[2] - 0.5) * 2.0 : (s[2] + 0.5) * 2.0;
        }
        else
        {
          // Throw an error
        }

        RefineablePlasticEquations<DIM>* child_pt =
          dynamic_cast<RefineablePlasticEquations<DIM>*>(
            this->tree_pt()->son_pt(child_num)->object_pt());

        // For the more general case where we may have a different time-stepper
        // for each of the plastic dofs we need a serialise function which
        // serialises previous time-step data too
        const unsigned ntstorage = this->plastic_data_pt(0, 0)->ntstorage();
        for (unsigned t = 0; t < ntstorage; t++)
        {
          Vector<double> child_data;
          child_pt->interpolate_plastic_data_serialised(child_data, s_child, t);
          this->assign_plastic_data_serialised(child_data, ipt, t);
        }
      }
    }

    void further_build() override
    {
      PlasticEquations<DIM>* cast_father_element_pt =
        dynamic_cast<PlasticEquations<DIM>*>(this->father_element_pt());

      this->Plastic_solve_use_fd =
        cast_father_element_pt->get_plastic_solve_by_fd();
      this->Plastic_consitutive_law_pt =
        cast_father_element_pt->plastic_constitutive_law_pt();
    }
  };

  template<unsigned DIM, unsigned NNODE>
  class RefineableQPlasticPVDElement
    : public virtual RefineableQPVDElement<DIM, NNODE>,
      public virtual RefineablePlasticEquations<DIM>
  {
  public:
    RefineableQPlasticPVDElement()
      : QPVDElement<DIM, NNODE>(),
        RefineableElement(),
        RefineableSolidElement(),
        RefineablePVDEquations<DIM>(),
        RefineableSolidQElement<DIM>(),
        RefineableQPVDElement<DIM, NNODE>(),
        PlasticEquations<DIM>()
    {
      this->compute_ipt_to_node_mapping();
    }

    void rebuild_from_sons(Mesh*& mesh_pt) override
    {
      RefineablePlasticEquations<DIM>::rebuild_from_sons(mesh_pt);
    }

    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // This is called at least twice per Newton solve but we only want one
      // So we only all when we compute residuals, NOT when we compute jacobian
      // since jacobian computation always follows a residual computation
      PlasticEquationsBase<DIM>::plastic_newton_solve();

      RefineablePVDEquations<DIM>::
        fill_in_generic_contribution_to_residuals_pvd(
          residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Fill in contribution to Jacobian (either by FD or analytically,
    /// control this via evaluate_jacobian_by_fd()
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      RefineablePVDEquations<
        DIM>::fill_in_generic_contribution_to_residuals_pvd(residuals,
                                                            jacobian,
                                                            1);
    }

    void further_build() override
    {
      RefineableQPVDElement<DIM, NNODE>::further_build();
      RefineablePlasticEquations<DIM>::further_build();

      RefineableQPlasticPVDElement<DIM, NNODE>* cast_father_element_pt =
        dynamic_cast<RefineableQPlasticPVDElement<DIM, NNODE>*>(
          this->father_element_pt());

      const unsigned n_ipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < n_ipt; ipt++)
      {
        // get the local coordinate of the integral point
        Vector<double> s(DIM);
        for (unsigned i = 0; i < DIM; i++)
        {
          s[i] = this->integral_pt()->knot(ipt, i);
        }

        Vector<double> s_father(DIM, 0.0);
        this->get_father_s(s, s_father);

        // For the more general case where we may have a different time-stepper
        // for each
        // of the plastic dofs we need a serialise function which serialises
        // previous time-step data too
        const unsigned ntstorage = this->plastic_data_pt(0, 0)->ntstorage();
        for (unsigned t = 0; t < ntstorage; t++)
        {
          Vector<double> father_data;
          cast_father_element_pt->interpolate_plastic_data_serialised(
            father_data, s_father, t);
          this->assign_plastic_data_serialised(father_data, ipt, t);
        }
      }
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      RefineableQPVDElement<DIM, NNODE>::output(outfile);
    }

    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      RefineableQPVDElement<DIM, NNODE>::output(outfile, n_plot);
    }

    /// C-style output function
    void output(FILE* file_pt)
    {
      RefineableQPVDElement<DIM, NNODE>::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      RefineableQPVDElement<DIM, NNODE>::output(file_pt, n_plot);
    }
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQPlasticPVDElement<2, NNODE_1D>>
    : public virtual SolidQElement<1, NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<1, NNODE_1D>() {}
  };

  template<unsigned NNODE_1D>
  class FaceGeometry<RefineableQPlasticPVDElement<3, NNODE_1D>>
    : public virtual SolidQElement<2, NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying solid element
    FaceGeometry() : SolidQElement<2, NNODE_1D>() {}
  };


} // namespace oomph

#endif