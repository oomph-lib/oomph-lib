#ifndef OOMPH_PLASTICITY_SOLID_ELEMENTS_HEADER
#define OOMPH_PLASTICITY_SOLID_ELEMENTS_HEADER

#include "solid_elements.h"
#include "constitutive/plastic_constitutive_laws.h"
#include "refineable_solid_elements.h"
#include "generic/interpolate_from_integral_points.h"

namespace oomph
{
  template<unsigned DIM>

  // ===========================================================================
  /// \short This class implements the extended subloading surface model as
  /// described in
  /// Hashiguchi, K. Multiplicative Hyperelastic-Based Plasticity for Finite
  /// Elastoplastic Deformation/Sliding: A Comprehensive Review. Arch Computat
  /// Methods Eng 26, 597–637 (2019). https://doi.org/10.1007/s11831-018-9256-5.
  // ===========================================================================
  class PlasticEquationsBase : public virtual InterpolateFromIntegralPointsBase,
                               public virtual PVDEquationsBase<DIM>
  {
  public:
    // =========================================================================
    /// \short Constructor: this initialises the unity matrix and the plastic
    /// data.
    // =========================================================================
    PlasticEquationsBase()
      : InterpolateFromIntegralPointsBase(), PVDEquationsBase<DIM>()
    {
      this->unity.resize(DIM, DIM, 0.0);
      for (unsigned int i = 0; i < DIM; i++) this->unity(i, i) = 1;

      // Construct all plastic data
      construct_plastic_data();
    }

    // =========================================================================
    /// \short Virtual default destructor
    // =========================================================================
    virtual ~PlasticEquationsBase() = default;

    // =========================================================================
    /// \short returns the Cauchy stress for an integration point
    /// \param[in] ipt The integration point.
    /// \param[out] sigma The computed Cauchy stress.
    // =========================================================================
    void get_cauchy_stress(const unsigned& ipt,
                           DenseMatrix<double>& sigma) const;

    // =========================================================================
    /// \short Prepares a string describing the plastic dofs.
    /// \returns The description.
    // =========================================================================
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

    // =========================================================================
    /// \short A function used for debugging which writes the initial condition
    /// to oomph_info.
    // =========================================================================
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

      oomph_info << "lambda(0) = " << get_lambda(0, ipt)
                 << " lambda(1) = " << get_lambda(1, ipt)
                 << " lambda(2) = " << get_lambda(2, ipt) << std::endl;

      oomph_info << "R(0) = " << get_r(0, ipt) << " R(1) = " << get_r(1, ipt)
                 << " R(2) = " << get_r(2, ipt) << std::endl;
    }

    // =========================================================================
    /// \short Enables the plastic solve through finite difference.
    /// Note: It is reccomended to rely on the analytical jacobian, whenever
    /// possible.
    // =========================================================================
    void enable_plastic_solve_by_fd()
    {
      Plastic_solve_use_fd = true;
    }

    // =========================================================================
    /// \short Disables the plastic solve through finite difference.
    // =========================================================================
    void disable_plastic_solve_by_fd()
    {
      Plastic_solve_use_fd = false;
    }

    // =========================================================================
    /// \returns if plastic solve is done by finite differencing.
    // =========================================================================
    const bool get_plastic_solve_by_fd() const
    {
      return Plastic_solve_use_fd;
    }

    // =========================================================================
    /// \short Access function to \ref Plastic_fd_jacobian_step_pt
    // =========================================================================
    double*& plastic_fd_jacobian_step_pt()
    {
      return Plastic_fd_jacobian_step_pt;
    }

    // =========================================================================
    /// \short Access function to the \ref Plastic_consitutive_law_pt
    // =========================================================================
    PlasticConstitutiveLaw*& plastic_constitutive_law_pt()
    {
      return Plastic_consitutive_law_pt;
    }

    // =========================================================================
    /// \short Assigns the fefault values of the plastic data based on the
    /// constitutive law. It is reccomended to call this function on every
    /// element, after the consitutive law has been set.
    ///
    /// This sets r to re.
    // =========================================================================
    void assign_default_values_based_on_constitutive_law()
    {
      double Re =
        this->Plastic_consitutive_law_pt->normal_yield_ratio_law_pt->get_re();
      for (unsigned int ipt = 0; ipt < this->integral_pt()->nweight(); ipt++)
      {
        set_r(ipt, Re);
      }
    }

    // =========================================================================
    /// \returns the plastic newton solver tolerance
    // =========================================================================
    double& plastic_newton_solver_tolerance()
    {
      return Plastic_Newton_Solver_Tolerance;
    }

    // =========================================================================
    /// \short Assigns the time-stepper for the plastic data only. It should be
    /// preffered to call the function set_internal_data_time_stepper when
    /// possible.
    ///
    /// \param[in] time_stepper_pt The time stepper to set.
    /// \param[in] preserve_existing_data Wether the internal data should be
    /// copied over to the new structure.
    // =========================================================================
    void assign_plastic_timestepper_pt(TimeStepper* time_stepper_pt,
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

  protected:
    ////////////////////////////////////////////////////////////////////////////
    // overriden functions
    ////////////////////////////////////////////////////////////////////////////

    // =========================================================================
    /// \short An Overridden version of the parent function \ref
    /// PVDEquationsBase<DIM>::set_internal_data_time_stepper. In addition to
    /// the parent function, this calls \ref assign_plastic_eqn_numbers.
    ///
    /// For a description of the arguments \see
    /// PVDEquationsBase<DIM>::set_internal_data_time_stepper
    // =========================================================================
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

    // =========================================================================
    /// \short An Overridden version of the parent function \ref
    /// PVDEquationsBase<DIM>::calculate_g. This function calculates g while
    /// taking the current plastic deformation into account.
    ///
    /// For a description of the arguments \see
    /// PVDEquationsBase<DIM>::calculate_g
    // =========================================================================
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

      // Calculate and pull back Fp^TFp to get g
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

    ////////////////////////////////////////////////////////////////////////////
    // Functions directly involved in the plastic solve
    ////////////////////////////////////////////////////////////////////////////

    // =========================================================================
    /// \short This function solves for the plastic variables of all integration
    /// points. For that, it first computes the global deformation and then
    /// calls the plastic_newton_solve on every data point.
    // =========================================================================
    void plastic_newton_solve()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      Vector<double> s(DIM);
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
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

        // THEN WE SOLVE THE PLASTIC EQUATIONS
        plastic_newton_solve(ipt, G);
      }
    }

    // =========================================================================
    /// \short This function solves for the plastic variables of one integration
    /// point. For that, it first checks, whether there is plastic deformation.
    /// If there is, it performs an internal newton solve. If there is not, only
    /// r is updated.
    ///
    /// \param[in] ipt The integration point
    /// \param[out] C The global deformation gradient in the undeformed frame of
    /// the plasticity model.
    // =========================================================================
    void plastic_newton_solve(const unsigned& ipt, const DenseMatrix<double>& C)
    {
      // Initialize all plastic variable such that their derivative is 0.
      // In any case, r is clamped to its validity range [Re, 1].
      initialise_solve(ipt);

      // Check if there is any plastic deformation. If not, r is automatically
      // updated to fullfill the yield surface equation.
      if (!is_there_plastic_deformation(ipt)) return;


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

      unsigned int nIter = 0;

      // Get the residuals only
      this->fill_in_residuals_plastic(residuals, ipt, C);
      double maxres = residuals.max();

      while (maxres > Plastic_Newton_Solver_Tolerance);
      {
        // Initialize and compute the jacobian.
        jacobian.initialise(0.0);
        fill_in_jacobian_plastic(residuals, jacobian, ipt, C);

        // Solve for the step
        DoubleVector resid(residuals);
        jacobian.solve(resid, residuals);

        // update values
        double* dx_pt = residuals.values_pt();
        for (unsigned i = 0; i < this->get_num_plastic_dofs(ipt); i++)
        {
          *Plastic_dof_data_pt[ipt][i] -= dx_pt[i];
        }

        // Compute the residual to determin the error
        this->fill_in_residuals_plastic(residuals, ipt, C);
        maxres = residuals.max();

        nIter++;
      }

      delete lin_alg_dist_pt;
      delete solver_pt;
      delete communicator_pt;
    }

    // =========================================================================
    /// \short This function does the actual computation of the time derivatives
    /// of the individual plastic parameters and computes a residual and its
    /// associated jacobian from that. The residuals are defined as
    ///     \dot{inv_fp}_time_stepper - \dot{inv_fp}_analytical
    ///     \dot{inv_bpks}_time_stepper - \dot{inv_bpks}_analytical
    ///     \dot{inv_bpcs}_time_stepper - \dot{inv_bpcs}_analytical
    ///     f(M) = R F(\lambda): The yield surface equation, Eq. (101)
    ///     R - R_predicted, where R is predicted by \ref NormalYieldRatioLaw
    ///
    /// \param[out] residuals The residuals computed.
    /// \param[out] jacobian The derivatives of the residuals w.r.t to the
    /// plastic quantities.
    /// \param[in] ipt The integration point.
    /// \param[in] C The deformation gradient tensor in the undeformed frame.
    /// \param[in] flag Wheter the jacobian shall be computed.
    // =========================================================================
    virtual void fill_in_generic_residual_and_jacobian_plastic(
      DoubleVector& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& ipt,
      const DenseMatrix<double>& C,
      const unsigned& flag);

    // =========================================================================
    /// \short A Convenience wrapper for the function \ref
    /// fill_in_generic_residual_and_jacobian_plastic which does not compute any
    /// derivative
    ///
    /// \param[out] residuals The residuals computed.
    /// \param[in] ipt The integration point.
    /// \param[in] C The deformation gradient tensor in the undeformed frame.
    // =========================================================================
    void fill_in_residuals_plastic(DoubleVector& residuals,
                                   const unsigned& ipt,
                                   const DenseMatrix<double>& C)
    {
      fill_in_generic_residual_and_jacobian_plastic(
        residuals, GeneralisedElement::Dummy_matrix, ipt, C, 0);
    }

    // =========================================================================
    /// \short A function to compute both the redidual and the jacobian of the
    /// plastic data.
    /// If \ref Plastic_solve_use_fd is true, this function calles \ref
    /// fill_in_residuals_plastic and to compute the plastic residual and it's
    /// jacobian through finite difference. If Plastic_solve_use_fd is false, it
    /// calls \ref fill_in_generic_residual_and_jacobian_plastic.
    ///
    /// \param[out] residuals The residuals computed.
    /// \param[out] jacobian The derivatives of the residuals w.r.t to the
    /// plastic quantities.
    /// \param[in] ipt The integration point.
    /// \param[in] C The deformation gradient tensor in the undeformed frame.
    // =========================================================================
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
            *Plastic_fd_jacobian_step_pt;

          fill_in_residuals_plastic(test_residuals, ipt, C);
          for (unsigned local_eqn = 0; local_eqn < get_num_plastic_dofs(ipt);
               local_eqn++)
          {
            jacobian(local_eqn, local_unknown) =
              (test_residuals[local_eqn] - residuals[local_eqn]) /
              (*Plastic_fd_jacobian_step_pt);
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

    // =========================================================================
    /// \short Initialised the plastic parameters after the global deformation
    /// has been updated
    /// This function initialises all plastic data, assuming there is no plastic
    /// deformation present. In the unsteady case this means, setting the
    /// derivatives of the plastic data to 0. In the steady case this results in
    /// either the last history value being copied, or, for zero history time
    /// steppers, the initial condition being applied, \see
    /// set_intial_condition.
    ///
    /// If the initialisation moves R outside its validity regime, it is moved
    /// back \see enforce_boundaries_of_r.
    ///
    /// \param[in] ipt The integrationpoint to operate on.
    // =========================================================================
    void initialise_solve(const unsigned int ipt);

    // =========================================================================
    /// \short Sets the values of the current plastic variables to their
    /// respective initial values
    ///
    /// The initial values are unit matrices for invFp, invBpks and invBpcs, Re
    /// for R, and 0 for \lambda.
    ///
    /// \param[in] ipt The integrationpoint to operate on.
    // =========================================================================
    void set_intial_condition(const unsigned int ipt);

    // =========================================================================
    /// \short Checks that R is within the range [Re, 1.0] and clamps it to the
    /// range, if it is not.
    ///
    /// \param[in] ipt The integrationpoint to operate on.
    // =========================================================================
    void enforce_boundaries_of_r(const unsigned& ipt)
    {
      const double r = get_r(ipt);
      const double re =
        this->Plastic_consitutive_law_pt->normal_yield_ratio_law_pt->get_re();
      if (r < re)
      {
        set_r(ipt, re);
      }
      else if (r > 1.0)
      {
        set_r(ipt, 1.0);
      }
    }

    /*!
     * \brief Determins, if there is plastic deformation and the plastic solve
     * routine has to be called.
     * \details Checks if the elastic stress is contained by the yield surface
     * and if the plastic multiplier is positive.
     */
    bool is_there_plastic_deformation(const unsigned int ipt);

    ////////////////////////////////////////////////////////////////////////////
    // Mathematical functions involved in computing the plastic residuals or
    // similar
    ////////////////////////////////////////////////////////////////////////////

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

    void compute_mandellike_kinematic_hardening(
      const DenseMatrix<double>& invBpks, DenseMatrix<double>& bar_Mk)
    {
      return compute_mandellike_kinematic_hardening(
        invBpks,
        bar_Mk,
        PlasticEquationsBase<DIM>::Dummy_rankfourtensor,
        false);
    }

    void compute_mandellike_elastic_core(
      const DenseMatrix<double>& invBpcs,
      DenseMatrix<double>& bar_Mc,
      RankFourTensor<double>& dbar_McdinvBpcs,
      bool computeDerivative);

    void compute_mandellike_elastic_core(const DenseMatrix<double>& invBpcs,
                                         DenseMatrix<double>& bar_Mc)
    {
      return compute_mandellike_elastic_core(
        invBpcs,
        bar_Mc,
        PlasticEquationsBase<DIM>::Dummy_rankfourtensor,
        false);
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

    ////////////////////////////////////////////////////////////////////////////
    // Functions for constructing the plastic data below                      //
    ////////////////////////////////////////////////////////////////////////////

    void construct_plastic_data()
    {
      resize_plastic_dof_numbers();

      construct_inv_fp_internal_data();
      construct_invBpks_internal_data();
      construct_invBpcs_internal_data();
      construct_r_internal_data();
      construct_lambda_internal_data();

      // We assign the equation numbers now because the user will likely want
      // all plastic data unpinned. If they pin any then they will need to call
      // assign_plastic_eqn_numbers again
      assign_plastic_eqn_numbers();
    }

    // Change to not have an argument - in fact if we're always building all
    // plastic data then we don't need thi Collapse all building plastic data
    // into a single function
    void resize_plastic_dof_numbers()
    {
      if (Plastic_dof_nunbers_has_been_resized) return;

      const unsigned nipt = this->integral_pt()->nweight();

      Plastic_data_pt.resize(nipt);
      Plastic_data_eqn_number.resize(nipt);
      Plastic_dof_data_pt.resize(nipt);
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Plastic_data_pt[ipt].resize(NUMBER_OF_PLASTIC_VARIABLE_TYPES);
        Plastic_data_eqn_number[ipt].resize(NUMBER_OF_PLASTIC_VARIABLE_TYPES);
      }
      Plastic_dof_nunbers_has_been_resized = true;
    }

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
            Plastic_data_eqn_number[ipt][data_type][i] = eqn_count++;
            Plastic_dof_data_pt[ipt].push_back(data_pt->value_pt(i));
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

    void construct_invBpcs_internal_data()
    {
      const unsigned nipt = this->integral_pt()->nweight();
      for (unsigned ipt = 0; ipt < nipt; ipt++)
      {
        Data* data_pt = new Data(DIM * DIM);
        Plastic_data_pt[ipt][invBpcs_INDEX] = data_pt;
        (void)this->add_internal_data(data_pt, false);
        for (unsigned i = 0; i < DIM * DIM; i++)
        {
          // Pin the plastic degree of freedom
          Plastic_data_pt[ipt][invBpcs_INDEX]->pin(i);
          // But we have to initialise the eqn number to something so it may as
          // well be a safe value
          Plastic_data_eqn_number[ipt][invBpcs_INDEX].push_back(-1);
          // We initialise the plastic deformation gradient tensors to the
          // identity
          if (i % (DIM + 1) == 0)
          {
            data_pt->set_value(i, 1.0);
          }
        }
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
        // But we have to initialise the eqn number to something so it may as
        // well be a safe value
        Plastic_data_eqn_number[ipt][Lambda_INDEX].push_back(-1);
        data_pt->set_value(0, 0.0);
      }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Functions used by the refineable elements
    ////////////////////////////////////////////////////////////////////////////

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
        serialise_plastic_data(ipt_data[ipt], ipt, t);
      }

      const unsigned n_data = ipt_data[0].size();
#ifdef PARANOID
      for (unsigned ipt = 1; ipt < n_ipt; ipt++)
      {
        if (ipt_data[ipt].size() != n_data)
        {
          throw OomphLibError("Size of data at ipts does not match",
                              OOMPH_EXCEPTION_LOCATION,
                              OOMPH_CURRENT_FUNCTION);
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

    ////////////////////////////////////////////////////////////////////////////
    // Access functions to plastic and associated data
    ////////////////////////////////////////////////////////////////////////////

    // Get the number of plastic data at the given ipt which is not pinned and
    // is to be solved for
    unsigned get_num_plastic_dofs(const unsigned& ipt)
    {
      return Plastic_dof_data_pt[ipt].size();
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

    unsigned plastic_invBpcs_eqn_number(const unsigned& ipt,
                                        const unsigned& i,
                                        const unsigned& j) const
    {
      return Plastic_data_eqn_number[ipt][invBpcs_INDEX][i * DIM + j];
    }

    unsigned plastic_r_eqn_number(const unsigned& ipt) const
    {
      return Plastic_data_eqn_number[ipt][R_INDEX][0];
    }

    unsigned plastic_lambda_eqn_number(const unsigned& ipt) const
    {
      return Plastic_data_eqn_number[ipt][Lambda_INDEX][0];
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
          invBpks(i, j) =
            Plastic_data_pt[ipt][invBpks_INDEX]->value(t, i * DIM + j);
        }
      }
    }

    /*!
     * \brief returns the time derivative of the kinematic hardening deformation
     * gradient as a matrix
     * \details Uses the data's time stepper to compute the derivative
     */
    void get_dot_invBpks_matrix(const unsigned& ipt,
                                DenseMatrix<double>& dotinvBpks)
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
            delta_invBpks(i, j) =
              get_invBpks(0, ipt, i, j) - get_invBpks(1, ipt, i, j);
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

    void get_dot_or_delta_invBpks_matrix(
      const unsigned& ipt, DenseMatrix<double>& dot_or_delta_invBpks)
    {
      const TimeStepper* invBpks_time_stepper_pt =
        Plastic_data_pt[ipt][invBpks_INDEX]->time_stepper_pt();

      if (invBpks_time_stepper_pt->is_steady())
      {
        return get_delta_invBpks_matrix(ipt, dot_or_delta_invBpks);
      }

      return get_dot_invBpks_matrix(ipt, dot_or_delta_invBpks);
    }

    double get_invBpcs(const unsigned& ipt,
                       const unsigned& i,
                       const unsigned& j) const
    {
      return get_invBpcs(0, ipt, i, j);
    }

    double get_invBpcs(const unsigned& t,
                       const unsigned& ipt,
                       const unsigned& i,
                       const unsigned& j) const
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      return Plastic_data_pt[ipt][invBpcs_INDEX]->value(t, i * DIM + j);
    }

    void get_invBpcs_matrix(const unsigned int& ipt,
                            DenseMatrix<double>& invBpcs)
    {
      get_invBpcs_matrix(0, ipt, invBpcs);
    }

    void get_invBpcs_matrix(const unsigned int t,
                            const unsigned int& ipt,
                            DenseMatrix<double>& invBpcs)
    {
      invBpcs.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          invBpcs(i, j) =
            Plastic_data_pt[ipt][invBpcs_INDEX]->value(t, i * DIM + j);
        }
      }
    }

    /*!
     * \brief returns the time derivative of the elastic core gradient as a
     * matrix
     * \details Uses the data's time stepper to compute the derivative
     */
    void get_dot_invBpcs_matrix(const unsigned& ipt,
                                DenseMatrix<double>& dotinvBpcs)
    {
      dotinvBpcs.resize(DIM);

      for (unsigned int i = 0; i < DIM; i++)
      {
        for (unsigned int j = 0; j < DIM; j++)
        {
          dotinvBpcs(i, j) = this->dinternal_data_dt(
            Plastic_data_pt[ipt][invBpcs_INDEX], i * DIM + j);
        }
      }
    }

    /*!
     * \details Returns the difference between the current and the previous
     * timestep. If the timestepper does not contain history, it will return
     * the difference between the current timestep and the initial value.
     */
    void get_delta_invBpcs_matrix(const unsigned& ipt,
                                  DenseMatrix<double>& delta_invBpcs)
    {
      const TimeStepper* invBpcs_time_stepper_pt =
        Plastic_data_pt[ipt][invBpcs_INDEX]->time_stepper_pt();

      delta_invBpcs.resize(DIM);

      if (invBpcs_time_stepper_pt->ntstorage() > 1)
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_invBpcs(i, j) =
              get_invBpcs(0, ipt, i, j) - get_invBpcs(1, ipt, i, j);
          }
        }
      }
      else
      {
        for (unsigned int i = 0; i < DIM; i++)
        {
          for (unsigned int j = 0; j < DIM; j++)
          {
            delta_invBpcs(i, j) = get_invBpcs(0, ipt, i, j);

            if (i == j) delta_invBpcs(i, j) -= 1;
          }
        }
      }
    }

    void get_dot_or_delta_invBpcs_matrix(
      const unsigned& ipt, DenseMatrix<double>& dot_or_delta_invBpcs)
    {
      const TimeStepper* invBpcs_time_stepper_pt =
        Plastic_data_pt[ipt][invBpcs_INDEX]->time_stepper_pt();

      if (invBpcs_time_stepper_pt->is_steady())
      {
        return get_delta_invBpcs_matrix(ipt, dot_or_delta_invBpcs);
      }

      return get_dot_invBpcs_matrix(ipt, dot_or_delta_invBpcs);
    }

    double get_r(const unsigned& ipt) const
    {
      return get_r(0, ipt);
    }

    double get_r(const unsigned t, const unsigned& ipt) const
    {
      if (Plastic_data_pt[ipt][R_INDEX]->time_stepper_pt()->ntstorage() < t + 1)
        return this->Plastic_consitutive_law_pt->normal_yield_ratio_law_pt
          ->get_re();
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
    void set_invBpcs(const unsigned& ipt,
                     const unsigned& i,
                     const unsigned& j,
                     const double& val)
    {
      MatrixHelpers::check_matrix_indices<DIM>(i, j);
      Plastic_data_pt[ipt][invBpcs_INDEX]->set_value(i * DIM + j, val);
    }
    void set_r(const unsigned& ipt, const double& val)
    {
      Plastic_data_pt[ipt][R_INDEX]->set_value(0, val);
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

    double compute_internal_data_from_time_derivative(Data* data_pt,
                                                      const unsigned value_idx,
                                                      const double dxdt) const
    {
      // Number of timsteps (past & present)
      const TimeStepper* data_time_stepper_pt = data_pt->time_stepper_pt();
      const unsigned n_time = data_time_stepper_pt->ntstorage();

      // Initialize the return value
      double x0 = 0;

      // If the timestepper is not steady
      if (!data_time_stepper_pt->is_steady())
      {
        double sum_remaining_steps = 0;
        for (unsigned t = 1; t < n_time; t++)
        {
          sum_remaining_steps +=
            data_time_stepper_pt->weight(1, t) * data_pt->value(t, value_idx);
        }

        x0 = (dxdt - sum_remaining_steps) / data_time_stepper_pt->weight(1, 0);
      }

      return x0;
    }

  public:
    Data* plastic_data_pt(const unsigned& ipt, const unsigned& data_type) const
    {
      return Plastic_data_pt[ipt][data_type];
    }

    double* plastic_dof_data_pt(const unsigned& ipt, const unsigned& ndof) const
    {
      return Plastic_dof_data_pt[ipt][ndof];
    }

  protected:
    /// Stores, if finite difference should be used for the plastic solve.
    bool Plastic_solve_use_fd = false;

    /// The step used to determine the jacobian, if plastic Plastic_solve_use_fd
    /// is set to true.
    double* Plastic_fd_jacobian_step_pt =
      &FiniteElement::Default_fd_jacobian_step;

    /// A pointer to the plastic constitutive law
    PlasticConstitutiveLaw* Plastic_consitutive_law_pt;

    // We need a dummy double matrix for the plastic residual fill in
    DenseMatrix<double> unity;

    /// We use an enum to define the indices in the vectors at which the
    /// different types of plastic data are stored. This simplifies things
    /// considerably when searching for plastic data.
    enum Plastic_Variables_Indexes
    {
      Lambda_INDEX,
      invFp_INDEX,
      invBpks_INDEX,
      R_INDEX,
      invBpcs_INDEX,
      NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE,
      NUMBER_OF_PLASTIC_VARIABLE_TYPES = NUMBER_OF_PLASTIC_VARIABLES_TO_SOLVE
    };

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
          // Divide by diag entry to push forward G
          G_test(i, j) = G(i, j) / diag_entry;
          G_test(j, i) = G(j, i) / diag_entry;
        }
      }

      for (unsigned i = 0; i < DIM; i++)
      {
        for (unsigned j = i; j < DIM; j++)
        {
          const double saved_value = G_test(i, j);
          // Perturb G in the reference state, hence the default step needs
          // to be pushed forward
          G_test(i, j) += FiniteElement::Default_fd_jacobian_step / diag_entry;

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
      RefineableSolidElement::further_build();

      PlasticEquations<DIM>* cast_father_element_pt =
        dynamic_cast<PlasticEquations<DIM>*>(this->father_element_pt());

      this->Plastic_solve_use_fd =
        cast_father_element_pt->get_plastic_solve_by_fd();
      this->Plastic_consitutive_law_pt =
        cast_father_element_pt->plastic_constitutive_law_pt();

      // Set the plastic timestepper
      this->assign_plastic_timestepper_pt(
        cast_father_element_pt->plastic_data_pt(0, 0)->time_stepper_pt(), true);
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