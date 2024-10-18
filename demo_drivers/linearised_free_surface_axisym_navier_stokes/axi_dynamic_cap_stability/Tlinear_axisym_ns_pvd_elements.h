#ifndef TLINEAR_AXISYM_NS_PVD_ELEMENTS_HEADER
#define TLINEAR_AXISYM_NS_PVD_ELEMENTS_HEADER

#include "Tfull_linearised_axisym_navier_stokes_elements.h"
#include "decomposed_pvd_elements.h"
#include "debug_jacobian_elements.h"

namespace oomph
{
  class TLinearisedAxisymNSPVDElement
    : public virtual TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations,
      public virtual DecomposedPVDEquations<2>,
      public virtual DebugJacobianFiniteElement
  {
  private:
    TimeStepper* Time_stepper_pt;

  public:
    /// Constructor: Call constructors for TElement and
    /// DecomposedLinearElasticity equations
    TLinearisedAxisymNSPVDElement()
      : TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations(),
        DecomposedPVDEquations<2>(),
        Time_stepper_pt(0)
    {
    }

    /// Broken copy constructor
    TLinearisedAxisymNSPVDElement(const TLinearisedAxisymNSPVDElement& dummy) = delete;

    void set_time_stepper_pt(TimeStepper* const& passed_time_stepper_pt)
    {
      Time_stepper_pt = passed_time_stepper_pt;
    }

    inline virtual unsigned required_nvalue(const unsigned& n) const
    {
      return TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::
        required_nvalue(n);
    }

    virtual inline unsigned u_index_pvd(const unsigned& n,
                                                      const unsigned& i,
                                                      const unsigned& j) const
    {
      return this->xhat_index_lin_axi_nst(n, i * 2 + j);
    }


    virtual inline unsigned u_index_lin_axi_nst(const unsigned& i) const
    {
      return i + 4;
    }

    virtual inline unsigned xhat_index_lin_axi_nst(const unsigned& n,
                                                   const unsigned& i) const
    {
      return i;
    }

    double p_lin_axi_nst(const unsigned& n_p, const unsigned& i) const
    {
      return this->nodal_value(Pconv[n_p], 10 + i);
    }

    virtual int p_local_eqn(const unsigned& n, const unsigned& i)
    {
      return this->nodal_local_eqn(Pconv[n], 10 + i);
    }

    virtual void pin_pressure(const unsigned& j)
    {
      for (unsigned n_p = 0; n_p < 3; n_p++)
      {
        this->node_pt(Pconv[n_p])->pin(10 + j);
      }
    }

    virtual void pin_pressure(const unsigned& n, const unsigned& j)
    {
      if (n < 3)
      {
        this->node_pt(Pconv[n])->pin(10 + j);
      }
      else
      {
        throw OomphLibError("Capillary number has not been set",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }

    void pin()
    {
      for (unsigned i = 0; i < 6; i++)
      {
        for (unsigned j = 0; j < 6; j++)
        {
          this->node_pt(i)->pin(4 + j);
        }
      }
      pin_pressure(0);
      pin_pressure(1);
    }

    /// Compute the element's residual Vector
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // and using a dummy matrix argument
      this->fill_in_generic_residual_contribution_lin_axi_nst(
        residuals,
        GeneralisedElement::Dummy_matrix,
        GeneralisedElement::Dummy_matrix,
        0);

      this->fill_in_generic_contribution_to_residuals_pvd(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// The jacobian is calculated by finite differences by default,
    /// We need only to take finite differences w.r.t. positional variables
    /// For this element
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      // Add the contribution to the residuals
      // this->fill_in_generic_contribution_to_residuals_pvd(
      //  residuals, jacobian, 1);
    }

    /// Add the element's contribution to its residuals vector,
    /// jacobian matrix and mass matrix
    // void fill_in_contribution_to_jacobian_and_mass_matrix(
    //  Vector<double>& residuals,
    //  DenseMatrix<double>& jacobian,
    //  DenseMatrix<double>& mass_matrix)
    //{
    //  FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
    //  const unsigned n_dofs = this->ndof();
    //  Vector<double> dummy_residuals(n_dofs);
    //  DenseMatrix<double> dummy_jacobian(n_dofs, n_dofs);
    //  // Call the generic routine with the flag set to 2
    //  fill_in_generic_residual_contribution_lin_axi_nst(
    //    dummy_residuals, dummy_jacobian, mass_matrix, 2);
    //}

    /// Add the elemental contribution to the jacobian and mass matrices
    /// and the residuals vector. Note that
    /// this function will NOT initialise the residuals vector or the jacobian
    /// matrix. It must be called after the residuals vector and
    /// jacobian matrix have been initialised to zero. The default
    /// is to use finite differences to calculate the jacobian
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Add the contribution to the residuals
      this->fill_in_contribution_to_residuals(residuals);

      // Allocate storage for the full residuals (residuals of entire
      // element)
      unsigned n_dof = this->ndof();
      Vector<double> full_residuals(n_dof);
      // Get the residuals for the entire element
      this->get_residuals(full_residuals);

      // Make the timestepper steady
      Time_stepper_pt->make_steady();

      // Use the fill in contribution to jacobian
      fill_in_contribution_to_jacobian(full_residuals, jacobian);

      Time_stepper_pt->undo_make_steady();

      DenseMatrix<double> unsteady_jacobian(n_dof, n_dof);
      fill_in_contribution_to_jacobian(full_residuals, unsteady_jacobian);

      const double dt = Time_stepper_pt->time_pt()->dt();
      for (unsigned i = 0; i < n_dof; i++)
      {
        for (unsigned j = 0; j < n_dof; j++)
        {
          /// The 2/3 is due to the BDF<2> scheme.
          mass_matrix(j, i) +=
            (2.0 / 3.0) * dt * (-unsteady_jacobian(j, i) + jacobian(j, i));
        }
      }
    }

    virtual void output(std::ostream& outfile)
    {
      TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::output(
        outfile);
    }

    virtual void output(std::ostream& outfile, const unsigned& n_plot)
    {
      TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::output(outfile,
                                                                        n_plot);
    }

    virtual void output(FILE* file_pt)
    {
      TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::output(
        file_pt);
    }

    virtual void output(FILE* file_pt, const unsigned& n_plot)
    {
      TaylorHoodFullLinearisedAxisymmetricNavierStokesEquations::output(file_pt,
                                                                        n_plot);
    }
  };

  //=======================================================================
  /// Face geometry for the TDecomposedLinearElasticityElement elements: The
  /// spatial dimension of the face elements is one lower than that of the bulk
  /// element but they have the same number of points along their 1D edges.
  //=======================================================================
  template<>
  class FaceGeometry<TLinearisedAxisymNSPVDElement> : public virtual TElement<1, 3>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : TElement<1, 3>() {}
  };

} // namespace oomph

#endif
