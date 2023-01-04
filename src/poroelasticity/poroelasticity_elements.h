// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_POROELASTICITY_ELEMENTS_HEADER
#define OOMPH_POROELASTICITY_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "generic/elements.h"
#include "generic/shape.h"

#include "elasticity_tensor.h"

namespace oomph
{
  /// Class implementing the generic maths of the poroelasticity
  /// equations: linear elasticity coupled with Darcy equations (using
  /// Raviart-Thomas elements with both edge and internal degrees of freedom)
  template<unsigned DIM>
  class PoroelasticityEquations : public virtual FiniteElement
  {
  public:
    /// Source function pointer typedef
    typedef void (*SourceFctPt)(const double& time,
                                const Vector<double>& x,
                                Vector<double>& f);

    /// Mass source function pointer typedef
    typedef void (*MassSourceFctPt)(const double& time,
                                    const Vector<double>& x,
                                    double& f);

    /// Constructor
    PoroelasticityEquations()
      : Elasticity_tensor_pt(0),
        Force_solid_fct_pt(0),
        Force_fluid_fct_pt(0),
        Mass_source_fct_pt(0),
        Lambda_sq_pt(&Default_lambda_sq_value),
        Density_ratio_pt(&Default_density_ratio_value),
        K_inv_pt(&Default_k_inv_value),
        Alpha_pt(&Default_alpha_value),
        Porosity_pt(&Default_porosity_value)
    {
    }

    /// Access function for timescale ratio (nondim density)
    const double& lambda_sq() const
    {
      return *Lambda_sq_pt;
    }

    /// Access function for pointer to timescale ratio (nondim density)
    double*& lambda_sq_pt()
    {
      return Lambda_sq_pt;
    }

    /// Access function for the density ratio
    const double& density_ratio() const
    {
      return *Density_ratio_pt;
    }

    /// Access function for pointer to the density ratio
    double*& density_ratio_pt()
    {
      return Density_ratio_pt;
    }

    /// Access function for the nondim inverse permeability
    const double& k_inv() const
    {
      return *K_inv_pt;
    }

    /// Access function for pointer to the nondim inverse permeability
    double*& k_inv_pt()
    {
      return K_inv_pt;
    }

    /// Access function for alpha
    const double& alpha() const
    {
      return *Alpha_pt;
    }

    /// Access function for pointer to alpha
    double*& alpha_pt()
    {
      return Alpha_pt;
    }

    /// Access function for the porosity
    const double& porosity() const
    {
      return *Porosity_pt;
    }

    /// Access function for pointer to the porosity
    double*& porosity_pt()
    {
      return Porosity_pt;
    }

    /// Access function: Pointer to solid force function
    SourceFctPt& force_solid_fct_pt()
    {
      return Force_solid_fct_pt;
    }

    /// Access function: Pointer to solid force function (const version)
    SourceFctPt force_solid_fct_pt() const
    {
      return Force_solid_fct_pt;
    }

    /// Access function: Pointer to fluid force function
    SourceFctPt& force_fluid_fct_pt()
    {
      return Force_fluid_fct_pt;
    }

    /// Access function: Pointer to fluid force function (const version)
    SourceFctPt force_fluid_fct_pt() const
    {
      return Force_fluid_fct_pt;
    }

    /// Access function: Pointer to mass source function
    MassSourceFctPt& mass_source_fct_pt()
    {
      return Mass_source_fct_pt;
    }

    /// Access function: Pointer to mass source function (const version)
    MassSourceFctPt mass_source_fct_pt() const
    {
      return Mass_source_fct_pt;
    }

    /// Indirect access to the solid force function - returns 0 if no
    /// forcing function has been set
    void force_solid(const double& time,
                     const Vector<double>& x,
                     Vector<double>& b) const
    {
      // If no function has been set, return zero vector
      if (Force_solid_fct_pt == 0)
      {
        // Get spatial dimension of element
        unsigned n = dim();
        for (unsigned i = 0; i < n; i++)
        {
          b[i] = 0.0;
        }
      }
      else
      {
        (*Force_solid_fct_pt)(time, x, b);
      }
    }

    /// Indirect access to the fluid forcing function - returns 0 if no
    /// forcing function has been set
    void force_fluid(const double& time,
                     const Vector<double>& x,
                     Vector<double>& b) const
    {
      // If no function has been set, return zero vector
      if (Force_fluid_fct_pt == 0)
      {
        // Get spatial dimension of element
        unsigned n = dim();
        for (unsigned i = 0; i < n; i++)
        {
          b[i] = 0.0;
        }
      }
      else
      {
        (*Force_fluid_fct_pt)(time, x, b);
      }
    }

    /// Indirect access to the mass source function - returns 0 if no
    /// mass source function has been set
    void mass_source(const double& time,
                     const Vector<double>& x,
                     double& b) const
    {
      // If no function has been set, return zero vector
      if (Mass_source_fct_pt == 0)
      {
        b = 0.0;
      }
      else
      {
        (*Mass_source_fct_pt)(time, x, b);
      }
    }

    /// Return the pointer to the elasticity_tensor
    ElasticityTensor*& elasticity_tensor_pt()
    {
      return Elasticity_tensor_pt;
    }

    /// Access function to the entries in the elasticity tensor
    const double E(const unsigned& i,
                   const unsigned& j,
                   const unsigned& k,
                   const unsigned& l) const
    {
      return (*Elasticity_tensor_pt)(i, j, k, l);
    }

    /// Return the Cauchy stress tensor, as calculated
    /// from the elasticity tensor at specified local coordinate
    void get_stress(const Vector<double>& s, DenseMatrix<double>& sigma) const;

    /// Return the strain tensor
    void get_strain(const Vector<double>& s, DenseMatrix<double>& strain) const;

    /// Number of values required at node n
    virtual unsigned required_nvalue(const unsigned& n) const = 0;

    /// Return the nodal index of the n-th solid displacement unknown
    virtual unsigned u_index(const unsigned& n) const = 0;

    /// Return the equation number of the n-th edge (flux) degree of freedom
    virtual int q_edge_local_eqn(const unsigned& n) const = 0;

    /// Return the equation number of the n-th internal (moment) degree of
    /// freedom
    virtual int q_internal_local_eqn(const unsigned& n) const = 0;

    /// Return the nodal index at which the nth edge unknown is stored
    virtual unsigned q_edge_index(const unsigned& n) const = 0;

    /// Return the index of the internal data where the q_internal
    /// degrees of freedom are stored
    virtual unsigned q_internal_index() const = 0;

    /// Return the number of the node where the nth edge unknown is stored
    virtual unsigned q_edge_node_number(const unsigned& n) const = 0;

    /// Return the values of the edge (flux) degrees of freedom
    virtual double q_edge(const unsigned& n) const = 0;

    /// Return the values of the edge (flux) degrees of freedom at time
    /// history level t
    virtual double q_edge(const unsigned& t, const unsigned& n) const = 0;

    /// Return the values of the internal (moment) degrees of freedom
    virtual double q_internal(const unsigned& n) const = 0;

    /// Return the values of the internal (moment) degrees of freedom at
    /// time history level t
    virtual double q_internal(const unsigned& t, const unsigned& n) const = 0;

    /// Return the total number of computational basis functions for q
    virtual unsigned nq_basis() const = 0;

    /// Return the number of edge basis functions for q
    virtual unsigned nq_basis_edge() const = 0;

    /// Returns the local form of the q basis at local coordinate s
    virtual void get_q_basis_local(const Vector<double>& s,
                                   Shape& q_basis) const = 0;

    /// Returns the local form of the q basis and dbasis/ds at local coordinate
    /// s
    virtual void get_div_q_basis_local(const Vector<double>& s,
                                       Shape& div_q_basis_ds) const = 0;

    /// Returns the transformed basis at local coordinate s
    void get_q_basis(const Vector<double>& s, Shape& q_basis) const
    {
      const unsigned n_node = this->nnode();
      Shape psi(n_node, DIM);
      const unsigned n_q_basis = this->nq_basis();
      Shape q_basis_local(n_q_basis, DIM);
      this->get_q_basis_local(s, q_basis_local);
      (void)this->transform_basis(s, q_basis_local, psi, q_basis);
    }

    /// Returns the number of gauss points along each edge of the element
    virtual unsigned nedge_gauss_point() const = 0;

    /// Returns the nth gauss point along an edge
    virtual double edge_gauss_point(const unsigned& edge,
                                    const unsigned& n) const = 0;

    /// Returns the global coordinates of the nth gauss point along an edge
    virtual void edge_gauss_point_global(const unsigned& edge,
                                         const unsigned& n,
                                         Vector<double>& x) const = 0;

    /// Pin the nth internal q value
    virtual void pin_q_internal_value(const unsigned& n) = 0;

    /// Return the equation number of the n-th pressure degree of freedom
    virtual int p_local_eqn(const unsigned& n) const = 0;

    /// Return the nth pressure value
    virtual double p_value(unsigned& n) const = 0;

    /// Return the total number of pressure basis functions
    virtual unsigned np_basis() const = 0;

    /// Return the pressure basis
    virtual void get_p_basis(const Vector<double>& s, Shape& p_basis) const = 0;

    /// Pin the nth pressure value
    virtual void pin_p_value(const unsigned& n, const double& p) = 0;

    /// Scale the edge basis to allow arbitrary edge mappings
    virtual void scale_basis(Shape& basis) const = 0;

    /// Performs a div-conserving transformation of the vector basis
    /// functions from the reference element to the actual element
    double transform_basis(const Vector<double>& s,
                           const Shape& q_basis_local,
                           Shape& psi,
                           DShape& dpsi,
                           Shape& q_basis) const;

    /// Performs a div-conserving transformation of the vector basis
    /// functions from the reference element to the actual element
    double transform_basis(const Vector<double>& s,
                           const Shape& q_basis_local,
                           Shape& psi,
                           Shape& q_basis) const
    {
      const unsigned n_node = this->nnode();
      DShape dpsi(n_node, DIM);
      return transform_basis(s, q_basis_local, psi, dpsi, q_basis);
    }

    /// Fill in contribution to residuals for the Darcy equations
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      this->fill_in_generic_residual_contribution(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Fill in the Jacobian matrix for the Newton method
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      this->fill_in_generic_residual_contribution(residuals, jacobian, 1);
    }

    /// Calculate the FE representation of u
    void interpolated_u(const Vector<double>& s, Vector<double>& disp) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      for (unsigned i = 0; i < DIM; i++)
      {
        // Index at which the nodal value is stored
        unsigned u_nodal_index = u_index(i);

        // Initialise value of u
        disp[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          disp[i] += nodal_value(l, u_nodal_index) * psi[l];
        }
      }
    }

    /// Calculate the FE representation of the i-th component of u
    double interpolated_u(const Vector<double>& s, const unsigned& i) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Get nodal index at which i-th velocity is stored
      unsigned u_nodal_index = u_index(i);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += nodal_value(l, u_nodal_index) * psi[l];
      }

      return (interpolated_u);
    }

    /// Calculate the FE representation of q
    void interpolated_q(const Vector<double>& s, Vector<double>& u) const
    {
      unsigned n_q_basis = nq_basis();
      unsigned n_q_basis_edge = nq_basis_edge();

      Shape q_basis(n_q_basis, DIM);

      get_q_basis(s, q_basis);
      for (unsigned i = 0; i < DIM; i++)
      {
        u[i] = 0.0;
        for (unsigned l = 0; l < n_q_basis_edge; l++)
        {
          u[i] += q_edge(l) * q_basis(l, i);
        }
        for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
        {
          u[i] += q_internal(l - n_q_basis_edge) * q_basis(l, i);
        }
      }
    }

    /// Calculate the FE representation of the i-th component of q
    double interpolated_q(const Vector<double>& s, const unsigned i) const
    {
      unsigned n_q_basis = nq_basis();
      unsigned n_q_basis_edge = nq_basis_edge();

      Shape q_basis(n_q_basis, DIM);

      get_q_basis(s, q_basis);
      double q_i = 0.0;
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        q_i += q_edge(l) * q_basis(l, i);
      }
      for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
      {
        q_i += q_internal(l - n_q_basis_edge) * q_basis(l, i);
      }

      return q_i;
    }

    /// Calculate the FE representation of div u
    void interpolated_div_q(const Vector<double>& s, double& div_q) const
    {
      // Zero the divergence
      div_q = 0;

      // Get the number of nodes, q basis function, and q edge basis functions
      unsigned n_node = nnode();
      const unsigned n_q_basis = nq_basis();
      const unsigned n_q_basis_edge = nq_basis_edge();

      // Storage for the divergence basis
      Shape div_q_basis_ds(n_q_basis);

      // Storage for the geometric basis and it's derivatives
      Shape psi(n_node);
      DShape dpsi(n_node, DIM);

      // Call the geometric shape functions and their derivatives
      this->dshape_local(s, psi, dpsi);

      // Storage for the inverse of the geometric jacobian (just so we can call
      // the local to eulerian mapping)
      DenseMatrix<double> inverse_jacobian(DIM);

      // Get the determinant of the geometric mapping
      double det = local_to_eulerian_mapping(dpsi, inverse_jacobian);

      // Get the divergence basis (wrt local coords) at local coords s
      get_div_q_basis_local(s, div_q_basis_ds);

      // Add the contribution to the divergence from the edge basis functions
      for (unsigned l = 0; l < n_q_basis_edge; l++)
      {
        div_q += 1.0 / det * div_q_basis_ds(l) * q_edge(l);
      }

      // Add the contribution to the divergence from the internal basis
      // functions
      for (unsigned l = n_q_basis_edge; l < n_q_basis; l++)
      {
        div_q += 1.0 / det * div_q_basis_ds(l) * q_internal(l - n_q_basis_edge);
      }
    }

    /// Calculate the FE representation of div q and return it
    double interpolated_div_q(const Vector<double>& s)
    {
      // Temporary storage for div u
      double div_q = 0;

      // Get the intepolated divergence
      interpolated_div_q(s, div_q);

      // Return it
      return div_q;
    }

    /// Calculate the FE representation of p
    void interpolated_p(const Vector<double>& s, double& p) const
    {
      // Get the number of p basis functions
      unsigned n_p_basis = np_basis();

      // Storage for the p basis
      Shape p_basis(n_p_basis);

      // Call the p basis
      get_p_basis(s, p_basis);

      // Zero the pressure
      p = 0;

      // Add the contribution to the pressure from each basis function
      for (unsigned l = 0; l < n_p_basis; l++)
      {
        p += p_value(l) * p_basis(l);
      }
    }

    /// Calculate the FE representation of p and return it
    double interpolated_p(const Vector<double>& s) const
    {
      // Temporary storage for p
      double p = 0;

      // Get the interpolated pressure
      interpolated_p(s, p);

      // Return it
      return p;
    }

    /// du/dt at local node n
    double du_dt(const unsigned& n, const unsigned& i) const
    {
      // Get the timestepper
      TimeStepper* time_stepper_pt = node_pt(n)->time_stepper_pt();

      // Storage for the derivative - initialise to 0
      double du_dt = 0.0;

      // If we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // Get the nodal index
        const unsigned u_nodal_index = u_index(i);

        // Get the number of values required to represent history
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Loop over history values
        for (unsigned t = 0; t < n_time; t++)
        {
          // Add the contribution to the derivative
          du_dt +=
            time_stepper_pt->weight(1, t) * nodal_value(t, n, u_nodal_index);
        }
      }

      return du_dt;
    }

    /// d^2u/dt^2 at local node n
    double d2u_dt2(const unsigned& n, const unsigned& i) const
    {
      // Get the timestepper
      TimeStepper* time_stepper_pt = node_pt(n)->time_stepper_pt();

      // Storage for the derivative - initialise to 0
      double d2u_dt2 = 0.0;

      // If we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // Get the nodal index
        const unsigned u_nodal_index = u_index(i);

        // Get the number of values required to represent history
        const unsigned n_time = time_stepper_pt->ntstorage();

        // Loop over history values
        for (unsigned t = 0; t < n_time; t++)
        {
          // Add the contribution to the derivative
          d2u_dt2 +=
            time_stepper_pt->weight(2, t) * nodal_value(t, n, u_nodal_index);
        }
      }

      return d2u_dt2;
    }

    /// dq_edge/dt for the n-th edge degree of freedom
    double dq_edge_dt(const unsigned& n) const
    {
      unsigned node_num = q_edge_node_number(n);

      // get the timestepper
      TimeStepper* time_stepper_pt = node_pt(node_num)->time_stepper_pt();

      // storage for the derivative - initialise to 0
      double dq_dt = 0.0;

      // if we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // get the number of values required to represent history
        const unsigned n_time = time_stepper_pt->ntstorage();

        // loop over history values
        for (unsigned t = 0; t < n_time; t++)
        {
          // add the contribution to the derivative
          dq_dt += time_stepper_pt->weight(1, t) * q_edge(t, n);
        }
      }

      return dq_dt;
    }

    /// dq_internal/dt for the n-th internal degree of freedom
    double dq_internal_dt(const unsigned& n) const
    {
      // get the internal data index for q
      unsigned internal_index = q_internal_index();

      // get the timestepper
      TimeStepper* time_stepper_pt =
        internal_data_pt(internal_index)->time_stepper_pt();

      // storage for the derivative - initialise to 0
      double dq_dt = 0.0;

      // if we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // get the number of values required to represent history
        const unsigned n_time = time_stepper_pt->ntstorage();

        // loop over history values
        for (unsigned t = 0; t < n_time; t++)
        {
          // add the contribution to the derivative
          dq_dt += time_stepper_pt->weight(1, t) * q_internal(t, n);
        }
      }

      return dq_dt;
    }

    /// Set the timestepper of the q internal data object
    void set_q_internal_timestepper(TimeStepper* const time_stepper_pt)
    {
      unsigned q_index = q_internal_index();

      this->internal_data_pt(q_index)->set_time_stepper(time_stepper_pt, false);
    }

    unsigned self_test()
    {
      return 0;
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    /// Output FE representation of soln: x,y,u1,u2,div_q,p at
    /// Nplot^DIM plot points
    void output(std::ostream& outfile, const unsigned& nplot);

    /// Output FE representation of exact soln: x,y,u1,u2,div_q,p at
    /// Nplot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// Output FE representation of exact soln: x,y,u1,u2,div_q,p at
    /// Nplot^DIM plot points. Unsteady version
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt);

    /// Compute the error between the FE solution and the exact solution
    /// using the H(div) norm for q and L^2 norm for p
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       Vector<double>& error,
                       Vector<double>& norm);

    /// Compute the error between the FE solution and the exact solution
    /// using the H(div) norm for q and L^2 norm for p. Unsteady version
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       Vector<double>& error,
                       Vector<double>& norm);

  protected:
    /// Returns the geometric basis, and the q, p and divergence basis
    /// functions and test functions at local coordinate s
    virtual double shape_basis_test_local(const Vector<double>& s,
                                          Shape& psi,
                                          DShape& dpsi,
                                          Shape& u_basis,
                                          Shape& u_test,
                                          DShape& du_basis_dx,
                                          DShape& du_test_dx,
                                          Shape& q_basis,
                                          Shape& q_test,
                                          Shape& p_basis,
                                          Shape& p_test,
                                          Shape& div_q_basis_ds,
                                          Shape& div_q_test_ds) const = 0;

    /// Returns the geometric basis, and the q, p and divergence basis
    /// functions and test functions at integration point ipt
    virtual double shape_basis_test_local_at_knot(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsi,
      Shape& u_basis,
      Shape& u_test,
      DShape& du_basis_dx,
      DShape& du_test_dx,
      Shape& q_basis,
      Shape& q_test,
      Shape& p_basis,
      Shape& p_test,
      Shape& div_q_basis_ds,
      Shape& div_q_test_ds) const = 0;

    // fill in residuals and, if flag==true, jacobian
    virtual void fill_in_generic_residual_contribution(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, bool flag);

    /// Pointer to the elasticity tensor
    ElasticityTensor* Elasticity_tensor_pt;

  private:
    /// Pointer to solid source function
    SourceFctPt Force_solid_fct_pt;

    /// Pointer to fluid source function
    SourceFctPt Force_fluid_fct_pt;

    /// Pointer to the mass source function
    MassSourceFctPt Mass_source_fct_pt;

    /// Timescale ratio (non-dim. density)
    double* Lambda_sq_pt;

    /// Density ratio
    double* Density_ratio_pt;

    /// 1/k
    double* K_inv_pt;

    /// Alpha
    double* Alpha_pt;

    /// Porosity
    double* Porosity_pt;

    /// Static default value for timescale ratio (1.0 -- for natural scaling)
    static double Default_lambda_sq_value;

    /// Static default value for the density ratio
    static double Default_density_ratio_value;

    /// Static default value for 1/k
    static double Default_k_inv_value;

    /// Static default value for alpha
    static double Default_alpha_value;

    /// Static default value for the porosity
    static double Default_porosity_value;
  };

} // namespace oomph

#endif
