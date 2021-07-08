// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision$
// LIC//
// LIC// $LastChangedDate$
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Header file for scalar advection elements

#ifndef OOMPH_EULER_ELEMENTS_HEADER
#define OOMPH_EULER_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#include "flux_transport_elements.h"
#include "generic/Qspectral_elements.h"
#include "generic/dg_elements.h"

namespace oomph
{
  //==============================================================
  /// Base class for Euler equations
  //=============================================================
  template<unsigned DIM>
  class EulerEquations : public FluxTransportEquations<DIM>
  {
    //\short Default value for gamma (initialised to 1.4)
    static double Default_Gamma_Value;

    // Pointer to the specific gas constant
    double* Gamma_pt;

    /// \short Storage for the average primitive values
    double* Average_prim_value;

    /// \short Storage for the approximated gradients
    double* Average_gradient;

  protected:
    /// \short DIM momentum-components, a density and an energy are transported
    inline unsigned nflux() const
    {
      return DIM + 2;
    }

    /// \short Return the flux as a function of the unknown
    void flux(const Vector<double>& u, DenseMatrix<double>& f);

    /// \short Return the flux derivatives as a function of the unknowns
    // void dflux_du(const Vector<double> &u, RankThreeTensor<double> &df_du);

  public:
    /// Constructor
    EulerEquations() :
      FluxTransportEquations<DIM>(), Average_prim_value(0), Average_gradient(0)
    {
      // Set the default value of gamma
      Gamma_pt = &Default_Gamma_Value;
    }

    /// Destructor
    virtual ~EulerEquations()
    {
      if (this->Average_prim_value != 0)
      {
        delete[] Average_prim_value;
        Average_prim_value = 0;
      }
      if (this->Average_gradient != 0)
      {
        delete[] Average_gradient;
        Average_gradient = 0;
      }
    }

    /// Calculate the pressure from the unknowns
    double pressure(const Vector<double>& u) const;

    /// \short Return the value of gamma
    double gamma() const
    {
      return *Gamma_pt;
    }

    /// \short Access function for the pointer to gamma
    double*& gamma_pt()
    {
      return Gamma_pt;
    }

    /// \short Access function for the pointer to gamma (const version)
    double* const& gamma_pt() const
    {
      return Gamma_pt;
    }

    ///\short The number of unknowns at each node is the number of flux
    /// components
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return DIM + 2;
    }

    /// Compute the error and norm of solution integrated over the element
    /// Does not plot the error in the outfile
    void compute_error(
      std::ostream& outfile,
      FiniteElement::UnsteadyExactSolutionFctPt exact_solution_pt,
      const double& t,
      Vector<double>& error,
      Vector<double>& norm)
    {
      // Find the number of fluxes
      const unsigned n_flux = this->nflux();
      // Find the number of nodes
      const unsigned n_node = this->nnode();
      // Storage for the shape function and derivatives of shape function
      Shape psi(n_node);
      DShape dpsidx(n_node, DIM);

      // Find the number of integration points
      unsigned n_intpt = this->integral_pt()->nweight();

      error.resize(n_flux);
      norm.resize(n_flux);
      for (unsigned i = 0; i < n_flux; i++)
      {
        error[i] = 0.0;
        norm[i] = 0.0;
      }

      Vector<double> s(DIM);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the shape functions at the knot
        double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

        // Get the integral weight
        double W = this->integral_pt()->weight(ipt) * J;

        // Storage for the local functions
        Vector<double> interpolated_x(DIM, 0.0);
        Vector<double> interpolated_u(n_flux, 0.0);

        // Loop over the shape functions
        for (unsigned l = 0; l < n_node; l++)
        {
          // Locally cache the shape function
          const double psi_ = psi(l);
          for (unsigned i = 0; i < DIM; i++)
          {
            interpolated_x[i] += this->nodal_position(l, i) * psi_;
          }

          for (unsigned i = 0; i < n_flux; i++)
          {
            // Calculate the velocity and tangent vector
            interpolated_u[i] += this->nodal_value(l, i) * psi_;
          }
        }

        // Now get the exact solution at this value of x and t
        Vector<double> exact_u(n_flux, 0.0);
        (*exact_solution_pt)(t, interpolated_x, exact_u);

        // Loop over the unknowns
        for (unsigned i = 0; i < n_flux; i++)
        {
          error[i] += pow((interpolated_u[i] - exact_u[i]), 2.0) * W;
          norm[i] += interpolated_u[i] * interpolated_u[i] * W;
        }
      }
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      unsigned nplot = 5;
      output(outfile, nplot);
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    void allocate_memory_for_averages()
    {
      // Find the number of fluxes
      const unsigned n_flux = this->nflux();
      // Resize the memory if necessary
      if (this->Average_prim_value == 0)
      {
        this->Average_prim_value = new double[n_flux];
      }
      // Will move this one day
      if (this->Average_gradient == 0)
      {
        this->Average_gradient = new double[n_flux * DIM];
      }
    }

    /// \short Return access to the average gradient
    double& average_gradient(const unsigned& i, const unsigned& j)
    {
      if (Average_gradient == 0)
      {
        throw OomphLibError("Averages not calculated yet",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      return Average_gradient[DIM * i + j];
    }

    /// \short Return the average values
    double& average_prim_value(const unsigned& i)
    {
      if (Average_prim_value == 0)
      {
        throw OomphLibError("Averages not calculated yet",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      return Average_prim_value[i];
    }
  };

  template<unsigned DIM, unsigned NNODE_1D>
  class QSpectralEulerElement :
    public virtual QSpectralElement<DIM, NNODE_1D>,
    public virtual EulerEquations<DIM>
  {
  public:
    ///\short  Constructor: Call constructors for QElement and
    /// Advection Diffusion equations
    QSpectralEulerElement() :
      QSpectralElement<DIM, NNODE_1D>(), EulerEquations<DIM>()
    {
    }

    /// Broken copy constructor
    QSpectralEulerElement(const QSpectralEulerElement<DIM, NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("QSpectralEulerElement");
    }

    /// Broken assignment operator
    // Commented out broken assignment operator because this can lead to a
    // conflict warning when used in the virtual inheritence hierarchy.
    // Essentially the compiler doesn't realise that two separate
    // implementations of the broken function are the same and so, quite
    // rightly, it shouts.
    /*void operator=(
     const QSpectralEulerElement<DIM,NNODE_1D>&)
     {
      BrokenCopy::broken_assign("QSpectralEulerElement");
      }*/

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      EulerEquations<DIM>::output(outfile);
    }

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      EulerEquations<DIM>::output(outfile, n_plot);
    }

    /*/// \short C-style output function:
   ///  x,y,u   or    x,y,z,u
   void output(FILE* file_pt)
    {
     EulerEquations<NFLUX,DIM>::output(file_pt);
    }

   ///  \short C-style output function:
   ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
   void output(FILE* file_pt, const unsigned &n_plot)
    {
     EulerEquations<NFLUX,DIM>::output(file_pt,n_plot);
    }

   /// \short Output function for an exact solution:
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   FiniteElement::SteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     EulerEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,exact_soln_pt);}


   /// \short Output function for a time-dependent exact solution.
   ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
   /// (Calls the steady version)
   void output_fct(std::ostream &outfile, const unsigned &n_plot,
                   const double& time,
                   FiniteElement::UnsteadyExactSolutionFctPt
                   exact_soln_pt)
    {
     EulerEquations<NFLUX,DIM>::
      output_fct(outfile,n_plot,time,exact_soln_pt);
      }*/

  protected:
    /// Shape, test functions & derivs. w.r.t. to global coords. Return
    /// Jacobian.
    inline double dshape_and_dtest_eulerian_flux_transport(
      const Vector<double>& s,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;

    /// \short Shape, test functions & derivs. w.r.t. to global coords. at
    /// integration point ipt. Return Jacobian.
    inline double dshape_and_dtest_eulerian_at_knot_flux_transport(
      const unsigned& ipt,
      Shape& psi,
      DShape& dpsidx,
      Shape& test,
      DShape& dtestdx) const;
  };

  // Inline functions:

  //======================================================================
  /// \short Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralEulerElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_flux_transport(const Vector<double>& s,
                                             Shape& psi,
                                             DShape& dpsidx,
                                             Shape& test,
                                             DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian(s, psi, dpsidx);

    // Loop over the test functions and derivatives and set them equal to the
    // shape functions
    for (unsigned i = 0; i < NNODE_1D; i++)
    {
      test[i] = psi[i];
      for (unsigned j = 0; j < DIM; j++)
      {
        dtestdx(i, j) = dpsidx(i, j);
      }
    }

    // Return the jacobian
    return J;
  }

  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  double QSpectralEulerElement<DIM, NNODE_1D>::
    dshape_and_dtest_eulerian_at_knot_flux_transport(const unsigned& ipt,
                                                     Shape& psi,
                                                     DShape& dpsidx,
                                                     Shape& test,
                                                     DShape& dtestdx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = this->dshape_eulerian_at_knot(ipt, psi, dpsidx);

    // Set the test functions equal to the shape functions (pointer copy)
    test = psi;
    dtestdx = dpsidx;

    // Return the jacobian
    return J;
  }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// \short Face geometry for the QEulerElement elements:
  /// The spatial dimension of the face elements is one lower than that
  /// of the bulk element but they have the same number of points along
  /// their 1D edges.
  //=======================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class FaceGeometry<QSpectralEulerElement<DIM, NNODE_1D>> :
    public virtual QSpectralElement<DIM - 1, NNODE_1D>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional QElement
    FaceGeometry() : QElement<DIM - 1, NNODE_1D>() {}
  };

  //======================================================================
  /// FaceElement for Discontinuous Galerkin Problems
  //======================================================================
  template<class ELEMENT>
  class DGEulerFaceElement :
    public virtual FaceGeometry<ELEMENT>,
    public virtual DGFaceElement
  {
    unsigned Nflux;

    bool Reflecting;

  public:
    /// Constructor
    DGEulerFaceElement(FiniteElement* const& element_pt,
                       const int& face_index) :
      FaceGeometry<ELEMENT>(), DGFaceElement()
    {
      // Attach geometric information to the element
      // N.B. This function also assigns nbulk_value from required_nvalue
      // of the bulk element.
      element_pt->build_face_element(face_index, this);
      // Set the value of the flux
      Nflux = 2 + element_pt->dim();
      // Use the Face integration scheme of the element for 2D elements
      if (Nflux > 3)
      {
        this->set_integration_scheme(
          dynamic_cast<ELEMENT*>(element_pt)->face_integration_pt());
      }
    }

    /// Specify the value of nodal zeta from the face geometry
    /// \short The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    // There is a single required n_flux
    unsigned required_nflux()
    {
      return Nflux;
    }

    /// Calculate the normal flux, so this is the dot product of the
    /// numerical flux with n_out
    void numerical_flux(const Vector<double>& n_out,
                        const Vector<double>& u_int,
                        const Vector<double>& u_ext,
                        Vector<double>& flux)
    {
      // Let's follow the yellow book and use local Lax-Friedrichs
      // This is almost certainly not the best flux to use ---
      // further investigation is required here.
      // Cache the bulk element
      ELEMENT* cast_bulk_element_pt =
        dynamic_cast<ELEMENT*>(this->bulk_element_pt());

      // Find the dimension of the problem
      const unsigned dim = cast_bulk_element_pt->dim();
      // Find the number of variables
      const unsigned n_flux = this->Nflux;

      const double gamma = cast_bulk_element_pt->gamma();

      // Now let's find the fluxes
      DenseMatrix<double> flux_int(n_flux, dim);
      DenseMatrix<double> flux_ext(n_flux, dim);

      cast_bulk_element_pt->flux(u_int, flux_int);
      cast_bulk_element_pt->flux(u_ext, flux_ext);

      DenseMatrix<double> flux_av(n_flux, dim);
      for (unsigned i = 0; i < n_flux; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          flux_av(i, j) = 0.5 * (flux_int(i, j) + flux_ext(i, j));
        }
      }

      flux.initialise(0.0);
      // Now set the first part of the numerical flux
      for (unsigned i = 0; i < n_flux; i++)
      {
        for (unsigned j = 0; j < dim; j++)
        {
          flux[i] += flux_av(i, j) * n_out[j];
        }
      }

      // Now let's find the normal jumps in the fluxes
      Vector<double> jump(n_flux);
      // Now we take the normal jump in the quantities
      // The first two are scalars
      for (unsigned i = 0; i < 2; i++)
      {
        jump[i] = u_int[i] - u_ext[i];
      }

      // The next are vectors so treat them accordingly ??
      // By taking the component dotted with the normal component
      double velocity_jump = 0.0;
      for (unsigned j = 0; j < dim; j++)
      {
        velocity_jump += (u_int[2 + j] - u_ext[2 + j]) * n_out[j];
      }
      // Now we multiply by the outer normal to get the vector flux
      for (unsigned j = 0; j < dim; j++)
      {
        jump[2 + j] = velocity_jump * n_out[j];
      }

      // Let's find the Roe average
      /*   Vector<double> u_average(n_flux);
      double sum = sqrt(u_int[0]) + sqrt(u_ext[0]);
      for(unsigned i=0;i<n_flux;i++)
       {
        u_average[i] = (sqrt(u_int[0])*u_int[i] + sqrt(u_ext[0])*u_ext[i])/sum;
       }

      //Find the average pressure
      double p = cast_bulk_element_pt->pressure(u_average);

      //Now do the normal velocity
      double vel = 0.0;
      for(unsigned j=0;j<dim;j++)
       {
        //vel += u_average[2+j]*u_average[2+j];
        vel += u_average[2+j]*n_out[j];
       }
      //vel = sqrt(vel)/u_average[0];
      vel /= u_average[0];

      double arg = std::fabs(gamma*p/u_average[0]);*/

      // Let's calculate the internal and external pressures and then enthalpies
      double p_int = cast_bulk_element_pt->pressure(u_int);
      double p_ext = cast_bulk_element_pt->pressure(u_ext);

      // Limit the pressures to zero if necessary, but keep the energy the
      // same
      if (p_int < 0)
      {
        oomph_info << "Negative int pressure" << p_int << "\n";
        p_int = 0.0;
      }
      if (p_ext < 0)
      {
        oomph_info << "Negative ext pressure " << p_ext << "\n";
        p_ext = 0.0;
      }

      double H_int = (u_int[1] + p_int) / u_int[0];
      double H_ext = (u_ext[1] + p_ext) / u_ext[0];

      // Also calculate the velocities
      Vector<double> vel_int(2), vel_ext(2);
      for (unsigned j = 0; j < dim; j++)
      {
        vel_int[j] = u_int[2 + j] / u_int[0];
        vel_ext[j] = u_ext[2 + j] / u_ext[0];
      }

      // Now we calculate the Roe averaged values
      Vector<double> vel_average(dim);
      double s_int = sqrt(u_int[0]);
      double s_ext = sqrt(u_ext[0]);
      double sum = s_int + s_ext;

      // Velocities
      for (unsigned j = 0; j < dim; j++)
      {
        vel_average[j] = (s_int * vel_int[j] + s_ext * vel_ext[j]) / sum;
      }

      // The enthalpy
      double H_average = (s_int * H_int + s_ext * H_ext) / sum;

      // Now calculate the square of the sound speed
      double arg = H_average;
      for (unsigned j = 0; j < dim; j++)
      {
        arg -= 0.5 * vel_average[j];
      }
      arg *= (gamma - 1.0);
      // Get the local sound speed
      if (arg < 0.0)
      {
        oomph_info << "Square of sound speed is negative!\n";
        arg = 0.0;
      }
      double a = sqrt(arg);

      // Calculate the normal average velocity
      // Now do the normal velocity
      double vel = 0.0;
      for (unsigned j = 0; j < dim; j++)
      {
        vel += vel_average[j] * n_out[j];
      }

      Vector<double> eigA(3, 0.0);

      eigA[0] = vel - a;
      eigA[1] = vel;
      eigA[2] = vel + a;

      double lambda = std::fabs(eigA[0]);
      for (unsigned i = 1; i < 3; i++)
      {
        if (std::fabs(eigA[i]) > lambda)
        {
          lambda = std::fabs(eigA[i]);
        }
      }

      // Find the magnitude of the external velocity
      /*double vel_mag = 0.0;
      for(unsigned j=0;j<dim;j++) {vel_mag += u_ext[2+j]*u_ext[2+j];}
      vel_mag = sqrt(vel_mag)/u_ext[0];

      //Get the pressure
      double p = cast_bulk_element_pt->pressure(u_ext);
      double lambda_ext = vel_mag + sqrt(std::fabs(gamma*p/u_ext[0]));

      //Let's do the same for the internal one
      vel_mag = 0.0;
      for(unsigned j=0;j<dim;j++) {vel_mag += u_int[2+j]*u_int[2+j];}
      vel_mag = sqrt(vel_mag)/u_int[0];

      //Get the pressure
      p = cast_bulk_element_pt->pressure(u_int);
      double lambda_int = vel_mag + sqrt(std::fabs(gamma*p/u_int[0]));

      //Now take the largest
      double lambda = lambda_int;
      if(lambda_ext > lambda) {lambda = lambda_ext;}*/

      for (unsigned i = 0; i < n_flux; i++)
      {
        flux[i] += 0.5 * lambda * jump[i];
      }
    }
  };

  //======================================================================
  /// \short FaceElement for Discontinuous Galerkin Problems with reflection
  /// boundary conditions
  //======================================================================
  template<class ELEMENT>
  class DGEulerFaceReflectionElement :
    public virtual FaceGeometry<ELEMENT>,
    public virtual DGFaceElement
  {
    unsigned Nflux;

  public:
    /// Constructor
    DGEulerFaceReflectionElement(FiniteElement* const& element_pt,
                                 const int& face_index) :
      FaceGeometry<ELEMENT>(), DGFaceElement()
    {
      // Attach geometric information to the element
      // N.B. This function also assigns nbulk_value from required_nvalue
      // of the bulk element.
      element_pt->build_face_element(face_index, this);
      // Set the value of the flux
      Nflux = 2 + element_pt->dim();
    }

    /// Specify the value of nodal zeta from the face geometry
    /// \short The "global" intrinsic coordinate of the element when
    /// viewed as part of a geometric object should be given by
    /// the FaceElement representation, by default (needed to break
    /// indeterminacy if bulk element is SolidElement)
    double zeta_nodal(const unsigned& n,
                      const unsigned& k,
                      const unsigned& i) const
    {
      return FaceElement::zeta_nodal(n, k, i);
    }

    // There is a single required n_flux
    unsigned required_nflux()
    {
      return Nflux;
    }

    /// We overload interpolated_u to reflect
    void interpolated_u(const Vector<double>& s, Vector<double>& u)
    {
      // Get the standard interpolated_u
      DGFaceElement::interpolated_u(s, u);

      // Now do the reflection condition for the velocities
      // Find dot product of normal and velocities
      const unsigned nodal_dim = this->nodal_dimension();
      Vector<double> n(nodal_dim);
      this->outer_unit_normal(s, n);

      double dot = 0.0;
      for (unsigned j = 0; j < nodal_dim; j++)
      {
        dot += n[j] * u[2 + j];
      }

      // Now subtract
      for (unsigned j = 0; j < nodal_dim; j++)
      {
        u[2 + j] -= 2.0 * dot * n[j];
      }
    }
  };

  //=================================================================
  /// General DGEulerClass. Establish the template parameters
  //===================================================================
  template<unsigned DIM, unsigned NNODE_1D>
  class DGSpectralEulerElement
  {
  };

  //==================================================================
  // Specialization for 1D DG Advection element
  //==================================================================
  template<unsigned NNODE_1D>
  class DGSpectralEulerElement<1, NNODE_1D> :
    public QSpectralEulerElement<1, NNODE_1D>,
    public DGElement
  {
    friend class DGEulerFaceElement<DGSpectralEulerElement<1, NNODE_1D>>;

    Vector<double> Inverse_mass_diagonal;

  public:
    /// Overload the required number of fluxes for the DGElement
    unsigned required_nflux()
    {
      return this->nflux();
    }

    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<1>::calculate_element_averages(average_value);
    }

    // Constructor
    DGSpectralEulerElement() : QSpectralEulerElement<1, NNODE_1D>(), DGElement()
    {
      // Need to up the order of integration for the accurate resolution
      // of the quadratic non-linearities
      this->set_integration_scheme(new Gauss<1, NNODE_1D>);
    }

    // Dummy
    Integral* face_integration_pt() const
    {
      return 0;
    }

    ~DGSpectralEulerElement() {}

    void build_all_faces()
    {
      // Make the two faces
      Face_element_pt.resize(2);
      // Make the face on the left
      Face_element_pt[0] =
        new DGEulerFaceElement<DGSpectralEulerElement<1, NNODE_1D>>(this, -1);
      // Make the face on the right
      Face_element_pt[1] =
        new DGEulerFaceElement<DGSpectralEulerElement<1, NNODE_1D>>(this, +1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QSpectralEulerElement<1, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }

    //============================================================================
    /// Function that returns the current value of the residuals
    /// multiplied by the inverse mass matrix (virtual so that it can be
    /// overloaded specific elements in which time and memory saving tricks can
    /// be applied)
    //============================================================================
    /* void get_inverse_mass_matrix_times_residuals(Vector<double> &minv_res)
    {
     //Now let's assemble stuff
     const unsigned n_dof = this->ndof();

     //Resize and initialise the vector that will holds the residuals
     minv_res.resize(n_dof);
     for(unsigned n=0;n<n_dof;n++) {minv_res[n] = 0.0;}

     //If we are recycling the mass matrix
     if(Mass_matrix_reuse_is_enabled && Mass_matrix_has_been_computed)
      {
       //Get the residuals
       this->fill_in_contribution_to_residuals(minv_res);
      }
     //Otherwise
     else
     {
      //Temporary mass matrix
      DenseDoubleMatrix M(n_dof,n_dof,0.0);

      //Get the local mass matrix and residuals
      this->fill_in_contribution_to_mass_matrix(minv_res,M);

      //Store the diagonal entries
      Inverse_mass_diagonal.clear();
      for(unsigned n=0;n<n_dof;n++)
    {Inverse_mass_diagonal.push_back(1.0/M(n,n));}

      //The mass matrix has been computed
      Mass_matrix_has_been_computed=true;
     }

     for(unsigned n=0;n<n_dof;n++) {minv_res[n] *= Inverse_mass_diagonal[n];}

    }*/
  };

  //=======================================================================
  /// Face geometry of the 1D  DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGSpectralEulerElement<1, NNODE_1D>> :
    public virtual PointElement
  {
  public:
    FaceGeometry() : PointElement() {}
  };

  //==================================================================
  /// Specialisation for 2D DG Elements
  //==================================================================
  template<unsigned NNODE_1D>
  class DGSpectralEulerElement<2, NNODE_1D> :
    public QSpectralEulerElement<2, NNODE_1D>,
    public DGElement
  {
    friend class DGEulerFaceElement<DGSpectralEulerElement<2, NNODE_1D>>;

    static Gauss<1, NNODE_1D> Default_face_integration_scheme;

  public:
    /// Overload the required number of fluxes for the DGElement
    unsigned required_nflux()
    {
      return this->nflux();
    }

    // Calculate averages
    void calculate_element_averages(double*& average_value)
    {
      FluxTransportEquations<2>::calculate_element_averages(average_value);
    }

    // Constructor
    DGSpectralEulerElement() : QSpectralEulerElement<2, NNODE_1D>(), DGElement()
    {
      // Need to up the order of integration for the accurate resolution
      // of the quadratic non-linearities
      this->set_integration_scheme(
        new GaussLobattoLegendre<2, 3 * NNODE_1D / 2>);
    }

    ~DGSpectralEulerElement() {}

    Integral* face_integration_pt() const
    {
      return &Default_face_integration_scheme;
    }

    void build_all_faces()
    {
      Face_element_pt.resize(4);
      Face_element_pt[0] =
        new DGEulerFaceElement<DGSpectralEulerElement<2, NNODE_1D>>(this, 2);
      Face_element_pt[1] =
        new DGEulerFaceElement<DGSpectralEulerElement<2, NNODE_1D>>(this, 1);
      Face_element_pt[2] =
        new DGEulerFaceElement<DGSpectralEulerElement<2, NNODE_1D>>(this, -2);
      Face_element_pt[3] =
        new DGEulerFaceElement<DGSpectralEulerElement<2, NNODE_1D>>(this, -1);
    }

    ///\short Compute the residuals for the Navier--Stokes equations;
    /// flag=1(or 0): do (or don't) compute the Jacobian as well.
    void fill_in_generic_residual_contribution_flux_transport(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix,
      unsigned flag)
    {
      QSpectralEulerElement<2, NNODE_1D>::
        fill_in_generic_residual_contribution_flux_transport(
          residuals, jacobian, mass_matrix, flag);

      this->add_flux_contributions_to_residuals(residuals, jacobian, flag);
    }
  };

  //=======================================================================
  /// Face geometry of the DG elements
  //=======================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<DGSpectralEulerElement<2, NNODE_1D>> :
    public virtual QSpectralElement<1, NNODE_1D>
  {
  public:
    FaceGeometry() : QSpectralElement<1, NNODE_1D>() {}
  };

} // namespace oomph

#endif
