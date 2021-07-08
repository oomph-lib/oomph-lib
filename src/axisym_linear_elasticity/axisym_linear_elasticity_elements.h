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

// Include guards to prevent multiple inclusion of the header
#ifndef OOMPH_AXISYMMETRIC_LINEAR_ELASTICITY_ELEMENTS_HEADER
#define OOMPH_AXISYMMETRIC_LINEAR_ELASTICITY_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

#ifdef OOMPH_HAS_MPI
#include "mpi.h"
#endif

// OOMPH-LIB headers
#include "generic/Qelements.h"
#include "generic/Telements.h"
#include "generic/projection.h"

namespace oomph
{
  //=======================================================================
  /// A base class for elements that solve the axisymmetric (in
  /// cylindrical polars) equations of linear elasticity.
  //=======================================================================
  class AxisymmetricLinearElasticityEquationsBase : public virtual FiniteElement
  {
  public:
    /// \short Return the index at which the i-th (i=0: r, i=1: z; i=2: theta)
    /// unknown displacement component is stored at the nodes.  The default
    /// assignment here (u_r, u_z, u_theta) is appropriate for single-physics
    /// problems.
    virtual inline unsigned u_index_axisymmetric_linear_elasticity(
      const unsigned& i) const
    {
      return i;
    }

    /// d^2u/dt^2 at local node n
    double d2u_dt2_axisymmetric_linear_elasticity(const unsigned& n,
                                                  const unsigned& i) const
    {
      // Get the timestepper
      TimeStepper* time_stepper_pt = node_pt(n)->time_stepper_pt();

      // Storage for the derivative - initialise to 0
      double d2u_dt2 = 0.0;

      // If we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // Get the nodal index
        const unsigned u_nodal_index =
          u_index_axisymmetric_linear_elasticity(i);

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

    /// du/dt at local node n
    double du_dt_axisymmetric_linear_elasticity(const unsigned& n,
                                                const unsigned& i) const
    {
      // Get the timestepper
      TimeStepper* time_stepper_pt = node_pt(n)->time_stepper_pt();

      // Storage for the derivative - initialise to 0
      double du_dt = 0.0;

      // If we are doing an unsteady solve then calculate the derivative
      if (!time_stepper_pt->is_steady())
      {
        // Get the nodal index
        const unsigned u_nodal_index =
          u_index_axisymmetric_linear_elasticity(i);

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

    /// Compute vector of FE interpolated displacement u at local coordinate s
    void interpolated_u_axisymmetric_linear_elasticity(
      const Vector<double>& s, Vector<double>& disp) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      for (unsigned i = 0; i < 3; i++)
      {
        // Index at which the nodal value is stored
        unsigned u_nodal_index = u_index_axisymmetric_linear_elasticity(i);

        // Initialise value of u
        disp[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          const double u_value = nodal_value(l, u_nodal_index);

          disp[i] += u_value * psi[l];
        }
      }
    }

    /// \short Return FE interpolated displacement u[i] (i=0: r, i=1: z; i=2:
    /// theta) at local coordinate s
    double interpolated_u_axisymmetric_linear_elasticity(
      const Vector<double>& s, const unsigned& i) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Get nodal index at which i-th velocity is stored
      unsigned u_nodal_index = u_index_axisymmetric_linear_elasticity(i);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Loop over the local nodes and sum
      for (unsigned l = 0; l < n_node; l++)
      {
        const double u_value = nodal_value(l, u_nodal_index);

        interpolated_u += u_value * psi[l];
      }

      return (interpolated_u);
    }

    /// Compute vector of FE interpolated velocity du/dt at local coordinate s
    void interpolated_du_dt_axisymmetric_linear_elasticity(
      const Vector<double>& s, Vector<double>& du_dt) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Loop over directions
      for (unsigned i = 0; i < 3; i++)
      {
        // Initialise value of u
        du_dt[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          du_dt[i] += du_dt_axisymmetric_linear_elasticity(l, i) * psi[l];
        }
      }
    }

    /// Compute vector of FE interpolated accel d2u/dt2 at local coordinate s
    void interpolated_d2u_dt2_axisymmetric_linear_elasticity(
      const Vector<double>& s, Vector<double>& d2u_dt2) const
    {
      // Find number of nodes
      unsigned n_node = nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      shape(s, psi);

      // Loop over directions
      for (unsigned i = 0; i < 3; i++)
      {
        // Initialise value of u
        d2u_dt2[i] = 0.0;

        // Loop over the local nodes and sum
        for (unsigned l = 0; l < n_node; l++)
        {
          d2u_dt2[i] += d2u_dt2_axisymmetric_linear_elasticity(l, i) * psi[l];
        }
      }
    }

    /// \short Function pointer to function that specifies the body force
    /// as a function of the Cartesian coordinates and time FCT(x,b) --
    /// x and b are  Vectors!
    typedef void (*BodyForceFctPt)(const double& time,
                                   const Vector<double>& x,
                                   Vector<double>& b);

    /// \short Constructor: Set null pointers for constitutive law.
    /// Set physical parameter values to
    /// default values, and set body force to zero.
    AxisymmetricLinearElasticityEquationsBase() :
      Youngs_modulus_pt(&Default_youngs_modulus_value),
      Nu_pt(0),
      Lambda_sq_pt(&Default_lambda_sq_value),
      Body_force_fct_pt(0)
    {
    }

    /// Return the pointer to Young's modulus
    double*& youngs_modulus_pt()
    {
      return Youngs_modulus_pt;
    }

    /// Access function to Young's modulus
    inline double youngs_modulus() const
    {
      return (*Youngs_modulus_pt);
    }

    /// Access function for Poisson's ratio
    double& nu() const
    {
#ifdef PARANOID
      if (Nu_pt == 0)
      {
        std::ostringstream error_message;
        error_message << "No pointer to Poisson's ratio set. Please set one!\n";
        throw OomphLibError(error_message.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return *Nu_pt;
    }

    /// Access function for pointer to Poisson's ratio
    double*& nu_pt()
    {
      return Nu_pt;
    }

    /// Access function for pointer to timescale ratio (nondim density)
    double*& lambda_sq_pt()
    {
      return Lambda_sq_pt;
    }

    /// Access function for timescale ratio (nondim density)
    const double& lambda_sq() const
    {
      return *Lambda_sq_pt;
    }

    /// Access function: Pointer to body force function
    BodyForceFctPt& body_force_fct_pt()
    {
      return Body_force_fct_pt;
    }

    /// Access function: Pointer to body force function (const version)
    BodyForceFctPt body_force_fct_pt() const
    {
      return Body_force_fct_pt;
    }

    /// \short Evaluate body force at Eulerian coordinate x at present time
    /// (returns zero vector if no body force function pointer has been set)
    inline void body_force(const double& time,
                           const Vector<double>& x,
                           Vector<double>& b) const
    {
      // If no function has been set, return zero vector
      if (Body_force_fct_pt == 0)
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
        (*Body_force_fct_pt)(time, x, b);
      }
    }

    /// \short The number of "DOF types" that degrees of freedom in this element
    /// are sub-divided into: for now lump them all into one DOF type.
    /// Can be adjusted later
    unsigned ndof_types() const
    {
      return 1;
    }

    /// \short Create a list of pairs for all unknowns in this element,
    /// so that the first entry in each pair contains the global equation
    /// number of the unknown, while the second one contains the number
    /// of the "DOF type" that this unknown is associated with.
    /// (Function can obviously only be called if the equation numbering
    /// scheme has been set up.)
    void get_dof_numbers_for_unknowns(
      std::list<std::pair<unsigned long, unsigned>>& dof_lookup_list) const
    {
      // temporary pair (used to store DOF lookup prior to being added
      // to list)
      std::pair<unsigned long, unsigned> dof_lookup;

      // number of nodes
      const unsigned n_node = this->nnode();

      // Integer storage for local unknown
      int local_unknown = 0;

      // Loop over the nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Loop over dimension
        for (unsigned i = 0; i < 3; i++)
        {
          // If the variable is free
          local_unknown = nodal_local_eqn(n, i);

          // ignore pinned values
          if (local_unknown >= 0)
          {
            // store DOF type lookup in temporary pair: First entry in pair
            // is global equation number; second entry is DOF type
            dof_lookup.first = this->eqn_number(local_unknown);
            dof_lookup.second = 0;

            // add to list
            dof_lookup_list.push_front(dof_lookup);
          }
        }
      }
    }

  protected:
    /// Pointer to the Young's modulus
    double* Youngs_modulus_pt;

    /// Pointer to Poisson's ratio
    double* Nu_pt;

    /// Timescale ratio (non-dim. density)
    double* Lambda_sq_pt;

    /// Pointer to body force function
    BodyForceFctPt Body_force_fct_pt;

    /// \short Static default value for Young's modulus (1.0 -- for natural
    /// scaling, i.e. all stresses have been non-dimensionalised by
    /// the same reference Young's modulus. Setting the "non-dimensional"
    /// Young's modulus (obtained by de-referencing Youngs_modulus_pt)
    /// to a number larger than one means that the material is stiffer
    /// than assumed in the non-dimensionalisation.
    static double Default_youngs_modulus_value;

    /// Static default value for timescale ratio (1.0 for natural scaling)
    static double Default_lambda_sq_value;
  };

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  //=======================================================================
  /// A class for elements that solve the axisymmetric (in cylindrical
  /// polars) equations of linear elasticity
  //=======================================================================
  class AxisymmetricLinearElasticityEquations :
    public AxisymmetricLinearElasticityEquationsBase
  {
  public:
    /// \short  Constructor
    AxisymmetricLinearElasticityEquations() {}

    /// Number of values required at node n.
    unsigned required_nvalue(const unsigned& n) const
    {
      return 3;
    }

    /// \short Return the residuals for the equations (the discretised
    /// principle of virtual displacements)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_generic_contribution_to_residuals_axisymmetric_linear_elasticity(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// The jacobian is calculated by finite differences by default,
    /// We need only to take finite differences w.r.t. positional variables
    /// For this element
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Add the contribution to the residuals
      this
        ->fill_in_generic_contribution_to_residuals_axisymmetric_linear_elasticity(
          residuals, jacobian, 1);
    }

    /// Get strain (3x3 entries; r, z, phi)
    void get_strain(const Vector<double>& s, DenseMatrix<double>& strain);

    /// Output exact solution: r,z, u_r, u_z, u_theta
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// Output exact solution: r,z, u_r, u_z, u_theta
    /// Time dependent version
    void output_fct(std::ostream& outfile,
                    const unsigned& nplot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt);

    /// Output: r,z, u_r, u_z, u_theta
    void output(std::ostream& outfile)
    {
      unsigned n_plot = 5;
      output(outfile, n_plot);
    }

    /// Output: r,z, u_r, u_z, u_theta
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// C-style output: r,z, u_r, u_z, u_theta
    void output(FILE* file_pt)
    {
      unsigned n_plot = 5;
      output(file_pt, n_plot);
    }

    /// Output:  r,z, u_r, u_z, u_theta
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Validate against exact solution.
    /// Solution is provided via function pointer.
    /// Plot at a given number of plot points and compute L2 error
    /// and L2 norm of displacement solution over element
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm);

    /// Validate against exact solution.
    /// Time-dependent version
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm);

  protected:
    /// \short Private helper function to compute residuals and (if requested
    /// via flag) also the Jacobian matrix.
    virtual void fill_in_generic_contribution_to_residuals_axisymmetric_linear_elasticity(
      Vector<double>& residuals, DenseMatrix<double>& jacobian, unsigned flag);
  };

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  //===========================================================================
  /// An Element that solves the equations of axisymmetric (in cylindrical
  /// polars) linear elasticity, using QElements for the geometry.
  //============================================================================
  template<unsigned NNODE_1D>
  class QAxisymmetricLinearElasticityElement :
    public virtual QElement<2, NNODE_1D>,
    public virtual AxisymmetricLinearElasticityEquations
  {
  public:
    /// Constructor
    QAxisymmetricLinearElasticityElement() :
      QElement<2, NNODE_1D>(), AxisymmetricLinearElasticityEquations()
    {
    }

    /// Output function
    void output(std::ostream& outfile)
    {
      AxisymmetricLinearElasticityEquations::output(outfile);
    }

    /// Output function
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      AxisymmetricLinearElasticityEquations::output(outfile, n_plot);
    }

    /// C-style output function
    void output(FILE* file_pt)
    {
      AxisymmetricLinearElasticityEquations::output(file_pt);
    }

    /// C-style output function
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      AxisymmetricLinearElasticityEquations::output(file_pt, n_plot);
    }
  };

  //============================================================================
  /// FaceGeometry of a linear
  /// QAxisymmetricLinearElasticityElement element
  //============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<QAxisymmetricLinearElasticityElement<NNODE_1D>> :
    public virtual QElement<1, NNODE_1D>
  {
  public:
    /// Constructor must call the constructor of the underlying element
    FaceGeometry() : QElement<1, NNODE_1D>() {}
  };

  /* //////////////////////////////////////////////////////////////////////// */
  /* //////////////////////////////////////////////////////////////////////// */
  /* //////////////////////////////////////////////////////////////////////// */

  /* //===========================================================================
   */
  /* /// An Element that solves the equations of axisymmetric (in cylindrical */
  /* /// polars) linear elasticity, using TElements for the geometry. */
  /* //============================================================================
   */
  /*   template<unsigned NNODE_1D> */
  /*    class TAxisymmetricLinearElasticityElement :  */
  /*   public virtual TElement<2,NNODE_1D>, */
  /*    public virtual AxisymmetricLinearElasticityEquations */
  /*    { */
  /*     public: */

  /*    /// Constructor */
  /*    TAxisymmetricLinearElasticityElement() :  */
  /*     TElement<2,NNODE_1D>(),  */
  /*     AxisymmetricLinearElasticityEquations() { } */

  /*     /// Output function */
  /*    void output(std::ostream &outfile)  */
  /*    {AxisymmetricLinearElasticityEquations::output(outfile);} */

  /*    /// Output function */
  /*    void output(std::ostream &outfile, const unsigned &n_plot) */
  /*    {AxisymmetricLinearElasticityEquations:: */
  /*      output(outfile,n_plot);} */

  /*    /// C-style output function */
  /*    void output(FILE* file_pt)  */
  /*    {AxisymmetricLinearElasticityEquations::output(file_pt);} */

  /*    /// C-style output function */
  /*    void output(FILE* file_pt, const unsigned &n_plot) */
  /*    {AxisymmetricLinearElasticityEquations:: */
  /*      output(file_pt,n_plot);} */

  /*   }; */

  /* //============================================================================
   */
  /* /// FaceGeometry of a linear  */
  /* /// TAxisymmetricLinearElasticityElement element */
  /* //============================================================================
   */
  /*  template<unsigned NNODE_1D> */
  /*   class FaceGeometry<TAxisymmetricLinearElasticityElement<NNODE_1D> > : */
  /*  public virtual TElement<1,NNODE_1D> */
  /*   { */
  /*     public: */
  /*    /// Constructor must call the constructor of the underlying element */
  /*     FaceGeometry() : TElement<1,NNODE_1D>() {} */
  /*   }; */

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  //==========================================================
  /// Axisym linear elasticity upgraded to become projectable
  //==========================================================
  template<class AXISYM_LINEAR_ELAST_ELEMENT>
  class ProjectableAxisymLinearElasticityElement :
    public virtual ProjectableElement<AXISYM_LINEAR_ELAST_ELEMENT>
  {
  public:
    /// \short Constructor [this was only required explicitly
    /// from gcc 4.5.2 onwards...]
    ProjectableAxisymLinearElasticityElement() {}

    /// \short Specify the values associated with field fld.
    /// The information is returned in a vector of pairs which comprise
    /// the Data object and the value within it, that correspond to field fld.
    /// In the underlying linear elasticity elements the
    /// the displacements are stored at the nodal values
    Vector<std::pair<Data*, unsigned>> data_values_of_field(const unsigned& fld)
    {
      // Create the vector
      Vector<std::pair<Data*, unsigned>> data_values;

      // Loop over all nodes and extract the fld-th nodal value
      unsigned nnod = this->nnode();
      for (unsigned j = 0; j < nnod; j++)
      {
        // Add the data value associated with the displacement components
        data_values.push_back(std::make_pair(this->node_pt(j), fld));
      }

      // Return the vector
      return data_values;
    }

    /// \short Number of fields to be projected: 3, corresponding to
    /// the displacement components
    unsigned nfields_for_projection()
    {
      return 3;
    }

    /// \short Number of history values to be stored for fld-th field.
    /// (includes present value!)
    unsigned nhistory_values_for_projection(const unsigned& fld)
    {
#ifdef PARANOID
      if (fld > 2)
      {
        std::stringstream error_stream;
        error_stream << "Elements only store two fields so fld can't be"
                     << " " << fld << std::endl;
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return this->node_pt(0)->ntstorage();
    }

    ///\short Number of positional history values: Read out from
    /// positional timestepper  (Note: count includes current value!)
    unsigned nhistory_values_for_coordinate_projection()
    {
      return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
    }

    /// \short Return Jacobian of mapping and shape functions of field fld
    /// at local coordinate s
    double jacobian_and_shape_of_field(const unsigned& fld,
                                       const Vector<double>& s,
                                       Shape& psi)
    {
      unsigned n_dim = this->dim();
      unsigned n_node = this->nnode();
      DShape dpsidx(n_node, n_dim);

      // Call the derivatives of the shape functions and return
      // the Jacobian
      return this->dshape_eulerian(s, psi, dpsidx);
    }

    /// \short Return interpolated field fld at local coordinate s, at time
    /// level t (t=0: present; t>0: history values)
    double get_field(const unsigned& t,
                     const unsigned& fld,
                     const Vector<double>& s)
    {
      unsigned n_node = this->nnode();

      // Local shape function
      Shape psi(n_node);

      // Find values of shape function
      this->shape(s, psi);

      // Initialise value of u
      double interpolated_u = 0.0;

      // Sum over the local nodes at that time
      for (unsigned l = 0; l < n_node; l++)
      {
        interpolated_u += this->nodal_value(t, l, fld) * psi[l];
      }
      return interpolated_u;
    }

    /// Return number of values in field fld
    unsigned nvalue_of_field(const unsigned& fld)
    {
      return this->nnode();
    }

    /// Return local equation number of value j in field fld.
    int local_equation(const unsigned& fld, const unsigned& j)
    {
      return this->nodal_local_eqn(j, fld);
    }
  };

  //=======================================================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<ProjectableAxisymLinearElasticityElement<ELEMENT>> :
    public virtual FaceGeometry<ELEMENT>
  {
  public:
    FaceGeometry() : FaceGeometry<ELEMENT>() {}
  };

  //=======================================================================
  /// Face geometry of the Face Geometry for element is the same as
  /// that for the underlying wrapped element
  //=======================================================================
  template<class ELEMENT>
  class FaceGeometry<
    FaceGeometry<ProjectableAxisymLinearElasticityElement<ELEMENT>>> :
    public virtual FaceGeometry<FaceGeometry<ELEMENT>>
  {
  public:
    FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT>>() {}
  };

} // namespace oomph

#endif
