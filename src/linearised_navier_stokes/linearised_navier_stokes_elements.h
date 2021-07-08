//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Header file for linearised axisymmetric Navier-Stokes elements

#ifndef OOMPH_LINEARISED_NAVIER_STOKES_ELEMENTS_HEADER
#define OOMPH_LINEARISED_NAVIER_STOKES_ELEMENTS_HEADER

// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif

// oomph-lib includes
#include "../generic/Qelements.h"
#include "../generic/fsi.h"

#include "./linearised_navier_stokes_eigenvalue_elements.h"

namespace oomph
{

#define DIM  2

 
 //=======================================================================
 /// \short A class for elements that solve the linearised version of the
 /// unsteady Navier--Stokes equations in cylindrical polar coordinates,
 /// where we have Fourier-decomposed in the azimuthal direction so that
 /// the theta-dependance is replaced by an azimuthal mode number.
 //=======================================================================
 class LinearisedNavierStokesEquations
  : public virtual FiniteElement
 {
   private:
  
  /// Static "magic" number that indicates that the pressure is not
  /// stored at a node
  static int Pressure_not_stored_at_node;
  
  /// Static default value for the physical constants
  /// (all initialised to zero)
  static double Default_Physical_Constant_Value;
  
  /// Static default value for the physical ratios (all initialised to one)
  static double Default_Physical_Ratio_Value;

   protected:
  
  // Physical constants
  // ------------------

  /// \short Pointer to the viscosity ratio (relative to the 
  /// viscosity used in the definition of the Reynolds number)
  double *Viscosity_Ratio_pt;
 
  /// \short Pointer to the density ratio (relative to the
  /// density used in the definition of the Reynolds number)
  double *Density_Ratio_pt;
 
  /// Pointer to global Reynolds number
  double *Re_pt;
 
  /// Pointer to global Reynolds number x Strouhal number (=Womersley)
  double *ReSt_pt;

  /// Pointer to eigenvalue
  double *Lambda_pt;

  /// Pointer to frequency
  double *Omega_pt;

  ///Pointer to the normalisation element
  LinearisedNavierStokesEigenfunctionNormalisationElement*
   Normalisation_element_pt;

  
  /// Index of datum where eigenvalue is stored
  unsigned Data_number_of_eigenvalue;

  unsigned Index_of_eigenvalue;
  
  /// Pointer to base flow solution (velocity components) function
  void (*Base_flow_u_fct_pt)(const double& time,
                             const Vector<double> &x, 
                             Vector<double> &result);

  /// \short Pointer to derivatives of base flow solution velocity
  /// components w.r.t. global coordinates (r and z) function
  void (*Base_flow_dudx_fct_pt)(const double& time,
                                const Vector<double> &x, 
                                DenseMatrix<double> &result);
  
  /// \short Boolean flag to indicate if ALE formulation is disabled when
  /// the time-derivatives are computed. Only set to true if you're sure
  /// that the mesh is stationary.
  bool ALE_is_disabled;

  /// \short Access function for the local equation number
  /// information for the i-th component of the pressure.
  /// p_local_eqn[n,i] = local equation number or < 0 if pinned.
  virtual int p_local_eqn(const unsigned &n, const unsigned &i)=0;

  /// \short Compute the shape functions and their derivatives 
  /// w.r.t. global coordinates at local coordinate s.
  /// Return Jacobian of mapping between local and global coordinates.
  virtual double dshape_and_dtest_eulerian_linearised_nst(
   const Vector<double> &s,
   Shape &psi, DShape &dpsidx,
   Shape &test, DShape &dtestdx) const=0;

  /// \short Compute the shape functions and their derivatives
  /// w.r.t. global coordinates at the ipt-th integration point.
  /// Return Jacobian of mapping between local and global coordinates.
  virtual double dshape_and_dtest_eulerian_at_knot_linearised_nst(
   const unsigned &ipt,
   Shape &psi, DShape &dpsidx,
   Shape &test, DShape &dtestdx) const=0;
 
  /// Compute the pressure shape functions at local coordinate s
  virtual void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi) const=0;

  /// Compute the pressure shape and test functions at local coordinate s
  virtual void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi, Shape &test) const=0;

  /// \short Calculate the velocity components of the base flow solution
  /// at a given time and Eulerian position
  virtual void get_base_flow_u(const double& time,
                               const unsigned& ipt,
                               const Vector<double>& x,
                               Vector<double>& result) const
  {
   // If the function pointer is zero return zero
   if(Base_flow_u_fct_pt==0)
    {
     // Loop over velocity components and set base flow solution to zero
     for(unsigned i=0;i<DIM;i++) { result[i] = 0.0; }
    }
   // Otherwise call the function
   else
    {
     (*Base_flow_u_fct_pt)(time,x,result);
    }
  }

  /// \short Calculate the derivatives of the velocity components of the 
  /// base flow solution w.r.t. global coordinates (r and z) at a given
  /// time and Eulerian position
  virtual void get_base_flow_dudx(const double& time,
                                  const unsigned& ipt,
                                  const Vector<double>& x,
                                  DenseMatrix<double>& result) const
  {
   // If the function pointer is zero return zero
   if(Base_flow_dudx_fct_pt==0)
    {
     // Loop over velocity components
     for(unsigned i=0;i<DIM;i++)
      {
       // Loop over coordinate directions and set to zero
       for(unsigned j=0;j<DIM;j++) { result(i,j) = 0.0; }
      }
    }
   // Otherwise call the function
   else
    {
     (*Base_flow_dudx_fct_pt)(time,x,result);
    }
  }


 inline int eigenvalue_local_eqn(const unsigned &i)
 {
  return this->external_local_eqn(this->Data_number_of_eigenvalue,
                                  this->Index_of_eigenvalue+i);
 }
 
  
  
  /// \short Compute the residuals for the Navier-Stokes equations; 
  /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
  virtual void fill_in_generic_residual_contribution_linearised_nst(
   Vector<double> &residuals,
   DenseMatrix<double> &jacobian, 
   DenseMatrix<double> &mass_matrix,
   unsigned flag);
  
   public:
  
  /// \short Constructor: NULL the base flow solution and the
  /// derivatives of the base flow function
  LinearisedNavierStokesEquations()
   : Base_flow_u_fct_pt(0), Base_flow_dudx_fct_pt(0), ALE_is_disabled(false)
   {
    // Set all the physical parameter pointers to the default value of zero
    Re_pt = &Default_Physical_Constant_Value;
    ReSt_pt = &Default_Physical_Constant_Value;

    Lambda_pt = &Default_Physical_Constant_Value;
    Omega_pt = &Default_Physical_Constant_Value;

    //Set to sensible defaults
    Data_number_of_eigenvalue=0;
    Index_of_eigenvalue=0;
    
    // Set the physical ratios to the default value of one
    Viscosity_Ratio_pt = &Default_Physical_Ratio_Value;
    Density_Ratio_pt = &Default_Physical_Ratio_Value;

    //Null out normalisation
    Normalisation_element_pt=0;
   }
   
   /// Vector to decide whether the stress-divergence form is used or not.
   //  N.B. This needs to be public so that the intel compiler gets things
   // correct. Somehow the access function messes things up when going to
   // refineable navier--stokes
   static Vector<double> Gamma;

   // Access functions for the physical constants
   // -------------------------------------------

   /// Reynolds number
   const double &re() const { return *Re_pt; }

   /// Product of Reynolds and Strouhal number (=Womersley number)
   const double &re_st() const { return *ReSt_pt; }

   const double &lambda() const {return *Lambda_pt;}

   const double &omega() const {return *Omega_pt;}
   
   /// Pointer to Reynolds number
   double* &re_pt() { return Re_pt; }
 
   /// Pointer to product of Reynolds and Strouhal number (=Womersley number)
   double* &re_st_pt() { return ReSt_pt; }

   /// Pointer to lambda
   double* &lambda_pt() { return Lambda_pt; }
 
   /// Pointer to frequency
   double* &omega_pt() { return Omega_pt; }

   ///Pointer to normalisation element
   LinearisedNavierStokesEigenfunctionNormalisationElement*
    normalisation_element_pt() {return Normalisation_element_pt;}
   
 /// the boolean flag check_nodal_data is set to false.
 void set_eigenfunction_normalisation_element(
  LinearisedNavierStokesEigenfunctionNormalisationElement* const 
  &normalisation_el_pt)
 {
  //Set the normalisation element
  Normalisation_element_pt = normalisation_el_pt;
  
  // Add eigenvalue unknown as external data to this element
  Data_number_of_eigenvalue=
   this->add_external_data(normalisation_el_pt->eigenvalue_data_pt());
  
  // Which value corresponds to the eigenvalue
  Index_of_eigenvalue=normalisation_el_pt->index_of_eigenvalue();

  //Now set the pointers to the eigenvalues
  Lambda_pt = normalisation_el_pt->eigenvalue_data_pt()
   ->value_pt(Index_of_eigenvalue);
  Omega_pt = normalisation_el_pt->eigenvalue_data_pt()
   ->value_pt(Index_of_eigenvalue+1);
 }

    
   /// \short Viscosity ratio for element: element's viscosity relative
   /// to the viscosity used in the definition of the Reynolds number
   const double &viscosity_ratio() const { return *Viscosity_Ratio_pt; }

   /// Pointer to the viscosity ratio
   double* &viscosity_ratio_pt() { return Viscosity_Ratio_pt; }

   /// \short Density ratio for element: element's density relative
   /// to the viscosity used in the definition of the Reynolds number
   const double &density_ratio() const { return *Density_Ratio_pt; }

   /// Pointer to the density ratio
   double* &density_ratio_pt() { return Density_Ratio_pt; }

   /// Access function for the base flow solution pointer
   void (* &base_flow_u_fct_pt())(const double& time,
                                  const Vector<double>& x, 
                                  Vector<double>& f) 
    {
     return Base_flow_u_fct_pt;
    }

   /// \short Access function for the derivatives of the base flow
   /// w.r.t. global coordinates solution pointer
   void (* &base_flow_dudx_fct_pt())(const double& time,
                                     const Vector<double>& x, 
                                     DenseMatrix<double>& f) 
    {
     return Base_flow_dudx_fct_pt;
    }

   /// \short Return the number of pressure degrees of freedom
   /// associated with a single pressure component in the element
   virtual unsigned npres_linearised_nst() const=0;
   
   /// \short Return the index at which the i-th unknown velocity
   /// component is stored. The default value, i, is appropriate for
   /// single-physics problems. In derived multi-physics elements, this
   /// function should be overloaded to reflect the chosen storage scheme.
   /// Note that these equations require that the unknowns are always
   /// stored at the same indices at each node.
   virtual inline unsigned u_index_linearised_nst(const unsigned &i)
    const { return i; }
 
   /// \short Return the i-th component of du/dt at local node n. 
   /// Uses suitably interpolated value for hanging nodes.
   double du_dt_linearised_nst(const unsigned &n, const unsigned &i) const
   {
    // Get the data's timestepper
    TimeStepper* time_stepper_pt = this->node_pt(n)->time_stepper_pt();

    // Initialise dudt
    double dudt = 0.0;

    // Loop over the timesteps, if there is a non-steady timestepper
    if (!time_stepper_pt->is_steady())
     {
      // Get the index at which the velocity is stored
      const unsigned u_nodal_index = u_index_linearised_nst(i);
      
      // Determine number of timsteps (past & present)
      const unsigned n_time = time_stepper_pt->ntstorage();
      
      // Add the contributions to the time derivative
      for(unsigned t=0;t<n_time;t++)
       {
        dudt += time_stepper_pt->weight(1,t)*nodal_value(t,n,u_nodal_index);
       }
     }
   
    return dudt;
   }
 
   /// \short Disable ALE, i.e. assert the mesh is not moving -- you do this
   /// at your own risk!
   void disable_ALE() { ALE_is_disabled = true; }

   /// \short (Re-)enable ALE, i.e. take possible mesh motion into account
   /// when evaluating the time-derivative. Note: By default, ALE is
   /// enabled, at the expense of possibly creating unnecessary work
   /// in problems where the mesh is, in fact, stationary.
   void enable_ALE() { ALE_is_disabled = false; }
   
   /// \short Return the i-th pressure value at local pressure "node" n_p.
   /// Uses suitably interpolated value for hanging nodes.
   virtual double p_linearised_nst(const unsigned &n_p,
                                       const unsigned &i) const=0; 

   /// \short Pin the real or imaginary part of the problem
   /// Input integer 0 for real 1 for imaginary
   virtual void pin_real_or_imag(const unsigned &real)=0;
   virtual void unpin_real_or_imag(const unsigned &real)=0;
   
   ///Pin the normalisation dofs
   virtual void pin_pressure_normalisation_dofs()=0;
   
   /// Which nodal value represents the pressure?
   //  N.B. This function has return type "int" (rather than "unsigned"
   //  as in the u_index case) so that we can return the "magic" number
   //  "Pressure_not_stored_at_node" ( = -100 )
   virtual inline int p_index_linearised_nst(const unsigned &i)
    const { return Pressure_not_stored_at_node; }

   /// \short Strain-rate tensor: \f$ e_{ij} \f$
   /// where \f$ i,j = r,z,\theta \f$ (in that order)
   void strain_rate(const Vector<double>& s, 
                    DenseMatrix<double>& strain_rate, const unsigned &real)
    const;
 
   /// \short Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
   /// in tecplot format. Default number of plot points
   void output(std::ostream &outfile)
   {
    const unsigned nplot = 5;
    output(outfile,nplot);
   }

   /// \short Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
   /// in tecplot format. nplot points in each coordinate direction
   void output(std::ostream &outfile, const unsigned &nplot);
   
   /// \short Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
   /// in tecplot format. Default number of plot points
   void output(FILE* file_pt)
   {
    const unsigned nplot = 5;
    output(file_pt,nplot);
   }

   /// \short Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S, P^C, P^S
   /// in tecplot format. nplot points in each coordinate direction
   void output(FILE* file_pt, const unsigned &nplot);
   
   /// \short Output function: r, z, U^C, U^S, V^C, V^S, W^C, W^S,
   /// in tecplot format. nplot points in each coordinate direction
   /// at timestep t (t=0: present; t>0: previous timestep)
   void output_veloc(std::ostream &outfile, const unsigned &nplot,
                     const unsigned& t);
   
   /// Compute the element's residual Vector
   void fill_in_contribution_to_residuals(Vector<double> &residuals)
   {
    // Call the generic residuals function with flag set to 0
    // and using a dummy matrix argument
    fill_in_generic_residual_contribution_linearised_nst(
     residuals,
     GeneralisedElement::Dummy_matrix,
     GeneralisedElement::Dummy_matrix,0);
   }

   /// \short Compute the element's residual Vector and the jacobian matrix.
   /// Virtual function can be overloaded by hanging-node version.
   /*void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                         DenseMatrix<double> &jacobian)
   {
    // Call the generic routine with the flag set to 1
    fill_in_generic_residual_contribution_linearised_nst(
     residuals,jacobian,GeneralisedElement::Dummy_matrix,1);
     }*/

   /// \short Add the element's contribution to its residuals vector,
   /// jacobian matrix and mass matrix
   /*void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double> &residuals, DenseMatrix<double> &jacobian, 
    DenseMatrix<double> &mass_matrix)
   {
    // Call the generic routine with the flag set to 2
    fill_in_generic_residual_contribution_linearised_nst(
     residuals,jacobian,mass_matrix,2);
     }*/
   
   /// \short Return the i-th component of the FE interpolated velocity
   /// u[i] at local coordinate s
   double interpolated_u_linearised_nst(const Vector<double> &s, 
                                            const unsigned &i) const
   {
    // Determine number of nodes in the element
    const unsigned n_node = nnode();
    
    // Provide storage for local shape functions
    Shape psi(n_node);

    // Find values of shape functions
    shape(s,psi);
   
    // Get the index at which the velocity is stored
    const unsigned u_nodal_index = u_index_linearised_nst(i);

    // Initialise value of u
    double interpolated_u = 0.0;
   
    // Loop over the local nodes and sum
    for(unsigned l=0;l<n_node;l++) 
     {
      interpolated_u += nodal_value(l,u_nodal_index)*psi[l];
     }
    
    return(interpolated_u);
   }

   /// \short Return the i-th component of the FE interpolated pressure
   /// p[i] at local coordinate s
   double interpolated_p_linearised_nst(const Vector<double> &s,
                                            const unsigned &i) const
   {
    // Determine number of pressure nodes in the element
    const unsigned n_pressure_nodes = npres_linearised_nst();
    
    // Provide storage for local shape functions
    Shape psi(n_pressure_nodes);
    
    // Find values of shape functions
    pshape_linearised_nst(s,psi);
    
    // Initialise value of p
    double interpolated_p = 0.0;
    
    // Loop over the local nodes and sum
    for(unsigned l=0;l<n_pressure_nodes;l++) 
     {
      // N.B. The pure virtual function p_linearised_nst(...)
      // automatically calculates the index at which the pressure value
      // is stored, so we don't need to worry about this here
      interpolated_p += p_linearised_nst(l,i)*psi[l];
     }
    
    return(interpolated_p);
   }
   
 }; // End of LinearisedNavierStokesEquations class definition



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


 
 //=======================================================================
 /// Crouzeix-Raviart elements are Navier-Stokes elements with quadratic
 /// interpolation for velocities and positions, but a discontinuous
 /// linear pressure interpolation
 //=======================================================================
 class LinearisedQCrouzeixRaviartElement
  : public virtual QElement<2,3>, 
  public virtual LinearisedNavierStokesEquations
  {
    private:
   
   /// Static array of ints to hold required number of variables at nodes
   static const unsigned Initial_Nvalue[];
   
    protected:
   
   /// \short Internal indices that indicate at which internal data the
   /// pressure values are stored. We note that there are two pressure
   /// values, corresponding to the functions P^C(r,z,t) and P^S(r,z,t)
   /// which multiply the cosine and sine terms respectively.
   Vector<unsigned> P_linearised_nst_internal_index;
   
   /// \short Velocity shape and test functions and their derivatives 
   /// w.r.t. global coordinates at local coordinate s (taken from geometry).
   /// Return Jacobian of mapping between local and global coordinates.
   inline double dshape_and_dtest_eulerian_linearised_nst(
    const Vector<double> &s, 
    Shape &psi, DShape &dpsidx,
    Shape &test, DShape &dtestdx) const;
   
   /// \short Velocity shape and test functions and their derivatives
   /// w.r.t. global coordinates at the ipt-th integation point
   /// (taken from geometry).
   /// Return Jacobian of mapping between local and global coordinates.
   inline double dshape_and_dtest_eulerian_at_knot_linearised_nst(
    const unsigned &ipt, 
    Shape &psi, DShape &dpsidx, 
    Shape &test, DShape &dtestdx) const;

   /// Compute the pressure shape functions at local coordinate s
   inline void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi) const;

   /// Compute the pressure shape and test functions at local coordinate s
   inline void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi, Shape &test) const;

    public:
   
   /// \short Constructor: there are three internal values for each
   /// of the two pressure components
   LinearisedQCrouzeixRaviartElement() : QElement<2,3>(),
    LinearisedNavierStokesEquations(),
    P_linearised_nst_internal_index(2)
    {
     // Loop over the two pressure components
     // and two normalisation constraints
     for(unsigned i=0;i<4;i++)
      {
       // Allocate and add one internal data object for each of the two
       // pressure components that store the three pressure values
       P_linearised_nst_internal_index[i]
        = this->add_internal_data(new Data(3));
      }
    }
    
    /// Return number of values (pinned or dofs) required at local node n
    virtual unsigned required_nvalue(const unsigned &n) const;
    
    /// \short Return the pressure value i at internal dof i_internal
    /// (Discontinous pressure interpolation -- no need to cater for
    /// hanging nodes)
    double p_linearised_nst(const unsigned &i_internal,
                                const unsigned &i) const
    {
     return internal_data_pt(P_linearised_nst_internal_index[i])
      ->value(i_internal);
    }

    //Pin the normalisation dofs
    void pin_pressure_normalisation_dofs()
    {
     for(unsigned i=2;i<4;i++)
      {
       this->internal_data_pt(P_linearised_nst_internal_index[i])->pin_all();
      }
    }

   void pin_real_or_imag(const unsigned &real_index)
   {
    unsigned n_node=this->nnode();
    for(unsigned n=0;n<n_node;n++)
     {
      Node* nod_pt = this->node_pt(n);
         
      for(unsigned i=0;i<DIM;++i)
       {
        //Provided it's not constrained then pin it
        if(!nod_pt->is_constrained(i))
         {
          this->node_pt(n)->pin(2*i +real_index);
         }
       }
     }

    //Similarly for the pressure
    this->internal_data_pt(P_linearised_nst_internal_index[real_index])
     ->pin_all();
   }
  
   void unpin_real_or_imag(const unsigned &real_index)
   {
    unsigned n_node=this->nnode();
    for(unsigned n=0;n<n_node;n++)
     {
      Node* nod_pt = this->node_pt(n);
         
      for(unsigned i=0;i<DIM;++i)
       {
        //Provided it's not constrained then unpin it
        if(!nod_pt->is_constrained(i))
         {
          nod_pt->unpin(2*i +real_index);
         }
       }
     }
    
    //Similarly for the pressure
    this->internal_data_pt(P_linearised_nst_internal_index[real_index])
     ->unpin_all();
   }

   void copy_efunction_to_normalisation()
   {
    unsigned n_node = this->nnode();
    for(unsigned n=0;n<n_node;n++)
     {
      Node* nod_pt = this->node_pt(n);

      //Transfer the eigenfunctions to the normalisation constraints
      for(unsigned i=0;i<DIM;++i)
       {
        for(unsigned j=0;j<2;++j)
         {
          nod_pt->set_value(2*(DIM+i) + j,nod_pt->value(2*i + j));
         }
       }
     }
    
    //Similarly for the pressure
    for(unsigned i=0;i<2;++i)
     {
      Data* local_data_pt =
       this->internal_data_pt(P_linearised_nst_internal_index[i]);
      Data* norm_local_data_pt =
       this->internal_data_pt(P_linearised_nst_internal_index[2+i]);
      for(unsigned j=0;j<3;j++)
       {
        norm_local_data_pt->set_value(j,local_data_pt->value(j));
       }
     }
   }

    
    /// \short Return number of pressure values corresponding to a
    /// single pressure component
    unsigned npres_linearised_nst() const { return 3; }

    /// \short Fix both components of the internal pressure degrees
    /// of freedom p_dof to pvalue
    void fix_pressure(const unsigned &p_dof, const double &pvalue)
    {
     // Loop over the two pressure components
     for(unsigned i=0;i<2;i++)
      {
       this->internal_data_pt(P_linearised_nst_internal_index[i])
        ->pin(p_dof);
       internal_data_pt(P_linearised_nst_internal_index[i])
        ->set_value(p_dof,pvalue);
      }
    }

    /// \short Overload the access function for the i-th component of the 
    /// pressure's local equation numbers
    inline int p_local_eqn(const unsigned &n, const unsigned &i) 
    {
     return internal_local_eqn(P_linearised_nst_internal_index[i],n);
    }
   
    /// Redirect output to NavierStokesEquations output
    void output(std::ostream &outfile) 
    { LinearisedNavierStokesEquations::output(outfile); }
    
    /// Redirect output to NavierStokesEquations output
    void output(std::ostream &outfile, const unsigned &n_plot)
    { LinearisedNavierStokesEquations::output(outfile,n_plot); }
    
    /// Redirect output to NavierStokesEquations output
    void output(FILE* file_pt) 
    { LinearisedNavierStokesEquations::output(file_pt); }
    
    /// Redirect output to NavierStokesEquations output
    void output(FILE* file_pt, const unsigned &n_plot)
    { LinearisedNavierStokesEquations::output(file_pt,n_plot); }
    
    /// \short The number of "dof-blocks" that degrees of freedom in this
    /// element are sub-divided into: Velocity and pressure.
    unsigned ndof_types() const { return 2*(DIM+1); }
        
  }; // End of LinearisedQCrouzeixRaviartElement class definition
 
 
 // Inline functions
 // ----------------

 //=======================================================================
 /// \short Derivatives of the shape functions and test functions w.r.t.
 /// global (Eulerian) coordinates at local coordinate s.
 /// Return Jacobian of mapping between local and global coordinates.
 //=======================================================================
 inline double LinearisedQCrouzeixRaviartElement::
  dshape_and_dtest_eulerian_linearised_nst(
   const Vector<double> &s,
   Shape &psi, DShape &dpsidx,
   Shape &test, DShape &dtestdx) const
  {
   // Call the geometrical shape functions and derivatives  
   const double J = this->dshape_eulerian(s,psi,dpsidx);
   //The test functions are equal to the shape functions
   test = psi;
   dtestdx = dpsidx;
   // Return the Jacobian
   return J;
  }

 //=======================================================================
 /// \short Derivatives of the shape functions and test functions w.r.t.
 /// global (Eulerian) coordinates at the ipt-th integration point.
 /// Return Jacobian of mapping between local and global coordinates.
 //=======================================================================
 inline double LinearisedQCrouzeixRaviartElement::
  dshape_and_dtest_eulerian_at_knot_linearised_nst(
   const unsigned &ipt, Shape &psi, 
   DShape &dpsidx, Shape &test, 
   DShape &dtestdx) const
  {
   
   // Call the geometrical shape functions and derivatives  
   const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

   // Loop over the test functions and derivatives and set them
   // equal to the shape functions
   test = psi;
   dtestdx = dpsidx;
   // Return the Jacobian
   return J;
  }

 //=======================================================================
 /// Pressure shape functions
 //=======================================================================
 inline void LinearisedQCrouzeixRaviartElement::
  pshape_linearised_nst(const Vector<double> &s, Shape &psi) const
  {
   psi[0] = 1.0;
   psi[1] = s[0];
   psi[2] = s[1];
  }
 
 //=======================================================================
 /// Define the pressure shape and test functions
 //=======================================================================
 inline void LinearisedQCrouzeixRaviartElement::
  pshape_linearised_nst(const Vector<double> &s,
                            Shape &psi, Shape &test) const
  {
   // Call the pressure shape functions
   pshape_linearised_nst(s,psi);
   
   // Loop over the test functions and set them equal to the shape functions
   for(unsigned i=0;i<3;i++) { test[i] = psi[i]; }
  }
 
 //=======================================================================
 /// Face geometry of the linearised axisym Crouzeix-Raviart elements
 //=======================================================================
 template<>
  class FaceGeometry<LinearisedQCrouzeixRaviartElement>
  : public virtual QElement<1,3>
  {
    public:
   FaceGeometry() : QElement<1,3>() {}
  };
 
 //=======================================================================
 /// \short Face geometry of face geometry of the linearised axisymmetric
 /// Crouzeix Raviart elements
 //=======================================================================
 template<>
  class FaceGeometry<FaceGeometry
  <LinearisedQCrouzeixRaviartElement> >
  : public virtual PointElement
 {
   public:
  FaceGeometry() : PointElement() {}
 };
 
 
 
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
 
 
 
 //=======================================================================
 /// Taylor--Hood elements are Navier--Stokes elements with quadratic
 /// interpolation for velocities and positions and continuous linear
 /// pressure interpolation
 //=======================================================================
 class LinearisedQTaylorHoodElement
  : public virtual QElement<2,3>, 
  public virtual LinearisedNavierStokesEquations
  {
    private:
 
   /// Static array of ints to hold number of variables at node
   static const unsigned Initial_Nvalue[];
   
    protected:
   
   /// \short Static array of ints to hold conversion from pressure 
   /// node numbers to actual node numbers
   static const unsigned Pconv[];
   
   /// \short Velocity shape and test functions and their derivatives
   /// w.r.t. global coordinates  at local coordinate s (taken from geometry).
   /// Return Jacobian of mapping between local and global coordinates.
   inline double dshape_and_dtest_eulerian_linearised_nst(
    const Vector<double> &s, 
    Shape &psi, DShape &dpsidx,
    Shape &test, DShape &dtestdx) const;
   
   /// \short Velocity shape and test functions and their derivatives
   /// w.r.t. global coordinates the ipt-th integation point
   /// (taken from geometry).
   /// Return Jacobian of mapping between local and global coordinates.
   inline double dshape_and_dtest_eulerian_at_knot_linearised_nst(
    const unsigned &ipt, 
    Shape &psi, DShape &dpsidx, 
    Shape &test, DShape &dtestdx) const;
 
   /// Compute the pressure shape functions at local coordinate s
   inline void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi) const;
 
   /// Compute the pressure shape and test functions at local coordinte s
   inline void pshape_linearised_nst(const Vector<double> &s,
                                         Shape &psi, Shape &test) const;
   
    public:
   
   /// Constructor, no internal data points
   LinearisedQTaylorHoodElement() : QElement<2,3>(),  
    LinearisedNavierStokesEquations() {}
   
    /// \short Number of values (pinned or dofs) required at node n. Can
    /// be overwritten for hanging node version
    inline virtual unsigned required_nvalue(const unsigned &n) const 
     { return Initial_Nvalue[n]; }

    /// \short Which nodal value represents the pressure? Overload version
    /// in base class which returns static int "Pressure_not_stored_at_node"
    virtual int p_index_linearised_nst(const unsigned &i) const
    {
     return (2*DIM+i);
    }
   
    /// \short Access function for the i-th component of pressure
    /// at local pressure node n_p (const version).
    double p_linearised_nst(const unsigned &n_p,
                                const unsigned &i) const
    {
     return nodal_value(Pconv[n_p],p_index_linearised_nst(i));
    }

    //Pin the normalisation dofs
    void pin_pressure_normalisation_dofs()
    {
     throw OomphLibError("This is not implemented yet\n",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }


    virtual void pin_real_or_imag(const unsigned &real) 
    {
     throw OomphLibError("This is not implemented yet\n",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
    
   virtual void unpin_real_or_imag(const unsigned &real)
    {
     throw OomphLibError("This is not implemented yet\n",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

    
    /// \short Return number of pressure values corresponding to a
    /// single pressure component
    unsigned npres_linearised_nst() const { return 4; }
    
    /// \short Fix both components of the pressure at local pressure
    /// node n_p to pvalue
    void fix_pressure(const unsigned &n_p, const double &pvalue)
    {
     // Loop over the two pressure components
     for(unsigned i=0;i<2;i++)
      {
       this->node_pt(Pconv[n_p])->pin(p_index_linearised_nst(i));
       this->node_pt(Pconv[n_p])
        ->set_value(p_index_linearised_nst(i),pvalue);
      }
    }

    /// \short Overload the access function for the i-th component of the
    /// pressure's local equation numbers
    inline int p_local_eqn(const unsigned &n, const unsigned &i)
     {
      return nodal_local_eqn(Pconv[n],p_index_linearised_nst(i));
     }
    
    /// Redirect output to NavierStokesEquations output
    void output(std::ostream &outfile) 
    { LinearisedNavierStokesEquations::output(outfile); }
    
    /// Redirect output to NavierStokesEquations output
    void output(std::ostream &outfile, const unsigned &n_plot)
    { LinearisedNavierStokesEquations::output(outfile,n_plot); }
    
    /// Redirect output to NavierStokesEquations output
    void output(FILE* file_pt) 
    { LinearisedNavierStokesEquations::output(file_pt); }
    
    /// Redirect output to NavierStokesEquations output
    void output(FILE* file_pt, const unsigned &n_plot)
    { LinearisedNavierStokesEquations::output(file_pt,n_plot); }
    
    /// \short Returns the number of "dof-blocks" that degrees of freedom
    /// in this element are sub-divided into: Velocity and pressure.
    unsigned ndof_types() const { return 8; }
    
  }; // End of LinearisedQTaylorHoodElement class definition
 

 // Inline functions
 // ----------------

 //=======================================================================
 /// \short Derivatives of the shape functions and test functions w.r.t
 /// global (Eulerian) coordinates at local coordinate s.
 /// Return Jacobian of mapping between local and global coordinates.
 //=======================================================================
 inline double LinearisedQTaylorHoodElement::
  dshape_and_dtest_eulerian_linearised_nst(
   const Vector<double> &s,
   Shape &psi, DShape &dpsidx,
   Shape &test, DShape &dtestdx) const
  {
   // Call the geometrical shape functions and derivatives  
   const double J = this->dshape_eulerian(s,psi,dpsidx);

   test = psi;
   dtestdx = dpsidx;
   
   // Return the Jacobian
   return J;
  }
 
 //=======================================================================
 /// \short Derivatives of the shape functions and test functions w.r.t
 /// global (Eulerian) coordinates at the ipt-th integration point.
 /// Return Jacobian of mapping between local and global coordinates.
 //=======================================================================
 inline double LinearisedQTaylorHoodElement::
  dshape_and_dtest_eulerian_at_knot_linearised_nst(
   const unsigned &ipt,
   Shape &psi, DShape &dpsidx, 
   Shape &test, DShape &dtestdx) const
  {
   // Call the geometrical shape functions and derivatives  
   const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

   test=psi;
   dtestdx = dpsidx;
   
   // Return the Jacobian
   return J;
  }
 
 //=======================================================================
 /// Pressure shape functions
 //=======================================================================
 inline void LinearisedQTaylorHoodElement::
  pshape_linearised_nst(const Vector<double> &s, Shape &psi) const
  {
   // Allocate local storage for the pressure shape functions
   double psi1[2], psi2[2];
   
   // Call the one-dimensional shape functions
   OneDimLagrange::shape<2>(s[0],psi1);
   OneDimLagrange::shape<2>(s[1],psi2);
   
   // Now let's loop over the nodal points in the element
   // s1 is the "r" coordinate, s2 the "z" 
   for(unsigned i=0;i<2;i++)
    {
     for(unsigned j=0;j<2;j++)
      {
       // Multiply the two 1D functions together to get the 2D function
       psi[2*i + j] = psi2[i]*psi1[j];
      }
    }
  }

 //=======================================================================
 /// Pressure shape and test functions
 //=======================================================================
 inline void LinearisedQTaylorHoodElement::
  pshape_linearised_nst(const Vector<double> &s,
                            Shape &psi, Shape &test) const
  {
   // Call the pressure shape functions
   pshape_linearised_nst(s,psi);
   
   // Loop over the test functions and set them equal to the shape functions
   for(unsigned i=0;i<4;i++) { test[i] = psi[i]; }
  }
 
 //=======================================================================
 /// Face geometry of the linearised axisymmetric Taylor Hood elements
 //=======================================================================
 template<>
  class FaceGeometry<LinearisedQTaylorHoodElement>
  : public virtual QElement<1,3>
  {
    public:
   FaceGeometry() : QElement<1,3>() {}
  };
 
 //=======================================================================
 /// \short Face geometry of the face geometry of the linearised
 /// axisymmetric Taylor Hood elements
 //=======================================================================
 template<>
  class FaceGeometry<FaceGeometry
  <LinearisedQTaylorHoodElement> >
  : public virtual PointElement
 {
   public:
  FaceGeometry() : PointElement() {}
 };
 
 
} // End of oomph namespace

#endif
