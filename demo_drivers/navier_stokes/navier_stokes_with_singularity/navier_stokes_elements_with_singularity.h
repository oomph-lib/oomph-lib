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
//Header file for Navier Stokes elements with singularity
#ifndef OOMPH_NAVIER_STOKES_ELEMENTS_WITH_SINGULARITY_HEADER
#define OOMPH_NAVIER_STOKES_ELEMENTS_WITH_SINGULARITY_HEADER


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
 
namespace oomph
{

//==================CLASS FOR THE ADDITIONAL UNKNOWN==================
/// We consider a singularity in the solution at the point O.
///
/// This class defines the singular functions.
/// 
/// velocity_singularity:
/// 
///    u_bar = C*velocity_singular_function
/// 
/// pressure_singularity:
///
///    p_bar = C*pressure_singular_function
///
/// and their gradients
///
/// The class also defines the function that computes the
/// residual associated to C which is:
///
/// R_C = \frac{\partial p_FE}{\partial x_i} (O)
///
/// and thus regularises the FE solution by setting the pressure gradient
/// in the coordinate direction x_i to zero. If the amplitude of the
/// singular solution is known, pin C.
//=====================================================================
template<class WRAPPED_NAVIER_STOKES_ELEMENT>
class SingularNavierStokesSolutionElement : public virtual GeneralisedElement
{

 public:

 /// Function pointer to the velocity singular function:
 typedef Vector<double> (*NavierStokesVelocitySingularFctPt)
  (const Vector<double>& x);

 /// Function pointer to the gradient of the velocity singular function:
 typedef Vector<Vector<double> > (*NavierStokesGradVelocitySingularFctPt)
  (const Vector<double>& x);

 /// Function pointer to the pressure singular function:
 typedef double (*NavierStokesPressureSingularFctPt)
  (const Vector<double>& x);

 /// Function pointer to the gradient of the pressure singular function:
 typedef Vector<double> (*NavierStokesGradPressureSingularFctPt)
  (const Vector<double>& x);

 /// Constructor
 SingularNavierStokesSolutionElement()
  {
   // Initialise Function pointer to velocity singular function to NULL
   Velocity_singular_fct_pt=0;

   // Initialise Function pointer to gradient of velocity singular
   // function to NULL
   Grad_velocity_singular_fct_pt=0;

   // Initialise Function pointer to pressure singular function to NULL
   Pressure_singular_fct_pt=0;

   // Initialise Function pointer to gradient of pressure singular
   // function to NULL
   Grad_pressure_singular_fct_pt=0;

   // Initalise pointer to the wrapped Navier-Stokes element which will be
   // used to compute the residual and which includes the point O
   Wrapped_navier_stokes_el_pt=0;

   // Initialise the pointer to the direction of the derivative used
   // for the residual of this element
   Direction_pt=0;

   // Safe assumption: Singular fct does not satisfy Stokes eqn
   Singular_function_satisfies_stokes_equation=false;

   // Create a single item of internal Data, storing one unknown which
   // represents the unknown C.
   add_internal_data(new Data(1));
  }

 /// Assert that singular function satisfies the Stokes equations by setting
 /// this to true or false.
 bool& singular_function_satisfies_stokes_equation()
  {
   return Singular_function_satisfies_stokes_equation;
  }

 /// Find the value of the unknown amplitude C
 double c() const 
 {
  return internal_data_pt(0)->value(0);
 }


 /// Find the value of the unknown amplitude C
 void set_c(const double& value) 
 {
  internal_data_pt(0)->set_value(0,value);
 }

 /// Pin the value of the unknown amplitude C
 void pin_c() 
 {
  return internal_data_pt(0)->pin(0);
 }

 /// Unpin the value of the unknown amplitude C
 void unpin_c() 
 {
  return internal_data_pt(0)->unpin(0);
 }
 
 ///  Set pointer to associated wrapped Navier-Stokes element which
 /// contains the singularity (at local coordinate s). Also specify the 
 /// direction in which the slope of the FE part of the pressure is 
 /// set to zero. (Could also set a velocity derivative to zero but this
 /// needs to be done with a separate function. Write it if you need it...)
 void set_wrapped_navier_stokes_element_pt
  (WRAPPED_NAVIER_STOKES_ELEMENT* wrapped_navier_stokes_el_pt,
   const Vector<double>& s,
   unsigned* direction_pt)
 {
  // Assign the pointer to the variable Wrapped_navier_stokes_el_pt
  Wrapped_navier_stokes_el_pt = wrapped_navier_stokes_el_pt;
  
  // Find number of nodes in the element
  unsigned nnod=wrapped_navier_stokes_el_pt->nnode();
  
  // Loop over the nodes of the element
  for (unsigned j=0;j<nnod;j++)
   {
    // Add the node as external data in the 
    // SingularNavierStokesSolutionElement class. Note that this
    // assumes that the pressure is stored at the nodes (Taylor Hood type
    // NSt elements, which is assumed elsewhere too...)
    add_external_data(Wrapped_navier_stokes_el_pt->node_pt(j));
   }
  
  // Assign the pointer to the local coordinate at which the residual
  // will be computed
  S_in_wrapped_navier_stokes_element = s;
  
  // Assign the pointer to the direction at which the derivative used
  // in the residual will be computed
  Direction_pt = direction_pt;
 }

 /// Access function to pointer to velocity singular function
 NavierStokesVelocitySingularFctPt& velocity_singular_fct_pt()
  {return Velocity_singular_fct_pt;}
 
 /// Access function to pointer to gradient of velocity singular function
 NavierStokesGradVelocitySingularFctPt& grad_velocity_singular_fct_pt()
  {return Grad_velocity_singular_fct_pt;}
 
 /// Access function to pointer to pressure singular function
 NavierStokesPressureSingularFctPt& pressure_singular_fct_pt()
  {return Pressure_singular_fct_pt;}
 
 /// Access function to pointer to gradient of pressure singular function
 NavierStokesGradPressureSingularFctPt& grad_pressure_singular_fct_pt()
  {return Grad_pressure_singular_fct_pt;}
 
 /// Evaluate velocity singular function at Eulerian position x
 Vector<double> velocity_singular_function(const Vector<double>& x) const
  {
#ifdef PARANOID
   if (Velocity_singular_fct_pt==0)
    {
     std::stringstream error_stream;
     error_stream 
      << "Pointer to velocity singular function hasn't been defined!"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // Evaluate velocity singular function
   return (*Velocity_singular_fct_pt)(x);
  }
 
 ///  Evaluate gradient of velocity singular function at Eulerian 
 /// position x. grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_velocity_singular_function(
  const Vector<double>& x) const
  {
#ifdef PARANOID
   if (Grad_velocity_singular_fct_pt==0)
    {
     std::stringstream error_stream;
     error_stream 
      << "Pointer to gradient of velocity singular function "
      << "hasn't been defined!"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // Evaluate gradient of velocity singular function
   return (*Grad_velocity_singular_fct_pt)(x);
  }
 
 /// Evaluate pressure singular function at Eulerian position x
 double pressure_singular_function(const Vector<double>& x) const
 {
#ifdef PARANOID
  if (Pressure_singular_fct_pt==0)
   {
    std::stringstream error_stream;
    error_stream 
     << "Pointer to pressure singular function hasn't been defined!"
     << std::endl;
    throw OomphLibError(
     error_stream.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Evaluate pressure singular function
  return (*Pressure_singular_fct_pt)(x);
 }
 
 /// Evaluate gradient of pressure singular function at Eulerian position x
 Vector<double> grad_pressure_singular_function(const Vector<double>& x) const 
  {
#ifdef PARANOID
   if (Grad_pressure_singular_fct_pt==0)
    {
     std::stringstream error_stream;
     error_stream 
      << "Pointer to gradient of pressure singular function "
      << "hasn't been defined!"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // Evaluate gradient of pressure singular function
   return (*Grad_pressure_singular_fct_pt)(x);
  }
 
 ///  Evaluate velocity singular function (including its amplitude):
 /// u_bar = C * velocity_singular
 Vector<double> u_bar(const Vector<double>& x)
  {
   // Find the value of C
   double c=internal_data_pt(0)->value(0);
   
   // Initialise the velocity vector
   Vector<double> u=velocity_singular_function(x);
   
   // Find the dimension of the problem
   unsigned cached_dim=Wrapped_navier_stokes_el_pt->dim();
   
   // Multiply the components of the velocity vector by the unknown C
   for (unsigned d=0;d<cached_dim;d++)
    {
     u[d] *= c;
    }
   
   // Value of u_bar at the position x
   return u;
  }
 
 ///  Evaluate gradient of velocity singular function
 /// (including its amplitude):
 /// grad_u_bar = C * grad_velocity_singular;
 /// grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_u_bar(const Vector<double>& x)
  {
   // Find the value of C
   double c=internal_data_pt(0)->value(0);
   
   // Initialise the gradient of velocity vector
   Vector<Vector<double> > grad_u=grad_velocity_singular_function(x);
   
   // Find the dimension of the problem 
   unsigned cached_dim=Wrapped_navier_stokes_el_pt->dim();
   
   // Multiply the components of the gradient of velocity by the unknown C
   for (unsigned d=0;d<cached_dim;d++)
    {
     for (unsigned i=0;i<cached_dim;i++)
      {
       grad_u[d][i] *= c;
      }
    }
   
   // Value of grad_u_bar at the position x
   return grad_u;
  }
 
 ///  Evaluate pressure singular function (including its amplitude):
 /// p_bar = C * pressure_singular
  double p_bar(const Vector<double>& x)
 {
  // Find the value of C
  double c=internal_data_pt(0)->value(0);
  
  // Value of p_bar at the position x
  return c*pressure_singular_function(x);
 }
 
 ///  Evaluate gradient of pressure singular function
 /// (including its amplitude):
 /// grad_p_bar = C * grad_pressure_singular
 Vector<double> grad_p_bar(const Vector<double>& x)
 {
  // Find the value of C
  double c=internal_data_pt(0)->value(0);

  // Initialise the gradient of pressure
  Vector<double> grad_p=grad_pressure_singular_function(x);

  // Find the dimension of the problem 
  unsigned cached_dim=Wrapped_navier_stokes_el_pt->dim();

  // Multiply the components of the gradient of pressure by the unknown C
  for (unsigned d=0;d<cached_dim;d++)
   {
    grad_p[d] *= c;
   }
  
  // Value of grad_p_bar at the position x
  return grad_p;
 }

 /// Compute residual
 void fill_in_contribution_to_residuals(Vector<double>& residual)
 {
  fill_in_generic_contribution_to_residuals
   (residual,GeneralisedElement::Dummy_matrix,0);
 }
 
 // Compute local residual and jacobian
 void fill_in_contribution_to_jacobian(Vector<double>& residual,
                                       DenseMatrix<double> &jacobian)
 {
  fill_in_generic_contribution_to_residuals(residual,jacobian,1);
 }
 
  private:
 
 /// Compute local residual, and, if flag=1, local jacobian matrix
 void fill_in_generic_contribution_to_residuals(Vector<double>& residual,
                                                DenseMatrix<double> &jacobian,
                                                const unsigned& flag)
 {
  // Get the local eqn number of our one-and-only
  // unknown
  int eqn_number=internal_local_eqn(0,0);
  if (eqn_number>=0)
   {
    residual[eqn_number] = Wrapped_navier_stokes_el_pt->dpdx_fe_only
     (S_in_wrapped_navier_stokes_element,Direction_pt);  
    
    // Do we want the Jacobian too?
    if (flag)
     {
      // Find the number of pressure dofs in the wrapped Navier-Stokes 
      // element pointed by
      // the SingularNavierStokesSolutionElement class
      unsigned n_pres = Wrapped_navier_stokes_el_pt->npres_nst();
      
      // Find the dimension of the problem
      unsigned cached_dim=Wrapped_navier_stokes_el_pt->dim();

      // Set up memory for the pressure shape functions and their derivatives
      Shape psip(n_pres), testp(n_pres);
      DShape dpsipdx(n_pres,cached_dim), dtestpdx(n_pres,cached_dim);
      
      // Compute the pressure shape functions and their derivatives
      // at the local coordinate S_in_wrapped_navier_stokes_element
      // (Test fcts not really needed but nobody's got around to writing
      // a fct that only picks out the basis fcts.
      Wrapped_navier_stokes_el_pt->dpshape_and_dptest_eulerian_nst
       (S_in_wrapped_navier_stokes_element,psip,dpsipdx,testp,dtestpdx);
      
      // Derivs
      for (unsigned j=0;j<n_pres;j++)
       {
        // Unknown
        int local_unknown = Wrapped_navier_stokes_el_pt->p_local_eqn(j);
        
        // If not pinned
        if (local_unknown >= 0)
         {
          // Add the contribution of the node to the local jacobian
          jacobian(eqn_number,local_unknown) = dpsipdx(j,*Direction_pt);
         }
       }
     }
   }
 }
 
 ///  Pointer to wrapped Navier-Stokes element
 WRAPPED_NAVIER_STOKES_ELEMENT* Wrapped_navier_stokes_el_pt;

 ///  Pointer to velocity singular function
 NavierStokesVelocitySingularFctPt Velocity_singular_fct_pt;

 ///  Pointer to gradient of velocity singular function;
 /// grad[i][j] = du_i/dx_j
 NavierStokesGradVelocitySingularFctPt Grad_velocity_singular_fct_pt;

 ///  Pointer to pressure singular function
 NavierStokesPressureSingularFctPt Pressure_singular_fct_pt;

 ///  Pointer to gradient of pressure singular function
 NavierStokesGradPressureSingularFctPt Grad_pressure_singular_fct_pt;

 /// Local coordinates of singulariity in wrapped Navier-Stokes element
 Vector<double> S_in_wrapped_navier_stokes_element;
 
 ///  Direction of the derivative used for the residual of the element
 unsigned* Direction_pt;

// Does singular fct satisfy Stokes eqn?
bool Singular_function_satisfies_stokes_equation;
 
}; // End of SingularNavierStokesSolutionElement class


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


//=======================================================================
/// Templated wrapper to add handling of singularities to
/// the underlying Navier-Stokes element (specified via the template
/// argument). Slightly inefficient because the integration loop
/// is repeated here. The only alternative is to add the additional
/// functionality into the Navier-Stokes elements. Not pretty either!
/// NOTE Element ia assumed to be of Taylor Hood type with pressures
/// stored at nodes.
///
/// Dirichlet BCs are applied in hijacking manner and must be imposed
/// from each element that shares a boundary node that is subject to
/// a Dirichlet condition.
//=======================================================================
template<class BASIC_NAVIER_STOKES_ELEMENT>
class NavierStokesElementWithSingularity :
public virtual BASIC_NAVIER_STOKES_ELEMENT
{

 public:

 /// Constructor
 NavierStokesElementWithSingularity()
  {

   // Find the number of nodes in the element
   unsigned n_node = this->nnode();

   // Find the dimension of the problem
   unsigned cached_dim = this->dim();

   // Initialise the vector indicating which node is subject to velocity
   // Dirichlet BCs. The size of the vector is equal to the number of nodes
   // in the element. Each component of the vector is a vector of booleans
   // indicating if the velocity components at the corresponding node are
   // subject to Dirichlet BC. By default, no node is subject to Dirichlet BC,
   // so the vector is full of false. 
   Node_is_subject_to_velocity_dirichlet_bcs.resize(n_node);
   for (unsigned j=0;j<n_node;j++)
    {
     Node_is_subject_to_velocity_dirichlet_bcs[j].resize(cached_dim);
     for (unsigned d=0;d<cached_dim;d++)
      {
       Node_is_subject_to_velocity_dirichlet_bcs[j][d] = false;
      }
    }

   // Initialise the vector of imposed velocity values on the nodes
   // subject to Dirichlet BC. The size of the vector is equal to the
   // number of nodes in the element. Each component of the vector is
   // a vector of length the dimension of the problem. This vector contains
   // the imposed values of the velocity vector at the corresponding node.
   // If a node is not subject to Dirichlet BC, its imposed values are zero.
   // By default, no node is subject to Dirichlet BC so the vector is full
   // of zeros
   Imposed_velocity_values_at_node.resize(n_node);
   for (unsigned j=0;j<n_node;j++)
    {
     Imposed_velocity_values_at_node[j].resize(cached_dim);
     for (unsigned d=0;d<cached_dim;d++)
      {
       Imposed_velocity_values_at_node[j][d] = 0.0;
      }
    }

   // Find the number of pressure dofs
   unsigned n_pres = this->npres_nst();

   // Initialise the vector indicating which pressure dof is subject to
   // Dirichlet BCs. The size of the vector is equal to the number of pressure
   // dofs in the element. Each component of the vector is a boolean
   // indicating if the corresponding pressure unknown is subject to Dirichlet
   // BC. By default, no pressure dof is subject to Dirichlet BC, so the vector
   // is full of false. 
   Pressure_dof_is_subject_to_dirichlet_bc.resize(n_pres);
   for (unsigned l=0;l<n_pres;l++)
    {
     Pressure_dof_is_subject_to_dirichlet_bc[l] = false;
    }

   // Initialise the vector of imposed values on the pressure dofs
   // subject to Dirichlet BC. The size of the vector is equal to the
   // number of pressure dofs in the element. Each component of 
   // the vector contains the imposed value at the corresponding pressure dof.
   // If a pressure dof is not subject to Dirichlet BC, its imposed value 
   // is zero. By default, no pressure dof is subject to Dirichlet BC 
   // so the vector is full of zeros
   Imposed_value_at_pressure_dof.resize(n_node);
   for (unsigned l=0;l<n_pres;l++)
    {
     Imposed_value_at_pressure_dof[l] = 0.0;
    }
   
  } // End of constructor



 ///  Impose Dirichlet BC on the d-th component of the velocity
 /// (including the singular contribution) at the j-th node
 void impose_velocity_dirichlet_bc_on_node(const unsigned& j, const unsigned& d)
 {
  Node_is_subject_to_velocity_dirichlet_bcs[j][d] = true;
 }
 
 ///  Undo Dirichlet BC on the d-th velocity component (including the
 /// singular contribution) of the jth node
 void undo_velocity_dirichlet_bc_on_node(const unsigned& j, const unsigned& d)
 {
  Node_is_subject_to_velocity_dirichlet_bcs[j][d] = false;
 }
 
 ///  Specify Dirichlet boundary value for the d-th velocity component
 /// (including the singular contribution) at the j-th local node
 void set_velocity_dirichlet_value_on_node(const unsigned& j,
                                           const unsigned& d,
                                           const double& value)
 {
  Imposed_velocity_values_at_node[j][d] = value;
 }
 

 /// Impose Dirichlet BC at the j-th pressure dof
 void impose_dirichlet_bc_on_pressure_dof(const unsigned& j)
 {
  Pressure_dof_is_subject_to_dirichlet_bc[j] = true;
 }

 /// Undo Dirichlet BC at the j-th pressure dof
 void undo_dirichlet_bc_on_pressure_dof(const unsigned& j)
 {
  Pressure_dof_is_subject_to_dirichlet_bc[j] = false;
 }

 ///  Specify Dirichlet boundary value for the j-th pressure dof
 void set_dirichlet_value_on_pressure_dof(const unsigned& j,
                                          const double& value)
 {
  Imposed_value_at_pressure_dof[j] = value;
 }
 
 /// Access function to vector of pointers to SingularNavierStokesSolutionElements
 Vector<SingularNavierStokesSolutionElement<NavierStokesElementWithSingularity
  <BASIC_NAVIER_STOKES_ELEMENT> >*> c_equation_elements_pt()
  {return C_equation_elements_pt;}
 
 ///  Add pointer to associated SingularNavierStokesSolutionElement that determines the
 /// value of the amplitude of the singular functions (and gives access
 /// to the singular functions). The unknown amplitude becomes
 /// external Data for this element so assign_eqn_numbers() must be
 /// called after this function has been called. 
 void add_c_equation_element_pt(
  SingularNavierStokesSolutionElement<NavierStokesElementWithSingularity
  <BASIC_NAVIER_STOKES_ELEMENT> >* c_pt)
 {
  // Add the element
  C_equation_elements_pt.push_back(c_pt);
  
  // Add the additional unknown of this object as external data in the
  // Navier-Stokes element
  this->add_external_data(c_pt->internal_data_pt(0));
 }


 /// Derivative of pressure in direction indicated by pointer to unsigned
 double dpdx_fe_only(Vector<double> s, const unsigned* direction_pt)
 {
  // Find the number of pressure dofs in the wrapped Navier-Stokes 
  // element pointed by
  // the SingularNavierStokesSolutionElement class
  unsigned n_pres = this->npres_nst();

  // Find the dimension of the problem
  unsigned cached_dim=this->dim();

  // Set up memory for the pressure shape functions and their derivatives
  Shape psip(n_pres), testp(n_pres);
  DShape dpsipdx(n_pres,cached_dim), dtestpdx(n_pres,cached_dim);

  // Compute the pressure shape functions and their derivatives
  // at the local coordinate S_in_wrapped_navier_stokes_element
  // (Test fcts not really needed but nobody's got around to writing
  // a fct that only picks out the basis fcts.
  this->dpshape_and_dptest_eulerian_nst
   (s,psip,dpsipdx,testp,dtestpdx);
  
  // Initialise the derivative used for the residual
  double interpolated_dpdx_fe_only=0.0;

  // Compute the derivative used for the residual. The direction of the
  // derivative is given by *Direction_pt
  for (unsigned j=0;j<n_pres;j++)
   {
    interpolated_dpdx_fe_only +=
     this->p_nst(j)*dpsipdx(j,*direction_pt);
   }

  return interpolated_dpdx_fe_only;
 }


 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  // Call the generic residuals function with flag set to 0
  // using a dummy matrix argument
  this->fill_in_generic_residual_contribution_wrapped_nst(
   residuals, 
   GeneralisedElement::Dummy_matrix,
   0);
 }
 
 ///  Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
 {  
  // Call the generic routine with the flag set to 1 and dummy mass matrix
  DenseMatrix<double> mass_matrix;
  this->fill_in_generic_residual_contribution_wrapped_nst(residuals,
                                                          jacobian,
                                                          1);

  /* // FD version if you want to mess around with anything...*/
  /* FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian); */
 }
 
 ///  Overload the output function 
 /// x, y, [z,] u, v, [w], p, u_fe, v_fe, [w_fe], p_fe,
 // u_sing, v_sing, [w_sing], p_sing
 void output(std::ostream &outfile,const unsigned &nplot)
 {
  // Find the dimension of the problem
  unsigned cached_dim = this->dim();

  //Vector of local coordinates
  Vector<double> s(cached_dim);
 
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
 
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {   
    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);

    Vector<double> x(cached_dim);
    for(unsigned i=0;i<cached_dim;i++) 
     {
      outfile << this->interpolated_x(s,i) << " ";
      x[i] = this->interpolated_x(s,i);
     }

    Vector<double> velocity = interpolated_u_nst(s);
    Vector<double> velocity_fe_only = interpolated_u_nst_fe_only(s);
    for (unsigned i=0;i<cached_dim;i++) 
     {
      outfile << velocity[i] << " ";
     }
    outfile << this->interpolated_p_nst(s) << " ";
    for (unsigned i=0;i<cached_dim;i++) 
     {
      outfile << velocity_fe_only[i] << " ";
     }
    outfile << this->interpolated_p_nst_fe_only(s) << " ";
    for (unsigned i=0;i<cached_dim;i++) 
     {
      outfile << velocity[i] - velocity_fe_only[i] << " ";
     }
    outfile 
     << this->interpolated_p_nst(s)-
     this->interpolated_p_nst_fe_only(s) 
     << " " << std::endl;
    
   }

  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);

 }


 /// Overloaded compute error function; uses FE+singular parts
void compute_error(std::ostream &outfile,
                   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                   const bool& include_pressure,
                   double& error, double& norm)
{

 unsigned cached_dim=this->dim();

 error=0.0;
 norm=0.0;

 //Vector of local coordinates
 Vector<double> s(cached_dim);

 // Vector for coordintes
 Vector<double> x(cached_dim);

 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
   

 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (u,v,[w],p)
 Vector<double> exact_soln(cached_dim+1);
 Vector<double> computed_soln(cached_dim+1);
   
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<cached_dim;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);

   // Get jacobian of mapping
   double J=this->J_eulerian(s);

   //Premultiply the weights and the Jacobian
   double W = w*J;

   // Get x position as Vector
   this->interpolated_x(s,x);

   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   Vector<double> u_comp=interpolated_u_nst(s);
   for(unsigned i=0;i<cached_dim;i++)
    {
     computed_soln[i]=u_comp[i];
    }
   unsigned hi_limit=cached_dim;
   if (include_pressure)
    {
     computed_soln[cached_dim]=interpolated_p_nst(s); 
     hi_limit=cached_dim+1;
    }


   // Velocity error
   for(unsigned i=0;i<hi_limit;i++)
    {
     norm+=exact_soln[i]*exact_soln[i]*W;
     error+=(exact_soln[i]-computed_soln[i])*
      (exact_soln[i]-computed_soln[i])*W;
    }

   //Output x,y,...,u_exact,...]
   for(unsigned i=0;i<cached_dim;i++)
    {
     outfile << x[i] << " ";
    }

   //Output x,y,[z],u_error,v_error,[w_error], [p_error]
   for(unsigned i=0;i<hi_limit;i++)
    {
     outfile << exact_soln[i]-computed_soln[i] << " ";
    }
   outfile << std::endl;
  }
}


 ///  Overloaded version of the interpolated velocity solution including
 /// the singular contributions
 inline Vector<double> interpolated_u_nst(const Vector<double> &s) 
  const
 {
  // Find the dimension of the problem
  unsigned cached_dim=this->dim();
  
  //Find number of nodes
  const unsigned n_node = this->nnode();
  
  //Initialise value of u
  Vector<double> interpolated_u(cached_dim,0.0);
  
  // Calculate the global coordinate associated with s
  Vector<double> x(cached_dim,0.0);
  
  //Local shape function
  Shape psif(n_node);
  
  //Find values of shape function
  this->shape(s,psif);
  
  //Loop over the spatial directions
  for (unsigned d=0;d<cached_dim;d++)
   {
    //Get the index at which the unknown is stored
    const unsigned u_nodal_index = this->u_index_nst(d);
    
    //Loop over the local nodes and sum
    for(unsigned j=0;j<n_node;j++) 
     {
      interpolated_u[d] += this->nodal_value(j,u_nodal_index)*psif[j];
      x[d]+=this->raw_nodal_position(j,d)*psif(j);
     }
   }
 
  // Add the contribution of the singularities (all of them, summed)
  Vector<double> u_bar_local=u_bar(x);
  for (unsigned d=0;d<cached_dim;d++)
   {
    interpolated_u[d] += u_bar_local[d]; 
   }
  return(interpolated_u);
 }
 
 /// Version of interpolated pressure including the singular contributions
 inline double interpolated_p_nst(const Vector<double>& s) const
 {
  // Initialise pressure value with fe part
  double interpolated_p = interpolated_p_nst_fe_only(s);

  // Calculate the global coordinate associated with s
  unsigned cached_dim = this->dim();
  Vector<double> x(cached_dim,0.0);

  //Find number of nodes
  const unsigned n_node = this->nnode();
  
  //Local shape function
  Shape psif(n_node);
  
  //Find values of shape function
  this->shape(s,psif);
  
  //Loop over the local nodes and sum
  for(unsigned j=0;j<n_node;j++) 
   {
    for (unsigned d=0;d<cached_dim;d++)
     {
      x[d]+=this->raw_nodal_position(j,d)*psif(j);
     }
   }

  // Add the singularities contribution (all of them, summed)
  interpolated_p+=p_bar(x);
  
  return interpolated_p;
}


  private:

 ///  Evaluate i-th u_bar (i-th velocity singular function 
 /// incl. the amplitudes) function at Eulerian position x
 Vector<double> u_bar(const unsigned& i,const Vector<double>& x) const
  {
   return C_equation_elements_pt[i]->u_bar(x);
  }
 
 ///  Evaluate sum of all velocity singular fcts
 /// (incl. the amplitude) at Eulerian position x
 Vector<double> u_bar(const Vector<double>& x) const
  {
   // Find the number of singularities
   unsigned n_sing = C_equation_elements_pt.size();
   
   // Find the dimension of the problem
   unsigned cached_dim=this->dim();  
   Vector<double> sum(cached_dim,0.0);
   for (unsigned s=0;s<n_sing;s++)
    {
     Vector<double> u_bar_local=
      C_equation_elements_pt[s]->u_bar(x);
     for (unsigned i=0;i<cached_dim;i++)
      {
       sum[i]+=u_bar_local[i];
      }
    }
   return sum;
 }

 ///  Evaluate i-th grad_u_bar (i-th gradient of velocity singular 
 /// fct incl. the amplitude) function at Eulerian position x;
 /// grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_u_bar(const unsigned& i,
                                    const Vector<double>& x) const
  {
   return C_equation_elements_pt[i]->grad_u_bar(x);
  }
 
 ///  Evaluate gradient of sum of all velocity singular fcts
 /// (incl. the amplitudes) at Eulerian position x: grad[i][j] = du_i/dx_j
 Vector<Vector<double> > grad_u_bar(const Vector<double>& x) const
  {
   // Find the number of singularities
   unsigned n_sing = C_equation_elements_pt.size();
   
   // Find the dimension of the problem
   unsigned cached_dim=this->dim();  
   Vector<Vector<double> > sum(cached_dim);
   for (unsigned i=0;i<cached_dim;i++)
    {
     sum[i].resize(cached_dim,0.0);
    }
   for (unsigned s=0;s<n_sing;s++)
    {
     Vector<Vector<double> > grad_u_bar_local=
      C_equation_elements_pt[s]->grad_u_bar(x);
     for (unsigned i=0;i<cached_dim;i++)
      {
       for (unsigned j=0;j<cached_dim;j++)
        {
         sum[i][j]+=grad_u_bar_local[i][j];
        }
      }
    }
   return sum;
 }


 /// Evaluate i-th "raw" velocity singular function at Eulerian position x
 Vector<double> velocity_singular_function(const unsigned& i, 
                                           const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->velocity_singular_function(x);
 }


 ///  Evaluate gradient of i-th "raw" velocity singular function at 
 /// Eulerian position x
 Vector<Vector<double> > grad_velocity_singular_function(const unsigned& i, 
                                           const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->grad_velocity_singular_function(x);
 }

 ///  Evaluate i-th pressure singular fct (without the 
 /// amplitude) at Eulerian position x
 double pressure_singular_function(const unsigned& i,
                                   const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->pressure_singular_function(x);
 }


 ///  Evaluate i-th p_bar (i-th pressure singular fct (incl. the 
 /// amplitude) at Eulerian position x
 double p_bar(const unsigned& i,const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->p_bar(x);
 }
 
 ///  Evaluate sum of all pressure singular fcts
 /// (incl. the amplitudes) at Eulerian position x
 double p_bar(const Vector<double>& x) const
 {
  // Find the number of singularities
  unsigned n_sing = C_equation_elements_pt.size();
  
  double sum=0.0;
  for (unsigned i=0;i<n_sing;i++)
   {
    sum+=C_equation_elements_pt[i]->p_bar(x);
   }
  return sum;
 }

 ///  Return FE representation of velocity solution WITHOUT singular
 /// contributions at local coordinate s
 Vector<double> interpolated_u_nst_fe_only
  (const Vector<double> &s) const
  {
   // Find number of nodes
   const unsigned n_node = this->nnode();
   
   // Find the dimension of the problem
   unsigned cached_dim = this->dim();
   
   //Initialise value of u
   Vector<double> interpolated_u(cached_dim,0.0);
   
   //Local shape function
   Shape psif(n_node);
   
   //Find values of shape function
   this->shape(s,psif);
   
   //Loop over the velocity components
   for (unsigned d=0;d<cached_dim;d++)
    {
     //Get the index at which the unknown is stored
     const unsigned u_nodal_index = this->u_index_nst(d);
     
     // Add the contribution of u_FE
     //Loop over the local nodes and sum
     for(unsigned j=0;j<n_node;j++) 
      {
       interpolated_u[d] += this->nodal_value(j,u_nodal_index)*psif[j];
      }
    }
   
   return interpolated_u;
  }
 

 ///  Return FE representation of pressure solution WITHOUT singular
 /// contributions at local coordinate s
 double interpolated_p_nst_fe_only(const Vector<double>& s) const
 {
  // Find number of pressure degrees of freedom
  unsigned n_pres = this->npres_nst();

  // Set up memory for pressure shape and test functions
  Shape psip(n_pres), testp(n_pres);

  // Compute the values of the pressure shape and test functions at the
  // local coordinate s
  this->pshape_nst(s,psip,testp);

  // Initialise pressure value
  double interpolated_p = 0.0;

  // Add the contribution of p_FE
  // Loop over the pressure dof and sum
  for (unsigned l=0;l<n_pres;l++)
   {
    interpolated_p += this->p_nst(l)*psip[l];
   }

  return interpolated_p;
 }


 /// Overloaded fill-in function 
void fill_in_generic_residual_contribution_wrapped_nst(
 Vector<double> &residuals, 
 DenseMatrix<double> &jacobian,
 const unsigned& flag)
{
 // Get the contribution from the underlying wrapped element first
 BASIC_NAVIER_STOKES_ELEMENT::fill_in_generic_residual_contribution_nst
  (residuals,jacobian,GeneralisedElement::Dummy_matrix,flag);
 
 // Find the dimension of the problem
 unsigned cached_dim = this->dim();
 
  // Do all the singular functions satisfy the Stokes eqn?
 bool all_singular_functions_satisfy_stokes_equation=true;
 
 // Find the number of singularities
 unsigned n_sing = C_equation_elements_pt.size();
 
 // Find the local equation number of the additional unknowns
 Vector<int> local_equation_number_C(n_sing);
 for (unsigned i=0;i<n_sing;i++)
  {
   local_equation_number_C[i] = this->external_local_eqn(i,0);
   if (!(C_equation_elements_pt[i]->
         singular_function_satisfies_stokes_equation()))
    {
     all_singular_functions_satisfy_stokes_equation=false;
    }
  }

 // Find out how many nodes there are
 unsigned n_node = this->nnode();
 
 // Find out how many pressure dofs there are
 unsigned n_pres = this->npres_nst();
 
 // Find the indices at which the local velocities are stored
 Vector<unsigned> u_nodal_index(cached_dim);
 for (unsigned d=0;d<cached_dim;d++)
  {
   u_nodal_index[d] = this->u_index_nst(d);
  }
 
 // integer to store the local equations
 int local_eqn=0, local_unknown=0;

 // Check that there's no time-dependence
 for (unsigned j=0;j<n_node;j++)
  {
   if (!(this->node_pt(j)->time_stepper_pt()->is_steady()))
    {
     std::stringstream error_stream;
     error_stream 
      << "Currently, the NavierStokesElementWithSingularity elements\n"
      << "only work for steady problems, but we have detected a\n"
      << "non-steady time-stepper for node " << j << ".\n"
      << "If your problem is time-dependent you're welcome to \n"
      << "volunteer and implement the required functionality in \n\n"
      << "   NavierStokesElementWithSingularity::fill_in_generic_residual_contribution_wrapped_nst()\n\n"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
  }



 // Check that ALE is disabled
 if (!this->ALE_is_disabled)
  {
   std::stringstream error_stream;
   error_stream 
    << "Currently, the NavierStokesElementWithSingularity elements\n"
    << "only work on fixed meshes, and to check this we require\n"
    << "that you assert that this is the case by calling \n\n"
    << "   NavierStokesElementWithSingularity::disable_ALE()"
    << "\n\n\n"
    << "for all bulk Navier-Stokes elements.\n"
    << "If your mesh is moving you're welcome to volunteer and implement\n"
    << "the required functionality in \n\n"
    << "   NavierStokesElementWithSingularity::fill_in_generic_residual_contribution_wrapped_nst()\n\n"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
  }


 // Throw an error for now...
 if ((!all_singular_functions_satisfy_stokes_equation)&&(flag==1))
  {
     std::stringstream error_stream;
     error_stream
      << "Currently, the analytical computation of the Jacobian for the\n"
      << "NavierStokesElementWithSingularity elements\n"
      << "only work with singular solutions that satisfy the Stokes eqns.\n"
      << "Of course, there's no way of checking this, so you'll have to\n"
      << "assert that this is the case by calling \n\n"
      << "  SingularNavierStokesSolutionElement::singular_function_satisfies_stokes_equation()"
      << "\n\n\n"
      << "for the SingularNavierStokesSolutionElement that specifies the singular functions.\n"
      << "If your singular functions do not satisfy the Stokes equations\n"
      << "you're welcome to volunteer and implement the required \n"
      << "functionality in \n\n"
      << "   NavierStokesElementWithSingularity::fill_in_generic_residual_contribution_wrapped_nst()\n\n"
      << "Fragments of code are already there...\n"
      << "Alternatively, use FD-based setup of Jacobian (still there, commented out, in code).\n"
      << "The reason this never got implemented is because this case typically only arises\n"
      << "when trying to \"blend\" solutions (to give them finite support) and this didn't\n"
      << "seem to work too well and was sort of abandoned. Next person to look at this\n"
      << "should consider only adding the non-integrated by parts (source-fct-like) terms\n"
      << "arising from the singular functions...\n"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
  }

 // Get Physical Variables from Element
 // Reynolds number must be multiplied by the density ratio
 double scaled_re = this->re()*this->density_ratio();
 
 // Do we need extra bulk terms?
 if ((scaled_re>0.0)||(!all_singular_functions_satisfy_stokes_equation))
  {
   
   // Set up memory for the velocity shape and test functions
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,cached_dim), dtestfdx(n_node,cached_dim);
   
   // Set up memory for pressure shape and test functions
   Shape psip(n_pres), testp(n_pres);
   
   // Number of integration points
   unsigned n_intpt = this->integral_pt()->nweight();
   
   // Set the vector to hold local coordinates
   Vector<double> s(cached_dim);
   
   // Cachec viscosity ratio
   double visc_ratio = this->viscosity_ratio();
   
   // Loop over the integration points
   for (unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     // Assign values of s
     for (unsigned d=0;d<cached_dim;d++)
      {
       s[d] = this->integral_pt()->knot(ipt,d);
      }
     
     // Get the integral weight
     double w = this->integral_pt()->weight(ipt);
     
     // Call the derivatives of the velocity shape and test functions
     double J =
      this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,
                                                  testf,dtestfdx);
     
     // Call the pressure shape and test functions
     this->pshape_nst(s,psip,testp);
     
     // Premultiply the weights and the Jacobian
     double W = w*J;
     
     // Initialise the global coordinate and velocity
     Vector<double> interpolated_x(cached_dim,0.0);
     Vector<double> interpolated_u(cached_dim,0.0);
     DenseMatrix<double> interpolated_dudx(cached_dim,cached_dim,0.0);    

     // Calculate the global coordinate associated with s
     // Loop over nodes
     for (unsigned l=0;l<n_node;l++)
      {
       // Loop over directions
       for (unsigned i=0;i<cached_dim;i++)
        {
         double u_value = this->raw_nodal_value(l,u_nodal_index[i]);
         interpolated_u[i] += u_value*psif[l];
         interpolated_x[i] += this->raw_nodal_position(l,i)*psif(l);
         for(unsigned j=0;j<cached_dim;j++)
          {                               
           interpolated_dudx(i,j) += u_value*dpsifdx(l,j);
          }
        }
      }
     
     // Get sum of singular functions
     Vector<double> u_bar_local=this->u_bar(interpolated_x);
     Vector<Vector<double> > grad_u_bar_local=this->grad_u_bar(interpolated_x);
     double p_bar_local=this->p_bar(interpolated_x);

     // Singular functions
     Vector<Vector<double> > u_hat_local(n_sing);
     for (unsigned s=0;s<n_sing;s++)
      {
       u_hat_local[s]=this->velocity_singular_function(s,interpolated_x);
      }
     Vector<Vector<Vector<double> > > grad_u_hat_local(n_sing);
     for (unsigned s=0;s<n_sing;s++)
      {
       grad_u_hat_local[s]=
        this->grad_velocity_singular_function(s,interpolated_x);
      }

     // MOMENTUM EQUATIONS
     //-------------------
     
     // Loop over the velocity test functions
     for (unsigned l=0;l<n_node;l++)
      {
       // Loop over the velocity components
       for (unsigned i=0;i<cached_dim;i++)
        {
         // Find its local equation number
         local_eqn = this->nodal_local_eqn(l,u_nodal_index[i]);
         
         // If it is not pinned
         if (local_eqn >= 0)
          {
           // If it is not a Dirichlet BC
           if (not(Node_is_subject_to_velocity_dirichlet_bcs[l][i]))
            {


             // Linear terms only needed if singular solution doesn't satisfy
             //--------------------------------------------------------------
             // Stokes eqn
             //-----------
             if (!all_singular_functions_satisfy_stokes_equation)
              {
               residuals[local_eqn]+=p_bar_local*dtestfdx(l,i)*W;
              for (unsigned k=0;k<cached_dim;k++)
               {
                residuals[local_eqn] -= visc_ratio*
                 (grad_u_bar_local[i][k]+this->Gamma[i]*grad_u_bar_local[k][i])
                 *dtestfdx(l,k)*W;
               }
              }
             
             
             // Nonlinear term. Always add (unless Re=0)
             //-----------------------------------------
             if (scaled_re>0.0)
              {
               // Add singular convective terms
               double sum = 0.0;
               for (unsigned k=0;k<cached_dim;k++)
                {
                 sum+=u_bar_local[k]*(grad_u_bar_local[i][k]+interpolated_dudx(i,k))
                  +interpolated_u[k]*grad_u_bar_local[i][k];
                }
               residuals[local_eqn] -= scaled_re*sum*testf[l]*W;
               
              }

             // Calculate the jacobian
             //-----------------------
             if(flag)
              {

               if (scaled_re>0.0)
                {
                 // Loop over the singularities and add the 
                 // contributions of the additional 
                 // unknowns associated with them to the 
                 // jacobian if they are not pinned
                 for (unsigned ss=0;ss<n_sing;ss++)
                  {
                   local_unknown=local_equation_number_C[ss];
                   if (local_unknown>=0)
                    {
                     double sum=0.0;
                     for (unsigned k=0;k<cached_dim;k++)
                      {
                       sum+=u_hat_local[ss][k]*(grad_u_bar_local[i][k]+
                                                interpolated_dudx(i,k))+
                        grad_u_hat_local[ss][i][k]*
                        (u_bar_local[k]+interpolated_u[k]);
                      }
                     jacobian(local_eqn,local_unknown)-=scaled_re*sum*testf[l]*W;
                    }
                  }
                 
                 
                 //Loop over the velocity shape functions again
                 for(unsigned l2=0;l2<n_node;l2++)
                  { 
                   //Loop over the velocity components again
                   for(unsigned i2=0;i2<cached_dim;i2++)
                    {
                     //If at a non-zero degree of freedom add in the entry
                     local_unknown = this->nodal_local_eqn(l2,u_nodal_index[i2]);
                     if(local_unknown >= 0)
                      {
                       double sum=0.0;
                       if (i==i2)
                        {
                         for (unsigned k=0;k<cached_dim;k++)
                          {
                           sum+=u_bar_local[k]*dpsifdx(l2,k);
                          }
                        }
                       sum+=psif(l2)*grad_u_bar_local[i][i2];
                       
                       //Add contribution to Elemental Matrix
                       jacobian(local_eqn,local_unknown) -=
                        scaled_re*sum*testf[l]*W;
                      }
                    }
                  }
                }
              }
            } // End of check of the Dirichlet status
          } // End of check of the pin status
        } // End of loop over velocity components
      } // End of loop over test functions
     

     
     // Linear terms only needed if singular solution doesn't satisfy
     //--------------------------------------------------------------
     // Stokes eqn
     //-----------
     if (!all_singular_functions_satisfy_stokes_equation)
     {
      // Loop over the pressure test functions
      for (unsigned l=0;l<n_pres;l++)
       {
        local_eqn = this->p_local_eqn(l);
        
        // If not pinned
        if (local_eqn >= 0)
         {
          // If not subject to Dirichlet BC
          if (not(Pressure_dof_is_subject_to_dirichlet_bc[l]))
           {
            double aux = 0.0;
            // Loop over velocity components
            for (unsigned k=0;k<cached_dim;k++)
             {
              aux += grad_u_bar_local[k][k];
             }
            residuals[local_eqn] += aux*testp[l]*W;
           }
         }
       }
     }
     
    } // End of loop over integration points
   
   
  }

  // VELOCITY DIRICHLET BCS
  //-----------------------
  Vector<double> u_bar_at_node(cached_dim);
  Vector<Vector<double> > u_hat_at_node(n_sing);
  for (unsigned i=0;i<n_sing;i++)
   {
    u_hat_at_node[i].resize(cached_dim);
   }

  // Loop over the nodes
  for (unsigned l=0;l<n_node;l++)
   {
    // Find the global coordinate of the node
    Vector<double> global_coordinate(cached_dim);
    for (unsigned d=0;d<cached_dim;d++)
     {
      global_coordinate[d] = this->raw_nodal_position(l,d);
     }
    
    // Get singular velocity at node
    u_bar_at_node=this->u_bar(global_coordinate); 
    for (unsigned i=0;i<n_sing;i++)
     {
      u_hat_at_node[i]=velocity_singular_function(i,global_coordinate);
     }

    // Loop over the velocity components
    for (unsigned d=0;d<cached_dim;d++)
     {
      // Find its local equation number
      local_eqn = this->nodal_local_eqn(l,u_nodal_index[d]);
      
      // If it is not pinned
      if (local_eqn >= 0)
       {
        // If it is a Dirichlet boundary condition
        if (Node_is_subject_to_velocity_dirichlet_bcs[l][d])
         {
          // Initialise the residual
          residuals[local_eqn] = 0.0;
          
          // Add the contribution of the nodal value
          residuals[local_eqn] += this->raw_nodal_value(l,u_nodal_index[d]);
          
          // Add the contribution of the singularities (all of them, summed)
          residuals[local_eqn] += u_bar_at_node[d];

          // Substract the imposed Dirichlet value
          residuals[local_eqn] -= Imposed_velocity_values_at_node[l][d];

          if (flag)
           {

            // Wipe the existing entries
            unsigned n_dof=this->ndof();
            for (unsigned j=0;j<n_dof;j++)
             {
              jacobian(local_eqn,j)=0.0;
             }

            // Add diagonal entry
            jacobian(local_eqn,local_eqn) += 1.0;


            // Add derivative w.r.t. the Cs
            for (unsigned i=0;i<n_sing;i++)
             {
              // Find the contribution of the additional unknowns to 
              // the jacobian
              local_unknown=local_equation_number_C[i];
              if (local_unknown>=0)
               {
                jacobian(local_eqn,local_unknown) += 
                 u_hat_at_node[i][d]; 
               }
             }
             
           }
         }
       }
     }
   }


  // PRESSURE DIRICHLET BCS
  //-----------------------
  
  // Loop over the pressure dofs
  for (unsigned l=0;l<n_pres;l++)
   {
    // Find its local equation number
    local_eqn = this->p_local_eqn(l);

    // If it is not pinned
    if (local_eqn >= 0)
     {
      // If it is subject to a Dirichlet BC
      if (Pressure_dof_is_subject_to_dirichlet_bc[l])
       {
        // Find its global coordinate
        // This conversionly works for Taylor Hood type elements
        // but there's not much point assigning pressure dofs
        Node* p_nod_pt=this->node_pt(this->Pconv[l]);
        
        oomph_info << "Constrained pressure node: " 
                   << this->Pconv[l] << " at: "; 
        
        Vector<double> global_coordinate(cached_dim,0.0);
        for (unsigned d=0;d<cached_dim;d++)
         {
          global_coordinate[d] = p_nod_pt->x(d);
          oomph_info << global_coordinate[d] << " ";
         }
        oomph_info << std::endl;

        // Initialise its residual component
        residuals[local_eqn] = 0.0;

        // Add the contribution of the pressure unknown
        residuals[local_eqn] += this->p_nst(l);

        // Add singular contributions
        residuals[local_eqn] += p_bar(global_coordinate);

        // Substract the imposed pressure value
        residuals[local_eqn] -= Imposed_value_at_pressure_dof[l];

        if (flag)
         {
          
          // Wipe the existing entries
          unsigned n_dof=this->ndof();
          for (unsigned j=0;j<n_dof;j++)
           {
            jacobian(local_eqn,j)=0.0;
           }
          
          // Add diagonal entry
          jacobian(local_eqn,local_eqn) += 1.0;
          
          // Add derivative w.r.t. the Cs
          for (unsigned i=0;i<n_sing;i++)
           {
            // Find the contribution of the additional unknowns to 
            // the jacobian
            local_unknown=local_equation_number_C[i];
            if (local_unknown>=0)
             {
              jacobian(local_eqn,local_unknown) += 
               pressure_singular_function(i,global_coordinate); 
             }
           }
         }
       }
     }
   }
  
  
 } // End of function



///  Vector of pointers to SingularNavierStokesSolutionElement objects
Vector<SingularNavierStokesSolutionElement<NavierStokesElementWithSingularity
<BASIC_NAVIER_STOKES_ELEMENT> >*> C_equation_elements_pt;

///  Vector indicating which velocity component of
/// which node is subject to Dirichlet BC
/// [size = number of nodes; initialised to false]
Vector<std::vector<bool> > Node_is_subject_to_velocity_dirichlet_bcs;

///  Imposed values of velocity component at nodes
/// that are subject to Dirichlet BC
/// [size = number of nodes; initialised to zero]
Vector<Vector<double> > Imposed_velocity_values_at_node;

///  Vector indicating which pressure dof is subject to Dirichlet BC
/// [size = number of pressure dofs; initialised to false]
std::vector<bool> Pressure_dof_is_subject_to_dirichlet_bc;

///  Imposed value at pressure dofs
/// that are subject to Dirichlet BC
/// [size = number of pressure dof; initialised to zero]
Vector<double> Imposed_value_at_pressure_dof;

}; // End of NavierStokesElementWithSingularity class

}

#endif
