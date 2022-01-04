//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
//Header file for Poisson elements with singularity
#ifndef OOMPH_POISSON_ELEMENTS_WITH_SINGULARITY_HEADER
#define OOMPH_POISSON_ELEMENTS_WITH_SINGULARITY_HEADER


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
 
//==================CLASS FOR THE ADDITIONAL UNKNOWN==================
/// We consider a singularity in the solution at the point O.
///
/// This class defines the functions u_sing and 
/// u_bar = C*u_sing and their gradients.
///
/// The class also defines the function that computes the
/// residual associated with C. 
/// R_C = \frac{\partial u_FE}{\partial x_d}(O)
/// where u_FE is : \sum_{i} U_i \Psi_i.
/// This sets the derivative of the finite-element part of the
/// solution to zero in the specified direction (d) and thus
/// imposes regularity on that part. 
//=====================================================================
template<class WRAPPED_POISSON_ELEMENT> 
class SingularPoissonSolutionElement : public virtual GeneralisedElement
{
  public :
 
 /// Function pointer to the singular function:
 typedef double (*PoissonSingularFctPt)(const Vector<double>& x); 
 
 /// Function pointer to the gradient of the singular function:
 typedef Vector<double> (*PoissonGradSingularFctPt)
  (const Vector<double>& x); 
 
 /// Constructor
 SingularPoissonSolutionElement()
  {
   // Initialise Function pointer to singular function 
   Singular_fct_pt=0;
   
   // Initialise Function pointer to gradient of the singular function
   Grad_singular_fct_pt=0;
   
   // Initalise pointer to the wrapped Poisson element which will be used to
   // compute the residual and which includes the point O
   Wrapped_poisson_el_pt=0;

   // Initialise the pointer to the direction of the derivative used
   // for the residual of this element
   Direction_pt=0;
   
   // Create a single item of internal Data, storing one unknown which
   // represents the unknown C.
   add_internal_data(new Data(1));
   
   // Safe assumption: Singular fct does not satisfy Laplace's eqn
   Singular_function_satisfies_laplace_equation=false;

  } // End of constructor


 /// Does the singular fct satisfy Laplace's eqn? Default assumption
 /// is that it doesn't; user can overwrite this here to avoid unnecessary
 /// computation of integrals in residual
 bool& singular_function_satisfies_laplace_equation()
  {
   return Singular_function_satisfies_laplace_equation;
  }
 
 /// Find the value of the unknown C
 double c() const 
 {
  return internal_data_pt(0)->value(0);
 }
 
 /// Set pointer to associated wrapped Poisson element which contains
 /// the singularity (at local coordinate s). Also specify the direction
 /// in which the slope of the FE part of the solution is set to zero 
 void set_wrapped_poisson_element_pt(WRAPPED_POISSON_ELEMENT* 
                                     wrapped_poisson_el_pt,
                                     const Vector<double>& s,
                                     unsigned* direction_pt)
 {
  // Assign the pointer to the variable Wrapped_poisson_el_pt
  Wrapped_poisson_el_pt = wrapped_poisson_el_pt;
  
  // Find number of nodes in the element
  unsigned nnod=wrapped_poisson_el_pt->nnode();
  
  // Loop over the nodes of the element
  for (unsigned j=0;j<nnod;j++)
   {
    // Add the node as external data in the SingularPoissonSolutionElement class
    // because they affect the slope that we set to zero to
    // determine the value of the amplitude C
    add_external_data(Wrapped_poisson_el_pt->node_pt(j));
   }
  
  // Assign the pointer to the local coordinate at which the residual
  // will be computed
  S_in_wrapped_poisson_element = s;
  
  // Assign the pointer to the direction at which the derivative used
  // in the residual will be computed
  Direction_pt = direction_pt;
 }
 
 
 /// Access function to pointer to singular function
 PoissonSingularFctPt& singular_fct_pt() {return Singular_fct_pt;}
 
 /// Evaluate singular function at Eulerian position x
 double singular_function(const Vector<double>& x) const
 {
#ifdef PARANOID
  if (Singular_fct_pt==0)
   {
    std::stringstream error_stream;
    error_stream 
     << "Pointer to singular function hasn't been defined!"
     << std::endl;
    throw OomphLibError(
     error_stream.str(),
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Evaluate singular function
  return (*Singular_fct_pt)(x);
 }
 
 /// The singular function (including its amplitude): ubar = C * sing
 /// where sing is defined via the function pointer.
 double u_bar(const Vector<double>& x)
 {
  // Find the value of C
  double c=internal_data_pt(0)->value(0);
  
  // Value of ubar at the position x
  return c*singular_function(x);
 } // End of function
 
 /// Access function to pointer to the gradient of sing function
 PoissonGradSingularFctPt& grad_singular_fct_pt()
  {return Grad_singular_fct_pt;}
 
 /// Evaluate the gradient of the sing function at Eulerian position x
 Vector<double> grad_singular_function(const Vector<double>& x) const
  {
   // Find the dimension of the problem
   unsigned cached_dim=Wrapped_poisson_el_pt->dim();
   
   // Declare the gradient of sing
   Vector<double> dsingular(cached_dim);
   
#ifdef PARANOID
   if (Grad_singular_fct_pt==0)
    {
     std::stringstream error_stream;
     error_stream 
      << "Pointer to gradient of singular function hasn't been defined!"
      << std::endl;
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
#endif
   
   // Evaluate gradient of the sing function
   return (*Grad_singular_fct_pt)(x);
  }
 
 /// Gradient of ubar (including the constant!)
 Vector<double> grad_u_bar(const Vector<double>& x)
  {
   // Compute the gradient of sing
   Vector<double> dsingular = grad_singular_function(x);
   
   // Initialise the gradient of ubar
   unsigned cached_dim=Wrapped_poisson_el_pt->dim();
   Vector<double> dubar(cached_dim);
   
   // Find the value of C
   double c=internal_data_pt(0)->value(0);
   for (unsigned i=0;i<cached_dim;i++)
    {
     dubar[i] = c*dsingular[i];
    }
   
   return dubar;
  } // End of function
 
 
 // Compute local residual
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
  
  // Get the derivative
  unsigned cached_dim=Wrapped_poisson_el_pt->dim();
  Vector<double> flux(cached_dim);
  Wrapped_poisson_el_pt->get_flux(S_in_wrapped_poisson_element,flux);
  double derivative = flux[*Direction_pt];
  
  // fill in the contribution to the residual
  residual[eqn_number] = derivative;
  
  if (flag)
   {
    // Find the number of nodes in the Poisson element associated to the
    // SingularPoissonSolutionElement class
    unsigned n_node = Wrapped_poisson_el_pt->nnode();
    
    // Compute the derivatives of the shape functions of the poisson
    // element associated to the SingularPoissonSolutionElement at the 
    // local coordinate S_in_wrapped_poisson_element_pt
    Shape psi(n_node);
    DShape dpsidx(n_node,cached_dim);
    Wrapped_poisson_el_pt->dshape_eulerian(
     S_in_wrapped_poisson_element,psi,dpsidx);
    
    // Loop over the nodes
    for (unsigned j=0;j<n_node;j++)
     {
      // Find the local equation number of the node
      int node_eqn_number = external_local_eqn(j,0);
      if (node_eqn_number>=0)
       {
       // Add the contribution of the node to the local jacobian
       jacobian(eqn_number,node_eqn_number) = dpsidx(j,*Direction_pt);
       }

     }
   }
  
 }
 
 /// Pointer to Poisson element that contains the singularity
 WRAPPED_POISSON_ELEMENT* Wrapped_poisson_el_pt; 
 
 /// Pointer to singular function
 PoissonSingularFctPt Singular_fct_pt;
 
 /// Pointer to the gradient of the singular function 
 PoissonGradSingularFctPt Grad_singular_fct_pt;
 
 /// Local coordinates of singularity in the wrapped Poisson element
 Vector<double> S_in_wrapped_poisson_element;
 
 /// Direction of the derivative used for the residual of the element
 unsigned* Direction_pt;
 
 /// Does singular fct satisfy Laplace's eqn?
 bool Singular_function_satisfies_laplace_equation;

}; // End of SingularPoissonSolutionElement Class


/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////


//=======================================================================
/// Templated wrapper to add handling of singularities to
/// the underlying Poisson element (specified via the template
/// argument). Slightly inefficient because the integration loop
/// is repeated here. The only alternative is to add the additional
/// functionality into the Poisson elements. Not pretty either!
///
/// Dirichlet BCs are applied in hijacking manner and must be imposed
/// from each element that shares a boundary node that is subject to
/// a Dirichlet condition.
//=======================================================================
template<class BASIC_POISSON_ELEMENT>
class PoissonElementWithSingularity : public virtual BASIC_POISSON_ELEMENT
{

  public:

 /// Constructor
 PoissonElementWithSingularity()
  {
   // Find the number of nodes in the element
   unsigned n_node = this->nnode();
   
   // Initialise the vector of booleans indicating which
   // node is subject to Dirichlet BC. The size of the
   // vector is equal to the number
   // of nodes in the element. By default, no node is subject
   // to Dirichlet BC, so the vector is full of false. 
   Node_is_subject_to_dirichlet_bc.resize(n_node);
   for (unsigned j=0;j<n_node;j++)
    {
     Node_is_subject_to_dirichlet_bc[j] = false;
    }
   
   // Initialise the vector of imposed values on the nodes subject
   // to Dirichlet BC. The size of the vector is equal to the number
   // of nodes in the element. If a node is not subject to Dirichlet
   // BC, its imposed value is zero. By default, no node is subject
   // to Dirichlet BC so the vector is full of zeros
   Imposed_value_at_node.resize(n_node);
   for (unsigned j=0;j<n_node;j++)
    {
     Imposed_value_at_node[j] = 0.0;
    }
   
  } // End of constructor
 
 
 /// Impose Dirichlet BC on sum of FE solution and u_bar on j-th local
 /// node
 void impose_dirichlet_bc_on_node(const unsigned& j)
 {
  Node_is_subject_to_dirichlet_bc[j] = true;
 }
 
 /// Undo Dirichlet BCs on jth node
 void undo_dirichlet_bc_on_node(const unsigned& j)
 {
  Node_is_subject_to_dirichlet_bc[j] = false;
 }
 
 /// Undo Dirichlet BCs on all nodes
 void undo_dirichlet_bc_on_all_nodes()
 {
  unsigned n_node = this->nnode();
  for (unsigned j=0;j<n_node;j++)
   {
    undo_dirichlet_bc_on_node(j);
   }
 }
 
 /// Specify Dirichlet boundary value for j-th local node
 void set_dirichlet_value_on_node(const unsigned& j,
                                  const double& value)
 {
  Imposed_value_at_node[j] = value;
 }
 
 
 /// Access function to pointer of vector of 
 /// SingularPoissonSolutionElements. These
 /// specify the singular functions that are added to the (regular)
 /// FE solution
 Vector<SingularPoissonSolutionElement<PoissonElementWithSingularity
  <BASIC_POISSON_ELEMENT> >*> c_equation_elements_pt()
  {
   return C_equation_elements_pt;
  }
 

 /// Add pointer to associated SingularPoissonSolutionElement that 
 /// determines the value of the amplitude of the singular function (and gives 
 /// access to the singular function). The unknown amplitude becomes
 /// external Data for this element so assign_eqn_numbers() must be
 /// called after this function has been called. 
 void add_c_equation_element_pt(
  SingularPoissonSolutionElement<PoissonElementWithSingularity
  <BASIC_POISSON_ELEMENT> >* c_pt)
 {
  // Add the element
  C_equation_elements_pt.push_back(c_pt);
  
  // Add the additional unknown of this object as external data in the
  // Poisson element
  this->add_external_data(c_pt->internal_data_pt(0));
 }
 
 /// Evaluate i-th u_bar function (i.e. i-th singular incl. the 
 /// amplitude) function at Eulerian position x
 double u_bar(const unsigned& i,const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->u_bar(x);
 }
 
 /// Evaluate i-th "raw" singular function at Eulerian position x
 double singular_function(const unsigned& i, const Vector<double>& x) const
 {
  return C_equation_elements_pt[i]->singular_function(x);
 }
 
 /// Evaluate the gradient the i-th singular (incl. the amplitude )
 /// at Eulerian position x
 Vector<double> grad_u_bar(const unsigned i,const Vector<double>& x) const
  {
   return C_equation_elements_pt[i]->grad_u_bar(x);
  }
 
 /// Evaluate the gradient of the i-th "raw" singular  at Eulerian position x
 Vector<double> grad_singular_function(const unsigned& i,
                                          const Vector<double>& x) const
  {
   return C_equation_elements_pt[i]->grad_singular_function(x);
  }
 
 
 /// Return FE representation of solution WITHOUT singular
 /// contributions at local coordinate s
 double interpolated_u_poisson_fe_only(const Vector<double> &s) const
 {
  //Find number of nodes
  const unsigned n_node = this->nnode();
  
  //Get the index at which the poisson unknown is stored
  const unsigned u_nodal_index = this->u_index_poisson();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  this->shape(s,psi);
  
  //Initialise value of u
  double interpolated_u = 0.0;
  
  // Add the contribution of u_FE
  //Loop over the local nodes and sum
  for(unsigned l=0;l<n_node;l++) 
   {
    interpolated_u += this->nodal_value(l,u_nodal_index)*psi[l];
   }
  
  return interpolated_u;
 }
 
 /// Overloaded version including the singular contributions
 inline double interpolated_u_poisson(const Vector<double> &s) const
 {
  //Find number of nodes
  const unsigned n_node = this->nnode();
  
  //Get the index at which the poisson unknown is stored
  const unsigned u_nodal_index = this->u_index_poisson();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  this->shape(s,psi);
  
  //Initialise value of u
  double interpolated_u = 0.0;
  
  // Calculate the global coordinate associated to s
  unsigned cached_dim=this->dim();
  Vector<double> x(cached_dim,0.0);

  // Add the contribution of u_FE
  //Loop over the local nodes and sum
  for(unsigned l=0;l<n_node;l++) 
   {
    interpolated_u += this->nodal_value(l,u_nodal_index)*psi[l];
    for (unsigned i=0;i<cached_dim;i++)
     {
      x[i]+=this->raw_nodal_position(l,i)*psi(l);
     }
   }
    
  // Loop over the singularities
  unsigned n_sing=C_equation_elements_pt.size();
  for (unsigned i=0;i<n_sing;i++)
   {
    // Add the contribution of ubari
    interpolated_u += u_bar(i,x);
   }
  
  return(interpolated_u);
 }
 
 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
 {
  //Call the generic residuals function with flag set to 0
  //using a dummy matrix argument
  this->fill_in_generic_residual_contribution_wrapped_poisson(
   residuals,GeneralisedElement::Dummy_matrix,0);
 }
 
 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
 {  
  //Call the generic routine with the flag set to 1
  fill_in_generic_residual_contribution_wrapped_poisson(residuals,jacobian,1);
 }


 /// Overloaded output fct: x, y [,z], u, u_fe_only, u-u_fe_only
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
    
    outfile << interpolated_u_poisson(s) << " "
            << interpolated_u_poisson_fe_only(s) << " "
            << interpolated_u_poisson(s)-interpolated_u_poisson_fe_only(s)
            << std::endl;
   
   }

  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);
 }

 
 /// Compute error
 void compute_error(std::ostream &outfile, 
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
 { 
 
  // Initialise
  error=0.0;
  norm=0.0;

  // Find the dimension of the problem
  unsigned cached_dim = this->dim();
 
  //Vector of local coordinates
  Vector<double> s(cached_dim);
 
  // Vector for coordintes
  Vector<double> x(cached_dim);
 
  //Find out how many nodes there are in the element
  unsigned n_node = this->nnode();
 
  Shape psi(n_node);
 
  //Set the value of n_intpt
  unsigned n_intpt = this->integral_pt()->nweight();
  
  // Tecplot 
  outfile << "ZONE" << std::endl;
 
  // Exact solution Vector (here a scalar)
  Vector<double> exact_soln(1);
 
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
   
    // Get FE function value
    double u_fe=interpolated_u_poisson(s);
   
    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);
   
    //Output x,y,...,error
    for(unsigned i=0;i<cached_dim;i++)
     {
      outfile << x[i] << " ";
     }
    outfile << exact_soln[0] << " " << exact_soln[0]-u_fe << std::endl;  
   
    // Add to error and norm
    norm+=exact_soln[0]*exact_soln[0]*W;
    error+=(exact_soln[0]-u_fe)*(exact_soln[0]-u_fe)*W;
   }
 }

 
  private:
 
 /// Overloaded fill-in function 
 void fill_in_generic_residual_contribution_wrapped_poisson(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian, 
  const unsigned& flag)
 {
  // Get the contribution from the underlying wrapped element first
  BASIC_POISSON_ELEMENT::fill_in_generic_residual_contribution_poisson
   (residuals,jacobian,flag);

  // Do all the singular functions satisfy laplace's eqn?
  bool all_singular_functions_satisfy_laplace_equation=true;
  
  // Find the number of singularities
  unsigned n_sing = C_equation_elements_pt.size();
  
  //Index at which the poisson unknown is stored
  const unsigned u_nodal_index = this->u_index_poisson();

  // Find the dimension of the problem
  unsigned cached_dim=this->dim();
  
  //Find out how many nodes there are
  const unsigned n_node = this->nnode();
  
  // Find the local equation number of the additional unknowns
  Vector<int> local_equation_number_C(n_sing);
  for (unsigned i=0;i<n_sing;i++)
   {
    local_equation_number_C[i] =
     this->external_local_eqn(i,0);
    if (!(C_equation_elements_pt[i]->
          singular_function_satisfies_laplace_equation()))
     {
      all_singular_functions_satisfy_laplace_equation=false;
     }
   }

  // Do we need to add the singular function's contribution to the
  // residuals or do they vanish by themselves (mathematically)
  if (!all_singular_functions_satisfy_laplace_equation)
   {
    
    //Set up memory for the shape and test functions
    Shape psi(n_node), test(n_node);
    DShape dpsidx(n_node,cached_dim), dtestdx(n_node,cached_dim);
    
    //Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();
    
    //Integers to store the local equation and unknown numbers
    int local_eqn=0;
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     
     {
      //Get the integral weight
      double w = this->integral_pt()->weight(ipt);
      
      //Call the derivatives of the shape and test functions
      double J = this->dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                                 test,dtestdx);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      //Calculate the global coordinate of the integration point
      Vector<double> interpolated_x(cached_dim,0.0);
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
       {
        // Loop over directions
        for(unsigned j=0;j<cached_dim;j++)
         {
          interpolated_x[j] += this->raw_nodal_position(l,j)*psi(l);
         }
       }

      // Precompute singular fcts
      Vector<Vector<double> > grad_u_bar_local(n_sing);
      for (unsigned i=0;i<n_sing;i++)
       {
        grad_u_bar_local[i]=grad_u_bar(i,interpolated_x);
       }
      Vector<Vector<double> > grad_singular_function_local(n_sing);
      for (unsigned i=0;i<n_sing;i++)
       {
        grad_singular_function_local[i]=grad_singular_function(i,interpolated_x);
       }
      
      // Assemble residuals and Jacobian
      //--------------------------------
      
      // The interior nodes
      //-------------------
      
      // Loop over the test functions
      for(unsigned l=0;l<n_node;l++)
       {
        //Get the local equation number
        local_eqn = this->nodal_local_eqn(l,u_nodal_index);
        
        // If it is not pinned
        if (local_eqn >= 0)
         {
          // IF it's not a boundary condition
          if(not(Node_is_subject_to_dirichlet_bc[l]))
           {
            // Add the contribution of the additional unknowns to the residual
            for(unsigned i=0;i<n_sing;i++)
             {
              if (!(C_equation_elements_pt[i]->
                    singular_function_satisfies_laplace_equation()))
               {
                for (unsigned k=0;k<cached_dim;k++)
                 {
                  residuals[local_eqn] += 
                   grad_u_bar_local[i][k]*dtestdx(l,k)*W;
                 }
               }
             }
            
            // Calculate the jacobian
            //-----------------------
            if(flag)
             {
              // Loop over the singularities and add the 
              // contributions of the additional 
              // unknowns associated with them to the 
              // jacobian if they are not pinned
              for (unsigned i=0;i<n_sing;i++)
               {
                if (!(C_equation_elements_pt[i]->
                      singular_function_satisfies_laplace_equation()))
                 {
                  if (local_equation_number_C[i]>=0)
                   {
                    for(unsigned d=0;d<cached_dim;d++)
                     {
                      jacobian(local_eqn,local_equation_number_C[i])
                       += grad_singular_function_local[i][d]*
                       dtestdx(l,d)*W;
                     }
                   }
                 }
               }
             }
           }
         }
       } // end of loop over shape functions
      
     } // End of loop over integration points
    
   }
  
  // The Dirichlet nodes
  //--------------------
  
  // Loop over the nodes to see if there is a node subject to
  // Dirichlet BC
  for (unsigned j=0;j<n_node;j++)
   {
    // If it is a dirichlet boundary condition
    if (Node_is_subject_to_dirichlet_bc[j])
     {
      // Find the global coordinate of the node
      Vector<double> global_coordinate_boundary_node(cached_dim);
      for (unsigned d=0;d<cached_dim;d++)
       {
        global_coordinate_boundary_node[d] = this->raw_nodal_position(j,d);
       }
      
      // If it is not pinned
      int local_eqn_number_boundary_node = 
       this->nodal_local_eqn(j,u_nodal_index);
      if (local_eqn_number_boundary_node >= 0)
       {
        // Add the contribution of the node unknown to the residual
        residuals[local_eqn_number_boundary_node] = 
        this->raw_nodal_value(j,u_nodal_index);
        
        // Add the contribution of the additional unknowns to the residual
        for (unsigned i=0;i<n_sing;i++)
         {
          residuals[local_eqn_number_boundary_node] += 
           u_bar(i,global_coordinate_boundary_node);
         }
        
        // Substract the value imposed by the Dirichlet BC from the residual
        residuals[local_eqn_number_boundary_node] -= Imposed_value_at_node[j];
        
        if (flag)
         {
          // Loop over the nodes
          for (unsigned l=0;l<n_node;l++)
           {
            // Find the equation number of the node
            int local_unknown = this->nodal_local_eqn(l,u_nodal_index);

            // If it is not pinned
            if (local_unknown >=0)
             {
              // Put to 0 the corresponding jacobian component
              jacobian(local_eqn_number_boundary_node,local_unknown)=0.0;
             }
           }
            
          // Find the contribution of the node to the local jacobian
          jacobian(local_eqn_number_boundary_node,
                   local_eqn_number_boundary_node) += 1.0;

          // Loop over the singularities
          for (unsigned i=0;i<n_sing;i++)
           {
            // Find the contribution of the additional unknowns to the jacobian
            int local_unknown=local_equation_number_C[i];
            if (local_unknown>=0)
             {
              jacobian(local_eqn_number_boundary_node,
                       local_unknown) += 
               singular_function(i,global_coordinate_boundary_node);
             }
           }
             
         }
       } // End of check of the pin status of the node
     } // End of check if the node is subject to Dirichlet BC-
   } // End of loop over the nodes
 } // End of function
   

   /// Vector of pointers to SingularPoissonSolutionElement objects
 Vector<SingularPoissonSolutionElement<PoissonElementWithSingularity
  <BASIC_POISSON_ELEMENT> >*> C_equation_elements_pt;

 /// Boolean indicating which node is subject to Dirichlet BC
 /// [size = number of nodes; initialised to false]
 vector<bool> Node_is_subject_to_dirichlet_bc;
 
 /// Imposed value at nodes that are subject to Dirichlet BC
 /// [size = number of nodes; initialised to zero]
 Vector<double> Imposed_value_at_node;
   
};



#endif

