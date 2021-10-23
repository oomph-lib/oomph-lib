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
//Header file for general soluble surfactant transport equations.
//These inherit from the insoluble surfactant transport equations
//but add a flux term between the surface concentration and a
//bulk concentration. 

//Oomph-lib headers, 
//We require the generic header
#include "generic.h"
//The fluid interface elements
#include "fluid_interface.h"

namespace oomph
{

 class SolubleSurfactantTransportInterfaceElement :
  public SurfactantTransportInterfaceElement
 
{
private:
  
 /// Pointer to a Biot number
 double *Bi_pt;

 /// Pointer to a Marangoni number
double *Ma_pt;

 /// Pointer to the reaction ratios
 double *K_pt;

 /// Pointer to beta
 double *Beta_b_pt;

  /// Index at which the bulk concentration is stored at the nodes
 unsigned C_bulk_index;

  /// Index at which the micelle concentration is stored at the nodes
 unsigned M_index;

protected:

 ///Get the micelle concentration
 double interpolated_M(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned m_index = M_index;

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   this->shape(s,psi);

   //Initialise value of m
   double M = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
      M += this->nodal_value(l,m_index)*psi(l);
    }

   return(M);
  }

   ///Get the bulk concentration
 double interpolated_C_bulk(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned c_index = C_bulk_index;

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   this->shape(s,psi);

   //Initialise value of C
   double C = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
      C += this->nodal_value(l,c_index)*psi(l);
    }

   return(C);
  }


 /// The surface tension function is linear in the
 /// temperature with constant of proportionality equal
 /// to the Marangoni number.
 double sigma(const Vector<double> &s)
  {
   //Get the value of the concentration
   const double C = this->interpolated_C(s);
   //Get the Elasticity numbers
   const double Beta = this->beta();
   //Return the variable surface tension
   return (1.0 + Beta*log(1.0-C));
  }

 //Overload the derivative function
 double dsigma_dC(const Vector<double> &s)
 {
  const double C = this->interpolated_C(s);
  const double Beta = this->beta();
  return (-Beta/(1.0 - C));
 }
   
 /// Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
  /*void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
    //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
  
   //Call the external stuff (but there isn't any!)
   //this->fill_in_jacobian_from_external_by_fd(jacobian);
   }*/

 ///Specify the flux from the bulk to the interface
 double flux_from_bulk(const double &C, const double &C_bulk)
 {
  const double Bi = this->bi();
  const double K = this->k();
  return Bi*(K*C_bulk*(1.0 - C) - C);
 }

  ///Specify the derivative of the flux from the bulk
  ///to the interface with respect to the bulk concentration
 double dflux_from_bulk_dC_bulk(const double &C, const double &C_bulk)
 {
  const double Bi = this->bi();
  const double K = this->k();
  return Bi*K*(1.0 - C);
 }

  ///Specify the derivative of the flux from the bulk
  ///to the interface with respect to the surface concentration
 double dflux_from_bulk_dC(const double &C, const double &C_bulk)
 {
  const double Bi = this->bi();
  const double K = this->k();
  return -Bi*(K*C_bulk + 1.0);
 }


 
 /// Overload the Helper function to calculate the residuals and 
 /// jacobian entries. This particular function ensures that the
 /// additional entries are calculated inside the integration loop
 void add_additional_residual_contributions_interface(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned &flag,
  const Shape &psif, const DShape &dpsifds,
  const DShape &dpsifdS, const DShape &dpsifdS_div,
  const Vector<double> &s,
  const Vector<double> &interpolated_x, 
  const Vector<double> &interpolated_n, 
  const double &W,
  const double &J)
  {
   //Call the underlying version from the Surfactant Transport equation
   SurfactantTransportInterfaceElement::
    add_additional_residual_contributions_interface(
     residuals,jacobian,flag,psif,dpsifds,dpsifdS, dpsifdS_div,s,
     interpolated_x,interpolated_n,W,J);

   //Now we need to add the bulk contribution
   
   //Find out how many nodes there are
   unsigned n_node = this->nnode();

   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;

   //Bulk flux condition
   //Find the index at which the bulk is stored
   unsigned c_bulk_index = this->C_bulk_index;

   
   //Now calculate the bulk concentration at this point
   //Assuming the same shape functions are used (which they are)
   double C_bulk = 0.0;
   double C = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     const double psi_ = psif(l);
     C_bulk += this->nodal_value(l,c_bulk_index)*psi_;
     C += this->nodal_value(l,this->C_index[l])*psi_;
    }

   //The transport between the two layers is given by the flux
   const double flux = this->flux_from_bulk(C,C_bulk);
   //Compute the derivatives if required for the Jacobian
   const double dflux_dC_bulk = this->dflux_from_bulk_dC_bulk(C,C_bulk);
   const double dflux_dC = this->dflux_from_bulk_dC(C,C_bulk);

   //Read out the beta (solubility) parameter
   const double Beta_b = this->beta_b();
   
   //Now we add the flux term to the bulk residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_bulk_index);

     //If not a boundary condition
     if(local_eqn >= 0)
      {
       //Add the flux out of the bulk
	residuals[local_eqn] -= Beta_b*flux*psif(l)*W*J;

       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the unknown
           local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
	      jacobian(local_eqn,local_unknown) -=
               Beta_b*dflux_dC_bulk*psif(l2)*psif(l)*W*J;
            }

           local_unknown = this->nodal_local_eqn(l2,this->C_index[l2]);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
	      jacobian(local_eqn,local_unknown) -= Beta_b*dflux_dC*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes


  //Now we add the flux term to the appropriate residuals
  for(unsigned l=0;l<n_node;l++)
   {
    //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,this->C_index[l]);
    
    //If not a boundary condition
    if(local_eqn >= 0)
     {
      //Add the flux term
      residuals[local_eqn] -= flux*psif(l)*W*J;
      
      //We also need to worry about the jacobian terms
      if(flag)
       {
        //Loop over the nodes again
        for(unsigned l2=0;l2<n_node;l2++)
         {
          //Get the unknown c_index
          local_unknown =this->nodal_local_eqn(l2,this->C_index[l2]);
          
          if(local_unknown >=0)
           {
	    //Add the flux term
	    jacobian(local_eqn,local_unknown) -=
             dflux_dC*psif(l2)*psif(l)*W*J;
	   }

	  //Local unknown the bulk concentration
	  local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
	  //If not a boundary condition
	  if(local_unknown >= 0)
	    {
             jacobian(local_eqn,local_unknown) -=
              dflux_dC_bulk*psif(l2)*psif(l)*W*J;
	    }
         }
       }
     }
  } //End of loop over the nodes

  }  

 
 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  DenseMatrix<double> &mass_matrix)
  {
   //Add the contribution to the jacobian
   this->fill_in_contribution_to_jacobian(residuals,jacobian);
   //No mass matrix terms, but should probably do kinematic bit here
   //for a stability analysis
  }

public:

 ///Constructor that passes the bulk element and face index
 ///down to the underlying 
  SolubleSurfactantTransportInterfaceElement() :
 SurfactantTransportInterfaceElement()
  {
   //Initialise the values
   Bi_pt = &Default_Physical_Constant_Value;
   Ma_pt = &Default_Physical_Constant_Value;
   K_pt = &Default_Physical_Constant_Value;
   Beta_b_pt = &Default_Physical_Constant_Value;
   //Big Hack
   M_index = 3;
  }

 //set the C bulk index
 inline void set_c_bulk_index(const unsigned &c_bulk_index)
 {C_bulk_index = c_bulk_index;}
  
  ///Return the Biot number
 double bi() {return *Bi_pt;}
 
 ///Return the Marangoni number
 double ma() {return *Ma_pt;}


 ///Return the reaction ratio
 double k() {return *K_pt;}

  //Return bulk beta
 double beta_b() {return *Beta_b_pt;}
  
 ///Access function for pointer to the Marangoni number
 double* &ma_pt() {return Ma_pt;}

 ///Access function for pointer to the Biot number
 double* &bi_pt() {return Bi_pt;}

 ///Access function for pointer to the reaction ratios
 double* &k_pt() {return K_pt;}

 ///Access function for pointer
 double* &beta_b_pt() {return Beta_b_pt;}


  //Calculate the mean square deviation from the
  //uniform state
  double l2_norm_of_height(const double &h0)
  {
    //Find the number of nodes
    const unsigned n_node = this->nnode();
    
    //Find out the number of surface coordinates
    const unsigned el_dim = this->dim();

    //Find the nodal dimension
    const unsigned n_dim = this->node_pt(0)->ndim();
    
    //Storage for the Shape functions
    Shape psif(n_node);

    //Storage for the local coordinate
    Vector<double> s(el_dim);
    
    //Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    //Storage for answer
    double norm = 0.0;

    //Work out x distance, assuming this is a line element
    //which it is
    double scale = (this->node_pt(n_node-1)->x(0) - this->node_pt(0)->x(0))/2.0;
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	//Get the value of the local coordiantes at the integration point
	for(unsigned i=0;i<el_dim;i++) {s[i] = this->integral_pt()->knot(ipt,i);}
	
	//Get the integral weight
	double W = this->integral_pt()->weight(ipt);

	//Call the derivatives of the shape function
	this->shape_at_knot(ipt,psif);
    
	//Define and zero the tangent Vectors and local velocities
	Vector<double> interpolated_x(n_dim,0.0);
	
	//Loop over the shape functions
	for(unsigned l=0;l<n_node;l++)
	  {
	    const double psi_ = psif(l);
	    //Loop over directional components 
	    for(unsigned i=0;i<n_dim;i++)
	      {
		// Coordinate
		interpolated_x[i] += this->nodal_position(l,i)*psi_;
	      }
	  }
	
	//Calculate the surface gradient and divergence
	norm += (interpolated_x[1] - h0)*(interpolated_x[1] - h0)*W*scale;
      }
    //Return the computed norm
    return norm;
  }
  

  //Calculate the total mass over the element
  double integrated_C()
  {
    //Find the number of nodes
    const unsigned n_node = this->nnode();
    
    //Find out the number of surface coordinates
    const unsigned el_dim = this->dim();

    //Find the nodal dimension
    //const unsigned n_dim = this->node_pt(0)->ndim();
    
    //Storage for the Shape functions
    Shape psif(n_node);

    //Storage for the local coordinate
    Vector<double> s(el_dim);
    
    //Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    //Storage for answer
    double integral = 0.0;
        
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
	//Get the value of the local coordiantes at the integration point
	for(unsigned i=0;i<el_dim;i++) {s[i] = this->integral_pt()->knot(ipt,i);}
	
	//Get the integral weight
	double W = this->integral_pt()->weight(ipt);

	//Calculate the Jacobian of the mapping
        double J = this->J_eulerian(s);
        
	//Call the derivatives of the shape function
	this->shape_at_knot(ipt,psif);
    
	//Define and zero the tangent Vectors and local velocities
	double interpolated_C=0.0;
	
	//Loop over the shape functions
	for(unsigned l=0;l<n_node;l++)
	  {
	    const double psi_ = psif(l);
	    //Loop over directional components 
            interpolated_C += this->nodal_value(l,this->C_index[l])*psi_;
	  }

        integral += interpolated_C*W*J;
      }
    //Return the computed norm
    return integral;
  }
  

 
void output(std::ostream &outfile, const unsigned &n_plot)
{
 //Set output Vector
 Vector<double> s(1);
 
 //Tecplot header info 
 outfile << "ZONE I=" << n_plot << std::endl;
 
 //Loop over plot points
 for(unsigned l=0;l<n_plot;l++)
  {
   s[0] = -1.0 + l*2.0/(n_plot-1);
   
   //Output the x,y,u,v 
   for(unsigned i=0;i<2;i++) outfile << this->interpolated_x(s,i) << " ";
   for(unsigned i=0;i<2;i++) outfile << this->interpolated_u(s,i) << " ";      
   //Output a dummy pressure
   outfile << 0.0 << " ";
   const double C = this->interpolated_C(s);
   const double C_bulk = this->interpolated_C_bulk(s);
   const double M = this->interpolated_M(s);
   
   //Output the concentrations
   outfile << C << " "
	   << C_bulk << " "
           << M << " "
	   << this->bi()*(this->k()*C_bulk*(1.0 - C) - C) << std::endl;
  }
 outfile << std::endl;
}


///Overload the output function
  void output(std::ostream &outfile) {FiniteElement::output(outfile);}
     
  ///Overload the C-style output function
  void output(FILE* file_pt) {FiniteElement::output(file_pt);}
  
  ///C-style Output function
  void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}


};

 


///=============================================================================
///This is the policy class for the surfactanttransport equations which require
///one additional value for the surface concentration
//=============================================================================
template<>
 class FluidInterfaceAdditionalValues<
 SolubleSurfactantTransportInterfaceElement>
 {
   public:
  FluidInterfaceAdditionalValues<
   SolubleSurfactantTransportInterfaceElement>() {}

  inline unsigned nadditional_values(const unsigned &n) {return 1;}

  inline void setup_equation_indices(
   SolubleSurfactantTransportInterfaceElement* const &element_pt,
   const unsigned &id)
  {
   const unsigned n_node = element_pt->nnode();
   Vector<unsigned> c_index(n_node);
   for(unsigned n=0;n<n_node;n++)
    {
     c_index[n] =
      dynamic_cast<BoundaryNodeBase*>(element_pt->node_pt(n))
      ->index_of_first_value_assigned_by_face_element(id);
    }
   
   element_pt->set_c_index(c_index);
   //This is a hack
   const unsigned C_bulk_index =
    dynamic_cast<AdvectionDiffusionReactionEquations<2,2>*>(
      element_pt->bulk_element_pt())->c_index_adv_diff_react(0);
    element_pt->set_c_bulk_index(C_bulk_index);
   
  }
  
 };


//-------------------------------GEOMETRIC SPECIALISATIONS----------------

///Specialise to the Line geometry
template<class ELEMENT>
 class SpineLineSolubleSurfactantTransportInterfaceElement :
 public SpineUpdateFluidInterfaceElement<
 SolubleSurfactantTransportInterfaceElement,
 LineDerivatives,ELEMENT>
 {
   public:
  
   SpineLineSolubleSurfactantTransportInterfaceElement(
  FiniteElement* const &element_pt, 
  const int &face_index) : 
  SpineUpdateFluidInterfaceElement<SolubleSurfactantTransportInterfaceElement,
				   LineDerivatives,ELEMENT>(element_pt,face_index) {}
};


//Define the bounding element type for the line elements 
//This will need to be updated
template<class ELEMENT>
 class BoundingElementType<SpineUpdateFluidInterfaceElement<
 SolubleSurfactantTransportInterfaceElement,LineDerivatives,ELEMENT> >:
 public SpinePointFluidInterfaceBoundingElement<ELEMENT>
 {
   public:
  
  BoundingElementType<SpineUpdateFluidInterfaceElement<
   SolubleSurfactantTransportInterfaceElement,LineDerivatives,ELEMENT> >() :
  SpinePointFluidInterfaceBoundingElement<ELEMENT>() { }
 };

 ///Specialise to the Line geometry
template<class ELEMENT>
 class ElasticLineSolubleSurfactantTransportInterfaceElement :
 public ElasticUpdateFluidInterfaceElement<
 SolubleSurfactantTransportInterfaceElement,
 LineDerivatives,ELEMENT>
 {
   public:
  
   ElasticLineSolubleSurfactantTransportInterfaceElement(
  FiniteElement* const &element_pt, 
  const int &face_index) : 
  ElasticUpdateFluidInterfaceElement<SolubleSurfactantTransportInterfaceElement,
   LineDerivatives,ELEMENT>(element_pt,face_index) {}
 };


//Define the bounding element type for the line elements 
//This will need to be updated
template<class ELEMENT>
 class BoundingElementType<ElasticUpdateFluidInterfaceElement<
 SolubleSurfactantTransportInterfaceElement,LineDerivatives,ELEMENT> >:
 public ElasticPointFluidInterfaceBoundingElement<ELEMENT>
 {
   public:
  
  BoundingElementType<ElasticUpdateFluidInterfaceElement<
   SolubleSurfactantTransportInterfaceElement,LineDerivatives,ELEMENT> >() :
  ElasticPointFluidInterfaceBoundingElement<ELEMENT>() { }
 };



 /// Old Implementation. Should ultimately be deleted
 

  template<class ELEMENT>
  class SpinePointMarangoniSurfactantFluidInterfaceBoundingElement : 
    public SpinePointFluidInterfaceBoundingElement<ELEMENT>
  
  {
    unsigned C_index;
    
    public:
 
    /// Constructor
    SpinePointMarangoniSurfactantFluidInterfaceBoundingElement() : 
      SpinePointFluidInterfaceBoundingElement<ELEMENT>(),
      C_index(0) {}
    
    /// Calculate the elemental residual vector and the Jacobian
    void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
					  DenseMatrix<double> &jacobian)
    {
      SpinePointFluidInterfaceBoundingElement<ELEMENT>::
	fill_in_contribution_to_jacobian(residuals,jacobian);
      
      //Make sure to include concentration variations
   {
    //Use finite differences to handle bulk concentration variations
    const unsigned n_node = this->nnode();
    //Find the number of dofs in the element
    const unsigned n_dof = this->ndof();
    //Create newres vector
    Vector<double> newres(n_dof);
    
    //Integer storage for local unknown
    int local_unknown=0;
    
    //Use the default finite difference step
    const double fd_step = this->Default_fd_jacobian_step;
    
    
    //Use finite differences to handle the interfacial concentration variations
    //Loop over the nodes again
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned c_index = this->C_index;

      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,c_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(c_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
       
        //Increment the value of the Nodal data
        *value_pt += fd_step;
       
        //Calculate the new residuals
        this->get_residuals(newres);
       
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
       
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
     }

    }

    //Allow the concentration index to be set
    unsigned &c_index_interface_boundary() {return C_index;}
    
 }; 

}









