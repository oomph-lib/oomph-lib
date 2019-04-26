///LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
   
 /// \short Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
    //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
  
   //Call the external stuff (but there isn't any!)
   //this->fill_in_jacobian_from_external_by_fd(jacobian);
  }

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


 
 /// \short Overload the Helper function to calculate the residuals and 
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
     dynamic_cast<DoubleBuoyantQCrouzeixRaviartElement<2>*>(
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
  
  
//==================================================================
///Spine-based Marangoni surface tension elements that add
///a linear dependence on temperature and concentration
///of a surface chemical to the surface tension, 
///which decreases with increasing temperature. In addition, 
///this element adds a flux contribution to the advection-diffusion
///equation to represent heat loss at the free surface. This
///introduces the Biot number.
//=================================================================
template<class ELEMENT>
class SpineLineMarangoniSurfactantFluidInterfaceElement :
public SpineLineFluidInterfaceElement<ELEMENT>
{
private:
  
 /// Pointer to a Biot number
 double *Bi_pt;

 /// Pointer to a Marangoni number
 double *Ma_pt;

 /// Pointer to an Elasticity number
 double *Beta_pt;

 /// Pointer to Surface Peclet number
 double *Peclet_S_pt;

 /// Pointer to the surface Peclect Strouhal number
 double *Peclet_Strouhal_S_pt;

 /// Pointer to the diffusion ratios
 double *D_pt;

 /// Pointer to the reaction ratios
 double *K_pt;

  /// Pointer to beta
  double *Beta_b_pt;

  /// Pointer to bulk peclet number
  double *Pe_b_pt;
  
 /// Index at which the bulk concentration is stored at the nodes
 unsigned C_bulk_index;

  /// Index at which the micelle concentration is stored at the nodes
 unsigned M_index;
  
 /// Index at which the surfactant concentration is stored at the
 /// nodes
 unsigned C_index;

 /// Default value of the physical constants
 static double Default_Physical_Constant_Value;

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


 ///Get the surfactant concentration
 double interpolated_C(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned c_index = C_index;

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


 /// The time derivative of the surface concentration
 double dcdt_surface(const unsigned &l) const
  {
   // Get the data's timestepper
   TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
   
   //Initialise dudt
   double dcdt=0.0;
   //Loop over the timesteps, if there is a non Steady timestepper
   if (time_stepper_pt->type()!="Steady")
    {
     //Find the index at which the variable is stored
     const unsigned c_index = C_index;

     // Number of timsteps (past & present)
     const unsigned n_time = time_stepper_pt->ntstorage();
     
     for(unsigned t=0;t<n_time;t++)
      {
       dcdt += time_stepper_pt->weight(1,t)*this->nodal_value(t,l,c_index);
      }
    }
   return dcdt;
  }

 /// The surface tension function is linear in the
 /// temperature with constant of proportionality equal
 /// to the Marangoni number.
 double sigma(const Vector<double> &s)
  {
   //Find the number of shape functions
   const unsigned n_node = this->nnode();
   //Now get the shape fuctions at the local coordinate
   Shape psi(n_node);
   this->shape(s,psi);

   //Now interpolate the surface surfactant concentration
   double C=0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     C += this->nodal_value(l,C_index)*psi(l);
    }
   
   //Get the Elasticity numbers
   double Beta = this->beta();
   //Return the variable surface tension
   return (1.0 + Beta*log(1.0-C));
  }

 /// \short Fill in the contribution to the residuals
  /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian)
  {
    //Call the generic routine with the flag set to 1
   this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
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


   //Call the external stuff
   this->fill_in_jacobian_from_external_by_fd(jacobian);

   //Call generic routine to handle spine variables
   this->fill_in_jacobian_from_geometric_data(jacobian);
   
   //Call the generic routine to handle the spine variables
   //SpineElement<FaceGeometry<ELEMENT> >::
   // fill_in_jacobian_from_geometric_data(jacobian);
  }
  
 
 /// \short Overload the Helper function to calculate the residuals and 
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
   //Find out how many nodes there are
   unsigned n_node = this->nnode();

   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;

   //Bulk flux condition
   //Find the index at which the bulk is stored
   unsigned c_bulk_index = this->C_bulk_index;

   //Find the index at which the surface surfactant
   unsigned c_index = this->C_index;

   
   //Now calculate the bulk concentration at this point
   //Assuming the same shape functions are used (which they are)
   double C_bulk = 0.0;
   double C = 0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     C_bulk += this->nodal_value(l,c_bulk_index)*psif(l);
     C += this->nodal_value(l,c_index)*psif(l);
    }

   //Get the reaction ratio
   const double K = this->k();
   //Read out the Bi number
   const double Bi = this->bi();

   //Read out the bulk peclet number and beta parameter
   const double Beta_b = this->beta_b();

   //The transport between the two layers is given by the flux
   double flux = Bi*(K*C_bulk*(1.0 - C) - C);
   
   //Now we add the flux term to the appropriate residuals
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
	      jacobian(local_eqn,local_unknown) -= Beta_b*Bi*K*(1.0-C)*psif(l2)*psif(l)*W*J;
            }
           
           local_unknown = this->nodal_local_eqn(l2,c_index);
           //If not a boundary condition
           if(local_unknown >= 0)
            {
	      jacobian(local_eqn,local_unknown) += Beta_b*Bi*(K*C_bulk + 1.0)*psif(l2)*psif(l)*W*J;
            }
          }
        }
      }
    } //End of loop over the nodes


   //Read out the surface peclect numbers
   const double Pe_s = this->peclet_s();
   //const double PeSt_s = this->peclet_strouhal_s();
   
    //Find the index at which the temperature is stored
   //unsigned c_index = this->C_index;
   Vector<unsigned> u_index = this->U_index_interface;

   //Now calculate the concentration at this point
   //Assuming the same shape functions are used (which they are)
   double interpolated_C = 0.0;
   //double interpolated_dCds = 0.0;
   double dCdt = 0.0;
   //The tangent vector
   const unsigned n_dim = this->node_pt(0)->ndim();
   Vector<double> interpolated_u(n_dim,0.0);
   Vector<double> mesh_velocity(n_dim,0.0);
   Vector<double> interpolated_grad_C(n_dim,0.0);
   double interpolated_div_u = 0.0;
  
  //Loop over the shape functions
  for(unsigned l=0;l<n_node;l++)
   {
    const double psi = psif(l);
    const double C_ = this->nodal_value(l,c_index);
    
    interpolated_C += C_*psi;
    dCdt += dcdt_surface(l)*psi;
    //Velocity and Mesh Velocity
    for(unsigned j=0;j<n_dim;j++)
     {
      const double u_ = this->nodal_value(l,u_index[j]);
      interpolated_u[j] += u_*psi;
      mesh_velocity[j] += this->dnodal_position_dt(l,j)*psi;
      interpolated_grad_C[j] += C_*dpsifdS(l,j);
      interpolated_div_u += u_*dpsifdS_div(l,j);
     }
   }
  
  //Pre-compute advection term
  double interpolated_advection_term = interpolated_C*interpolated_div_u;
  for(unsigned i=0;i<n_dim;i++)
   {
    interpolated_advection_term += (interpolated_u[i] - mesh_velocity[i])*interpolated_grad_C[i];
   }
  
  
  //Now we add the flux term to the appropriate residuals
  for(unsigned l=0;l<n_node;l++)
   {
    //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_index);
    
    //If not a boundary condition
    if(local_eqn >= 0)
     {
      //Time derivative term
      residuals[local_eqn] += dCdt*psif(l)*W*J;
      
      //First Advection term
      residuals[local_eqn] += interpolated_advection_term*psif(l)*W*J;
      
      //Diffusion term
      double diffusion_term = 0.0;
      for(unsigned i=0;i<n_dim;i++)
       {
        diffusion_term += interpolated_grad_C[i]*dpsifdS(l,i);
       }
      residuals[local_eqn] += (1.0/Pe_s)*diffusion_term*W*J;

      //Add the flux term
      residuals[local_eqn] -= flux*psif(l)*W*J;
      
      //We also need to worry about the jacobian terms
      if(flag)
       {
        //Loop over the nodes again
        for(unsigned l2=0;l2<n_node;l2++)
         {
          //Get the time stepper
          TimeStepper* time_stepper_pt=this->node_pt(l2)->time_stepper_pt();
          
          //Get the unknown c_index
          local_unknown =this->nodal_local_eqn(l2,c_index);
          
          if(local_unknown >=0)
           {
            jacobian(local_eqn,local_unknown) += 
             time_stepper_pt->weight(1,0)*psif(l2)*psif(l)*W*J;
            
            jacobian(local_eqn,local_unknown) +=
             psif(l2)*interpolated_div_u*psif(l)*W*J;
            
            for(unsigned i=0;i<n_dim;i++)
             {
              jacobian(local_eqn,local_unknown) +=
               (interpolated_u[i] - mesh_velocity[i])*dpsifdS(l2,i)*psif(l)*W*J;
             }
            
            for(unsigned i=0;i<n_dim;i++)
             {
              jacobian(local_eqn,local_unknown) += (1.0/Pe_s)*dpsifdS(l2,i)*dpsifdS(l,i)*W*J;
             }

	    //Add the flux term
	    jacobian(local_eqn,local_unknown) += Bi*(K*C_bulk + 1.0)*psif(l2)*psif(l)*W*J;
	   }

	  //Local unknown the bulk concentration
	  local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
	  //If not a boundary condition
	  if(local_unknown >= 0)
	    {
	      jacobian(local_eqn,local_unknown) -= Bi*K*(1.0 - C)*psif(l2)*psif(l)*W*J;
	    }
	            
          //Loop over the velocity components
          for(unsigned i2=0;i2<n_dim;i2++)
           {
            
            //Get the unknown
            local_unknown = this->nodal_local_eqn(l2,u_index[i2]);
            
            
            //If not a boundary condition
            if(local_unknown >= 0)
             {
              //Bits from the advection term
              jacobian(local_eqn,local_unknown) += (interpolated_C*dpsifdS_div(l2,i2) +
                                                    psif(l2)*interpolated_grad_C[i2])
               *psif(l)*W*J;
             }
           }
         }
       }
     }

    //Analytic derivative of surface tension relationship
    /*  if(flag)
     {
      const double dsigma = this->dsigma_dC(s);
      const double Ca = this->ca();
      for(unsigned l2=0;l2<n_node;l2++)
       {
        local_unknown = this->nodal_local_eqn(l2,c_index);
        if(local_unknown >= 0)
         {
          const double psi_ = psif(l2);
          for(unsigned i=0;i<n_dim;i++)
           {
            //Add the Jacobian contribution from the surface tension
            local_eqn = this->nodal_local_eqn(l,u_index[i]);
            if(local_eqn >= 0)
             {
              jacobian(local_eqn,local_unknown) -= (dsigma/Ca)*psi_*dpsifdS_div(l,i)*J*W;
             }
           }
         }
       }
       }*/
    
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
 SpineLineMarangoniSurfactantFluidInterfaceElement(
 FiniteElement* const &element_pt, const int &face_index) : 
  SpineLineFluidInterfaceElement<ELEMENT>
  (element_pt,face_index)
  {
   //Initialise the values
   Bi_pt = &Default_Physical_Constant_Value;
   Ma_pt = &Default_Physical_Constant_Value;
   Beta_pt = &Default_Physical_Constant_Value;
   Peclet_S_pt = &Default_Physical_Constant_Value;
   Peclet_Strouhal_S_pt = &Default_Physical_Constant_Value;
   K_pt = &Default_Physical_Constant_Value;
   D_pt = &Default_Physical_Constant_Value;
   Beta_b_pt = &Default_Physical_Constant_Value;
   Pe_b_pt = &Default_Physical_Constant_Value;
   
   //Cast the bulk element 
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   //Now find the index at which the micelle concentration is stored from the 
   //advection-diffusion part of the bulk element
   M_index = cast_element_pt->c_index_adv_diff_react(1);
   
   //Find the index at which the bulk surfactant concentration is stored
   C_bulk_index = cast_element_pt->c_index_adv_diff_react(0);

   //Add the additional surfactant terms to these surface elements
   
   //Read out the number of nodes on the face
   //For some reason I need to specify the this pointer here(!)
   unsigned n_node_face = this->nnode();
   //Set the additional data values in the face
   //There is one additional values at each node --- the
   //surface concentration
   Vector<unsigned> additional_data_values(n_node_face);
   for(unsigned i=0;i<n_node_face;i++) additional_data_values[i] = 1;
   //Resize the data arrays accordingly 
   this->resize_nodes(additional_data_values);

   //The C_index is the new final value
   //HACK HERE
   C_index = this->node_pt(0)->nvalue()-1;
  }



  //Make the bounding element
  FluidInterfaceBoundingElement* make_bounding_element(
   const int &face_index)
  {
   //Create a temporary pointer to the appropriate FaceElement read our from
   //our policy class
    SpinePointMarangoniSurfactantFluidInterfaceBoundingElement<ELEMENT> *face_el_pt =
      new  SpinePointMarangoniSurfactantFluidInterfaceBoundingElement<ELEMENT>;
    
    //Attach the geometrical information to the new element
    this->build_face_element(face_index,face_el_pt);
    
    //Set the index at which the velocity nodes are stored
    face_el_pt->u_index_interface_boundary() = this->U_index_interface;
    face_el_pt->c_index_interface_boundary() = this->C_index;
    
    //Set the value of the nbulk_value, the node is not resized
    //in this bounding element,
    //so it will just be the actual nvalue here
    // There is some ambiguity about what this means (however)
    // We are interpreting it to mean the number of
    // values in this FaceElement before creating the new
    // bounding element.
    const unsigned n_node_bounding = face_el_pt->nnode();
    for(unsigned n=0;n<n_node_bounding;n++)
      {
	face_el_pt->nbulk_value(n) =
	  face_el_pt->node_pt(n)->nvalue();
      }
    
    //Set of unique geometric data that is used to update the bulk,
    //but is not used to update the face
    std::set<Data*> unique_additional_geom_data;
    
    //Get all the geometric data for this (bulk) element
    this->assemble_set_of_all_geometric_data(unique_additional_geom_data);
    
    //Now assemble the set of geometric data for the face element
    std::set<Data*> unique_face_geom_data_pt;
    face_el_pt->assemble_set_of_all_geometric_data(unique_face_geom_data_pt);
    
    //Erase the face geometric data from the additional data
    for(std::set<Data*>::iterator it=unique_face_geom_data_pt.begin();
	it!=unique_face_geom_data_pt.end();++it)
      {unique_additional_geom_data.erase(*it);}
    
    //Finally add all unique additional data as external data
    for(std::set<Data*>::iterator it = unique_additional_geom_data.begin();
	it!= unique_additional_geom_data.end();++it)
      {
	face_el_pt->add_external_data(*it);
      }
    
    //Return the value of the pointer
    return face_el_pt;
  }


  ///Return the Biot number
 double bi() {return *Bi_pt;}
 
 ///Return the Marangoni number
 double ma() {return *Ma_pt;}

 ///Return the Elasticity number
 double beta() {return *Beta_pt;}

 ///Return the surface peclect number
 double peclet_s() {return *Peclet_S_pt;}

 ///Return the surface peclect strouhal number
 double peclet_strouhal_s() {return *Peclet_Strouhal_S_pt;}

 ///Return the diffusion ratio
 double d() {return *D_pt;}

 ///Return the reaction ratio
 double k() {return *K_pt;}

  //Return bulk beta
 double beta_b() {return *Beta_b_pt;}

 ///Return the reaction ratio
 double pe_b() {return *Pe_b_pt;}

  
 ///Access function for pointer to the Marangoni number
 double* &ma_pt() {return Ma_pt;}

 ///Access function for pointer to the Biot number
 double* &bi_pt() {return Bi_pt;}

 ///Access function for pointer to the Elasticity number
 double* &beta_pt() {return Beta_pt;}

 ///Access function for pointer to the surface Peclet number
 double* &peclet_s_pt() {return Peclet_S_pt;}

 ///Access function for pointer to the surface Peclet x Strouhal number
 double* &peclet_strouhal_s_pt() {return Peclet_Strouhal_S_pt;}

 ///Access function for pointer to the diffusion ratios
 double* &d_pt() {return D_pt;}

 ///Access function for pointer to the reaction ratios
 double* &k_pt() {return K_pt;}

 ///Access function for pointer
 double* &beta_b_pt() {return Beta_b_pt;}

 ///Access function for pointer to the reaction ratios
 double* &pe_b_pt() {return Pe_b_pt;}


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

    //Storage for the index at which surface concentration is stored
    const unsigned c_index = this->C_index;
        
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
            interpolated_C += this->nodal_value(l,c_index)*psi_;
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



};


//Define the default physical value to be one
template<class ELEMENT>
double SpineLineMarangoniSurfactantFluidInterfaceElement<ELEMENT>::
Default_Physical_Constant_Value = 1.0;


}









