//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
//Header file for elements that are used to provide streamfunction traction
//at inlet and outlet

#ifndef OOMPH_POLAR_STREAMFUNCTION_TRACTION_ELEMENTS
#define OOMPH_POLAR_STREAMFUNCTION_TRACTION_ELEMENTS

//OOMPH-LIB headers
#include "generic/Qelements.h"

namespace oomph
{

//======================================================================
/// A class for elements that allow the imposition of an applied traction
/// to the Navier--Stokes equations 
//======================================================================
template <class ELEMENT>
class PolarStreamfunctionTractionElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
private:

 /// The highest dimension of the problem 
 unsigned Dim;

protected:


 /// Access function that returns the local equation number
 /// for the streamfunction components.
 /// The default is to asssume that n is the local node number
 /// and the streamfunction component is the 1st unknown stored at the node.
 virtual inline int s_local_eqn(const unsigned &n)
  {return nodal_local_eqn(n,0);}
 
 /// Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping 
 inline double shape_and_test_at_knot(const unsigned &ipt, 
                                      Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();
   //Calculate the shape functions
   shape_at_knot(ipt,psi);
   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}
   //Return the value of the jacobian
   return J_eulerian_at_knot(ipt);
  }

 /// This function returns the residuals for the 
 /// traction function.
 /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
 void fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                            DenseMatrix<double> &jacobian,
                                            DenseMatrix<double> &mass_matrix,
                                            unsigned flag);
 /// Pointer to the angle alpha
 double *Alpha_pt;

 //Traction elements need to know whether they're at the inlet or outlet
 //as the unit outward normal has a differing sign dependent on
 //the boundary
 // -1=inlet, 1=outlet
 int Boundary;

public:

 /// Alpha
 const double &alpha() const {return *Alpha_pt;}

 /// Pointer to Alpha
 double* &alpha_pt() {return Alpha_pt;}

 /// Boundary
 const int boundary() const {return Boundary;}

 /// Function to set boundary
 void set_boundary(int bound) {Boundary=bound;}

 /// Constructor, which takes a "bulk" element and the value of the index
 /// and its limit
 PolarStreamfunctionTractionElement(FiniteElement* const &element_pt, 
                                    const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 
#ifdef PARANOID
 {
  //Check that the element is not a refineable 3d element
  ELEMENT* elem_pt = new ELEMENT;
  //If it's three-d
  if(elem_pt->dim()==3)
   {
    //Is it refineable
    if(dynamic_cast<RefineableElement*>(elem_pt))
     {
      //Issue a warning
      OomphLibWarning(
       "This flux element will not work correctly if nodes are hanging\n",
       "PolarNavierStokesTractionElement::Constructor",
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 }
#endif

   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
 element_pt->build_face_element(face_index,this);

   //Set the dimension from the dimension of the first node
   Dim = this->node_pt(0)->ndim();
 }

 /// Destructor should not delete anything
 ~PolarStreamfunctionTractionElement() { }
 
 /// This function returns just the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution(
    residuals,GeneralisedElement::Dummy_matrix,GeneralisedElement::Dummy_matrix,0);
  }
 
 /// This function returns the residuals and the jacobian
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,jacobian,GeneralisedElement::Dummy_matrix,1);
  }

 /// Compute the element's residual Vector and the jacobian matrix
 /// Plus the mass matrix especially for eigenvalue problems
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {   
   //Call the generic routine with the flag set to 2
   fill_in_generic_residual_contribution(residuals,jacobian,GeneralisedElement::Dummy_matrix,2);
  }
 
 /// Overload the output function
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

/// Output function: x,y,[z],u,v,[w],p in tecplot format
void output(std::ostream &outfile, const unsigned &nplot)
 {FiniteElement::output(outfile,nplot);}

 /// local streamfunction
 double s(const unsigned &l )
  {return nodal_value(l,0);}

 /// local position
 double x(const unsigned &l, const unsigned &i )
  {return nodal_position(l,i);}

}; 



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//============================================================================
/// Function that returns the residuals for the imposed traction Navier_Stokes
/// equations
//============================================================================
template<class ELEMENT>
void PolarStreamfunctionTractionElement<ELEMENT>::
fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                  DenseMatrix<double> &jacobian,
				  DenseMatrix<double> &mass_matrix,
                                  unsigned flag)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();
 
 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = integral_pt()->nweight();

 //Get Alpha
 const double Alpha = alpha();

 //Get boundary multiplier
 //This is necessary because the sign of the traction is 
 //dependent on the boundary
 const int multiplier = boundary();
 
 //Integers to store local equation numbers
 int local_eqn=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test_at_knot(ipt,psif,testf);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Need to find position to feed into Traction function
   Vector<double> interpolated_x(Dim);
   double interpolated_v;   

   //Initialise to zero
   interpolated_v = 0.0;
   for(unsigned i=0;i<Dim;i++) interpolated_x[i] = 0.0;
   
   //Calculate velocities and derivatives:
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value
     //Assume streamfunction stored as first nodal value and v as third.
     interpolated_v += this->nodal_value(l,2)*psif[l];
     //Loop over directions
     for(unsigned i=0;i<Dim;i++) interpolated_x[i] += this->nodal_position(l,i)*psif[l];
    }
   
   //Now add to the appropriate equations
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
      local_eqn = s_local_eqn(l);
      /*IF it's not a boundary condition*/
      if(local_eqn >= 0)
      {
       //Add the user-defined traction terms
       //Essentially (v*boundary)
       residuals[local_eqn] += multiplier*interpolated_v*testf[l]*interpolated_x[0]*Alpha*W;  
	  
               //CALCULATE THE JACOBIAN
                if(flag)
                 {
		   /// No Jacobian terms
    
		 } /*End of Jacobian calculation*/

        } //end of if not boundary condition statement

    } //End of loop over shape functions 

  } //End of loop over integration points
 
} //End of fill_in_generic_residual_contribution

} //End of namespace oomph

#endif
