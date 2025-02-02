//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Header file for Refineable Streamfunction elements
#ifndef OOMPH_REFINEABLE_POLAR_STREAMFUNCTION
#define OOMPH_REFINEABLE_POLAR_STREAMFUNCTION

//Oomph-lib headers
//Should already be looking in build/include/ for generic.h
#include "generic/refineable_quad_element.h"
#include "generic/error_estimator.h"
#include "streamfunction_elements.h"
#include "polar_streamfunction_traction_elements.h"

namespace oomph
{

/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////

//======================================================================
/// Refineable version of Streamfunction equations
///
///
//======================================================================
class RefineablePolarStreamfunctionEquations : 
public virtual PolarStreamfunctionEquations,
public virtual RefineableElement,
public virtual ElementWithZ2ErrorEstimator
{
  public:

 /// Constructor, simply call other constructors
 RefineablePolarStreamfunctionEquations() : 
 PolarStreamfunctionEquations(),
 RefineableElement(), 
 ElementWithZ2ErrorEstimator() {} 

 /// Broken copy constructor
 RefineablePolarStreamfunctionEquations(const RefineablePolarStreamfunctionEquations& dummy) 
  { 
   BrokenCopy::broken_copy("RefineablePolarStreamfunctionEquations");
  } 
 
 /// Number of 'flux' terms for Z2 error estimation 
 unsigned num_Z2_flux_terms() {return 2;}

 /// Get 'flux' for Z2 error recovery:  Standard flux.from Streamfunction equations
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {this->get_flux(s,flux);}


/// Get the function value u in Vector.
/// Note: Given the generality of the interface (this function
/// is usually called from black-box documentation or interpolation routines),
/// the values Vector sets its own size in here.
void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
 {
  // Set size of Vector: streamfunction, u and v
  values.resize(3);
  
  //Find number of nodes
  unsigned n_node = nnode();
  
  //Local shape function
  Shape psi(n_node);
  
  //Find values of shape function
  shape(s,psi);
  
  //Initialise values 
  values[0] = 0.0;
  values[1] = 0.0;
  values[2] = 0.0;

  //Find the index at which the unknowns are stored
  unsigned s_nodal_index = this->u_index_streamfunction();
  //Find the indices at which the local velocities are stored
  unsigned u_nodal_index[2];
  for(unsigned i=0;i<2;i++) {u_nodal_index[i] = u_index_velocity(i);}

  
  //Loop over the local nodes and sum up the values
  for(unsigned l=0;l<n_node;l++)
   {
    values[0] += this->nodal_value(l,s_nodal_index)*psi[l];
    values[1] += this->nodal_value(l,u_nodal_index[0])*psi[l];
    values[2] += this->nodal_value(l,u_nodal_index[1])*psi[l];
   }
 }

 /// Get the function value u in Vector.
 /// Note: Given the generality of the interface (this function
 /// is usually called from black-box documentation or interpolation routines),
 /// the values Vector sets its own size in here.
void get_interpolated_values(const unsigned& t,
                             const Vector<double>&s, 
                              Vector<double>& values)
  {
   if (t!=0)
    {
     std::string error_message =
      "Time-dependent version of get_interpolated_values() ";
     error_message += "not implemented for this element \n";
     throw 
      OomphLibError(error_message,
                    "RefineablePolarStreamfunctionEquations::get_interpolated_values()",
                    OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     //Make sure that we call this particular object's steady 
     //get_interpolated_values (it could get overloaded lower down)
     RefineablePolarStreamfunctionEquations::get_interpolated_values(s,values);
    }
  }

 
 ///  Further build: Copy source function pointer from father element
 void further_build()
  {
   //Find the father element
   RefineablePolarStreamfunctionEquations* cast_father_element_pt
    = dynamic_cast<RefineablePolarStreamfunctionEquations*>
    (this->father_element_pt());
   
   //Set pointer to alpha
   this->Alpha_pt = cast_father_element_pt->alpha_pt();
  }

  private:


/// Add element's contribution to elemental residual vector and/or 
/// Jacobian matrix 
/// flag=1: compute both
/// flag=0: compute only residual vector
void fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                       DenseMatrix<double> &jacobian, 
                                       unsigned flag); 

};


//======================================================================
/// Refineable version of 2D StreamfunctionElement elements
///
///
//======================================================================
class RefineablePolarStreamfunctionElement : 
public PolarStreamfunctionElement,
public virtual RefineablePolarStreamfunctionEquations,
public virtual RefineableQElement<2>
{
  public:

 /// Constructor, simply call the other constructors 
 RefineablePolarStreamfunctionElement() : 
  RefineableElement(),
  RefineablePolarStreamfunctionEquations(),
  RefineableQElement<2>(),
  PolarStreamfunctionElement() {} 


 /// Broken copy constructor
 RefineablePolarStreamfunctionElement(const RefineablePolarStreamfunctionElement& 
                           dummy) 
  { 
   BrokenCopy::broken_copy("RefineableQuadStreamfunctionElement");
  } 
 
 /// Number of continuously interpolated values: 3
 unsigned ncont_interpolated_values() const {return 3;}

 /// Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return PolarStreamfunctionElement::nvertex_node();}

 /// Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return PolarStreamfunctionElement::vertex_node_pt(j);}

 /// Rebuild from sons: empty
 void rebuild_from_sons(Mesh* &mesh_pt) {}

 /// Order of recovery shape functions for Z2 error estimation:
 /// Same order as shape functions.
 unsigned nrecovery_order() {return (2);}

 ///  Perform additional hanging node procedures for variables
 /// that are not interpolated by all nodes. Empty.
 void further_setup_hanging_nodes(){}

};

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//=======================================================================
/// Face geometry for the RefineableQuadStreamfunctionElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<>
class FaceGeometry<RefineablePolarStreamfunctionElement > : 
public virtual FaceGeometry<PolarStreamfunctionElement >
{
  public:
 FaceGeometry() : FaceGeometry<PolarStreamfunctionElement >() {}
};

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//========================================================================
/// Add element's contribution to the elemental 
/// residual vector and/or Jacobian matrix.
/// flag=1: compute both
/// flag=0: compute only residual vector
//========================================================================
void RefineablePolarStreamfunctionEquations::
fill_in_generic_residual_contribution(Vector<double> &residuals, 
                                  DenseMatrix<double> &jacobian, 
                                  unsigned flag)
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Get Alpha
 const double Alpha = alpha();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,2), dtestdx(n_node,2);

 //Indicies at which the unknowns are stored
 const unsigned s_nodal_index = u_index_streamfunction();

 //Find the indices at which the local velocities are stored
 unsigned u_nodal_index[2];
 for(unsigned i=0;i<2;i++) {u_nodal_index[i] = u_index_velocity(i);}
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 // Local storage for pointers to hang_info objects
 HangInfo *hang_info_pt=0, *hang_info2_pt=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {    
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                        test,dtestdx);
       
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Allocate and initialise to zero
   Vector<double> interpolated_x(2,0.0);
   double interpolated_s=0.0;
   Vector<double> interpolated_dsdx(2,0.0);
   Vector<double> interpolated_u(2);
   DenseMatrix<double> interpolated_dudx(2,2);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the poisson unknown
     double s_value = nodal_value(l,s_nodal_index);
     interpolated_s += s_value*psi(l);
     // Loop over directions1
     for(unsigned i=0;i<2;i++)
      {
       interpolated_x[i] += nodal_position(l,i)*psi(l);
       interpolated_dsdx[i] += s_value*dpsidx(l,i);
       double u_value = this->nodal_value(l,u_nodal_index[i]);
       interpolated_u[i] += u_value*psi(l);
       // Loop over directions2
       for(unsigned j=0;j<2;j++)
        {                               
         interpolated_dudx(i,j) += u_value*dpsidx(l,j);
        }
      }
    }

   // Assemble residuals and Jacobian
   //--------------------------------
 
 // Loop over the nodes for the test functions 
 for(unsigned l=0;l<n_node;l++)
  {
   //Local variables used to store the number of master nodes and the
   //weight associated with the shape function if the node is hanging
   unsigned n_master=1; double hang_weight=1.0;
   //Local bool (is the node hanging)
   bool is_node_hanging = this->node_pt(l)->is_hanging();

   //If the node is hanging, get the number of master nodes
   if(is_node_hanging)
    {
     hang_info_pt = this->node_pt(l)->hanging_pt();
     n_master = hang_info_pt->nmaster(); 
    }
   //Otherwise there is just one master node, the node itself
   else
    {
     n_master = 1; 
    }
   
   //Loop over the master nodes
   for(unsigned m=0;m<n_master;m++)
    {
     //Get the local equation number and hang_weight
     //If the node is hanging
     if(is_node_hanging)
      {
       //Read out the local equation number from the m-th master node
       local_eqn =  this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                   s_nodal_index);
       //Read out the weight from the master node
       hang_weight = hang_info_pt->master_weight(m); 
      }
     //If the node is not hanging
     else
      {
       //The local equation number comes from the node itself
       local_eqn = this->nodal_local_eqn(l,s_nodal_index);
       //The hang weight is one
       hang_weight = 1.0;
      }
     
     //If the nodal equation is not a boundary condition
     if(local_eqn >= 0)
      {       
	// Laplacian of the streamfunction (integrated by parts)
	residuals[local_eqn] += interpolated_dsdx[0]*dtestdx(l,0)*interpolated_x[0]*Alpha*W*hang_weight;
	residuals[local_eqn] += (1./(interpolated_x[0]*Alpha))*interpolated_dsdx[1]
	                       *(1./(interpolated_x[0]*Alpha))*dtestdx(l,1)
	  *interpolated_x[0]*Alpha*W*hang_weight;

	// Should equal the vorticity
	residuals[local_eqn] -= (interpolated_dudx(1,0)+(interpolated_u[1]/interpolated_x[0]))*test(l)*interpolated_x[0]*Alpha*W*hang_weight;
	residuals[local_eqn] += (1./(interpolated_x[0]*Alpha))*interpolated_dudx(0,1)*test(l)*interpolated_x[0]*Alpha*W*hang_weight;

       // Calculate the Jacobian
	if(flag)
        {
         //Local variables to store the number of master nodes
         //and the weights associated with each hanging node
         unsigned n_master2=1; double hang_weight2=1.0;
         //Loop over the nodes for the variables
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           //Local bool (is the node hanging)
           bool is_node2_hanging = this->node_pt(l2)->is_hanging();
           //If the node is hanging, get the number of master nodes
           if(is_node2_hanging)
            {
             hang_info2_pt = this->node_pt(l2)->hanging_pt();
             n_master2 = hang_info2_pt->nmaster();
            }
           //Otherwise there is one master node, the node itself
           else
            {
             n_master2 = 1; 
            }
           
           //Loop over the master nodes
           for(unsigned m2=0;m2<n_master2;m2++)
            {
              //Get the local unknown and weight
              //If the node is hanging
              if(is_node2_hanging)
               {
                //Read out the local unknown from the master node
                local_unknown = 
		  this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
				       s_nodal_index);
		//Read out the hanging weight from the master node
		hang_weight2 = hang_info2_pt->master_weight(m2); 
	       }
	      //If the node is not hanging
	      else
	       {
                //The local unknown number comes from the node
                local_unknown = this->nodal_local_eqn(l2,s_nodal_index);
		//The hang weight is one
		hang_weight2 = 1.0;
	       }

	      //If the unknown is not pinned
	      if(local_unknown >= 0)
               {
               //Add contribution to Elemental Matrix
                 jacobian(local_eqn,local_unknown) 
		   += dpsidx(l2,0)*dtestdx(l,0)*interpolated_x[0]*Alpha*W*hang_weight*hang_weight2;
                 jacobian(local_eqn,local_unknown) 
		   += (1./(interpolated_x[0]*Alpha))*dpsidx(l2,1)*(1./(interpolated_x[0]*Alpha))*dtestdx(l,1)
			*interpolated_x[0]*Alpha*W*hang_weight*hang_weight2;
	       }

	     } //End of loop over master nodes

          } //End of loop over nodes
        } //End of Jacobian calculation
       
      } //End of case when residual equation is not pinned

    } //End of loop over master nodes for residuals
  } //End of loop over nodes
 
} // End of loop over integration points
}

}

#endif

