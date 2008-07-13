//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
//Non-inline member functions for discontinuous galerkin elements

//oomph-lib includes
#include "dg_elements.h"
#include "shape.h"
#include <iomanip>

namespace oomph
{

//===================================================================
///Find pointers to neighbouring faces and the local coordinates
///in those faces that correspond to the integration points in the
///present face.
///This is achieved by moving up to the bulk element and thence
///the mesh which MUST have implemented a neighbour finding scheme
//==================================================================
void DGFaceElement::setup_neighbour_info()
{
 //Cache the pointer to the bulk element
 DGElement* const bulk_element_pt = 
  dynamic_cast<DGElement*>(this->bulk_element_pt());
 
 //Find the number of points in the integration scheme
 const unsigned n_intpt = integral_pt()->nweight();
 //Resize the storage in the element
 Neighbour_face_pt.resize(n_intpt);
 Neighbour_local_coordinate.resize(n_intpt);
 
 //Get the dimension of the present element
 const unsigned el_dim = this->dim();
 //Local coordinate in the face element
 Vector<double> s(el_dim);
 
 //Get the dimension of the bulk element
 const unsigned n_dim = bulk_element_pt->dim(); 
 //Local coordinate in the bulk element
 Vector<double> s_bulk(n_dim);
 
 //Now loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the local coordinate of the integration points
   for(unsigned i=0;i<el_dim;i++)
    {s[i] = integral_pt()->knot(ipt,i);}
   
   //Now get the bulk coordinate
   this->get_local_coordinate_in_bulk(s,s_bulk);

   //Now we find the neighbouring face via the bulk element's member
   //function which calls the Mesh's member function
   bulk_element_pt->
    get_neighbouring_face_and_local_coordinate(this->face_index(),
     s_bulk,Neighbour_face_pt[ipt],Neighbour_local_coordinate[ipt]);
  }
}


//======================================================================
//Report the global coordinates corresponding to the integration points
//and the corresponding coordinates in the neighbouring faces
//======================================================================
void DGFaceElement::report_info()
{
 //Find the number of nodes
 const unsigned n_node = this->nnode();
 //Allocate storage for the shape functions
 Shape psi(n_node);
 
 //Find the dimension of the problem
 const unsigned dim = this->nodal_dimension();
 //Storage for local coordinates in this face and its neighbour
 Vector<double> x(dim), face_x(dim);
 
 //Find the dimension of the element
 const unsigned el_dim = this->dim();
 //Storage for local coordinate
 Vector<double> s(el_dim);
 
 //Calculate the number of integration points from the array
 const unsigned n_intpt = this->integral_pt()->nweight();
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Find the local coordinate at the knot point
   for(unsigned i=0;i<el_dim;i++) 
    {s[i] = this->integral_pt()->knot(ipt,i);}
   //Get the value of the global coordinate in the present face
   this->interpolated_x(s,x);
   
   //Get the value of x in the neighbouring face
   Neighbour_face_pt[ipt]->interpolated_x(Neighbour_local_coordinate[ipt],
                                          face_x);
   
   //Let's compare
   std::cout <<  "In Face                   In Neighbour\n";
   for(unsigned i=0;i<dim;i++)
    {
     if(i==0) {std::cout << "(";}
     else {std::cout << ", ";}
     std::cout << std::setw(5) << std::left <<  x[i];
    }
   std::cout << ")";
   
   std::cout << "                   ";
   
   for(unsigned i=0;i<dim;i++)
    {
     if(i==0) {std::cout << "(";}
     else {std::cout << ", ";}
     std::cout << std::setw(5) << std::left << face_x[i];
    }
   std::cout << ")";
   std::cout << std::endl;
  }
}


//=====================================================================
///Return the interpolated values of the unknown fluxes
//=====================================================================
void DGFaceElement::interpolated_u(const Vector<double> &s, 
                                   Vector<double> &u)
{
 //Find the number of nodes
 const unsigned n_node = nnode();
 //If there are no nodes then return immediately
 if(n_node==0) {return;}
 
 //Get the shape functions at the local coordinate
 Shape psi(n_node);
 this->shape(s,psi);
 
 //Find the number of fluxes
 const unsigned n_flux = this->required_nflux();
 
 //Find the indices at which the local fluxes are stored
 Vector<unsigned> flux_nodal_index(n_flux);
 for(unsigned i=0;i<n_flux;i++)
  {
   flux_nodal_index[i] = this->flux_index(i);
  }
 //Initialise the fluxes to zero
 for(unsigned i=0;i<n_flux;i++) {u[i] = 0.0;}
 
 //Loop over the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   const double psi_ = psi[n];
   for(unsigned i=0;i<n_flux;i++)
    {
     u[i] += this->node_pt(n)->value(flux_nodal_index[i])*psi_;
    }
  }
}

//====================================================================
///Calculate the numerical flux at the knot point ipt. This is the 
///most general interface than can be overloaded if desired. The shape
///functions at the knot point will be passed into this function.
//====================================================================
void DGFaceElement::numerical_flux_at_knot(const unsigned &ipt,
                                           const Shape &psi,
                                           Vector<double> &flux)
{
 //Find the number of nodes
 const unsigned n_node = this->nnode();
 //Find the nodal dimension
 const unsigned nodal_dim = this->nodal_dimension();
 //Number of fluxes
 const unsigned n_flux = this->required_nflux();
 //Find the indices at which the local fluxes are stored
 Vector<unsigned> flux_nodal_index(n_flux);
 for(unsigned i=0;i<n_flux;i++)
  {
   flux_nodal_index[i] = this->flux_index(i);
  }
 
 //Calculate the local unknowns
 Vector<double> interpolated_u(n_flux,0.0);
 
 //Loop over the shape functions
 for(unsigned l=0;l<n_node;l++)
  {
   //Cache the shape functions
   const double psi_ = psi(l);
   //Loop over the fluxes
   for(unsigned i=0;i<n_flux;i++)
    {
     //Calculate the velocity from the most recent nodal values
     interpolated_u[i] += this->nodal_value(l,flux_nodal_index[i])*psi_;
    }
  }
 
 //Now calculate the outer unit normal Vector
 Vector<double> interpolated_n(nodal_dim);
 this->outer_unit_normal(ipt,interpolated_n);
 
 //Get the pointer to the neighbour
 DGFaceElement* neighbour_element_pt =
  dynamic_cast<DGFaceElement*>(Neighbour_face_pt[ipt]);
 
 //Get the neighbour's fluxes 
 Vector<double> interpolated_u_neigh(n_flux); 
 
 neighbour_element_pt->
  interpolated_u(Neighbour_local_coordinate[ipt],
                 interpolated_u_neigh);

 //Call the "standard" numerical flux function
 this->numerical_flux(interpolated_n,interpolated_u,
                      interpolated_u_neigh,flux);
 
}


//===================================================================
///Calculate the integrated (numerical) flux out of the face and add
///it to the residuals vector
//===================================================================
void DGFaceElement::add_flux_contributions(Vector<double> &residuals)
{
 //Find the number of nodes
 const unsigned n_node = nnode();
 //Dimension of the face
 const unsigned el_dim = dim();
 //Storage for the shape functions
 Shape psi(n_node);

 //Number of integration points
 const unsigned n_intpt = this->integral_pt()->nweight();
 //Number of fluxes
 const unsigned n_flux = this->required_nflux();
 //Find the indices at which the local fluxes are stored
 Vector<unsigned> flux_nodal_index(n_flux);
 for(unsigned i=0;i<n_flux;i++)
  {
   flux_nodal_index[i] = this->flux_index(i);
  }

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double W = this->integral_pt()->weight(ipt);
   //Get the shape functions at the knot
   this->shape_at_knot(ipt,psi);
     
   //Calculate the Jacobian 
   //For a point element, it's one
   double J=W;
   //Otherwise calculate for the element
   if(el_dim != 0) {J *= this->J_eulerian_at_knot(ipt);}
   
   //Now calculate the numerical flux
   Vector<double> F(n_flux);
   this->numerical_flux_at_knot(ipt,psi,F);
   
   //Limit if desired here
   
   //Finally  we need to assemble the appropriate contributions
   //to the residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Loop over the fluxes
     for(unsigned i=0;i<n_flux;i++)
      {
       //Get the local equation number in the bulk
       int local_eqn = dynamic_cast<DGElement*>(this->bulk_element_pt())
        ->nodal_local_eqn(bulk_node_number(l),flux_nodal_index[i]);
       
       //If it's not a boundary condition
       if(local_eqn >= 0)
        {
         //Add the flux multiplied by the shape function and the jacobian
         residuals[local_eqn] -= psi(l)*F[i]*J;
        }
      }
    }
  }
}

//========================================================================
///\short Function that computes and stores the (inverse) mass matrix
//========================================================================
void DGElement::pre_compute_mass_matrix()
{
 //Now let's assemble stuff
 const unsigned n_dof = this->ndof();
 //Allocate storage for the local mass matrix (if required)
 if(M_pt==0) {M_pt = new DenseDoubleMatrix;}
 
 //Resize and initialise the vector that will holds the residuals
 Vector<double> dummy(n_dof,0.0);
 
 //Resize the mass matrix
 M_pt->resize(n_dof,n_dof);
 //Initialise the entries to zero
 M_pt->initialise(0.0);
 //Get the local mass matrix and residuals
 this->fill_in_contribution_to_mass_matrix(dummy,*M_pt);
 
 //Now invert the mass matrix it will always be small
 //This can possibly be streamlined (for example in spectral
 //elements the mass matrix is diagonal)
 M_pt->ludecompose();
 
 //The mass matrix has been computed
 Mass_matrix_has_been_computed=true;
}



//============================================================================
///Function that returns the current value of the residuals 
///multiplied by the inverse mass matrix (virtual so that it can be overloaded
///specific elements in which time and memory saving tricks can be applied)
//============================================================================
void DGElement::
get_inverse_mass_matrix_times_residuals(Vector<double> &minv_res)
{
 //Now let's assemble stuff
 const unsigned n_dof = this->ndof();
 //Allocate storage for the local mass matrix (if required)
 if(M_pt==0) {M_pt = new DenseDoubleMatrix;}
 
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
  //Resize the mass matrix
  M_pt->resize(n_dof,n_dof);
  //Initialise the entries to zero
  M_pt->initialise(0.0);
  //Get the local mass matrix and residuals
  this->fill_in_contribution_to_mass_matrix(minv_res,*M_pt);
 
  //Now invert the mass matrix it will always be small
  //This can possibly be streamlined (for example in spectral
  //elements the mass matrix is diagonal)
  M_pt->ludecompose();

  //The mass matrix has been computed
  Mass_matrix_has_been_computed=true;
 }
 
 //Always do the backsubstitution
 M_pt->lubksub(minv_res);
}


void DGElement::get_neighbouring_face_and_local_coordinate(
 const int &face_index,
 const Vector<double> &s, FaceElement* &face_element_pt,
 Vector<double> &s_face)
{
 DG_mesh_pt->neighbour_finder(this,face_index,s,face_element_pt,s_face);
}


double DGMesh::FaceTolerance = 1.0e-10;

}
