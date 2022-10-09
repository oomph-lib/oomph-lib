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
// Extra elements for the ST problem

// C++ includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio> 
#include <complex>
#include <set>
 
// The oomphlib headers   
#include "generic.h"
#include "navier_stokes.h"
 
using namespace oomph;

//=======================================================================
/// Point element that is used to set the bubble pressure
//=======================================================================
class FixSpineHeightElement : public virtual SpineElement<PointElement>
{
private:
 
 /// Pointer to the desired value of the spine height
 double *Height_pt;
 
 /// The local eqn number for the pressure that has been traded for
 /// the volume constraint
 int Ptraded_local_eqn;

 /// The Data that contains the traded pressure
 Data *Ptraded_data_pt;
 

 protected:
 
 /// Calculate the geometric shape functions at local coordinate s. 
 void shape(const Vector<double> &s, Shape &psi) const {psi[0] = 1.0;}

 /// Calculate the geometric shape functions and
 /// derivatives w.r.t. local coordinates at local coordinate s
 void dshape_local(const Vector<double> &s, Shape &psi,
                           DShape &dpsids) const
  { 
   psi[0] = 1.0;
   dpsids(0,0) = 0.0;
  }


public:
 
 /// Constructor, there are no internal values. The pointer to the 
 /// element's (single) spine is set on construction
 FixSpineHeightElement(SpineNode* const &spine_node_pt) : 
  SpineElement<PointElement>() 
  {
   // Initialise pointer to prescribed spine height
   Height_pt=0;
   // Initialise pointer to "traded" pressure Data.
   Ptraded_data_pt=0;
   //Add the node to the element's node pointer
   //this->set_n_node(1);
   node_pt(0) = spine_node_pt;
   //Add the node's spine to the element
  }

 /// Access function to the prescribed spine height
 double*& height_pt() {return Height_pt;}

 /// Return the spatial (Eulerian) dimension of the element
 unsigned dim() const {return 0;} 

 /// Return the spatial dimension of local node n
 unsigned required_ndim() const {return 1;}

 /// Return the 'linear order' = number of nodal points along the edge
 unsigned nnode_1d() const {return 1;}

 /// Calculate the residuals
 void fill_in_contribution_to_residuals(Vector<double> &residuals);
 
 /// Calculate the residuals and the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                                   DenseMatrix<double> &jacobian);

 /// Assign the local equation numbers and their coincidence with
 /// the global values
 void assign_additional_local_eqn_numbers();

 /// Set the Data that contains the single pressure value
 /// that is "traded" for the volume constraint.
 /// The Data item must only contain a single value!
 void set_traded_pressure_data(Data* traded_pressure_data_pt)
  {
#ifdef PARANOID
   if (traded_pressure_data_pt->nvalue()!=1)
    {
     std::ostringstream error_stream;
     error_stream
      << "Traded pressure Data must only contain a single value!\n"
      << "This one contains " << traded_pressure_data_pt->nvalue()<< "\n";

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
     // Store pointer explicitly
     Ptraded_data_pt=traded_pressure_data_pt;
     // Add to the element's external data so it gets included
     // in the black-box local equation numbering scheme
     add_external_data(traded_pressure_data_pt);
  }
 
 /// Overload the output function
 void output(std::ostream &outfile) { }

 /// Output function: x,y,[z],u,v,[w],p in tecplot format
 void output(std::ostream &outfile, const unsigned &Np) { }

 /// Overload the self test
 //unsigned self_test() {return GeneralisedElement::self_test();}
}; 

//==========================================================================
/// Residuals for the spine-based volumetric constraint (point) element
//==========================================================================
void FixSpineHeightElement::
fill_in_contribution_to_residuals(Vector<double>&residuals)
{
#ifdef PARANOID
 if (Height_pt==0)
  {
   throw OomphLibError(
    "Please set the pointer to the prescribed height",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Integer for the local equation number
 int local_eqn = Ptraded_local_eqn;

 //Set up the contribution to the global conservation equation,
 //if we have a degree of freedom to trade for it
 if(local_eqn >= 0)
  {
   residuals[local_eqn] =
    static_cast<SpineNode*>(node_pt(0))->spine_pt()->height() - *Height_pt;
  }
}

//=========================================================================
/// Calculate the residual vector and the Jacobian matrix
/// for the spine-based volumetric constraint element
//=========================================================================
void FixSpineHeightElement::
fill_in_contribution_to_jacobian(Vector<double> &residuals, 
                             DenseMatrix<double> &jacobian)
{
 //Get the residuals
 fill_in_contribution_to_residuals(residuals);

 //Call the generic finite difference routine to handle the spine variables
 //contributions
 fill_in_jacobian_from_geometric_data(jacobian);
}

//==========================================================================
/// Assign the local equation numbers and their coincidence with the
/// global ones.
//==========================================================================
void FixSpineHeightElement::assign_additional_local_eqn_numbers()
{
 //If there is no "traded" pressure data assigned, set 
 //Ptraded_local_equation to -1
 if(Ptraded_data_pt==0)
  {
   Ptraded_local_eqn = -1;
  }
 //Otherwise, copy its local equation number across from the generic
 //equation numbering scheme -- here we're relying on the fact that 
 //the relevant external Data item only stores a single value.
 else
  {
#ifdef PARANOID
   if (external_data_pt(0)->nvalue()!=1)
    {
     std::cout<< "The external Data item that stores the traded pressure in "
         << std::endl;
     std::cout<<"SpineVolumeConstraintPointElement should only have a single"<<std::endl;
     std::cout<<"value but is has " 
         << external_data_pt(0)->nvalue() 
         << std::endl;
     exit(1);//assert(false);
    }
#endif
   Ptraded_local_eqn = external_local_eqn(0,0);
  }
 
 //Local equation numbering for the single spine is done automatically
 //in the underlying SpineElement.

}



//======================================================================
/// Special Face element used to calculate the additional inlet 
/// velocities and tractions when gravity has a component in the axial 
/// direction
//=======================================================================
template <class ELEMENT>
class SpineGravityTractionElement : 
 public virtual SpineElement<FaceGeometry<ELEMENT> >, 
 public virtual FaceElement
{
 /// The highest dimension of the problem 
 unsigned Dim;

 /// Pointer to the global Reynold number divided by the Froude number
 double *ReInvFr_pt;

 /// Pointer to global gravity Vector
 Vector<double> *G_pt;

 /// Pointer to the viscosity ratio (relative to the 
 /// viscosity used in the definition of the Reynolds number)
 double *Viscosity_Ratio_pt;
 
 /// Pointer to the density ratio (relative to the density used in the 
 /// definition of the Reynolds number)
 double *Density_Ratio_pt;

protected:

 /// Array to hold local eqn number information for veloc: 
 /// U_local_eqn(jnod,i) = local equation number or < 0 if pinned
 DenseMatrix<int> U_local_eqn;

 /// Array to hold the local eqn number information for the
 /// external data (other nodes in the bulk element)
 DenseMatrix<int> External_u_local_eqn;

 /// Vector to keep track of the external data associated
 /// with each bulk node
 Vector<unsigned> External_node;

//BEGIN, bit fo rconsidering Ca as an unknown

 /// The Data that contains the capillary number is stored
 /// as external Data for the element. Which external Data item is it?
 unsigned External_data_number_of_invca;

 /// Pointer to the Data item that stores the capillary number
 Data* invca_data_pt;

// Pointer to the BOnd number. This is neccesary because in the equations
// the bond number always appears as ReSt
 double* bond_pt;

// Return the equation number corresponding to 1ovCa (which is the fixing of the flow rate)
int invca_local_eqn()
{
 return external_local_eqn(External_data_number_of_invca,0);

}
//END, bit fo rconsidering Ca as an unknown


public:


 /// Constructor, which takes a "bulk" element and the value of the index
 /// and its limit
 SpineGravityTractionElement(FiniteElement *element_pt, 
                             int face_index) :
  SpineElement<FaceGeometry<ELEMENT> >(), FaceElement()
  { 
   
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_data from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
   //Set the dimension from the dimension of the first node
   Dim = node_pt(0)->ndim();

   //Initializa pointer 
   invca_data_pt = 0;

   //Set the Physical values from the bulk elemenet
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   this->ReInvFr_pt = cast_element_pt->re_invfr_pt();
   this->G_pt = cast_element_pt->g_pt();
   this->Viscosity_Ratio_pt = cast_element_pt->viscosity_ratio_pt();
   this->Density_Ratio_pt = cast_element_pt->density_ratio_pt();

 
   //Hijack the nodes in the bulk element in the axial coordinate
   unsigned n_node = this->nnode();
   
   for(unsigned m=0;m<n_node;m++)
    {
     delete cast_element_pt->hijack_nodal_value(bulk_node_number(m),1);
    }
   
   //The other nodes of the bulk element must be external data, because
   //they can affect the derivatives that are used in this element
   //Loop over the nodes of the parent element
   unsigned n_node_parent = cast_element_pt->nnode();
   for(unsigned n=0;n<n_node_parent;n++)
    {
     bool external=true;
     //Loop over the face nodes
     for(unsigned m=0;m<n_node;m++)
      {
       //If the parent's node is  one of the face nodes continue
       if(n == this->bulk_node_number(m)) {external=false;}
      }

     //If it's external data add it
     if(external) 
      {
       this->add_external_data(cast_element_pt->node_pt(n));
       External_node.push_back(n);
      }
    }

   //Now add the spines of the bulk elemnt as external data, which 
   //we will finite difference
   
   //Set of unique geometric data that is used to update the bulk,
   //but is not used to update the face
   std::set<Data*> unique_additional_geom_data;
   //Get all the geometric data for this (bulk) element
   cast_element_pt->assemble_set_of_all_geometric_data(
    unique_additional_geom_data);
   
   //Now assemble the set of geometric data for the face element
   std::set<Data*> unique_face_geom_data_pt;
   this->assemble_set_of_all_geometric_data(unique_face_geom_data_pt);
   //Erase the face geometric data from the additional data
   for(std::set<Data*>::iterator it=unique_face_geom_data_pt.begin();
       it!=unique_face_geom_data_pt.end();++it)
    {unique_additional_geom_data.erase(*it);}

   //Finally add all unique additional data as geometric data
   /*for(std::set<Data*>::iterator it = unique_additional_geom_data.begin();
       it!= unique_additional_geom_data.end();++it)
    {
     this->add_external_data(*it);
     }*/
   this->identify_geometric_data(unique_additional_geom_data);
  }

 
 //Overload the function Node Update 
 //so that it all the nodes of the parent element are updated.
 void node_update()
  {
   this->bulk_element_pt()->node_update();
  }


 /// Return the number of external velocity data
 unsigned nexternal_u_data() {return External_node.size();}
 
 /// Viscosity ratio for element: Element's viscosity relative to the 
 /// viscosity used in the definition of the Reynolds number
 const double &viscosity_ratio() const {return *Viscosity_Ratio_pt;}

 /// Pointer to Viscosity Ratio
 double* &viscosity_ratio_pt() {return Viscosity_Ratio_pt;}

 /// Density ratio for element: Element's density relative to the 
 ///  viscosity used in the definition of the Reynolds number
 const double &density_ratio() const {return *Density_Ratio_pt;}

 /// Pointer to Density ratio
 double* &density_ratio_pt() {return Density_Ratio_pt;}

 /// Pointer to Reynolds number divided by Froude number
 double* &re_invfr_pt() {return ReInvFr_pt;}

 /// Return the value of the Re/Fr number
 const double &re_invfr() const  {return *ReInvFr_pt;}

 /// Vector of gravitational components
 const Vector<double> &g() const {return *G_pt;}

 /// Pointer to Vector of gravitational components
 Vector<double>* &g_pt() {return G_pt;}

 /// Calculate the flow across the element
 double flow()
  {
   //Storage for the normal flux
   double flux = 0.0;
   
   //Set the value of n_intpt
   const unsigned n_intpt = integral_pt()->nweight();
 
   //Set the Vector to hold local coordinates
   Vector<double> s(Dim-1);
   
   //Set the Vector to hold the local coordinates of the parent element
   Vector<double> s_parent(Dim);

   //Get a pointer to the parent element
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s and s_parent
     for(unsigned i=0;i<(Dim-1);i++)  {s[i] = integral_pt()->knot(ipt,i);}
     
     //Get the local coordinate in the bulk
     this->get_local_coordinate_in_bulk(s,s_parent);
     
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Get the local jacobian for the FaceElement
     double J = J_eulerian(s);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Get the velocity from the parent
     Vector<double> interpolated_u(Dim,0.0);
     bulk_el_pt->interpolated_u_nst(s_parent,interpolated_u);

     //Get the outer unit normal
     Vector<double> normal(Dim,0);
     outer_unit_normal(s,normal);

     //Find the normal flux
     double normal_flux = 0.0;
     for(unsigned i=0;i<Dim;i++) {normal_flux += interpolated_u[i]*normal[i];}

     //Add to the flux
     flux += normal_flux*W;
     
    } //End of loop over integration points
   //Return the flux
   return flux;
  }

 /// Access function for the velocity. N. B. HEAVY ASSUMPTIONS HERE
 double u(const unsigned &l,const unsigned &i)
  {return this->nodal_value(l,i);}
 
 /// /// Velocity i at local node l at timestep t (t=0: present; 
 /// t>0: previous). SIMILAR HEAVY ASSUMPTIONS
 double u(const unsigned &t, const unsigned &l, 
          const unsigned &i) const
  {return this->nodal_value(t,l,i);}
 
 /// i-th component of du/dt at local node l. 
 double du_dt(const unsigned &l, const unsigned &i) const
  {
   // Get the data's timestepper
	  TimeStepper* time_stepper_pt=node_pt(l)->time_stepper_pt(); // cgj: previously Node_pt[l] -- but (interpreting as FiniteElement::Node_pt as clang does) this is private!

   // Number of timsteps (past & present)
   unsigned n_time = time_stepper_pt->ntstorage();
   
   double dudt=0.0;
   
   //Loop over the timesteps
   if (time_stepper_pt->type()!="Steady")
    {
     for(unsigned t=0;t<n_time;t++)
      {
       dudt+=time_stepper_pt->weight(1,t)*u(t,l,i);
      }
    }
   
   return dudt;
  }

 /// Add the contribution to the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Create a dummy matrix
   DenseMatrix<double> dummy(1);
   //Call the generic residuals function with flag set to 0
   add_generic_residual_contribution(residuals,dummy,0);
  }

 /// This function returns the residuals and the jacobian
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   add_generic_residual_contribution(residuals,jacobian,1);
   //Add the external data
   this->fill_in_jacobian_from_external_by_fd(jacobian);
   //Add the spine contributions
   this->fill_in_jacobian_from_geometric_data(jacobian);
  }


 //----------------------------------------------------------------------
 /// This function returns the residuals for the Navier--Stokes equations; 
 /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
 //----------------------------------------------------------------------
 void add_generic_residual_contribution(Vector<double> &residuals, 
                                        DenseMatrix<double> &jacobian, 
                                        unsigned flag)
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
   
   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
 
   //Set the Vector to hold local coordinates
   Vector<double> s(Dim-1);
   
   //Set the Vector to hold the local coordinates of the parent element
   Vector<double> s_parent(Dim);

   //Get a pointer to the parent element
   ELEMENT* bulk_el_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());

   //Find the number of nodes in the parent element
   unsigned n_node_parent = bulk_el_pt->nnode();
   //Set up memory for shape functions and their derivatives
   Shape psif_parent(n_node_parent), testf_parent(n_node_parent);
   DShape dpsifdx_parent(n_node_parent,Dim), 
    dtestfdx_parent(n_node_parent,Dim);

   //Storage for the local equation number
   int local_eqn=0, local_unknown=0;

   //Get the Physical Variable
   double ReInvFr = re_invfr()*density_ratio();
   double Viscosity_Ratio = viscosity_ratio();
   Vector<double> G = g();

   
   unsigned n_external =  this->nexternal_u_data();
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s and s_parent
     for(unsigned i=0;i<(Dim-1);i++)  {s[i] = integral_pt()->knot(ipt,i);}

     //Get the local coordinate in the bulk
     this->get_local_coordinate_in_bulk(s,s_parent);
     
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
   
     //Find the shape functions and derivatives of the parent
     (void)bulk_el_pt->
      dshape_eulerian(s_parent,psif_parent,dpsifdx_parent);
     
     //Get the local jacobian for the FaceElement
     double J = J_eulerian(s);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Need only to find the velocity derivatives
     DenseMatrix<double> interpolated_dudx(Dim,Dim,0.0);
     
     Vector<double> interpolated_u(Dim,0.0);

     //Calculate velocities and derivatives
     for(unsigned l=0;l<n_node_parent;l++) 
      {
       //Loop over velocity components
       for(unsigned i=0;i<Dim;i++)
        {
         for(unsigned j=0;j<Dim;j++)
          {
           interpolated_dudx(i,j) += 
            bulk_el_pt->u_nst(l,i)*dpsifdx_parent(l,j);
          }
         dtestfdx_parent(l,i) = dpsifdx_parent(l,i);
         interpolated_u[i] += bulk_el_pt->u_nst(l,i)*psif_parent(l);
        }
       //Set the test functions to be the same as the shape functions
       testf_parent[l] = psif_parent[l];
      }

     //Storage for the outer unit normal
     Vector<double> normal(Dim);
     outer_unit_normal(s,normal);
     
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //Do the first (axial) velocity component -- Poisson equation
       {
        unsigned i=1;

        local_eqn = U_local_eqn(l,i);
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          //Add the gravitational body force term
          residuals[local_eqn] += 
           ReInvFr*testf_parent[bulk_node_number(l)]*G[i]*W;
          
          //Add in the equation term (equation constant traction in the axial direction d(Tyy+p)/dy =0)  
          residuals[local_eqn] -= Viscosity_Ratio*W* (
           ( interpolated_dudx(i,0) /*+ interpolated_dudx(0,i)*/ )
             *dtestfdx_parent(bulk_node_number(l),0) +           
           ( interpolated_dudx(i,2) /*+ interpolated_dudx(2,i) */)
            *dtestfdx_parent(bulk_node_number(l),2) );

          //Now add the jacobian terms
          if(flag)
           {
          
            //Loop over all nodes again
            for(unsigned l2=0;l2<n_node;l2++)
             {

              /*   {
               unsigned i2=0;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(bulk_node_number(l2),i)*
                                   dtestfdx_parent(bulk_node_number(l),0))*W;
                }
                } */
              
              //Poisson equation term
              {
               unsigned i2=1;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(bulk_node_number(l2),0)*
                                   dtestfdx_parent(bulk_node_number(l),0)  +
                                   dpsifdx_parent(bulk_node_number(l2),2)*
                                   dtestfdx_parent(bulk_node_number(l),2) )*W;
                }
              }

              /* {
               unsigned i2=2;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(bulk_node_number(l2),i)*
                                   dtestfdx_parent(bulk_node_number(l),2))*W;
                }
                }*/


             }
            
            //Loop over external data
            for(unsigned l2=0;l2<n_external;l2++)
             {
              /* {
               unsigned i2=0;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(External_node[l2],i)*
                                   dtestfdx_parent(bulk_node_number(l),0))*W;
                }
                }*/
           
              //Poisson equation term 
              {
               unsigned i2=1;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(External_node[l2],0)*
                                   dtestfdx_parent(bulk_node_number(l),0)  +
                                   dpsifdx_parent(External_node[l2],2)*
                                   dtestfdx_parent(bulk_node_number(l),2)  )*W;
                }
              }
   
              /*{
               unsigned i2=2;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(External_node[l2],i)*
                                   dtestfdx_parent(bulk_node_number(l),2))*W;
                }
                }*/
             }
           } //End of if flag
         }
        
       }

       //Now do the other (traction) components
       /* for(unsigned i =0;i<3;i+=2)
       {
        local_eqn = U_local_eqn(l,i);
        //If it's not a boundary condition
        if(local_eqn >= 0)
         {
          for(unsigned k=0;k<Dim;k++)
           {
            residuals[local_eqn] += Viscosity_Ratio*
             (interpolated_dudx(i,k) + interpolated_dudx(k,i))*normal[k]
             *testf_parent(bulk_node_number(l))*W;
           }

          //Now add the jacobian terms
          if(flag)
           {
            //Loop over all nodes again
            for(unsigned l2=0;l2<n_node;l2++)
             {
              for(unsigned i2=0;i2<Dim;i2++)     
              {      
               if(i2!=i) 
               {
                 // unsigned i2=0;
                 local_unknown = U_local_eqn(l2,i2);
                 if(local_unknown >= 0)
                  {
                    jacobian(local_eqn,local_unknown) +=
                    Viscosity_Ratio*dpsifdx_parent(bulk_node_number(l2),i)*
                    normal[i2]*testf_parent(bulk_node_number(l))*W;
                  }
                }
                 else
                 {
                 // unsigned i2=1;
                  local_unknown = U_local_eqn(l2,i2);
                  if(local_unknown >= 0)
                  {
                     for(unsigned k=0;k<Dim;k++)
                     {
                       jacobian(local_eqn,local_unknown) +=
                       Viscosity_Ratio*dpsifdx_parent(bulk_node_number(l2),k)*
                       normal[k]*testf_parent(bulk_node_number(l))*W;
                     }
                 
                     jacobian(local_eqn,local_unknown) +=
                     Viscosity_Ratio*dpsifdx_parent(bulk_node_number(l2),i)*
                     normal[i2]*testf_parent(bulk_node_number(l))*W;
                  }
                }
              } 
             }


            //Loop over all external data
            for(unsigned l2=0;l2<n_external;l2++)
             {
              for(unsigned i2=0;i2<Dim;i2++)
              {
               if (i2!=i)
                {
                 //unsigned i2=0;
                 local_unknown = External_u_local_eqn(l2,i2);
                 if(local_unknown >= 0)
                  {
                   jacobian(local_eqn,local_unknown) +=
                    Viscosity_Ratio*dpsifdx_parent(External_node[l2],i)*
                    normal[i2]*testf_parent(bulk_node_number(l))*W;
                  }
                }
                else
                {
                   //unsigned i2=1;
                 local_unknown = External_u_local_eqn(l2,i2);
                 if(local_unknown >= 0)
                  {
                   for(unsigned k=0;k<Dim;k++)
                    {
                     jacobian(local_eqn,local_unknown) +=
                      Viscosity_Ratio*dpsifdx_parent(External_node[l2],k)*
                      normal[k]*testf_parent(bulk_node_number(l))*W;
                    }
                 
                   jacobian(local_eqn,local_unknown) +=
                    Viscosity_Ratio*dpsifdx_parent(External_node[l2],i)*
                    normal[i2]*testf_parent(bulk_node_number(l))*W;
                  }
                }
               }
               }
           }
         }
         } */ 
      }
    }

  }
 

/// Define the local equation numbering schemes
 void assign_additional_local_eqn_numbers()
  {
   //Get number of nodes
   unsigned n_node = nnode();
   //Resize the equation counters
   U_local_eqn.resize(n_node,Dim);
   
   //Loop over the nodes
   for(unsigned i=0;i<n_node;i++)
    {
     //Loop over the nodal values
     for(unsigned j=0;j<Dim;j++)
      {
       U_local_eqn(i,j) = this->nodal_local_eqn(i,j);
      }
    }

   //Now sort out the external degrees of freedom
   //Find the number of external degrees of freedom
  //ACHTUNG !!!!!!
//   unsigned n_external = nexternal_data();
 unsigned n_external =  External_node.size();

   //Resize the external degree of freedom counter
   External_u_local_eqn.resize(n_external,Dim);
   //Loop over the external degrees of freedom
   for(unsigned e=0;e<n_external;e++)
    {
     //Loop over the velocity values
     for(unsigned i=0;i<Dim;i++)
      {
       External_u_local_eqn(e,i) = external_local_eqn(e,i);
      }
    }
  }

 /// Overload the output function
 void output(std::ostream &outfile) { }

/// Output function: x,y,[z],u,v,[w],p in tecplot format
void output(std::ostream &outfile, const unsigned &Np) 
 {
  SpineElement<FaceGeometry<ELEMENT> >::output(outfile,Np);
 }

}; 
