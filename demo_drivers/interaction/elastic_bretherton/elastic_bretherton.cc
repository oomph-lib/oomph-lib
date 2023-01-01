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
//Driver code for the elastic-walled Bretherton (airway-reopening) 
//problem that can also include the effects of transverse gravity.
//For full details see Hazel & Heil Physics of Fluids 20, 092102 (2008)
//The initial setup is rather involved:
// 1) Start from the Bretherton problem with stiff springs (almost rigid walls)
// 2) Change the outlet boundary conditions to impose a given flux
//    and decrease the flux until it is the correct value (-2.0) for the
//   airway-reopening problem
// 3) Crank up the Capillary number to 0.1
// 4) Crank down the wall stiffness until we reach the production value
//    0.5e-7

// Once we have got to the appropriate initial guess then we restart
// from the dumped solution, by re-compiling and running again
// with the restart_flag set to be true.

//In order to get relatively quick convergence for self-test,
//the mesh is relatively coarse compared to the proper production runs

// Note that this is/was production code, so is a little bit messy and
// may be tidied up at some point, by it's equally likely that it won't
// The main reason for inclusion
// is to check that "upgrades" to the library don't break it....
// It also serves as a prototype for complex FSI interaction problems with free
// surfaces.

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
#include "fluid_interface.h"
#include "beam.h"

// The mesh
#include "meshes/bretherton_spine_mesh.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;
 
using namespace oomph;

//======================================================================
/// Namepspace for global parameters
//======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re;

 /// Womersley = Reynolds times Strouhal
 double ReSt; 
 
 /// Product of Reynolds and Froude number
 double ReInvFr;

 /// Reynolds divided by Capillary number
 double ReCa=0.0;

 /// Capillary number
 double Ca;  

 /// Bond number
 double Bo;

 double Pext = 0.0;
 /// Direction of gravity
 Vector<double> G(2);
 
 /// The Poisson ratio
 double Nu = 0.49;

 // The wall thickness
 double H = 5.0e-4;

 // The ratio
 double Gamma = 1.0e-7; 

 // The axial tension
 double T = 0.0;

 /// Spring stiffness
 double Kstiff = 0.5*50.0e-7; //50 times Matthias's Value

 /// Huge stiffness for hitting
 double Kseparation = H;

 /// Huge stiffness for the table
 double Ktable = 0.0001;

 /// Natural width of the open tube
 double Tube_width = 10.0;//2.0;

 /// Position of the table
 double Table_position = 10.78;

 /// Ratio of scales
 double Q;

 /// Rest length of the linear springs
 double Rest_length_linear=1.0;

 /// Do we include the non-linear springs
 bool Non_linear_springs = false;

 /// Pointer to the upper wall
 GeomObject* Upper_wall_pt;

 /// Pointer to the lower wall
 GeomObject* Lower_wall_pt;

 /// Upper wall map
 map<double,pair<GeomObject*,Vector<double> > > upper_map;

 /// Lower wall map
 map<double,pair<GeomObject*,Vector<double> > > lower_map;

 /// Function that prescribes the hydrostatic pressure field at the outlet
 void hydrostatic_pressure(const double &time, const Vector<double> &x, 
                           const Vector<double>& n, Vector<double> &traction)
 {
  traction[0] = -ReInvFr*G[1]*x[1];
  traction[1] = 0.0;
 }

 /// Load function for the wall elements
 void spring_load(const Vector<double> &xi, const Vector<double> &x,
                  const Vector<double> &N, Vector<double> &load)
 {
  load[0] = -Pext*N[0]; 
  load[1] = -Pext*N[1];// + Kstiff*(Rest_length_linear - x[1]);
  
  //Now we wish to penalise the tube if it goes above the natural width
  //Find the position of the lower wall
  pair<GeomObject*,Vector<double> > paired_lower
   = lower_map[xi[0]];
  Vector<double> lower_x(2);
  paired_lower.first->position(paired_lower.second,lower_x);
  
  //Now work out the distance between the walls
  Vector<double> d_vec(2);
  for(unsigned i=0;i<2;++i) {d_vec[i] = lower_x[i] - x[i];}
  const double distance = sqrt(d_vec[0]*d_vec[0] + d_vec[1]*d_vec[1]);

  const double width = 2*Rest_length_linear;
  
  const double strain = (distance - width)/width;
  
  for(unsigned i=0;i<2;i++)
   {
    load[i] += Kstiff*strain*d_vec[i]/distance;
   }

  //Note, it turns out the a linear stiffness of 0.0002 
  //is exactly correct for
  //a tube of wall thickness, whatever, so let's stick with it.
    
  if(Non_linear_springs)
   {
    //Now add on the penalty
    if(strain > (-3.0*H*H/12.0))
     {
      for(unsigned i=0;i<2;i++) 
       {
        load[i] += H*d_vec[i]/distance*strain;
       }
     }
    //Now the non-linear part of the tube law can be approximated
    //by 0.7*sqrt(P-3), see my other notes
    else
     {
      for(unsigned i=0;i<2;i++)
       {
        load[i] -= (1.0/12.0)*H*H*H*(3.0 + 1.98*strain*strain)*
         d_vec[i]/distance;
       }
     }
   }
 }


 /// Load function for the wall elements
 void spring_load_lower(const Vector<double> &xi, const Vector<double> &x,
                  const Vector<double> &N, Vector<double> &load)
 {
  load[0] = Pext*N[0]; 
  load[1] = Pext*N[1];// - Kstiff*(Rest_length_linear + x[1]);
  //Add a table
  double table_strain = Table_position + x[1];
  if(table_strain < 0.0)
   {
    //Add a huge restoring force
    load[1] -= Ktable*table_strain;
   }

    //Now we wish to penalise the tube if it goes above the natural width
    //Find the position of the upper wall
    pair<GeomObject*,Vector<double> > paired_upper
     = upper_map[xi[0]];
    Vector<double> upper_x(2);
    paired_upper.first->position(paired_upper.second,upper_x);
    
    //Now work out the distance between the walls
    Vector<double> d_vec(2);
    for(unsigned i=0;i<2;++i) {d_vec[i] = upper_x[i] - x[i];}
    const double distance = sqrt(d_vec[0]*d_vec[0] + d_vec[1]*d_vec[1]);

    const double width = 2*Rest_length_linear;
    
    const double strain = (distance - width)/width;

    for(unsigned i=0;i<2;i++)
     {
      load[i] += Kstiff*strain*d_vec[i]/distance;
     }

    if(Non_linear_springs)
     {
      //Now add on the penalty
    if(strain > (-3.0*H*H/12.0))
     {
      for(unsigned i=0;i<2;i++) 
       {
        load[i] += H*d_vec[i]/distance*strain;
       }
     }
    else
     {
      for(unsigned i=0;i<2;i++)
       {
        load[i] -= (1.0/12.0)*H*H*H*(3.0 + 1.98*strain*strain)*
         d_vec[i]/distance;
       }
     }
   }
 }

}


namespace No_Slip
{
 /// Function that is used to set and update the no-slip boundary condition
 void no_slip_condition_first(Node* node_pt)
 {
  //Cast the node to a spine node
  SpineNode* spine_node_pt = static_cast<SpineNode*>(node_pt);
  //Get the wall coordinate
  Vector<double> s(1);
  s[0] = spine_node_pt->spine_pt()->geom_parameter(0);
  // Prepare the storage for the derivative
  DenseMatrix<double> drdxi(1,2);
  //Get the derivative
  dynamic_cast<FSIHermiteBeamElement*>
   (spine_node_pt->spine_pt()->geom_object_pt(0))->
   dposition_dlagrangian_at_local_coordinate(s,drdxi);
  
  for(unsigned i=0;i<2;i++) {*node_pt->value_pt(i) = -drdxi(0,i);}
 }

 void no_slip_condition_second(Node* node_pt)
 {
  //Cast the node to a spine node
  SpineNode* spine_node_pt = static_cast<SpineNode*>(node_pt);
  //Get the wall coordinate
  Vector<double> s(1);
  s[0] = spine_node_pt->spine_pt()->geom_parameter(1);
  // Prepare the storage for the derivative
  DenseMatrix<double> drdxi(1,2);
  //Get the derivative
  dynamic_cast<FSIHermiteBeamElement*>
   (spine_node_pt->spine_pt()->geom_object_pt(1))->
   dposition_dlagrangian_at_local_coordinate(s,drdxi);
  
  for(unsigned i=0;i<2;i++) {*node_pt->value_pt(i) = -drdxi(0,i);}
 }
}


//================================================================
/// Function-type-object to perform comparison of elements
//================================================================
class ElementCmp
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x, GeneralisedElement* const &y) 
  const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);

   if((cast_x ==0) || (cast_y==0)) {return 0;}
   else
    {return cast_x->node_pt(0)->x(0) < cast_y->node_pt(0)->x(0);}
  }
};



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

 /// Pointer to an External Data object that represents the
 /// an unknown pressure gradient
 Data* Delta_P_pt;

protected:

 /// Array to hold local eqn number information for veloc: 
 /// U_local_eqn(jnod,i) = local equation number or < 0 if pinned
 DenseMatrix<int> U_local_eqn;

 /// The local equation number for the external data associated
 /// with the unknown pressure gradient
 int Delta_P_local_eqn;

 /// The index in the external data at which the Delta_p data is
 /// stored
 unsigned External_Delta_P_index;

 /// Array to hold the local eqn number information for the
 /// external data (other nodes in the bulk element)
 DenseMatrix<int> External_u_local_eqn;

 /// Vector to keep track of the external data associated
 /// with each bulk node
 Vector<unsigned> External_node;

public:
 
 /// Constructor, which takes a "bulk" element and the value of the index
 /// and its limit
 SpineGravityTractionElement(FiniteElement* const &element_pt, 
                             const int &face_index) :
  SpineElement<FaceGeometry<ELEMENT> >(), FaceElement()
  { 
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_data from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);
   
   //Set the dimension from the dimension of the first node
   Dim = node_pt(0)->ndim();

   //Set the Physical values from the bulk elemenet
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(element_pt);
   this->ReInvFr_pt = cast_element_pt->re_invfr_pt();
   this->G_pt = cast_element_pt->g_pt();
   this->Viscosity_Ratio_pt = cast_element_pt->viscosity_ratio_pt();
   this->Density_Ratio_pt = cast_element_pt->density_ratio_pt();

   //Initialise the external pressure gradient pointer to zero
   Delta_P_pt=0;

   //Hijack the nodes in the bulk element in the axial coordinate
   unsigned n_node = this->nnode();
   for(unsigned m=0;m<n_node;m++)
    {
     delete cast_element_pt->hijack_nodal_value(bulk_node_number(m),0);
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

     //If it's external data add it, but do not finite difference
     if(external) 
      {
       this->add_external_data(cast_element_pt->node_pt(n),false);
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

   //Finally add all unique additional data as external data
   for(std::set<Data*>::iterator it = unique_additional_geom_data.begin();
       it!= unique_additional_geom_data.end();++it)
    {
     this->add_external_data(*it);
    }

  }

 void hijack_all_nodes()
  {
   ELEMENT* cast_element_pt = dynamic_cast<ELEMENT*>(bulk_element_pt());
   //Hijack the nodes in the bulk element in the axial coordinate
   unsigned n_node = this->nnode();
   for(unsigned m=0;m<n_node;m++)
    {
     delete cast_element_pt->hijack_nodal_value(bulk_node_number(m),0);
    }
  }

 /// Access function to the pointer to pressure gradient data
 void set_delta_p_pt(Data* const &delta_p_pt) 
  {
   Delta_P_pt = delta_p_pt;
   //Add the external data and do not finite difference it!
   External_Delta_P_index = add_external_data(delta_p_pt,false);
  }
 
 /// Return the value of the pressure gradient
 const double delta_p() const 
  {
   if(Delta_P_pt==0) {return 0.0;}
   else {return Delta_P_pt->value(0);}
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

 /// Reynolds number divided by Froude number
 const double &re_invfr() {return *ReInvFr_pt;}

 /// Vector of gravitational components
 const Vector<double> &g() const {return *G_pt;}

 /// Pointer to Vector of gravitational components
 Vector<double>* &g_pt() {return G_pt;}

 /// Access function for the velocity. N. B. HEAVY ASSUMPTIONS HERE
 double u(const unsigned &l,const unsigned &i)
  {return nodal_value(l,i);}
 
 /// /// Velocity i at local node l at timestep t (t=0: present; 
 /// t>0: previous). SIMILAR HEAVY ASSUMPTIONS
 double u(const unsigned &t, const unsigned &l, 
          const unsigned &i) const
  {return nodal_value(t,l,i);}
 
 /// i-th component of du/dt at local node l. 
 double du_dt(const unsigned &l, const unsigned &i) const
  {
   // Get the data's timestepper
   TimeStepper* time_stepper_pt=node_pt(l)->time_stepper_pt();

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
   //Add in external data contributions
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
   double Delta_p = delta_p();
   
   unsigned n_external_u = nexternal_u_data();
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s
     for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}

     //Get the values of s_parent
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
        unsigned i=0;

        local_eqn = U_local_eqn(l,i);
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          //Add the pressure gradient
          residuals[local_eqn] += 
           Delta_p*testf_parent[bulk_node_number(l)]*W;
          
          //Add the gravitational body force term
          residuals[local_eqn] += 
           ReInvFr*testf_parent[bulk_node_number(l)]*G[i]*W;

          //Add in the Poisson terms (Only in y(1) direction)
          residuals[local_eqn] -= Viscosity_Ratio*
           (interpolated_dudx(i,1) /*+ interpolated_dudx(1,i)*/)
           *dtestfdx_parent(bulk_node_number(l),1)*W;
          
          //Now add the jacobian terms
          if(flag)
           {
            //Loop over all nodes again
            for(unsigned l2=0;l2<n_node;l2++)
             {
              {
               unsigned i2=0;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(bulk_node_number(l2),1)*
                                   dtestfdx_parent(bulk_node_number(l),1))*W;
                }
              }
              
              /*{
               unsigned i2=1;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(bulk_node_number(l2),i)*
                                   dtestfdx_parent(bulk_node_number(l),1))*W;
                }
                }*/
             }
            
            //Loop over external data
            for(unsigned l2=0;l2<n_external_u;l2++)             {
              {
               unsigned i2=0;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(External_node[l2],1)*
                                   dtestfdx_parent(bulk_node_number(l),1))*W;
                }
              }
              
              /*{
               unsigned i2=1;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) -=
                  Viscosity_Ratio*(dpsifdx_parent(External_node[l2],i)*
                                   dtestfdx_parent(bulk_node_number(l),1))*W;
                }
                }*/
             }

            //Add in the pressure gradient term
            local_unknown = Delta_P_local_eqn;
            if(local_unknown >= 0)
             {
              jacobian(local_eqn,local_unknown) += 
               testf_parent[bulk_node_number(l)]*W;
             }
           } //End of if flag
         }
        
       }

       //Now do the second (traction) component
       /*     {
        unsigned i=1;

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
              {
               unsigned i2=0;
               local_unknown = U_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) +=
                  Viscosity_Ratio*dpsifdx_parent(bulk_node_number(l2),i)*
                  normal[i2]*testf_parent(bulk_node_number(l))*W;
                }
              }
              
              {
               unsigned i2=1;
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


            //Loop over all external data
            for(unsigned l2=0;l2<n_external_u;l2++)
             {
              {
               unsigned i2=0;
               local_unknown = External_u_local_eqn(l2,i2);
               if(local_unknown >= 0)
                {
                 jacobian(local_eqn,local_unknown) +=
                  Viscosity_Ratio*dpsifdx_parent(External_node[l2],i)*
                  normal[i2]*testf_parent(bulk_node_number(l))*W;
                }
              }
              
              {
               unsigned i2=1;
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
         }*/
      }

     //Add the contribution to the pressure gradient constraint,
     //which states that the flux must be 2.0
     local_eqn = Delta_P_local_eqn;
     if(local_eqn >= 0)
      {
       residuals[local_eqn] += interpolated_u[0]*W;

       //Add the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Only the x-velocities at the nodes on the boundary  
           //affect the residuals
           local_unknown = U_local_eqn(l2,0);
           if(local_unknown >= 0)
            {
             //Add the appropriate jacobian term
             jacobian(local_eqn,local_unknown) += 
              psif_parent(bulk_node_number(l2))*W;
            }
          }
        }
      }
    }
  }


 /// Calculate the flux 
 double get_flux()
  {
   //Initialise flux to zero
   double flux = 0.0;

   //Find out how many nodes there are
   unsigned n_node = nnode();
   
   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
 
   //Set the Vector to hold local coordinates
   Vector<double> s(Dim-1);
   
   //Set up memory for shape functions and their derivatives
   Shape psif(n_node);
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Assign values of s and s_parent
     for(unsigned i=0;i<(Dim-1);i++) 
      {s[i] = integral_pt()->knot(ipt,i);}
     
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
   
     //Find the shape functions and derivatives of the parent
     shape(s,psif);
     
     //Get the local jacobian for the FaceElement
     double J = J_eulerian(s);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     Vector<double> interpolated_u(Dim,0.0);

     //Calculate velocities and derivatives
     for(unsigned l=0;l<n_node;l++) 
      {
       //Loop over velocity components
       for(unsigned i=0;i<Dim;i++)
        {
         interpolated_u[i] += u(l,i)*psif(l);
        }
      }


     flux += interpolated_u[0]*W;
    }

   return flux;
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
   unsigned n_external_u = nexternal_u_data();
   //Resize the external degree of freedom counter
   External_u_local_eqn.resize(n_external_u,Dim);
   //Loop over the external degrees of freedom
   for(unsigned e=0;e<n_external_u;e++)
    {
     //Loop over the velocity values
     for(unsigned i=0;i<Dim;i++)
      {
       External_u_local_eqn(e,i) = external_local_eqn(e,i);
      }
    }
   
   //Add the final additional data object for the external pressure
   //gradient
   if(Delta_P_pt!=0)
    {
     //This is the final external data object
     Delta_P_local_eqn = external_local_eqn(External_Delta_P_index,0);
    }
   else
    {
     Delta_P_local_eqn = -1;
    }
  }

 /// Overload the output function
 void output(ostream &outfile) { }

/// Output function: x,y,[z],u,v,[w],p in tecplot format
void output(ostream &outfile, const unsigned &Np) { }

}; 



//=======================================================================
/// Point element that is used to set the bubble pressure
//=======================================================================
class FixSpineHeightElement : public SpineElement<PointElement>
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
      << "This one contains " << traded_pressure_data_pt->nvalue()<<"\n";;
     
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
 void output(ostream &outfile) { }

 /// Output function: x,y,[z],u,v,[w],p in tecplot format
 void output(ostream &outfile, const unsigned &Np) { }

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

 //Call the generic routine to handle the spine variables
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
     cout<< "The external Data item that stores the traded pressure in "
         << std::endl;
     cout<<"SpineVolumeConstraintPointElement should only have a single"<<std::endl;
     cout<<"value but is has " 
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

//=======================================================================
/// Generalised element used to set the pressure gradient
//=======================================================================
class FluxConstraint : public GeneralisedElement
{
 double Flux;
public:
 /// Constructor, there is one bit of internal data, the flux
 FluxConstraint() {add_internal_data(new Data(1));}
 
 void set_flux(const double &flux_value)
  {Flux = flux_value;}

 const double &read_flux() const {return Flux;}

 //Add the contribution to the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   int local_eqn = internal_local_eqn(0,0);
   if(local_eqn >= 0) 
    {
     residuals[local_eqn] = - Flux;
    }
  }
 
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the residuals
   fill_in_contribution_to_residuals(residuals);
   //No Jacobian terms :)
  }

};

//======================================================================
/// Bretherton problem
//======================================================================
template<class ELEMENT>
class AirwayReopeningProblem : public Problem
{
public:

 /// Constructor: 
 AirwayReopeningProblem();

 /// Destructor, clean up all allocated memory
 ~AirwayReopeningProblem()
  {
   //Delete objects created in constructor in reverse order
   delete Constraint_mesh_pt;
   delete Bulk_mesh_pt;
   //Delete objects associated with lower wall mesh
   delete Global_Physical_Variables::Lower_wall_pt;
   delete Lower_wall_mesh_pt;
   delete Undeformed_lower_wall_geom_pt;
   //Delete objects associated with upper wall mesh
   delete Global_Physical_Variables::Upper_wall_pt;
   delete Upper_wall_mesh_pt;
   delete Undeformed_upper_wall_geom_pt;
   
   //Delete global data
   delete Mesh_fraction_at_transition_pt;
   delete Bubble_pressure_data_pt;
  //Delete the linear solver, if allocated
   if(Frontal_solver) {delete linear_solver_pt();}
  }

 /// Overload the continuation actions because we're 
 /// continuing in Ca which does not affect the mesh
	void actions_after_change_in_global_parameter(double* 
						      const &parameter_pt)
  {
   //Check that the multipliers have worked
   using namespace Global_Physical_Variables;
   ReInvFr = Bo/Ca;
   Re = ReCa*Ca;
   Q = Ca*Gamma;
   T = 100.0*Gamma/H;
  }

 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check()
  {
    actions_after_change_in_global_parameter(&Global_Physical_Variables::Ca);

   // Mesh Update
   Bulk_mesh_pt->node_update();

   // This driver code cannot be allowed to use the analytical form of
   // get_dresidual_dnodal_coordinates(...) that is implemented in the
   // NavierStokesEquations class, since the elemental residuals have
   // contributions from external data which is not taken into account
   // by that routine. We therefore force the bulk elements to use the
   // fully-finite differenced version.
   const unsigned n_bulk_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_bulk_element;e++)
    {
     ElementWithMovingNodes* el_pt =
      dynamic_cast<ElementWithMovingNodes*>(Bulk_mesh_pt->element_pt(e));
     el_pt->evaluate_shape_derivs_by_direct_fd();
    }
  }
 
 /// Update before solve: empty
 void actions_before_newton_solve() {}

 /// Update after solve can remain empty, because the update 
 /// is performed automatically after every Newton step.
 void actions_after_newton_solve() {}

 /// Fix pressure value l in element e to value p_value
 void fix_pressure(const unsigned &e, const unsigned &l, 
                   const double &pvalue)
  {
   //Fix the pressure at that element
   dynamic_cast<ELEMENT *>(Bulk_mesh_pt->element_pt(e))->
    fix_pressure(l,pvalue);
  }

 /// Run a parameter study; perform specified number of steps
 void parameter_study(const unsigned& nsteps,
                      const bool& restart); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Set up the actual outlet elements
 void construct_poisson_outlet_elements()
  {
   //Loop over the elements on the outlet boundary
   unsigned b=1;
   unsigned nbound_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nbound_element;e++)
    {
     SpineGravityTractionElement<ELEMENT>* outlet_element
      = new SpineGravityTractionElement<ELEMENT>(
       Bulk_mesh_pt->boundary_element_pt(b,e),
       Bulk_mesh_pt->face_index_at_boundary(b,e));
     
     //Add the flux constraint
     outlet_element->set_delta_p_pt(Flux_constraint_pt->internal_data_pt(0));
     
     //Add the elements to the local storage
     Outlet_traction_element_pt.push_back(outlet_element);
     
     //Add the elements to the mesh
     Bulk_mesh_pt->add_element_pt(outlet_element);
    }
  }

 /// Delete the outlet elements
 void delete_outlet_elements()
  {
   unsigned n_outlet = Outlet_traction_element_pt.size();
   for(unsigned e=n_outlet;e>0;e--)
    {
     //Remove and delete the outlet elements
     Bulk_mesh_pt->element_pt().pop_back();
     delete Outlet_traction_element_pt[e-1];
     Outlet_traction_element_pt[e-1] = 0;
    }
   Outlet_traction_element_pt.clear();
  }
 
 /// Set up the inital outlet elements
 void construct_pressure_gradient_outlet_elements()
  {
   //Loop over the elements next to the inflow boundarys
   unsigned b=1;
   unsigned nbound_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nbound_element;e++)
    {
     SpineElement<NavierStokesTractionElement<ELEMENT> >* outlet_element
      = new SpineElement<NavierStokesTractionElement<ELEMENT> >(
       Bulk_mesh_pt->boundary_element_pt(b,e),
       Bulk_mesh_pt->face_index_at_boundary(b,e));
     
     //Add the elements to the local storage
     Outlet_traction_element_pt.push_back(outlet_element);
     
     //Add the traction function
     outlet_element->traction_fct_pt() = 
      &Global_Physical_Variables::hydrostatic_pressure;
     
     //Add the elements to the mesh
     Bulk_mesh_pt->add_element_pt(outlet_element);
    }
  }

 double get_outlet_flux()
  {
   if(dynamic_cast<SpineGravityTractionElement<ELEMENT>*>
      (Outlet_traction_element_pt[0])==0)
    {return 0.0;}

   double flux = 0.0;
   //Loop over the outlet and calculate the flux
   for(unsigned e=0;e<Outlet_traction_element_pt.size();e++)
    {
     flux += dynamic_cast<SpineGravityTractionElement<ELEMENT>*>
      (Outlet_traction_element_pt[e])->get_flux();
    }
   return flux;
  }
   

 void change_boundary_conditions()
  {
   //Delete the traction elements and set what we really want
   delete_outlet_elements();

   //Construct the new elements
   construct_poisson_outlet_elements();
   
   //Rebuild the global mesh
   rebuild_global_mesh();
   
   //Now set a consistent flux
   Flux_constraint_pt->set_flux(get_outlet_flux());
   //Unpin the flux constraint
   Flux_constraint_pt->internal_data_pt(0)->unpin(0);
   
   //Hijack the nodes again
   for(unsigned i=0;i<Outlet_traction_element_pt.size();i++)
    {
     dynamic_cast<SpineGravityTractionElement<ELEMENT>*>(
      Outlet_traction_element_pt[i])->hijack_all_nodes();
    }
   
   //Hijack the nodes again
   for(unsigned i=0;i<Gravity_traction_element_pt.size();i++)
    {
     dynamic_cast<SpineGravityTractionElement<ELEMENT>*>(
      Gravity_traction_element_pt[i])->hijack_all_nodes();
    }

   
   //Do equation numbering
   cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 


   //Sort the elements if using the frontal solver
   if(Frontal_solver)
    {
     std::sort(mesh_pt()->element_pt().begin(),
               mesh_pt()->element_pt().end(),
               ElementCmp());
    }

  }

 /// Dump the entire problem
 void dump(ofstream &dump_file) const
  {
   using namespace Global_Physical_Variables;
   dump_file << Re << std::endl;
   dump_file << ReSt << std::endl;
   dump_file << ReInvFr << std::endl;
   dump_file << Ca << std::endl;
   dump_file << Bo << std::endl;
   dump_file << Pext << std::endl;
   dump_file << G[0] << " " << G[1] << std::endl;
   dump_file << Nu << std::endl;
   dump_file << H << std::endl;
   dump_file << Gamma << std::endl;
   dump_file << T << std::endl;
   dump_file << Kstiff << std::endl;
   dump_file << Q << std::endl;
   dump_file << Kseparation << std::endl;
   dump_file << Tube_width << std::endl;
   dump_file << Ktable << std::endl;

   dump_file << Fixed_spine_height << std::endl;
   //Standard problem dumper
   Problem::dump(dump_file);
  }

 /// Read the entire problem
 void read(ifstream &restart_file)
  {
   using namespace Global_Physical_Variables;
   restart_file >> Re;
   restart_file >> ReSt;
   restart_file >> ReInvFr;
   restart_file >> Ca;
   restart_file >> Bo;
   restart_file >> Pext;
   restart_file >> G[0];
   restart_file >> G[1];
   restart_file >> Nu;
   restart_file >> H;
   restart_file >> Gamma;
   restart_file >> T;
   restart_file >> Kstiff;
   restart_file >> Q;
   restart_file >> Kseparation;
   restart_file >> Tube_width;
   restart_file >> Ktable;

   restart_file >> Fixed_spine_height;

   //Standard problem reader
   Problem::read(restart_file);
  }

 /// Setup the coupling between the upper and lower walls
 void connect_walls(Mesh* const &wall_mesh_pt, 
                    GeomObject* const &geom_object_pt,
                    map<double, pair<GeomObject*,Vector<double> > >& globalmap)
  {
   //Loop over the elements in the wall mesh
   unsigned n_wall_element = wall_mesh_pt->nelement();
   for(unsigned e=0;e<n_wall_element;e++)
    {
     //A set of geometric data to be added to the element
     set<Data*> extra_data;
     
     //Cast each element to an FSIWallElement
     FSIHermiteBeamElement *solid_element_pt = 
      dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt->element_pt(e));
     
     //Find the number of Gauss points of the element
     unsigned nintpt = solid_element_pt->integral_pt()->nweight();
     //Find the dimension of the element
     unsigned el_dim = solid_element_pt->dim();
     //Set storage for the local coordinates in the solid and bulk elements
     Vector<double> s(el_dim), s_bulk(el_dim);
     //Set storage for the zeta (intrinsic coordinate) at the Gauss point 
     Vector<double> zeta(el_dim);
     //Define storage for the geometric object that contains the Gauss point
     GeomObject *sub_obj_pt; 
     
     //Loop over the integration points
     for(unsigned ipt=0;ipt<nintpt;ipt++)
      {
       //Loop over the dimension of the solid element and find the local
       //coordinates of the Gauss points
       for(unsigned i=0;i<el_dim;i++) 
        {s[i] = solid_element_pt->integral_pt()->knot(ipt,i);}
       //Get the value of zeta at the integration point
       //Note that zeta coincides with xi (the Lagrangian coordinate)
       //in solid elements
       solid_element_pt->interpolated_xi(s,zeta);
       
       //Find the geometric object and local coordinate within that
       //object for the given value of the intrinsic coordinate, zeta.
       geom_object_pt->locate_zeta(zeta, sub_obj_pt, s_bulk);

       //Stick the result in the map
       globalmap[zeta[0]] = make_pair(sub_obj_pt,s_bulk);

       //Now sort add the external data to the set
       unsigned n_sub_geom_data = sub_obj_pt->ngeom_data();
       for(unsigned i=0;i<n_sub_geom_data;i++)
        {
         extra_data.insert(sub_obj_pt->geom_data_pt(i));
        }
      }
     
     //Add the data to the external data of the object
     for(set<Data*>::iterator it=extra_data.begin();
         it!=extra_data.end();++it)
      {
       solid_element_pt->add_external_data(*it);
      }
    }
  }

private:

 /// Pointer to control element
 ELEMENT* Control_element_pt;

 /// Trace file
 ofstream Trace_file;

 /// Pointer to bulk mesh
 BrethertonSpineMesh<ELEMENT,
                     SpineLineFluidInterfaceElement<ELEMENT> >* 
 Bulk_mesh_pt;

 /// Pointer to Wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Upper_wall_mesh_pt;

 /// Pointer to lower wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Lower_wall_mesh_pt;

 /// Pointer to the constraint mesh
 Mesh* Constraint_mesh_pt;

 /// Pointer to the Undeformed wall
 GeomObject* Undeformed_upper_wall_geom_pt;

 /// Pointer to the Undeformed wall
 GeomObject* Undeformed_lower_wall_geom_pt;

 /// Fixed height value
 double Fixed_spine_height;

 /// Data value that represents the bubble pressure
 Data* Bubble_pressure_data_pt;

 /// Data value that represents the vertical fraction 
 /// between the upper and lower walls of the centre of the spines
 /// at the point of transition of the spine mesh
 Data* Mesh_fraction_at_transition_pt;
 
 /// Vector of pointers to the traction elements
 Vector<FiniteElement*> Gravity_traction_element_pt; 

 /// Vector of pointers to the outlet traciton elements
 Vector<FiniteElement*> Outlet_traction_element_pt;

 /// Flux constraint element
 FluxConstraint* Flux_constraint_pt;

 /// Boolean flag  to specify the use of the frontal solver
 bool Frontal_solver;
};


//====================================================================
/// Problem constructor
//====================================================================
template<class ELEMENT>
AirwayReopeningProblem<ELEMENT>::AirwayReopeningProblem()
{
 Frontal_solver = false;

 // Number of elements in the deposited film region
 unsigned nx1=100;//300;

 // Number of elements in the bottom part of the transition region
 unsigned nx2=6;//12;

 // Number of elements in channel region
 unsigned nx3=100;//300;

 // Number of elements in the vertical part of the transition region
 // (=half the number of elements in the liquid filled region ahead
 // of the finger tip)
 unsigned nhalf=2;//4;

 // Number of elements through thickness of deposited film
 unsigned nh=2;//3;
 
 // Number of elements in each wall
 unsigned nwall = 100;// 450;

 // Thickness of deposited film
 double h = 0.128;

 // Start coordinate on wall
 double xi0= -300.0;

 double start_transition = -1.0;

 // End of transition region on wall
 double xi1 = 1.5;

 // End of liquid filled region (inflow) on wall
 double xi2 = 150.0;//300.0

 //Set the frontal solver, if desired
 if(Frontal_solver)
  {
   linear_solver_pt() = new HSL_MA42;
  }

 //Create a bubble presure data
 Bubble_pressure_data_pt = new Data(1);
 //This will be global data
 add_global_data(Bubble_pressure_data_pt);

 //Create a mesh fraction data
 Mesh_fraction_at_transition_pt = new Data(1);
 add_global_data(Mesh_fraction_at_transition_pt);
 //Pin the value
 Mesh_fraction_at_transition_pt->set_value(0,0.5);
 Mesh_fraction_at_transition_pt->pin(0);

 // The underformed wall is a straight line at y = 1.0
 Undeformed_upper_wall_geom_pt = new StraightLine(1.0);

 // Create the Lagrangian Mesh
 Upper_wall_mesh_pt = 
  new OneDLagrangianMesh<FSIHermiteBeamElement>(nwall,xi0,xi2, 
                                                Undeformed_upper_wall_geom_pt);

 //Loop over the wall elements
 unsigned n_upper_wall_element = Upper_wall_mesh_pt->nelement();
 for(unsigned e=0;e<n_upper_wall_element;e++)
  {
   FSIHermiteBeamElement *element_pt = dynamic_cast<FSIHermiteBeamElement*>(
    Upper_wall_mesh_pt->element_pt(e));
   
   //Assign the undeformed beam shape
   element_pt->undeformed_beam_pt() = Undeformed_upper_wall_geom_pt;
 
   //Set the load vector
   element_pt->load_vector_fct_pt() = &Global_Physical_Variables::spring_load;

   //Set the wall thickness
   element_pt->h_pt() = &Global_Physical_Variables::H;

   //Set the axial pre-stress
   element_pt->sigma0_pt() = &Global_Physical_Variables::T;

   //Function that specifies the load ratios
   element_pt->q_pt() = &Global_Physical_Variables::Q;

   //In this case the normal is directed out of the fluid
   element_pt->set_normal_pointing_out_of_fluid();
 }

 //Create a geometric object that represents the wall geometry
 MeshAsGeomObject* upper_wall_element_pt = 
  new MeshAsGeomObject(Upper_wall_mesh_pt);
 
 // The underformed lower wall is a straight line at y = -1.0
 Undeformed_lower_wall_geom_pt = new StraightLine(-1.0);
 
 // Create the Lagrangian Mesh
 Lower_wall_mesh_pt = 
  new OneDLagrangianMesh<FSIHermiteBeamElement>(nwall,xi0,xi2, 
                                                Undeformed_lower_wall_geom_pt);
 
 //Loop over the wall elements
 unsigned n_lower_wall_element = Lower_wall_mesh_pt->nelement();
 for(unsigned e=0;e<n_lower_wall_element;e++)
  {
   FSIHermiteBeamElement *element_pt = dynamic_cast<FSIHermiteBeamElement*>(
    Lower_wall_mesh_pt->element_pt(e));
   
   //Assign the undeformed beam shape
   element_pt->undeformed_beam_pt() = Undeformed_lower_wall_geom_pt;
 
   //Set the load vector
   element_pt->load_vector_fct_pt() = 
    &Global_Physical_Variables::spring_load_lower;
   
   //Set the wall thickness
   element_pt->h_pt() = &Global_Physical_Variables::H;

   //Set the axial pre-stress
   element_pt->sigma0_pt() = &Global_Physical_Variables::T;

   //Function that specifies the load ratios
   element_pt->q_pt() = &Global_Physical_Variables::Q;

   //In this case the normal is OK
   //element_pt->set_normal_pointing_into_fluid();
 }

 //Create a geometric object that represents the wall geometry
 MeshAsGeomObject* lower_wall_element_pt = 
  new MeshAsGeomObject(Lower_wall_mesh_pt);
 
 // Create wall geom objects
 //GeomObject* lower_wall_pt = new StraightLine(-1.0);
 GeomObject* lower_wall_pt = lower_wall_element_pt;
 Global_Physical_Variables::Lower_wall_pt  = lower_wall_pt;


 //GeomObject* upper_wall_pt=new StraightLine( 1.0);
 GeomObject* upper_wall_pt = upper_wall_element_pt;
 Global_Physical_Variables::Upper_wall_pt = upper_wall_pt;

 //Now create the mesh
 Bulk_mesh_pt = new  BrethertonSpineMesh<ELEMENT,
  SpineLineFluidInterfaceElement<ELEMENT> >
  (nx1,nx2,nx3,nh,nhalf,h,lower_wall_pt,upper_wall_pt,xi0,start_transition,
   xi1,xi2);

 //Set the vertical fraction
 Bulk_mesh_pt->set_spine_centre_fraction_pt(
  Mesh_fraction_at_transition_pt->value_pt(0));

 // Store the control element
 Control_element_pt=Bulk_mesh_pt->control_element_pt();

 //Create a fixed element using the central spine
 FixSpineHeightElement* fix_spine_element_pt = 
  new FixSpineHeightElement(
   static_cast<SpineNode*>(Control_element_pt->node_pt(8)));

 //Set the fixed spine height
 Fixed_spine_height = 
  static_cast<SpineNode*>(Control_element_pt
                          ->node_pt(8))->spine_pt()->height();

 //Add the fixed height to it
 fix_spine_element_pt->height_pt() = &Fixed_spine_height;
 
 //Set the pressure data
 fix_spine_element_pt->set_traded_pressure_data(Bubble_pressure_data_pt);

 //Add the Fixed element to the mesh
 Bulk_mesh_pt->add_element_pt(fix_spine_element_pt);

 //Loop over the elements next to the inflow boundarys
 for(unsigned b=3;b<=5;b+=2)
  {
   unsigned nbound_element = Bulk_mesh_pt->nboundary_element(b);
   for(unsigned e=0;e<nbound_element;e++)
    {
     SpineGravityTractionElement<ELEMENT>* inlet_element
      = new SpineGravityTractionElement<ELEMENT>(
       Bulk_mesh_pt->boundary_element_pt(b,e),
       Bulk_mesh_pt->face_index_at_boundary(b,e));
     
     //Add the elements to the local storage
     Gravity_traction_element_pt.push_back(inlet_element);
     
     //Add the elements to the mesh
     Bulk_mesh_pt->add_element_pt(inlet_element);
    }
  }

 //Set a prescribed flux element
 Flux_constraint_pt = new FluxConstraint;
 
 Constraint_mesh_pt = new Mesh;
 
 Constraint_mesh_pt->add_element_pt(Flux_constraint_pt);
 
 Flux_constraint_pt->internal_data_pt(0)->pin(0);

 construct_pressure_gradient_outlet_elements();
  
 Bulk_mesh_pt->node_update();

 {
  //Set the boundary coordinates on the upper wall
  unsigned nbound = Bulk_mesh_pt->nboundary_node(2);
  Vector<double> zeta(1);
  for(unsigned n=0;n<nbound;n++)
   {
    Node* node_pt = Bulk_mesh_pt->boundary_node_pt(2,n);
    zeta[0] = node_pt->x(0);
    node_pt->set_coordinates_on_boundary(2,zeta);
   }
 }
 
 //Sort out the interaction the upper wall
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,2,Bulk_mesh_pt,Upper_wall_mesh_pt);

 {
  //Set the boundary coordinates on the lower wall
  unsigned nbound = Bulk_mesh_pt->nboundary_node(0);
  Vector<double> zeta(1);
  for(unsigned n=0;n<nbound;n++)
   {
    Node* node_pt = Bulk_mesh_pt->boundary_node_pt(0,n);
    zeta[0] = node_pt->x(0);
    node_pt->set_coordinates_on_boundary(0,zeta);
   }
 }

 //Sort out the interaction the lower wall
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,0,Bulk_mesh_pt,Lower_wall_mesh_pt);

 add_sub_mesh(Bulk_mesh_pt);

 add_sub_mesh(Upper_wall_mesh_pt);

 add_sub_mesh(Lower_wall_mesh_pt);

 add_sub_mesh(Constraint_mesh_pt);

 build_global_mesh();


 // Connect the upper mesh
 connect_walls(Lower_wall_mesh_pt,upper_wall_pt,
               Global_Physical_Variables::upper_map);

 connect_walls(Upper_wall_mesh_pt,lower_wall_pt,
               Global_Physical_Variables::lower_map);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here
 
 //Pin one end of the mesh (prevent transvere slide)
 Upper_wall_mesh_pt->boundary_node_pt(1,0)->pin_position(0,0);
 Upper_wall_mesh_pt->boundary_node_pt(0,0)->pin_position(0,0);

 //Pin both ends of the lower mesh
 Lower_wall_mesh_pt->boundary_node_pt(1,0)->pin_position(0,0);
 Lower_wall_mesh_pt->boundary_node_pt(0,0)->pin_position(0,0);

 //Pin the lower mesh vertically
 Lower_wall_mesh_pt->boundary_node_pt(0,0)->pin_position(0,1);

 // No slip on boundaries 0 1 (inflow) and 2
 for(unsigned long ibound=0;ibound<=2;ibound++)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Leave the inflow traction free in the both directions
     if(ibound!=1)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
      }
    }
  }

 // Set the values of the velocities on the wall boundaries
 for (unsigned ibound=0;ibound<=2;ibound+=2)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parallel flow
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,-1.0);
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1, 0.0);
    }
  }
 
 //Loop over the upper wall
 {
  unsigned b=2;
  unsigned num_nod=Bulk_mesh_pt->nboundary_node(b);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    //cast to a spine node
    SpineNode* spine_node_pt = static_cast<SpineNode*>
     (Bulk_mesh_pt->boundary_node_pt(b,inod));

    //Find the id of the node
    unsigned id = spine_node_pt->node_update_fct_id();
    switch(id)
     {
      
     case 4:
     case 5:
      //Add the auxilliary function, so that the no-slip condition is
      //applied properly
      spine_node_pt->set_auxiliary_node_update_fct_pt(
       No_Slip::no_slip_condition_first);
      break;
      
     case 6:
      spine_node_pt->set_auxiliary_node_update_fct_pt(
       No_Slip::no_slip_condition_second);
      break;

     default:
      cout << "Should not get here" << std::endl;
      break;
     }
   }
 }
  
 //Loop over the lower wall
 {
  unsigned b=0;
  unsigned num_nod=Bulk_mesh_pt->nboundary_node(b);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    //cast to a spine node
    SpineNode* spine_node_pt = static_cast<SpineNode*>
     (Bulk_mesh_pt->boundary_node_pt(b,inod));

    //Find the id of the node
    unsigned id = spine_node_pt->node_update_fct_id();
    switch(id)
     {
      
     case 0:
     case 1:
     case 6:
      //Add the auxilliary function, so that the no-slip condition is
      //applied properly
      spine_node_pt->set_auxiliary_node_update_fct_pt(
       No_Slip::no_slip_condition_first);
      break;

     default:
      cout << "Should not get here" << std::endl;
      break;
     }
   }
 }

 //Parallel, uniform outflow on boundaries 3 and 5
 for (unsigned ibound=3;ibound<=5;ibound+=2)
  {
   unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Parallel inflow
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0, -1.0);
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1, 0.0);
    }
  }

 //Loop over the elements in the layer
 unsigned long n_bulk=Bulk_mesh_pt->nbulk();
 for(unsigned long i=0;i<n_bulk;i++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           bulk_element_pt(i));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Loop over 1D interface elements and set capillary number
 unsigned interface_element_pt_range = Bulk_mesh_pt->ninterface_element();
 for(unsigned i=0;i<interface_element_pt_range;i++)
  {
   //Cast to a interface element
   SpineLineFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>
    (Bulk_mesh_pt->interface_element_pt(i));

   //Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;
   //Set the external pressure data
   el_pt->set_external_pressure_data(Bubble_pressure_data_pt);

   //We need to make sure that we hijack the nodes on the boundaries
   if(i==0) {el_pt->hijack_nodal_value(0,0);}
   if(i==interface_element_pt_range-1) 
    {el_pt->hijack_nodal_value(el_pt->nnode()-1,0);}
  }

 //Loop over the traction elements
 unsigned n_traction = Gravity_traction_element_pt.size();
 for(unsigned e=0;e<n_traction;e++)
  {
   //Cast to the traction element
   SpineGravityTractionElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineGravityTractionElement<ELEMENT>*>
    (Gravity_traction_element_pt[e]);
   
   //Set the Reynolds number divided by the Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Do equation numbering
 cout << "Number of unknowns: " << assign_eqn_numbers() << std::endl; 

 //Sort the elements if using the frontal solver
 //If this is not done, it's very slow indeed
 //If it is done, then it's very fast compared to SuperLU
 if(Frontal_solver)
  {
   std::sort(mesh_pt()->element_pt().begin(),
             mesh_pt()->element_pt().end(),
             ElementCmp());
  }
}

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void AirwayReopeningProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Number of interface elements
// unsigned ninterface=Bulk_mesh_pt->ninterface_element();
 
 // Number of spines
 //unsigned nspine=Bulk_mesh_pt->nspine();

 // Control coordinate: At bubble tip
 Vector<double> s(2);
 s[0]=1.0;
 s[1]=1.0;
 
 // Last proper spine
 unsigned last_spine=Bulk_mesh_pt->nfree_surface_spines()-1;

 // Doc
 Trace_file << " " << Global_Physical_Variables::Ca;
 Trace_file << " " << Global_Physical_Variables::Bo;
 Trace_file << " " << Global_Physical_Variables::Re;
 Trace_file << " " << Global_Physical_Variables::Tube_width;
 Trace_file << " " << get_outlet_flux();
 Trace_file << " " <<  
  Upper_wall_mesh_pt->boundary_node_pt(0,0)->x(1);
 Trace_file << " " << 
  Upper_wall_mesh_pt->boundary_node_pt(1,0)->x(1);
 Trace_file << " " <<  
  Lower_wall_mesh_pt->boundary_node_pt(0,0)->x(1);
 Trace_file << " " << 
  Lower_wall_mesh_pt->boundary_node_pt(1,0)->x(1);
 Trace_file << " " << Bubble_pressure_data_pt->value(0);
  Trace_file << " " << Bulk_mesh_pt->spine_pt(0)->height();
 Trace_file << " " << Bulk_mesh_pt->spine_pt(last_spine)->height();
 Trace_file << " " << 1.3375*pow(Global_Physical_Variables::Ca,2.0/3.0);
 Trace_file << " " << 
  (Bubble_pressure_data_pt->value(0)-Control_element_pt->interpolated_p_nst(s))*
                      Global_Physical_Variables::Ca;
 Trace_file << " " << 1.0+3.8*pow(Global_Physical_Variables::Ca,2.0/3.0);
 Trace_file << " " << Control_element_pt->interpolated_u_nst(s,0);
 Trace_file << " " << Control_element_pt->interpolated_u_nst(s,1);
 Trace_file << " "  
//             << abs(dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>(
//                Bulk_mesh_pt->interface_element_pt(0))->
//                actual_contact_angle_left())*
//                180.0/MathematicalConstants::Pi << " " ;
//  Trace_file << " "  
//             << abs(dynamic_cast<SpineLineFluidInterfaceElement<ELEMENT>*>(
//                Bulk_mesh_pt->interface_element_pt(ninterface-1))->
//                actual_contact_angle_right())*
//                180.0/MathematicalConstants::Pi 
            << " ";
 //Trace_file << get_outlet_flux() << " ";
 Trace_file << std::endl;


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned n_bulk=Bulk_mesh_pt->nbulk();
 for(unsigned i=0;i<n_bulk;i++)
  {
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->
                                           bulk_element_pt(i));
   el_pt->output(some_file,npts);
  }
 some_file.close();

}




//=============================================================================
/// Parameter study
//=============================================================================
template<class ELEMENT>
void AirwayReopeningProblem<ELEMENT>::parameter_study(const unsigned& nsteps,
                                                 const bool &restart)
{
 double flux = 0.0;
 // Increase maximum residual
 Problem::Max_residuals=100.0;
 Problem::Max_newton_iterations=100;
 //Problem::Newton_solver_tolerance=1.0e-6;

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 doc_info.number()=0;

 /*double zeta_trans_start_lo = -1.0;
 double zeta_trans_end_lo = 1.5;
 double fraction = 0.45;*/

 // Open trace file
 char filename[100], dumpfile[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 Trace_file << "VARIABLES=\"Ca\",";
 Trace_file << "\"h<sub>bottom</sub>\",\"h<sub>too</sub>\",";
 Trace_file << "\"h<sub>Bretherton</sub>\",\"p<sub>tip</sub>\",";
 Trace_file << "\"p<sub>tip (Bretherton)</sub>\",\"u<sub>stag</sub>\",";
 Trace_file << "\"v<sub>stag</sub>\"";
 Trace_file << "\"<greek>a</greek><sub>bottom</sub>\",";
 Trace_file << "\"<greek>a</greek><sub>top</sub>\"";
 Trace_file << std::endl;

 if(restart) 
  {
   ifstream dumpfile("dump.dat");
   read(dumpfile); 
   dumpfile.close();
  }

 // Initial scaling factor for Ca reduction
 //double factor=1.5;

 // Update so the initial domain shape is plotted properly
 Bulk_mesh_pt->node_update();

 //Doc initial solution
 doc_solution(doc_info);

 //Loop over the steps
 for(unsigned step=1;step<=nsteps;step++)
  {
   cout << std::endl << "STEP " << step << std::endl;

   //Solve
   cout << "Solving for capillary number: " 
        << Global_Physical_Variables::Ca 
        << " Bo = " 
        << Global_Physical_Variables::Bo << std::endl;

   if((!restart) || (step > 1))
    {
     using namespace Global_Physical_Variables;
     newton_solve();
    }
   else
    {
     cout << "Restart skipping" << std::endl;
    }
   cout << "Bubble pressure is " << Bubble_pressure_data_pt->value(0) << std::endl;

   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info);
   
   //Dump if we're solving a real problem
   if(step > 1)
    {
     using namespace Global_Physical_Variables;
     sprintf(dumpfile,
             "%s/dump.Ca_%g.Bo_%g.Re_%g.H_%g.TW_%g_Frac_%g.Kstiff_%g.dat",
             doc_info.directory().c_str(),
             Ca,Bo,Re,H,Tube_width,
             Mesh_fraction_at_transition_pt->value(0),Kstiff);
     ofstream dumpstream(dumpfile);
     dump(dumpstream);
     dumpstream.close();
     cout << "Actual flux is " << get_outlet_flux() << std::endl;
    }
   
   //First time round, unpin things
   if(step==1)
    {
     change_boundary_conditions();
     flux = Flux_constraint_pt->read_flux();
    }
   //Otherwise do our parameter study
   else
    {
     using namespace Global_Physical_Variables;
     //Now start to increase the flux
     //The computed flux is -0.389 or something so
     //if we start from -0.35 and subtract 0.05 we get 0.4
     //and then continue from the first time
     if(step==2) {flux = -0.35;}
     
     //Then worry about the flux
     if(flux > -2.0) {flux -= 0.05; Flux_constraint_pt->set_flux(flux);}
     else
      {
       if(flux < -2.0) {flux = -2.0; Flux_constraint_pt->set_flux(flux);}
       
       //Now increase the capillary number
       if(Ca < 0.1) {Ca += 0.01;}
       else
        {
         Ca = 0.1;
         //Finally reduce the spring stiffness
         Kstiff -= 0.5*1.0e-7;
        }
      }
    }
  }
}

       
//======================================================================
/// Driver code for unsteady two-layer fluid problem. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only a single step.
//======================================================================
int main(int argc, char* argv[]) 
{
  // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set physical parameters:
 using namespace Global_Physical_Variables;
 
 //Set direction of gravity: Vertically downwards
 G[0] = 0.0;
 G[1] = -1.0;

 // Womersley number = Reynolds number (St = 1)
 ReSt = 0.0;
 Re = Global_Physical_Variables::ReSt;

 // The Capillary number
 Ca = 0.05;

 // The Bond number
 Bo = 0.0;

 // Re/Fr -- a measure of gravity...
 ReInvFr = Bo/Ca;

 // The ratio of scales
 Q = Ca*Gamma;

 // The value of sigma0
 T = 100.0*Gamma/H;

 //Set up the problem
 AirwayReopeningProblem<Hijacked<SpineElement<QCrouzeixRaviartElement<2> > > > 
  problem;

 // Self test:
 problem.self_test();

 // Number of steps: 
 unsigned nstep;
 if (CommandLineArgs::Argc>1)
  {
   // Validation run: Just one step
   nstep=2;
  }
 else
  {
   // Full run otherwise
   nstep=200;
  }

 // Run the parameter study: Perform nstep steps
 problem.parameter_study(nstep,false);
}

