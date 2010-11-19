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
#include <fenv.h> 

//Generic routines
#include "generic.h"


// The equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "rigid_body.h"


// The mesh
//#include "my_triangle_mesh.template.h"
//#include "my_triangle_mesh.template.cc"
#include "meshes/triangle_mesh.h"

#include <algorithm>
using namespace std;
using namespace oomph;

void imposed_torque(const double &time, 
            double &external_torque)
{
 external_torque = time*time*time;
}


namespace oomph
{
//==start_of_namespace==============================
/// Namespace for Problem Parameter
//==================================================
 namespace Problem_Parameter
 {    
  /// Doc info
  DocInfo Doc_info;
  
  /// Reynolds number
  double Re=0.0;

  /// Pseudo-solid Poisson ratio
  double Nu=0.3;

  /// \short Pseudo solid "density" -- set to zero because we don't want
  /// inertia in the node update!
  double Lambda_sq=0.0;

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt=0;

  /// Trace file
  ofstream Trace_file;

  /// File to document the norm of the solution (for validation purposes)
  ofstream Norm_file;

  /// File to document the motion of the centre of gravity
  ofstream Cog_file;

  /// Direction of gravity
  Vector<double> G(2);

  /// Magnitude of gravitational body force in Navier--Stokes
  /// much, much smaller than the rigid body
  double ReInvFr = 0.0;

  /// Magnitude of body force acting on rigid body
  double Grav = 1.0;

  /// Parameter that determines whether the flow is on
  double Alpha = 0.0;

  /// Axis in x-direction
  double A = 0.25;
  
  /// Axis in y-direction
  double B = 0.5;

 } // end_of_namespace


//=======================================================================
///Exact solution for the rotation of an ellipse in unbounded shear flow
///as computed by Jeffery (1922)
//=======================================================================
namespace Jeffery_Solution
{
 //Null function
 double null(const double &t) {return 0.0;}

 //Angular position as a function of time t
 double angle(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::B;
  
  return atan((b/a)*tan((a*b*t)/(b*b + a*a)));
 }

 //Angular velocity as function of time t
 double velocity(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::B;

  //Get the angle
  double chi = angle(t);

  //Now return the velocity
  return (a*a*sin(chi)*sin(chi) + b*b*cos(chi)*cos(chi))/(a*a + b*b);
 }


 //Angular acceleration as a function of time t
 double acceleration(const double &t)
 {
  const double a = Problem_Parameter::A;
  const double b = Problem_Parameter::A;

  //Get the angle and velocity
  double chi = angle(t);
  double chi_dot = velocity(t);

  //Now return the acceleration
  return 2.0*(a*a - b*b)*(sin(chi)*cos(chi))*chi_dot/(a*a + b*b);
 }
}

//My own Ellipse class
class GeneralEllipse : public GeomObject
{
private:
 //Internal data to store the centre and semi-axes
 double *centre_x_pt, *centre_y_pt, *a_pt, *b_pt;

public:
 
 //Constructor
 GeneralEllipse(const double &centre_x, const double &centre_y,
                const double &a, const double &b)
  : GeomObject(1,2), centre_x_pt(0), centre_y_pt(0), a_pt(0), b_pt(0)
  {
   centre_x_pt = new double(centre_x);
   centre_y_pt = new double(centre_y);
   a_pt = new double(a);
   b_pt = new double(b);
  }

 //Destructor
 ~GeneralEllipse()
  {
   delete centre_x_pt;
   delete centre_y_pt;
   delete a_pt;
   delete b_pt;
  }

 //Return the position
 void position(const Vector<double> &xi, Vector<double> &r) const
  {
   r[0] = *centre_x_pt + *a_pt*cos(xi[0]);
   r[1] = *centre_y_pt + *b_pt*sin(xi[0]);
  }

 //Return the position which is always fixed
 void position(const unsigned &t,
               const Vector<double> &xi, Vector<double> &r) const
  {
   return position(xi,r);
  }

};
 

//==============================================================
/// Overload TaylorHood element to modify output
//==============================================================
 class MyTaylorHoodElement : 
  public virtual PseudoSolidNodeUpdateElement<TTaylorHoodElement<2>, 
  TPVDElement<2,3> >
 {
  
 private:
  
  /// Storage for elemental error estimate -- used for post-processing
  double Error;

 public:

  /// Constructor initialise error
  MyTaylorHoodElement()
   {
    Error=0.0;
   }

  /// Set error value for post-processing
  void set_error(const double& error){Error=error;}
  
  /// Return variable identifier
  std::string variable_identifier()
   {
    std::string txt="VARIABLES=";
    txt+="\"x\",";
    txt+="\"y\",";
    txt+="\"u\",";
    txt+="\"v\",";
    txt+="\"p\",";   
    txt+="\"du/dt\",";
    txt+="\"dv/dt\",";
    txt+="\"u_m\",";   
    txt+="\"v_m\",";
    txt+="\"x_h1\",";
    txt+="\"y_h1\",";   
    txt+="\"x_h2\",";
    txt+="\"y_h2\",";   
    txt+="\"u_h1\",";
    txt+="\"v_h1\",";   
    txt+="\"u_h2\",";
    txt+="\"v_h2\",";   
    txt+="\"error\",";   
    txt+="\"size\",";   
    txt+="\n";
    return txt;
   }

  
  /// Overload output function
  void output(std::ostream &outfile, 
              const unsigned &nplot)
   {
    
    // Assign dimension 
    unsigned el_dim=2;
    
    // Vector of local coordinates
    Vector<double> s(el_dim);
    
    // Acceleration
    Vector<double> dudt(el_dim);
    
    // Mesh elocity
    Vector<double> mesh_veloc(el_dim,0.0);
   
    // Tecplot header info
    outfile << tecplot_zone_string(nplot);
   
    // Find out how many nodes there are
    unsigned n_node = nnode();
   
    //Set up memory for the shape functions
    Shape psif(n_node);
    DShape dpsifdx(n_node,el_dim);
   
    // Loop over plot points
    unsigned num_plot_points=nplot_points(nplot);
    for (unsigned iplot=0;iplot<num_plot_points;iplot++)
     {
     
      // Get local coordinates of plot point
      get_s_plot(iplot,nplot,s);
     
      //Call the derivatives of the shape and test functions
      dshape_eulerian(s,psif,dpsifdx);
     
      //Allocate storage
      Vector<double> mesh_veloc(el_dim);
      Vector<double> dudt(el_dim);
      Vector<double> dudt_ALE(el_dim);
      DenseMatrix<double> interpolated_dudx(el_dim,el_dim);
     
      //Initialise everything to zero
      for(unsigned i=0;i<el_dim;i++)
       {
        mesh_veloc[i]=0.0;
        dudt[i]=0.0;
        dudt_ALE[i]=0.0;
        for(unsigned j=0;j<el_dim;j++)
         {
          interpolated_dudx(i,j) = 0.0;
         }
       }
     
      //Calculate velocities and derivatives

      //Loop over directions
      for(unsigned i=0;i<el_dim;i++)
       {
        //Get the index at which velocity i is stored
        unsigned u_nodal_index = u_index_nst(i);
        // Loop over nodes
        for(unsigned l=0;l<n_node;l++) 
         {
          dudt[i]+=du_dt_nst(l,u_nodal_index)*psif[l];
          mesh_veloc[i]+=dnodal_position_dt(l,i)*psif[l];
          
          //Loop over derivative directions for velocity gradients
          for(unsigned j=0;j<el_dim;j++)
           {                               
            interpolated_dudx(i,j) += nodal_value(l,u_nodal_index)*
             dpsifdx(l,j);
           }
         }
       }
     
     
      // Get dudt in ALE form (incl mesh veloc)
      for(unsigned i=0;i<el_dim;i++)
       {
        dudt_ALE[i]=dudt[i];
        for (unsigned k=0;k<el_dim;k++)
         {
          dudt_ALE[i]-=mesh_veloc[k]*interpolated_dudx(i,k);
         }
       }
     
     
      // Coordinates
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_x(s,i) << " ";
       }
     
      // Velocities
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << interpolated_u_nst(s,i) << " ";
       }
     
      // Pressure
      outfile << interpolated_p_nst(s)  << " ";
     
      // Accelerations
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << dudt_ALE[i] << " ";
       }
     
      // Mesh velocity
      for(unsigned i=0;i<el_dim;i++) 
       {
        outfile << mesh_veloc[i] << " ";
       }
     
      // History values of coordinates
      unsigned n_prev=node_pt(0)->position_time_stepper_pt()->ntstorage();
      for (unsigned t=1;t<n_prev;t++)
       {
        for(unsigned i=0;i<el_dim;i++) 
         {
          outfile << interpolated_x(t,s,i) << " ";
         }
       }
     
      // History values of velocities
      n_prev=node_pt(0)->time_stepper_pt()->ntstorage();
      for (unsigned t=1;t<n_prev;t++)
       {
        for(unsigned i=0;i<el_dim;i++) 
         {
          outfile << interpolated_u_nst(t,s,i) << " ";
         }
       }

      outfile << Error << " " 
              << size() << std::endl;        
     }
    
    // Write tecplot footer (e.g. FE connectivity lists)
    write_tecplot_zone_footer(outfile,nplot); 
    }


  /// Get square of L2 norm of velocity 
  double square_of_l2_norm()
   {

    // Assign dimension 
    unsigned el_dim=2; 
    // Initalise
    double sum=0.0;
    
    //Find out how many nodes there are
    unsigned n_node = nnode();
    
    //Find the indices at which the local velocities are stored
    unsigned u_nodal_index[el_dim];
    for(unsigned i=0;i<el_dim;i++) {u_nodal_index[i] = u_index_nst(i);}
    
    //Set up memory for the velocity shape fcts
    Shape psif(n_node);
    DShape dpsidx(n_node,el_dim);
    
    //Number of integration points
    unsigned n_intpt = integral_pt()->nweight();
    
    //Set the Vector to hold local coordinates
    Vector<double> s(el_dim);
    
    //Loop over the integration points
    for(unsigned ipt=0;ipt<n_intpt;ipt++)
     {
      //Assign values of s
      for(unsigned i=0;i<el_dim;i++) s[i] = integral_pt()->knot(ipt,i);
      
      //Get the integral weight
      double w = integral_pt()->weight(ipt);
      
      // Call the derivatives of the veloc shape functions
      // (Derivs not needed but they are free)
      double J = this->dshape_eulerian_at_knot(ipt,psif,dpsidx);
      
      //Premultiply the weights and the Jacobian
      double W = w*J;
      
      //Calculate velocities 
      Vector<double> interpolated_u(el_dim,0.0);      
      
      // Loop over nodes
      for(unsigned l=0;l<n_node;l++) 
       {
        //Loop over directions
        for(unsigned i=0;i<el_dim;i++)
         {
          //Get the nodal value
          double u_value = raw_nodal_value(l,u_nodal_index[i]);
          interpolated_u[i] += u_value*psif[l];
         }
       }

      //Assemble square of L2 norm
      for(unsigned i=0;i<el_dim;i++)
       {
        sum+=interpolated_u[i]*interpolated_u[i]*W;
       }           
     }

    return sum;

   }

 };





//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<>
 class FaceGeometry<MyTaylorHoodElement>
  : public virtual SolidTElement<1,3> 
 {
 public:
  FaceGeometry() : SolidTElement<1,3>() {}
 };



} //End of namespace extension



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Unstructured Navier-Stokes ALE Problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem :
 public virtual ProjectionProblem<ELEMENT>

{

public:

 /// Constructor
 UnstructuredFluidProblem();
 
 /// Destructor
 ~UnstructuredFluidProblem()
  {
   // Fluid timestepper
   delete this->time_stepper_pt(0);

   // Kill data associated with outer boundary
   unsigned n=Outer_boundary_polyline_pt->npolyline();
   for (unsigned j=0;j<n;j++)
    {
     delete Outer_boundary_polyline_pt->polyline_pt(j);
    }
   delete Outer_boundary_polyline_pt;
   
   // Flush Lagrange multiplier mesh
   delete_lagrange_multiplier_elements();
   delete Lagrange_multiplier_mesh_pt;

   // Delete error estimator
   delete Fluid_mesh_pt->spatial_error_estimator_pt();

   // Delete fluid mesh
   delete Fluid_mesh_pt;

   // Kill const eqn
   delete Problem_Parameter::Constitutive_law_pt;

  }

 /// Reset the boundary conditions when timestepping
 void actions_before_implicit_timestep()
  {
   this->set_boundary_velocity(true);
  }
 
 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt()
  {
   // Kill the  elements and wipe surface mesh
   delete_lagrange_multiplier_elements();
   
   //Kill the drag element
   delete_drag_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
  
  }// end of actions_before_adapt

 
 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt()
  {
   // Create the elements that impose the displacement constraint 
   create_lagrange_multiplier_elements();
   
   // Create the drag elements anew
   create_drag_elements();

   // Rebuild the Problem's global mesh from its various sub-meshes
   this->rebuild_global_mesh();
   
   // Setup the problem again -- remember that fluid mesh has been
   // completely rebuilt and its element's don't have any
   // pointers to Re etc. yet
   complete_problem_setup();

   // Output solution after adaptation/projection
   bool doc_projection=true;
   doc_solution("new mesh with projected solution",doc_projection);
  }// end of actions_after_adapt

 
 /// \short Re-apply the no slip condition (imposed indirectly via enslaved
 /// velocities)
 void actions_before_newton_convergence_check()
  {
   // Update mesh -- this applies the auxiliary node update function
   Fluid_mesh_pt->node_update();
  }
 
 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
 
 
 /// \short Set boundary condition, assign auxiliary node update fct.
 /// Complete the build of all elements, attach power elements that allow
 /// computation of drag vector
 void complete_problem_setup()
  {   
   // Set the boundary conditions for fluid problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   unsigned nbound=Fluid_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<nbound;ibound++)
    {
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Pin everywhere apart from outflow and inflow (boundaries 0, 2)
       // where we only impose parallel flow)
       if((ibound!=2) && (ibound!=0))
        {
         nod_pt->pin(0);
        }
       nod_pt->pin(1);
       
       // Pin pseudo-solid positions apart from hole boundary we want to move
       SolidNode* solid_node_pt = dynamic_cast<SolidNode*>(nod_pt);
       
       // Unpin the position of all the nodes on hole boundaries
       // since they will be moved using Lagrange Multiplier
       if(ibound > 3)
        {
         solid_node_pt->unpin_position(0);
         solid_node_pt->unpin_position(1);
         
         // Assign auxiliary node update fct if we're dealing with a 
         // hole boundary
         // A more accurate version may be obtained by using velocity
         // based on the actual position of the geometric object
         nod_pt->set_auxiliary_node_update_fct_pt(
          FSI_functions::apply_no_slip_on_moving_wall); 
        }
       else
        {
         solid_node_pt->pin_position(0);
         solid_node_pt->pin_position(1);
        }
      }

    } // end loop over boundaries
   
   // Complete the build of all elements so they are fully functional
   unsigned n_element = Fluid_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));

     // Set pointer to continous time
     el_pt->time_pt()=this->time_pt();
   
     // Set the Reynolds number
     el_pt->re_pt() = &Problem_Parameter::Re;

     // Set the Wormesley number (same as Re since St=1)
     el_pt->re_st_pt() = &Problem_Parameter::Re;

     // Set the constitutive law for pseudo-elastic mesh deformation
     el_pt->constitutive_law_pt()=Problem_Parameter::Constitutive_law_pt;

     // Set the "density" for pseudo-elastic mesh deformation
     el_pt->lambda_sq_pt()=&Problem_Parameter::Lambda_sq;
    }
   
   // Re-apply Dirichlet boundary conditions (projection ignores
   // boundary conditions!)
   
   // Zero y - velocity and history values of velocity on walls 
   // (boundaries 0, 1 and 3)
   nbound=this->Fluid_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<4;ibound++)
    {
     unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       
       // Get number of previous (history) values
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       
       // Zero all current and previous veloc values  on all boundaries
       // apart from x-component on boundary 2
       for (unsigned t=0;t<=n_prev;t++)
        {
         if((ibound!=2) && (ibound!=0))
          {
           nod_pt->set_value(t,0,0.0);
          }
         nod_pt->set_value(t,1,0.0);
        }
      }
    }
  

   //Turn on a shear flow
   this->set_boundary_velocity();

   //Reset the solid boundary conditions
   //Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  }


 ///Set the boundary velocity
 void set_boundary_velocity(const bool &current_value_only=false)
  {
   //Loop over top and bottom walls and inlet
   for(unsigned ibound=0;ibound<4;ibound++)
    {
     //if not the outlet or inlet
     if((ibound!=2) && (ibound != 0))
      {
       unsigned num_nod=this->Fluid_mesh_pt->nboundary_node(ibound);
       for (unsigned inod=0;inod<num_nod;inod++)
        {
         // Get node
         Node* nod_pt=this->Fluid_mesh_pt->boundary_node_pt(ibound,inod);
         
         // Get number of previous (history) values
         unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
         
         //If only doing the current step set n_prev = 0
         if(current_value_only) {n_prev=0;}

         //Now set the boundary velocity
         double y = nod_pt->x(1);
         //Get the previous time
         for(unsigned t=0;t<=n_prev;t++)
          {
           //Get the time
           double time_ = this->time_pt()->time(t);
           
           //Get the velocity ramp
           //Initially zero (nothing at all is going on)
           double ramp = 0.0;
           double delta = 5.0;
           
           
           double e1 = exp(-delta);
           double a1 = 1.0 - (1.0 + delta + 0.5*delta*delta)*e1;
           double b1 = (3.0 + 3.0*delta + delta*delta)*e1 - 3.0;
           double c1 = 3.0 - (3.0 + 2.0*delta + 0.5*delta*delta)*e1;
           //Smooth start
           if((time_ >= 0.0) & (time_ <= 1.0)) 
            { 
             ramp = a1*time_*time_*time_
              + b1*time_*time_
              + c1*time_;
            }
           //Coupled to exponential levelling
           else if (time_ > 1.0)
            {
             ramp = 1.0 - exp(-delta*time_);
            }
           
           nod_pt->set_value(t,0,-y*ramp);
          }
    
         // Zero all current and previous veloc values
         // for the v-velocity
         for (unsigned t=0;t<=n_prev;t++)
          {
           nod_pt->set_value(t,1,0.0);
          }
        }
      }
    }
  }

 /// Set the history of the rigid body from the Jeffery solution
 /*void set_consistent_rigid_body_history()
  {
   //Get the newmark timestepper
   Newmark<2>* newmark_pt = 
    dynamic_cast<Newmark<2>*>(
     dynamic_cast<GeneralisedElement*>(this->Rigid_body_pt[0])
     ->internal_data_pt(0)->time_stepper_pt());

   //Prepare the vectors of initial conditions, velocities and accelerations
   Vector<Newmark<2>::
    InitialConditionFctPt> initial_value_fct(3,Jeffery_Solution::null);
   Vector<Newmark<2>::
    InitialConditionFctPt> initial_veloc_fct(3,Jeffery_Solution::null);
   Vector<Newmark<2>::
    InitialConditionFctPt> initial_accel_fct(3,Jeffery_Solution::null);

   //Explicity set the rotation components
   initial_value_fct[2] = Jeffery_Solution::angle;
   initial_veloc_fct[2] = Jeffery_Solution::velocity;
   initial_accel_fct[2] = Jeffery_Solution::acceleration;

   //Now set the data
   newmark_pt->assign_initial_data_values(
    dynamic_cast<GeneralisedElement*>(
     this->Rigid_body_pt[0])->internal_data_pt(0),
    initial_value_fct, initial_veloc_fct,initial_accel_fct);

   //Tell me the answer
   Vector<double> veloc(3);
   Vector<double> accel(3);

   newmark_pt->time_derivative(1,  dynamic_cast<GeneralisedElement*>(
                                this->Rigid_body_pt[0])->internal_data_pt(0), 
                               veloc);

   newmark_pt->time_derivative(2,  dynamic_cast<GeneralisedElement*>(
                                this->Rigid_body_pt[0])->internal_data_pt(0), 
                               accel);
   
   //Should also set the position histories of the nodes on the boundary
   //consistently
   //Loop over all the ellipse boundaries
   unsigned n_boundary = Fluid_mesh_pt->nboundary();
   for(unsigned b=4;b<n_boundary;b++)
    {
     //Find the number of nodes on the boundary 
     unsigned n_node = Fluid_mesh_pt->nboundary_node(b);
     //Loop over the nodes
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(b,n);
       //Get the boundary coordinate
       Vector<double> zeta(1);
       nod_pt->get_coordinates_on_boundary(b,zeta);
       
       //Loop over the time histories and positions
       unsigned n_prev=nod_pt->time_stepper_pt()->nprev_values();
       for(unsigned t=0;t<=n_prev;t++)
        {
         //Get the position history from the rigid body
         Vector<double> r(2), drdt(2);
         this->Rigid_body_pt[0]->position(t,zeta,r);

         //Set the position history
         for(unsigned i=0;i<2;i++) {nod_pt->x(t,i) = r[i];}

         //Set the velocity histories
         this->Rigid_body_pt[0]->dposition_dt(zeta,1,drdt);
         for(unsigned i=0;i<2;i++)
          {nod_pt->set_value(t,i,drdt[i]);}
        }
       
       //Reset the auxilliary node update function
       nod_pt->set_auxiliary_node_update_fct_pt(
        FSI_functions::apply_no_slip_on_moving_wall); 
      }
    }
    }*/

 /// Set the boundary velocity on the rigid body
 /*void set_jeffery_velocity()
  {
   //Get the angular velocity from the jeffery solution
   double chi_dot = Jeffery_Solution::velocity(this->time());

   //Loop over all the ellipse boundaries
   unsigned n_boundary = Fluid_mesh_pt->nboundary();
   for(unsigned b=4;b<n_boundary;b++)
    {
     //Find the number of nodes on the boundary 
     unsigned n_node = Fluid_mesh_pt->nboundary_node(b);
     //Loop over the nodes
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(b,n);
       //Get the boundary coordinate
       Vector<double> zeta(1);
       nod_pt->get_coordinates_on_boundary(b,zeta);
       //So now I can work out the velocities
       nod_pt->set_value(0,-Problem_Parameter::A*sin(zeta[0])*chi_dot);
       nod_pt->set_value(1,Problem_Parameter::B*cos(zeta[0])*chi_dot);
       //Unset the auxilliary node update function
       nod_pt->set_auxiliary_node_update_fct_pt(0);
      }
    }
    }*/


 ///Pin the degrees of freedom associated with the solid bodies
 void pin_rigid_body()
  {
   unsigned n_rigid_body = Rigid_body_pt.size();
   for(unsigned i=0;i<n_rigid_body;++i)
    {
     unsigned n_geom_data = Rigid_body_pt[i]->ngeom_data();
     for(unsigned j=0;j<n_geom_data;j++)
      {
       Rigid_body_pt[i]->geom_data_pt(j)->pin_all();
      }
    }
  }

 ///Pin the degrees of freedom associated with the solid bodies
 void pin_rigid_body_position()
  {
   unsigned n_rigid_body = Rigid_body_pt.size();
   for(unsigned i=0;i<n_rigid_body;++i)
    {
     for(unsigned j=0;j<2;j++)
      {
       Rigid_body_pt[i]->geom_data_pt(0)->pin(j);
      }
    }
  }



 ///Unpin the degrees of freedom associated with the solid bodies
 void unpin_rigid_body()
  {
   unsigned n_rigid_body = Rigid_body_pt.size();
   for(unsigned i=0;i<n_rigid_body;++i)
    {
     unsigned n_geom_data = Rigid_body_pt[i]->ngeom_data();
     for(unsigned j=0;j<n_geom_data;j++)
      {
       Rigid_body_pt[i]->geom_data_pt(j)->unpin_all();
      }
    }
  }


 /// Doc the solution
 void doc_solution(const std::string& comment="", const bool& project=false);
 
 /// Compute the error estimates and assign to elements for plotting
 void compute_error_estimate(double& max_err,
                             double& min_err);
  
 /// Sanity check: Doc boundary coordinates from mesh and GeomObject
 void doc_boundary_coordinates();
  

private:
 

 /// \short Create elements that enforce prescribed boundary motion
 /// for the pseudo-solid fluid mesh by Lagrange multipliers
 void create_lagrange_multiplier_elements();

 /// \short Delete elements that impose the prescribed boundary displacement
 /// and wipe the associated mesh
 void delete_lagrange_multiplier_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Lagrange_multiplier_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();
   
  } // end of delete_lagrange_multiplier_elements

 /// \short Create elements that calculate the drag and torque on
 /// the boundaries
 void create_drag_elements();
 
 /// \short Delete elements that calculate the drag and torque on the 
 /// boundaries
 void delete_drag_elements()
  {
   unsigned n_bodies = Drag_mesh_pt.size();
   for(unsigned n=0;n<n_bodies;n++)
    {
     // How many surface elements are in the surface mesh
     unsigned n_element = Drag_mesh_pt[n]->nelement();
     
     // Loop over the surface elements
     for(unsigned e=0;e<n_element;e++)
      {
       // Kill surface element
       delete Drag_mesh_pt[n]->element_pt(e);
      }
     
     // Wipe the mesh
     Drag_mesh_pt[n]->flush_element_and_node_storage();
    }
  } // end of delete_drag_elements

 
 /// Pointers to mesh of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

public: 
 /// Pointers to Fluid_mesh
 RefineableSolidTriangleMesh<ELEMENT>* Fluid_mesh_pt;
 
 /// Triangle mesh polygon for outer boundary 
 TriangleMeshPolygon* Outer_boundary_polyline_pt; 

 /// Mesh of drag elements
 Vector<Mesh*> Drag_mesh_pt;

 /// Mesh of the generalised elements for the rigid bodies
 Mesh* Rigid_body_mesh_pt;

 /// Storage for the geom object
 Vector<GeomObject*> Rigid_body_pt;

}; // end_of_problem_class


//==start_constructor=====================================================
/// Constructor: build the first mesh with TriangleMeshPolygon and
///              TriangleMeshHolePolygon object
//========================================================================
template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
 // Allow for rough startup
 this->Problem::Max_residuals=1000.0;

 // Output directory
 Problem_Parameter::Doc_info.set_directory("RESLT");

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 this->add_time_stepper_pt(new BDF<2>);

 // Allocate a timestepper for the rigid body
 this->add_time_stepper_pt(new Newmark<2>);

 // Define the boundaries: Polyline with 4 different
 // boundaries for the outer boundary and 2 internal holes, 
 // egg shaped, with 2 boundaries each
 
 // Build the boundary segments for outer boundary, consisting of
 //--------------------------------------------------------------
 // four separeate polyline segments
 //---------------------------------
 Vector<TriangleMeshPolyLine*> boundary_segment_pt(4);
 
 // Initialize boundary segment
 Vector<Vector<double> > bound_seg(2);
 for(unsigned i=0;i<2;i++)
  {
   bound_seg[i].resize(2);
  }
 
 //Set the length of the channel
 double half_length = 5.0;
 double half_height = 5.0;

 // First boundary segment
 bound_seg[0][0]=-half_length;
 bound_seg[0][1]=-half_height;
 bound_seg[1][0]=-half_length;
 bound_seg[1][1]=half_height;
 
 // Specify 1st boundary id
 unsigned bound_id = 1;

 // Build the 1st boundary segment
 boundary_segment_pt[0] = new TriangleMeshPolyLine(bound_seg,bound_id);
 
 // Second boundary segment
 bound_seg[0][0]=-half_length;
 bound_seg[0][1]=half_height;
 bound_seg[1][0]=half_length;
 bound_seg[1][1]=half_height;

 // Specify 2nd boundary id
 bound_id = 2;

 // Build the 2nd boundary segment
 boundary_segment_pt[1] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Third boundary segment
 bound_seg[0][0]=half_length;
 bound_seg[0][1]=half_height;
 bound_seg[1][0]=half_length;
 bound_seg[1][1]=-half_height;

 // Specify 3rd boundary id
 bound_id = 3;

 // Build the 3rd boundary segment
 boundary_segment_pt[2] = new TriangleMeshPolyLine(bound_seg,bound_id);

 // Fourth boundary segment
 bound_seg[0][0]=half_length;
 bound_seg[0][1]=-half_height;
 bound_seg[1][0]=-half_length;
 bound_seg[1][1]=-half_height;

 // Specify 4th boundary id
 bound_id = 4;

 // Build the 4th boundary segment
 boundary_segment_pt[3] = new TriangleMeshPolyLine(bound_seg,bound_id);
  
 // Create the triangle mesh polygon for outer boundary using boundary segment
 Outer_boundary_polyline_pt = new TriangleMeshPolygon(boundary_segment_pt);


 // Now deal with the moving holes
 //-------------------------------

 // We have one rigid body
 Rigid_body_pt.resize(1);

 // Build first hole
 //-----------------
 double x_center = 0.0;
 double y_center = 0.0;
 double A = Problem_Parameter::A;//0.25;
 double B = Problem_Parameter::B; //0.5
 GeomObject* temp_hole_pt = new GeneralEllipse(x_center,y_center,A,B);
 Rigid_body_pt[0] = new RigidBodyElement(temp_hole_pt,
                                         this->time_stepper_pt(1));
   
 //Now set the split coordinates
 Vector<Vector<double> > split_coord(2);
 split_coord[0].resize(2);
 split_coord[0][0] = 0.0;
 split_coord[0][1] = MathematicalConstants::Pi;
 split_coord[1].resize(2);
 split_coord[1][0] = MathematicalConstants::Pi;
 split_coord[1][1] = 2.0*MathematicalConstants::Pi;
 
 // Now build the mesh, based on the boundaries specified by
 //---------------------------------------------------------
 // polygons just created
 //----------------------
 double uniform_element_area=0.2;
 Fluid_mesh_pt = 
  new RefineableSolidTriangleMesh<ELEMENT>(Outer_boundary_polyline_pt, 
                                           Rigid_body_pt,
                                           split_coord,
                                           uniform_element_area,
                                           this->time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;


 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=0.005;
 Fluid_mesh_pt->min_permitted_error()=0.001; 
 Fluid_mesh_pt->max_element_size()=0.2;
 Fluid_mesh_pt->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   Fluid_mesh_pt->min_element_size()=0.01; 
  }

 // Set the problem pointer
 Fluid_mesh_pt->problem_pt()=this;
   
 // Output boundary and mesh
 this->Fluid_mesh_pt->output_boundaries("boundaries.dat");
 this->Fluid_mesh_pt->output("mesh.dat");
   

 // Set boundary condition, assign auxiliary node update fct,
 // complete the build of all elements, attach power elements that allow
 // computation of drag vector
 complete_problem_setup();
 
 // Create Lagrange multiplier mesh for boundary motion
 //----------------------------------------------------
 // Construct the mesh of elements that enforce prescribed boundary motion
 // of pseudo-solid fluid mesh by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();
 
 //Set the parameters of the rigid body elements
 RigidBodyElement* rigid_element1_pt = 
  dynamic_cast<RigidBodyElement*>(Rigid_body_pt[0]);
 rigid_element1_pt->initial_centre_of_mass(0) = x_center;
 rigid_element1_pt->initial_centre_of_mass(1) = y_center; 
 rigid_element1_pt->mass() = MathematicalConstants::Pi*A*B;
 rigid_element1_pt->moment_of_inertia() = 
  0.25*MathematicalConstants::Pi*A*B*(A*A + B*B);
 rigid_element1_pt->g_pt() = &Problem_Parameter::G;

 // Create the drag mesh for the rigid bodies
 Drag_mesh_pt.resize(1);
 for(unsigned m=0;m<1;m++) {Drag_mesh_pt[m] = new Mesh;}
 this->create_drag_elements();

 //Add the drag meshes to the appropriate rigid bodies
 rigid_element1_pt->drag_mesh_pt() = Drag_mesh_pt[0];

 //rigid_element1_pt->external_torque_fct_pt() = imposed_torque;

 // Create the mesh for the rigid bodies
 Rigid_body_mesh_pt = new Mesh;
 Rigid_body_mesh_pt->add_element_pt(rigid_element1_pt);

 // Combine meshes
 //---------------
 
 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Lagrange_multiplier sub meshes
 this->add_sub_mesh(this->Lagrange_multiplier_mesh_pt);

 this->add_sub_mesh(this->Rigid_body_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();
  
 // Sanity check: Doc boundary coordinates from mesh and GeomObject
 doc_boundary_coordinates();
  
 // Setup equation numbering scheme
 cout <<"Number of equations: " << this->assign_eqn_numbers() << std::endl;
 
} // end_of_constructor




//============start_doc_solid_zeta=======================================
/// Doc boundary coordinates in mesh and plot GeomObject representation
/// of inner boundary.
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_boundary_coordinates()
{

 ofstream some_file;
 char filename[100];
 
 // Number of plot points (should this be increased)
 unsigned npoints = 5;
 
 // Output solution and projection files
 sprintf(filename,"RESLT/inner_hole_boundary_from_geom_obj.dat");
 some_file.open(filename);

 //Initialize zeta and r
 Vector<double> zeta(1);
 zeta[0]=0;
 
 Vector<double> r(2);
 r[0]=0;
 r[1]=0;
   
 // Get the boundary geometric objects associated with the fluid mesh
 unsigned n_boundary = Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary;b++)
  {
   //The geometric object associated with the b-th boundary
   GeomObject* boundary_geom_obj_pt = 
    Fluid_mesh_pt->boundary_geom_object_pt(b);
   
   //Only bother to do anything if there is a geometric object
   if(boundary_geom_obj_pt!=0)
    {
     // Zone label 
     some_file <<"ZONE T=boundary"<<b<<std::endl;
     
     //Get the coordinate limits
     Vector<double> zeta_limits = 
      Fluid_mesh_pt->boundary_coordinate_limits(b);
     
     //Set the increment
     double zeta_inc = 
      (zeta_limits[1] - zeta_limits[0])/(double)(npoints-1);
     
     for(unsigned i=0;i<npoints;i++)
      {
       // Get coordinate
       zeta[0] = zeta_limits[0] + i*zeta_inc;
       //Get the position
       boundary_geom_obj_pt->position(zeta,r);
       
       // Print it
       some_file <<r[0]<<" "<<r[1]<<" "<<zeta[0]<<std::endl;  
      }
    }
  }
 some_file.close();

// Doc boundary coordinates using Lagrange_multiplier_mesh_pt
std::ofstream the_file("RESLT/inner_hole_boundary_from_mesh.dat");

 // Initialise max/min boundary coordinate
 double zeta_min= DBL_MAX;
 double zeta_max=-DBL_MAX;

 // Loop over Lagrange_multiplier elements
 unsigned n_face_element = this->Lagrange_multiplier_mesh_pt->nelement();
 
 for(unsigned e=0;e<n_face_element;e++)
  {
   
   //Cast the element pointer
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt=
    dynamic_cast< ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(e));

   // Doc boundary coordinate
   Vector<double> s(1);
   Vector<double> zeta(1);
   Vector<double> x(2);
   unsigned n_plot=5;

   the_file << el_pt->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {         
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,n_plot,s);         
     el_pt->interpolated_zeta(s,zeta);
     el_pt->interpolated_x(s,x);
     for (unsigned i=0;i<2;i++)
      {
       the_file << x[i] << " ";
      }
     the_file << zeta[0] << " ";

     // Update max/min boundary coordinate
     if (zeta[0]<zeta_min) zeta_min=zeta[0];
     if (zeta[0]>zeta_max) zeta_max=zeta[0];

     the_file << std::endl;
    }
  }
 // Close doc file
 the_file.close();
 
  
} //end doc_solid_zeta

//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
/// for the pseudo-solid fluid mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_lagrange_multiplier_elements()
{ 
 // The idea is to apply Lagrange multipliers to the boundaries in 
 // the mesh that have associated geometric objects

 //Find the number of boundaries
 unsigned n_boundary = Fluid_mesh_pt->nboundary();

 // Loop over the boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   //Get the geometric object associated with the boundary
   GeomObject* boundary_geom_obj_pt = 
    Fluid_mesh_pt->boundary_geom_object_pt(b);

   //Only bother to do anything if there is a geometric object
   if(boundary_geom_obj_pt!=0)
    {
     // How many bulk fluid elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk fluid elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk fluid element that is 
       // adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element. Note that we use different Lagrange
       // multiplier fields for each distinct boundary (here indicated
       // by b.
       ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt =
        new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
         bulk_elem_pt,face_index,b);   
       
       // Add it to the mesh
       Lagrange_multiplier_mesh_pt->add_element_pt(el_pt);
       
       // Set the GeomObject that defines the boundary shape and set
       // which bulk boundary we are attached to (needed to extract
       // the boundary coordinate from the bulk nodes)
       el_pt->set_boundary_shape_geom_object_pt(
        boundary_geom_obj_pt,b);
       
       // Loop over the nodes to pin Lagrange multiplier
       unsigned nnod=el_pt->nnode();
       for(unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt = el_pt->node_pt(j);
         
         // How many nodal values were used by the "bulk" element
         // that originally created this node?
         unsigned n_bulk_value=el_pt->nbulk_value(j);
         
         // Pin two of the four Lagrange multipliers at vertices
         // This is not totally robust, but will work in this application
         unsigned nval=nod_pt->nvalue();
         if (nval==7)
          {
           for (unsigned i=0;i<2;i++) 
            { 
             // Pin lagrangian multipliers
             nod_pt->pin(n_bulk_value+2+i);
            }
          }
        }
      } // end loop over the element
    } //End of  case if there is a geometric object
  } //End of loop over boundaries
}
// end of create_lagrange_multiplier_elements



//============start_of_create_drag_elements===============
/// Create elements that calculate the drag and torque on
/// the obstacles in the fluid mesh
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_drag_elements()
{ 
 //Once again the idea is only to attach drag elements to those
 //boundaries with associated geometric objects (and therefore holes)
 
 // Get the number of rigid bodies
 unsigned n_rigid = Rigid_body_pt.size();

 // Get the number of boundaries in the mesh
 unsigned n_boundary = Fluid_mesh_pt->nboundary();

 //Loop over the rigid bodies
 for(unsigned r=0;r<n_rigid;r++)
  {
   //Allocate storage for all geometric data
   std::set<Data*> bulk_geometric_data_pt;
   //Allocate storage for all load data
   std::set<std::pair<Data*,unsigned> > bulk_load_data_pt;
   
   //Get the rigid boundary geometric object
   RigidBodyElement* rigid_el_pt = 
    dynamic_cast<RigidBodyElement*>(this->Rigid_body_pt[r]);

   //Flush any exisiting external data of the rigid body data
   rigid_el_pt->flush_external_data();
   
   // Loop over all boundaries
   for(unsigned b=0;b<n_boundary;b++)
    {
     //Does the boundary correspond to the current rigid body
     if(dynamic_cast<RigidBodyElement*>
        (Fluid_mesh_pt->boundary_geom_object_pt(b)) == rigid_el_pt)
      {
       // How many bulk fluid elements are adjacent to boundary b?
       unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
       
       // Loop over the bulk fluid elements adjacent to boundary b?
       for(unsigned e=0;e<n_element;e++)
        {
         // Get pointer to the bulk fluid element that is 
         // adjacent to boundary b
         ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
          Fluid_mesh_pt->boundary_element_pt(b,e));
       
         //Find the index of the face of element e along boundary b
         int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
         
         // Create new element. Note that we use different Lagrange
         // multiplier fields for each distinct boundary (here indicated
         // by b.
         NavierStokesSurfaceDragTorqueElement<ELEMENT>* el_pt =
          new NavierStokesSurfaceDragTorqueElement<ELEMENT>(
           bulk_elem_pt,face_index);   
         
         //Add the geometric and load data
         bulk_elem_pt->identify_geometric_data(bulk_geometric_data_pt);
         bulk_elem_pt->identify_load_data(bulk_load_data_pt);
         
         // Add it to the mesh
         Drag_mesh_pt[r]->add_element_pt(el_pt);
         
       //Set the original centre of rotation and the geometric data
         for(unsigned i=0;i<2;i++) 
          {el_pt->centre_of_rotation(i) = 
            rigid_el_pt->initial_centre_of_mass(i);}
         
         //Set the translation as well
         el_pt->set_translation_and_rotation(rigid_el_pt->geom_data_pt(0));
        } // end loop over the element
      }
    } //End of loop over boundaries
 
   //How much data do we have
   std::cout << bulk_geometric_data_pt.size() << " geom data\n";
   std::cout << bulk_load_data_pt.size() << " load data\n";
   
   //Need to add all these data as external data to the object as external data
   for(std::set<Data*>::iterator it = bulk_geometric_data_pt.begin();
       it!=bulk_geometric_data_pt.end();++it)
    {
     rigid_el_pt->add_external_data(*it);
    }
   
   //Now do the same but make custom data for the load data
   for(std::set<std::pair<Data*,unsigned> >::iterator it = 
        bulk_load_data_pt.begin();
       it!=bulk_load_data_pt.end();++it)
    {
     Data* temp_data_pt = new HijackedData(it->second,it->first);
     rigid_el_pt->add_external_data(temp_data_pt);
    }
  } // end loop over the rigid bodies
}
// end of create_drag_elements


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_solution(
 const std::string& comment,
 const bool& project)
{ 

 oomph_info << "Docing step: " << Problem_Parameter::Doc_info.number()
            << std::endl;

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 


 // Compute errors and assign to each element for plotting
 double max_err;
 double min_err;
 compute_error_estimate(max_err,min_err);
 
 // Output solution and projection files
 if(!project)
  {
   sprintf(filename,"%s/soln%i.dat",
           Problem_Parameter::Doc_info.directory().c_str(),
           Problem_Parameter::Doc_info.number());
  }
 else
  {
   sprintf(filename,"%s/proj%i.dat",
           Problem_Parameter::Doc_info.directory().c_str(),
           Problem_Parameter::Doc_info.number()-1);
  }


 // Assemble square of L2 norm 
 double square_of_l2_norm=0.0;
 unsigned nel=Fluid_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   square_of_l2_norm+=
    dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(e))->
    square_of_l2_norm();
  }
 std::cout << "Checking " << sqrt(square_of_l2_norm) << "\n";
 Problem_Parameter::Norm_file << sqrt(square_of_l2_norm) << "\n";
 

 some_file.open(filename);
 some_file << dynamic_cast<ELEMENT*>(this->Fluid_mesh_pt->element_pt(0))
  ->variable_identifier();
 this->Fluid_mesh_pt->output(some_file,npts);   
 some_file << "TEXT X = 25, Y = 78, CS=FRAME T = \"Global Step " 
           << Problem_Parameter::Doc_info.number() << "  " 
           << comment << "\"\n";
 some_file.close();

 // No trace file writing after projection
 if(project) return;

 // Get max/min area
 double max_area;
 double min_area;
 Fluid_mesh_pt->max_and_min_area(max_area, min_area);

 // Write trace file
 Problem_Parameter::Trace_file 
  << this->time_pt()->time() << " " 
  << Fluid_mesh_pt->nelement() << " "
  << max_area << " "
  << min_area << " "
  << max_err << " "
  << min_err << " "
  << sqrt(square_of_l2_norm) << " "
  << std::endl;

 //Output the motion of the centre of gravity
 dynamic_cast<RigidBodyElement*>(this->Rigid_body_pt[0])->
  output_centre_of_gravity(Problem_Parameter::Cog_file);

 // Increment the doc_info number
 Problem_Parameter::Doc_info.number()++;



}

//========================================================================
/// Compute error estimates and assign to elements for plotting
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::compute_error_estimate(double& max_err,
                                                               double& min_err)
{ 
 // Get error estimator
 ErrorEstimator* err_est_pt=Fluid_mesh_pt->spatial_error_estimator_pt();
 
 // Get/output error estimates
 unsigned nel=Fluid_mesh_pt->nelement();
 Vector<double> elemental_error(nel);
 
 // We need a dynamic cast, get_element_errors from the Fluid_mesh_pt
 // Dynamic cast is used because get_element_errors require a Mesh* ans
 // not a SolidMesh*
 Mesh* fluid_mesh_pt=dynamic_cast<Mesh*>(Fluid_mesh_pt);
 err_est_pt->get_element_errors(this->communicator_pt(),
                                fluid_mesh_pt,
                                elemental_error);

 // Set errors for post-processing and find extrema
 max_err=0.0;
 min_err=DBL_MAX;
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<MyTaylorHoodElement*>(Fluid_mesh_pt->element_pt(e))->
    set_error(elemental_error[e]);

   max_err=std::max(max_err,elemental_error[e]);
   min_err=std::min(min_err,elemental_error[e]);
  }
  
}


//============================================================
///Driver code for moving block problem
//============================================================
int main(int argc, char **argv)
{
 // feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 Problem_Parameter::G.resize(2);
 Problem_Parameter::G[0]= 0.0;
 Problem_Parameter::G[1] = 0.0;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation?
 CommandLineArgs::specify_command_line_flag("--validation");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 Problem_Parameter::Constitutive_law_pt = 
  new GeneralisedHookean(&Problem_Parameter::Nu);
 
 // Open trace file
 Problem_Parameter::Trace_file.open("RESLT/trace.dat");
 
 // Open norm file
 Problem_Parameter::Norm_file.open("RESLT/norm.dat");

 // Open file to trace the centre of gravity
 Problem_Parameter::Cog_file.open("RESLT/cog_trace.dat");

 // Create problem in initial configuration
 UnstructuredFluidProblem<ProjectableTaylorHoodElement<MyTaylorHoodElement> > 
  problem;  
 
 //Let's pin the rigid body dofs initially
 problem.pin_rigid_body();
 problem.assign_eqn_numbers();
 //Do a steady solve to map the nodes to the boundary of the ellipse
 problem.steady_newton_solve();

 //Now unpin the rigid body
 problem.unpin_rigid_body();
 //but repin the position of the centre of mass
 problem.pin_rigid_body_position();
 problem.assign_eqn_numbers();
 
 // Initialise timestepper
 double dt=0.05;//0.05;
 problem.initialise_dt(dt);
 
 // Perform impulsive start (this could be a very bad idea)
 problem.assign_initial_values_impulsive();

 //Now set the values for the rigid body and boundary nodes 
 //from the Jefferey solution
// problem.set_consistent_rigid_body_history();

 // Output initial conditions
 problem.doc_solution();

 //Set the velocity
 Problem_Parameter::Alpha=1.0;
 problem.set_boundary_velocity();

 // Solve problem a few times on given mesh
 unsigned nstep=3;
 for (unsigned i=0;i<nstep;i++)
  {
   // Solve the problem
   problem.unsteady_newton_solve(dt);    
   problem.doc_solution();
  }

 // Now do a couple of adaptations
 unsigned ncycle=200;
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   ncycle=1;
   oomph_info << "Only doing one cycle during validation\n";
  }

 for (unsigned j=0;j<ncycle;j++)
  {       
   // Adapt
   problem.adapt();

   //problem.actions_before_implicit_timestep();
   //problem.newton_solve(dt);

   //problem.doc_solution("hmm");
   //exit(1);

   //Solve problem a few times
   for (unsigned i=0;i<nstep;i++)
    {     
     // Solve the problem
     problem.unsteady_newton_solve(dt);
     // Build the label for doc
     std::stringstream label;
     label << "Cycle " <<j << " Step "<< i;
     problem.doc_solution(label.str());
    }
  }

 // Close norm and trace files
 Problem_Parameter::Cog_file.close();
 Problem_Parameter::Norm_file.close();
 Problem_Parameter::Trace_file.close();

} //End of main
