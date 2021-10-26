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
/// Driver code for the 3D Bretherton problem in a square tube

// Standard includes
#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>

//Generic includes
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "merge_meshes.h"
#include "st_mesh.h"

//Inlet and fixed_spine elements
#include "extra_elements.h"

using namespace std;

using namespace oomph;
 

//==================================================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re= 0.0;
 
 /// Capillary number
 double Ca= 0.1;
 
 //Bond Nunber
 double Bo = 0.0;
 
 // Product of reynolds and Froude number
 double ReInvFr = Bo/Ca;
 
 // Direction of gravity
 Vector<double> G(3); 
 
 //Aspect ratio of the tube
 double Alpha = 1.0 ;

 //Height of the tube
 double Height = 1.0;
 
 double Length_tip = 0.50;
 
 double Length_liq = 1.0;

 double Length_can = 1.5;
 
 unsigned Ncan = 2;
 
 unsigned Ntip = 2;
 unsigned Nliq = 2;

 double Radius = 0.923;
 
 // Rated pressure
 double Rat_press = 2.87;
 
 /// Function that prescribes the hydrostatic pressure field at the outlet
 void hydrostatic_pressure(const double &time, const Vector<double> &x, 
                           const Vector<double> &n, Vector<double> &traction)
 {
  traction[0] = 0.0;
  traction[1] = ReInvFr*G[2]*x[2]; //Perpendicular to tube axis
  traction[2] = 0.0;
 }

}



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//====================================================================
/// Brethreton  problem in square tube domain
//====================================================================
template<class ELEMENT>
class ThreeDimBethertonProblem : public Problem
{

public:
 
 /// Constructor
 ThreeDimBethertonProblem();

 /// Destructor to clean up memory
 ~ThreeDimBethertonProblem();

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned int &e, const unsigned int &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }

 /// Remesh before any convergence checke
 void actions_before_newton_convergence_check()
  {
   mesh_pt()->node_update();
   
   // This driver code cannot be allowed to use the analytical form of
   // get_dresidual_dnodal_coordinates(...) that is implemented in the
   // NavierStokesEquations class, since the elemental residuals have
   // contributions from external data which is not taken into account
   // by that routine. We therefore force the bulk elements to use the
   // fully-finite differenced version.
   const unsigned n_element = mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
   {
    if(dynamic_cast<ElementWithMovingNodes*>(mesh_pt()->element_pt(e)))
  {
     ElementWithMovingNodes* el_pt =
       dynamic_cast<ElementWithMovingNodes*>(mesh_pt()->element_pt(e));
      el_pt->evaluate_shape_derivs_by_direct_fd();
     }
   }
  }


// Calulate the flow across the Outlet
 double get_outflow()
  {
   double flow = 0;
   unsigned int nel = Outlet_traction_element_pt.size();
   for(unsigned int i =0;i<nel;i++)
    {
     flow += dynamic_cast<SpineGravityTractionElement<Hijacked<ELEMENT> >* >
      (Outlet_traction_element_pt[i])->flow();
    }
   return flow;
  }

 // Calculate the finger width at the end of the cell
 double get_lambda()
  { 
   //Initial max
   double max_x = 0.0;
   //Loop over all in interface elements
   unsigned n_line_interface = Interface_line_element_pt.size();
   for(unsigned e=0;e<n_line_interface;e++)
    {
     //Cache the element
     FiniteElement* const element_pt = Interface_line_element_pt[e];
     //Loop over all nodes
     const unsigned n_node = element_pt->nnode();
     for(unsigned n=0;n<n_node;n++)
      {
       if(element_pt->node_pt(n)->x(0) > max_x) 
        {max_x = element_pt->node_pt(n)->x(0);}
      }
    }
   //Now scale by the half-width
   return max_x/(Global_Physical_Variables::Alpha/2);
  }
 
 
 //If mult_width = 0, it multiplu only the channel length
 void multiply_aspect_ratio(const double factor, const bool& mult_width)
  {
   // until now we have just worked with a square channel
   // lets change the aspect ratio
   // we will change also the y aspect ratio 
   // so that the tip is also streched in the y direction can be changed

   //Loop over all nodes in the mesh
   const unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     //Multiply the x-coordinate and y-coordinate by the specified factor
     Node* const node_pt =  mesh_pt()->node_pt(n);
     if(mult_width){node_pt->x(0) =  node_pt->x(0) * factor ;}
     node_pt->x(1) =  node_pt->x(1) * factor;
    }
   
   // Update also the origin of the spines to the new position
   // of the sidewalls
   const unsigned n_spine = mesh_pt()->nspine();
   for(unsigned s = 0;s<n_spine;s++)
    {
     Spine* spine_pt = mesh_pt()->spine_pt(s);
     const double x_spine = spine_pt->geom_parameter(0);
     const double y_spine = spine_pt->geom_parameter(1);
     if(mult_width){spine_pt->geom_parameter(0) = x_spine * factor;}
     spine_pt->geom_parameter(1) = y_spine * factor;
    }
   
   
   // Update the spines heights according to the new geometry
   // This is done by looping over all nodes on the free surface
   // and calculating the distance from the new side wall
   const unsigned n_boundary_node = mesh_pt()->nboundary_node(5);
   
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     //Get the root of the spine
     SpineNode* spine_node_pt =  
      dynamic_cast<SpineNode*>(mesh_pt()->boundary_node_pt(5,n));
     Spine* spine_pt = spine_node_pt->spine_pt();
     double x_spine = spine_pt->geom_parameter(0);
     double y_spine = spine_pt->geom_parameter(1);
     double z_spine = spine_pt->geom_parameter(2);
     
     //Find the distance of the node on the free surface from the root
     //of the spine
     double  distance = sqrt( 
      (spine_node_pt->x(0) -  x_spine )*(spine_node_pt->x(0) -  x_spine )+
      (spine_node_pt->x(1) -  y_spine)*(spine_node_pt->x(1) -  y_spine)+ 
      (spine_node_pt->x(2) -  z_spine)*(spine_node_pt->x(2) -  z_spine) );
     //Set the spine length to be the new distance
     spine_node_pt->h() = distance;
    }
   
   //Change also the distance of the fixed spine by the factor
   Fixed_spine_height =  Fixed_spine_height *factor;
   // At the end we change the value of Alpha
   if(mult_width){Global_Physical_Variables::Alpha =  
                   Global_Physical_Variables::Alpha * factor;}
   if(!mult_width){ Global_Physical_Variables::Radius =  
                     Global_Physical_Variables::Radius *factor;}
   Global_Physical_Variables::Length_tip =  
    Global_Physical_Variables::Length_tip * factor;
   Global_Physical_Variables::Length_can =  
    Global_Physical_Variables::Length_can * factor;
   Global_Physical_Variables::Length_liq =  
    Global_Physical_Variables::Length_liq * factor; 
  }
 
 //Access function for the specific mesh
 STSpineMesh<Hijacked<ELEMENT>, 
             SpineSurfaceFluidInterfaceElement<ELEMENT> >* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast< STSpineMesh<Hijacked<ELEMENT>, 
    SpineSurfaceFluidInterfaceElement<ELEMENT> >* >(Problem::mesh_pt());
  }


 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:
 
 //Vector of pointers to the inlet traction elements
 Vector<FiniteElement*> Inlet_traction_element_pt;
 
 //Vector of pointers to the outlet traction elements
 Vector<FiniteElement*> Outlet_traction_element_pt;
 
//Vector of pointers to the interface line elements
 Vector<FiniteElement*> Interface_line_element_pt;
 
 //Pointer to the element that fixes the spine height
 FiniteElement* Fix_spine_height_element_pt;

//Data value that represents the bubble pressure
 Data* Bubble_pressure_data_pt;

public:
 
//Fixed height value
 double Fixed_spine_height;

 private:

/// Trace file
 ofstream Trace_file;

// Trace Spine poniters

 Spine* Trace_spine_r_h_top;

 Spine* Trace_spine_r_h_middle;

 Spine* Trace_spine_r_h_bottom;
 
 Spine* Trace_spine_r_d_top;

 Spine* Trace_spine_r_d_bottom;

};


//========================================================================
/// Constructor for RectangularDrivenCavity problem 
//========================================================================
template<class ELEMENT>
ThreeDimBethertonProblem<ELEMENT>::ThreeDimBethertonProblem()
{ 
 
 double alpha = Global_Physical_Variables::Alpha;
 double height = Global_Physical_Variables::Height;
 double length_tip =  Global_Physical_Variables::Length_tip;
 double length_liq =  Global_Physical_Variables::Length_liq;
 double length_can =  Global_Physical_Variables::Length_can;
 double radius =  Global_Physical_Variables::Radius;
 double rat_press =  Global_Physical_Variables::Rat_press;
 this->Max_residuals = 10.0;
 this->Max_newton_iterations = 15;

 unsigned  ncan = Global_Physical_Variables::Ncan;
 unsigned  ntip = Global_Physical_Variables::Ntip;
 unsigned  nliq = Global_Physical_Variables::Nliq;

// Build and assign mesh
 Problem::mesh_pt() = new  
  STSpineMesh<Hijacked<ELEMENT>, SpineSurfaceFluidInterfaceElement<ELEMENT>  >
  (ncan,ntip,nliq, alpha,length_can,length_tip,length_liq,height,radius);
 
  cout<<"Nnodes boundary 5: "<<mesh_pt()->nboundary_node(5)<<endl;
  cout<<"Nnodes spines : "<<mesh_pt()->nspine()<<endl;
  
  // Complete the build of all elements so they are fully functional
  
  // Set the fixed spine element
  unsigned num_nodes = mesh_pt()->nnode();
  SpineNode* spinenode_fixed=0;
  
  //Look fot the spine
  for(unsigned inod = 0;inod< num_nodes;inod++)
   {
    SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
    double height = Global_Physical_Variables::Height;
    double length_tip =  Global_Physical_Variables::Length_tip;
    double distance =
     sqrt(  (spnd_pt->x(0))* (spnd_pt->x(0)) + 
            (spnd_pt->x(1) + (length_tip * radius))* 
            (spnd_pt->x(1) + (length_tip*radius)) + 
            (spnd_pt->x(2) - height/2)* (spnd_pt->x(2) - height/2) );
    
    if(distance < 10E-8)
     {
      cout<<"Spine to be pinned found"<<endl;   
      spinenode_fixed = spnd_pt;
     }
   }
  
  //Create a fixed element using the central spine
  FixSpineHeightElement* fix_spine_element_pt = 
   new FixSpineHeightElement(spinenode_fixed);
  
  //Set the fixed spine height
  Fixed_spine_height = spinenode_fixed->spine_pt()->height();
  
  //Add the fixed height to it
  fix_spine_element_pt->height_pt() = &Fixed_spine_height;
  
  //Store the element in the problem
  Fix_spine_height_element_pt = fix_spine_element_pt;

  //Create a bubble presure data
  Bubble_pressure_data_pt = new Data(1);
  
  //Set values
  Bubble_pressure_data_pt->set_value(0,rat_press);
  
  //This will be global data
  add_global_data(Bubble_pressure_data_pt);
  
  //Set the pressure data
  fix_spine_element_pt->set_traded_pressure_data(Bubble_pressure_data_pt);
  
  //Add the Fixed element to the mesh
  mesh_pt()->add_element_pt(fix_spine_element_pt);
  
  cout<<"Bubble pressure "<<Bubble_pressure_data_pt->value(0)<<endl;
  
  //Find number of bulk elements in mesh
  unsigned n_bulk = mesh_pt()->nbulk();
  
  // Loop over the elements to set up element-specific 
  // things that cannot be handled by constructor
  for(unsigned e=0;e<n_bulk;e++)
   {
    // Upcast from GeneralisedElement to the present element
    Hijacked<ELEMENT>* el_pt = 
     dynamic_cast<Hijacked<ELEMENT>* >(mesh_pt()->bulk_element_pt(e));
    
    //Set the Reynolds number
    el_pt->re_pt() = &Global_Physical_Variables::Re;
    
    //Set the  ReInvFr number
    el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
    
    //Set the gravity
    el_pt->g_pt() = &Global_Physical_Variables::G;

    //Insist that shape derivatives are calculated by "full" finite
    //differences (which is the silent default)
    //el_pt->evaluate_shape_derivs_by_direct_fd();
   }

  
  //Find number of interface elements in mesh
  unsigned n_interface = mesh_pt()->ninterface_element();
  unsigned nhi = 0;
  
  for(unsigned e=0;e<n_interface;e++)
   {
    // Upcast from GeneralisedElement to the present element
    SpineSurfaceFluidInterfaceElement<ELEMENT>* el_pt = 
     dynamic_cast<SpineSurfaceFluidInterfaceElement<ELEMENT>*>
     (mesh_pt()->interface_element_pt(e));
    el_pt->ca_pt() =  &Global_Physical_Variables::Ca;
    
    // Set the external pressure data
    el_pt->set_external_pressure_data(Bubble_pressure_data_pt);
    
    // Quite inefficient way of Hijack the nodes
    unsigned nnodes = el_pt->nnode();
    for(unsigned i =0; i<nnodes;i++)
     {
      if( ( (el_pt->node_pt(i)->x(1)- length_can)*
            (el_pt->node_pt(i)->x(1)- length_can) )< 10E-10 )
       {
        el_pt->hijack_nodal_value(i,1);
        nhi++;
       }
     }

    //Insist that shape derivatives are calculated by "full" finite
    //differences (which is the silent default)
    //el_pt->evaluate_shape_derivs_by_direct_fd();
   }
  
  cout<<nhi<<" hijacked nodes at the surface elements."<<endl;
  
  
  //Create and set the inlet elements 
  unsigned long ninlet = mesh_pt()->nbulkinlet();
  for(unsigned long i = 0;i<ninlet;i++)
   {
    SpineElement<NavierStokesTractionElement<Hijacked<ELEMENT> > >* 
     inlet_element_pt = 
     new  SpineElement<NavierStokesTractionElement<Hijacked<ELEMENT> > >(
      mesh_pt()->bulk_inlet_element_pt(i),
      mesh_pt()->face_index_inlet());
    
    //Add Elements to the local storage
    Inlet_traction_element_pt.push_back(inlet_element_pt);
    
    //Set the traction function
    inlet_element_pt->traction_fct_pt() = 
     &Global_Physical_Variables::hydrostatic_pressure;

     //Insist that shape derivatives are calculated by "full" finite
    //differences (which is the silent default)
    //inlet_element_pt->evaluate_shape_derivs_by_direct_fd();

    
    //Add element to the mesh
    mesh_pt()->add_element_pt(inlet_element_pt);
   }
  
  //Create and set the outlet elements 
  unsigned long noutlet = mesh_pt()->nbulkoutlet();
  for(unsigned long i = 0;i<noutlet;i++)
   {
    SpineGravityTractionElement<Hijacked<ELEMENT> >* outlet_element_pt = 
     new  SpineGravityTractionElement<Hijacked<ELEMENT> >(
      mesh_pt()->bulk_outlet_element_pt(i),
      mesh_pt()->face_index_outlet());   
    
    //Set the  ReInvFr number
    outlet_element_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
    
    //Set the gravity
    outlet_element_pt->g_pt() = &Global_Physical_Variables::G;
    
    //Add Elements to the local storage
    Outlet_traction_element_pt.push_back(outlet_element_pt);

    //Insist that shape derivatives are calculated by "full" finite
    //differences (which is the silent default)
    //outlet_element_pt->evaluate_shape_derivs_by_direct_fd();
    
    //Add element to the mesh
    mesh_pt()->add_element_pt(outlet_element_pt);
   }
  

  //Create and set the line elements elements 
  //(the coordinates s constants are the same as in the outlet elements
  unsigned long nlineel = mesh_pt()->ninterfaceline();
  for(unsigned long i = 0;i<nlineel;i++)
   {
    //Cast the return type to  the edge element and also 
    //cast the interface element to the specific interface element
    SpineLineFluidInterfaceBoundingElement<ELEMENT>* line_element_pt = 
     dynamic_cast<SpineLineFluidInterfaceBoundingElement<ELEMENT>*>(
     dynamic_cast<SpineSurfaceFluidInterfaceElement<ELEMENT>*>
     ( mesh_pt()->interface_line_element_pt(i))->make_bounding_element(
      mesh_pt()->face_index_outlet()));
     
    //Set the capillary number
    line_element_pt->ca_pt() =  &Global_Physical_Variables::Ca;
    
    //Add Elements to the local storagea
    Interface_line_element_pt.push_back(line_element_pt);
    
    //Add element to the mesh
    mesh_pt()->add_element_pt(line_element_pt);
    
    //Hijack the nodes
    for(unsigned e=0;e<line_element_pt->nnode();e++)
     line_element_pt ->hijack_nodal_value(e,1);
   }

  unsigned n_element = mesh_pt()->nelement();
  cout<<"We have "<<n_bulk<<" bulk elemens, "<<n_interface<<" interface elements,  "<<n_element<<" generic elements in this mesh"<<endl;
  
  mesh_pt()->node_update();
  
  //Now set the boundary conditions by pinning those variables for 
  //which we apply Dirichlet boundary conditions
  
// Just pin the boundary Dirichlet conditions
  unsigned ibound = 0;
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);

  //Base of the tube (pin all values)
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Loop over values (ux uy and uz velocities)
    for (unsigned i=0;i<3;i++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
     }
   }
  
  //Inlet boundary (pin the transverse velocity components)
  ibound = 1;
  num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Loop over values (ux uy and uz velocities)
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
    mesh_pt()->boundary_node_pt(ibound,inod)->pin(2);  
   }

  //Outer wall of the tube (pin all velocitites)
  ibound = 2;
  num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Loop over values (ux uy and uz velocities)
    for (unsigned i=0;i<3;i++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
     }
   }
  
 //Outlet (leave free)
 //ibound = 3;
 //num_nod= mesh_pt()->nboundary_node(ibound);
 //for (unsigned inod=0;inod<num_nod;inod++)
 // {
 // }

 //Symmetry boundary (apply symmetry conditions in x)
 ibound = 4;   
 num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);  
  }

 //Free surface (leave free
 //ibound = 5;
 //num_nod= mesh_pt()->nboundary_node(ibound);
 //for (unsigned inod=0;inod<num_nod;inod++)
 // {
 // }
     
 //Top of the tube (pin all velocities)
 ibound = 6;
 num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Loop over values (ux uy and uz velocities)
   for (unsigned i=0;i<3;i++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
    }
  }
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;    
 
//Set initial guess that all velocities are uniform in the y direction
 const unsigned n_node = mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   mesh_pt()->node_pt(n)->set_value(0,0.0);
   mesh_pt()->node_pt(n)->set_value(1,1.0);
   mesh_pt()->node_pt(n)->set_value(2,0.0);
  }
   
 //Set boundary conditions
 //Set the no-slip conditions on all tube boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //If we are not on the symmetry boundary (4)
   //Outlet (3), Inlet(1) or free surface (5), set the no-slip conditions
   if((ibound != 4) && (ibound != 3) && (ibound != 1) && (ibound !=5) )
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {     
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,1.0);
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(2,0.0);
      }
    }
  }

 //Prepare output information

//Open trace file
 char filename[100];
 sprintf(filename,"RESLT/trace.dat");
 Trace_file.open(filename);
 Trace_file << "VARIABLES=";
 Trace_file << "\"Ca\",\t \"Bo\", \t \"Alpha\",  \t \"m\", \t \"p<sub>bubble</sub>\", \t \"r<sub>h_top</sub>\", \t \"r<sub>d_top</sub>\", \t \"r<sub>h_middle</sub>\", \t \"r<sub>d_bottom</sub>\", \t \"r<sub>h_bottom</sub>\",\t \"lambda\"  ";
 
 Trace_file<<endl;

//Find the trace spines
 
// Coordinates of the node which indicate the spine to trace
 double x_trace; 
 double y_trace; 
 double z_trace;

 // All the nodes to trace are in the last slice
 y_trace = length_can;

//Find the spine of the top (horizontal spine)
 x_trace = 0.0;
 z_trace = height; 


 //Loop over all nodes
 for(unsigned inod = 0;inod<n_node;inod++)
  {
   SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
   
   double distance =sqrt( 
    (spnd_pt->x(0) - x_trace)* (spnd_pt->x(0) - x_trace) + 
    (spnd_pt->x(1) - y_trace)* (spnd_pt->x(1) - y_trace) + 
    (spnd_pt->x(2) - z_trace)* (spnd_pt->x(2) - z_trace) );
   if(distance < 10E-8)
    {
     cout<<"Top horizontal Spine found"<<endl;   
     Trace_spine_r_h_top = spnd_pt->spine_pt();
    }
  }
 
 
//Find the spine of the top (diagonal spine)
 x_trace = alpha/2;
 z_trace = height; 

 for(unsigned inod = 0;inod<n_node;inod++)
  {
   SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
   
   double distance =sqrt( 
    (spnd_pt->x(0) - x_trace)* (spnd_pt->x(0) - x_trace) + 
    (spnd_pt->x(1) - y_trace)* (spnd_pt->x(1) - y_trace) + 
    (spnd_pt->x(2) - z_trace)* (spnd_pt->x(2) - z_trace) );
   if(distance < 10E-8)
    {
     cout<<"Top diagonal Spine found"<<endl;   
     Trace_spine_r_d_top = spnd_pt->spine_pt();
    }
  }
 
 
//Find the spine of the middle (horizontal spine)
 x_trace = alpha/2;
 z_trace = height/2; 
 
 for(unsigned inod = 0;inod<n_node;inod++)
  {
   SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
   
   double distance =sqrt( 
    (spnd_pt->x(0) - x_trace)* (spnd_pt->x(0) - x_trace) + 
    (spnd_pt->x(1) - y_trace)* (spnd_pt->x(1) - y_trace) + 
    (spnd_pt->x(2) - z_trace)* (spnd_pt->x(2) - z_trace) );
   if(distance < 10E-8)
    {
     cout<<"Middle horizontal Spine found"<<endl;   
     Trace_spine_r_h_middle = spnd_pt->spine_pt();
    }
  }
 
 
//Find the spine of the bottom (horizontal spine)
 x_trace = 0.0;
 z_trace = 0.0; 



 for(unsigned inod = 0;inod<n_node;inod++)
  {
   SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
   
   double distance =sqrt(  
    (spnd_pt->x(0) - x_trace)* (spnd_pt->x(0) - x_trace) + 
    (spnd_pt->x(1) - y_trace)* (spnd_pt->x(1) - y_trace) + 
    (spnd_pt->x(2) - z_trace)* (spnd_pt->x(2) - z_trace) );
   if(distance < 10E-8)
    {
     cout<<"Bottom horizontal Spine found"<<endl;   
     Trace_spine_r_h_bottom = spnd_pt->spine_pt();
    }
  }
 

//Find the spine of the bottom (diagonal spine)
 x_trace = alpha/2;
 z_trace = 0.0; 



 for(unsigned inod = 0;inod<n_node;inod++)
  {
    SpineNode* spnd_pt = dynamic_cast<SpineNode*>(mesh_pt()->node_pt(inod));
   
    double distance =sqrt(
     (spnd_pt->x(0) - x_trace)* (spnd_pt->x(0) - x_trace) + 
     (spnd_pt->x(1) - y_trace)* (spnd_pt->x(1) - y_trace) + 
     (spnd_pt->x(2) - z_trace)* (spnd_pt->x(2) - z_trace) );
    if(distance < 10E-8)
     {
      cout<<"Bottom diagonal Spine found"<<endl;   
      Trace_spine_r_d_bottom = spnd_pt->spine_pt();
     }
    
  }
}


//========================================================================
/// Destructor for RectangularDrivenCavity problem 
//========================================================================
template<class ELEMENT>
ThreeDimBethertonProblem<ELEMENT>::~ThreeDimBethertonProblem()
{ 

 // Timestepper gets killed in general problem destructor

 // Mesh gets killed in general problem destructor

}




//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ThreeDimBethertonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];
 
 double height = Global_Physical_Variables::Height/2 ;  
 double width =  Global_Physical_Variables::Alpha/2; 
 double diagonal = sqrt(height*height + width*width ) ;
 double scale = 1/height;  //scale for the radius 

 double outflow =  this->get_outflow();

 double capillary_number = Global_Physical_Variables::Ca;

 double wet_fraction = outflow/width/(height*2);  
//The inflow is calculated in half channel

 double  p_drop = Bubble_pressure_data_pt->value(0);

// Get finger width at end of the computational domain
 double lambda = get_lambda();
 
 Trace_file<<capillary_number;

 Trace_file<<" \t"<<Global_Physical_Variables::Bo;

 Trace_file<<" \t"<<Global_Physical_Variables::Alpha;

 Trace_file<<" \t"<<wet_fraction;

 Trace_file<<" \t"<<p_drop;;
 
 Trace_file<<" \t"<<(height- Trace_spine_r_h_top->height() )*scale ;
 Trace_file<<" \t"<<(diagonal- Trace_spine_r_d_top->height() )*scale ;
 Trace_file<<" \t"<<(width - Trace_spine_r_h_middle->height() )*scale ;
 Trace_file<<" \t"<<(diagonal - Trace_spine_r_d_bottom->height() )*scale ;
 Trace_file<<" \t"<<(height - Trace_spine_r_h_bottom->height() )*scale ;
 
 Trace_file<<" \t"<<lambda ;
 
 Trace_file<<endl;
 
 // Number of plot points
 unsigned npts;
 npts=5; 
 
 mesh_pt()->node_update();
  
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",
         doc_info.directory().c_str(),doc_info.number());
 some_file.open(filename);
 unsigned n_bulk = mesh_pt()->nbulk();
 for(unsigned e=0;e<n_bulk;e++)
  {
   Hijacked<ELEMENT>* el_pt = 
    dynamic_cast<Hijacked<ELEMENT>* >(mesh_pt()->bulk_element_pt(e));
   el_pt->output(some_file,npts);
  }
 some_file.close();
}



//=====================================================================
/// Driver for RectangularDrivenCavity test problem -- test drive
/// with two different types of element.
//=====================================================================
int main(int argc,char *argv[])
{
 // Doc info object
 DocInfo doc_info;

 // Set output directory
 char name_dir[20];
 sprintf(name_dir,"RESLT");
 
 doc_info.set_directory(name_dir);
 doc_info.number()=0;

 // Direction of gravity
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = 0.0;
 Global_Physical_Variables::G[2] = -1.0;
 
 //Construct the problem
 ThreeDimBethertonProblem<SpineElement<QCrouzeixRaviartElement<3> >  >
  problem;

 //Solve and document without gravity
 problem.newton_solve();
 problem.doc_solution(doc_info);
 
 //Crank up the Bond number
 for(unsigned i=0;i<2;i++) 
  {
   doc_info.number()++;
   using namespace Global_Physical_Variables;
   Bo += 0.1;
   ReInvFr = Bo/Ca;
   problem.newton_solve();
   problem.doc_solution(doc_info);
  }

}
