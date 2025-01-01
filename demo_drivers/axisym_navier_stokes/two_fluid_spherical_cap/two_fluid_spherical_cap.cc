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
// Driver code for an axisymmetric interface hydrostatics problem.
// The system consists of equal volumes of two fluids with 
// the same physical properties, but a non-zero surface tension at 
// their interface, in a cylinder of height 2 and radius 0.5.
// The program solves for the interface position as the contact angle
// at the wall, alpha, decreases from pi/2. The resulting shapes should all be
// spherical caps and the pressure jump across the interface should be
// 2 cos(alpha)/0.5 = 4 cos(alpha).

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/two_layer_spine_mesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables
{
 /// Direction of the wall normal vector
 Vector<double> Wall_normal;

 /// Function that specifies the wall unit normal
 void wall_unit_normal_fct(const Vector<double> &x, 
                      Vector<double> &normal)
 {
  normal=Wall_normal;
 }
}


//============================================================================
/// A Problem class that solves the Navier--Stokes equations
/// in an axisymmetric geometry
/// N.B. we usually template problems by element-type and a timestepper in
/// any problems where different elements or timesteppers MAY be used
//============================================================================
template<class ELEMENT>
class CapProblem : public Problem
{
private:

 /// The Capillary number 
 double Ca;

 /// The volume of the fluid
 double Volume;

 /// The contact angle
 double Angle;

public:
 //Constructor:
 //Nr: Number of elements in the r (radial) direction
 //Nz1: Number of elements in the z (axial) direction in the lower fluid
 //Nz2: Number of elements in the z (axial) direction in the upper fluid
 CapProblem(const unsigned &Nr, const unsigned &Nz1, 
            const unsigned &Nz2, const double & width);

 /// Set boundary conditions on the wall. 
 void set_boundary_conditions();

 /// Solve the problem
 void solve_system();

 /// Create the volume constraint elements
 void create_volume_constraint_elements();

 /// Pointer to the mesh
 TwoLayerSpineMesh<SpineElement<ELEMENT> >* Bulk_mesh_pt;

 /// Pointer to the volume constraint mesh
 Mesh* Volume_constraint_mesh_pt;

 /// Pointer to the surface mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to the point mesh
 Mesh* Point_mesh_pt;

 /// Data that is traded for the volume constraint
 Data* Traded_pressure_data_pt;
 
 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {Bulk_mesh_pt->node_update();}

};

/// Constructor 
template<class ELEMENT>
CapProblem<ELEMENT>::CapProblem
(const unsigned &Nr, const unsigned &Nz1, 
 const unsigned &Nz2, const double &width) :
 Ca(1.0),  //Initialise value of Ca to one
 Volume(0.125),  //Initialise the value of the volume to 0.125
 Angle(1.57)     //Initialise the value of the contact angle
{
 //Set the wall normal (in the r-direction)
 Global_Physical_Variables::Wall_normal.resize(2); 
 Global_Physical_Variables::Wall_normal[0] = 1.0; 
 Global_Physical_Variables::Wall_normal[1] = 0.0;

 //Construct our mesh
 Bulk_mesh_pt = 
  new TwoLayerSpineMesh<SpineElement<ELEMENT> >(Nr,Nz1,Nz2,width,1.0,1.0);

 Surface_mesh_pt = new Mesh;

 //Loop over the horizontal elements
 unsigned n_interface_lower = Bulk_mesh_pt->ninterface_lower();
 for(unsigned i=0;i<n_interface_lower;i++)
  {
   //Construct a new 1D line element adjacent to the interface
   FiniteElement *interface_element_pt =  new 
    SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >(
     Bulk_mesh_pt->interface_lower_boundary_element_pt(i),
     Bulk_mesh_pt->interface_lower_face_index_at_boundary(i));
   
    this->Surface_mesh_pt->add_element_pt(interface_element_pt);
   }
 
//Create the point mesh
 Point_mesh_pt = new Mesh;
 {
  //Make the point (contact) element from the last surface element
  FiniteElement* point_element_pt = 
   dynamic_cast<SpineAxisymmetricFluidInterfaceElement<
    SpineElement<ELEMENT> >*>(Surface_mesh_pt->element_pt(n_interface_lower-1))
   ->make_bounding_element(1);
  
  //Add it to the mesh
  this->Point_mesh_pt->add_element_pt(point_element_pt);
 }
 
 
 //Set the linear solver (frontal solver) (function defined in Problem class)
 //linear_solver_pt() = new HSL_MA42;

 //We are going to hijack one of the pressure values and make a copy of it
 //in global data
 flush_global_data();


 //Hijack one of the pressure values in the upper fluid
 add_global_data(dynamic_cast<ELEMENT*>(Bulk_mesh_pt->upper_layer_element_pt(0))
  ->hijack_internal_value(0,0));

 Traded_pressure_data_pt = this->global_data_pt(0);

 //Loop over the elements on the free surface
 unsigned Ninterface = Surface_mesh_pt->nelement();
 for(unsigned e=0;e<Ninterface;e++)
  {
   //Cast to a 1D element
   SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> > 
    *el_pt = 
    dynamic_cast<
    SpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >*>
    (Surface_mesh_pt->element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;
  }

 //Set the boundary conditions
 //Pin all velocity components on the bottom, wall and top
 for(unsigned b=0;b<4;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     for(unsigned i=0;i<3;i++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(i);
      }
    }
  }

 //Axis boundary, pin u (component 0) and v (component 2), but leave w free
 for(unsigned b=4;b<6;b++)
  {
   //Find the number of nodes on the boundary
   unsigned Nboundary_node = Bulk_mesh_pt->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<Nboundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
    }
  }
 
 //Set the contact angle on the RHS
 dynamic_cast<FluidInterfaceBoundingElement*>
  (Point_mesh_pt->element_pt(0))->set_contact_angle(&Angle);
 
 dynamic_cast<FluidInterfaceBoundingElement*>
  (Point_mesh_pt->element_pt(0))->ca_pt() = &Ca;
 
 //Set the wall normal on the rhs
 dynamic_cast<FluidInterfaceBoundingElement*>
  (Point_mesh_pt->element_pt(0))->wall_unit_normal_fct_pt() = 
  &Global_Physical_Variables::wall_unit_normal_fct;

 //Pin a single pressure value
 dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0))->fix_pressure(0,0.0);

 this->create_volume_constraint_elements();

 this->add_sub_mesh(Bulk_mesh_pt);
 this->add_sub_mesh(Surface_mesh_pt);
 this->add_sub_mesh(Point_mesh_pt);
 this->add_sub_mesh(Volume_constraint_mesh_pt);
 
 this->build_global_mesh();

 //Setup all the equation numbering and look-up schemes 
 //(function defined in Problem class)
 cout << assign_eqn_numbers() << std::endl; 
}


//======================================================================
/// Create the volume constraint elements
//========================================================================
template<class ELEMENT>
void CapProblem<ELEMENT>::create_volume_constraint_elements()
{
 //The single volume constraint element
 Volume_constraint_mesh_pt = new Mesh;
 VolumeConstraintElement* vol_constraint_element = 
  new VolumeConstraintElement(&Volume,Traded_pressure_data_pt,0);
 Volume_constraint_mesh_pt->add_element_pt(vol_constraint_element);

 //Loop over lower boundaries (or just 1, why?)
 unsigned lower_boundaries[4]={0,1,5};
 for(unsigned i=0;i<3;i++)
  {
   //Read out the actual boundaries
   unsigned b = lower_boundaries[i];

   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));

     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Create new element
     SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }

 //Loop over one side of the interface
  {
   // How many bulk fluid elements are adjacent to boundary b?
   unsigned n_element = Bulk_mesh_pt->ninterface_lower();
   
   // Loop over the bulk fluid elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->interface_lower_boundary_element_pt(e));

     //Find the index of the face of element e along boundary b
     int face_index = Bulk_mesh_pt->interface_lower_face_index_at_boundary(e);
     
     // Create new element
     SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>* el_pt =
      new SpineAxisymmetricVolumeConstraintBoundingElement<ELEMENT>(
       bulk_elem_pt,face_index);   

     //Set the "master" volume control element
     el_pt->set_volume_constraint_element(vol_constraint_element);

     // Add it to the mesh
     Volume_constraint_mesh_pt->add_element_pt(el_pt);     
    }
  }

}


//Set the boundary conditions
template<class ELEMENT>
void CapProblem<ELEMENT>::set_boundary_conditions()
{
 //NOTE: The default value of all parameters is zero, so we don't need 
 //to set anything here
}


template<class ELEMENT>
void CapProblem<ELEMENT>::solve_system()
{
 //Define a string that we can set to be the name of the output file
 char filename[100];
 
 //Set the boundary conditions (only need to do this once)
 //If the boundary conditions depend up on time, or Re, we need to
 //reset them every time something changes. This is most easily
 //achieved using the actions_before_newton_solve() {} function.
 set_boundary_conditions();

 //Solve the problem (function defined in Problem class)
 //The default tolerance is that the maximum residual must be less that 1.0e-8
 steady_newton_solve();

 //Output first solution
 ofstream file("step0_spine.dat");
 Bulk_mesh_pt->output(file,5);
 Surface_mesh_pt->output(file,5);
 file.close();

 //Gradually increase the contact angle
 for(unsigned i=1;i<6;i++)
  {
   //Decrease the contact angle
   Angle -= 0.1;

   //Solve the problem
   steady_newton_solve();

   //Create the filename, including the array index
   sprintf(filename,"step%i_spine.dat",i);
   //Actually, write the data
   file.open(filename);
   Bulk_mesh_pt->output(file,5);
   Surface_mesh_pt->output(file,5);
   file.close();
  }

} 

//Main driver loop
int main()
{
 //Construct the problem, you can use either element types
 CapProblem<Hijacked<AxisymmetricQCrouzeixRaviartElement>  > 
 problem(4,4,4,0.5);
 //Solve the problem :)
 problem.solve_system();
}








