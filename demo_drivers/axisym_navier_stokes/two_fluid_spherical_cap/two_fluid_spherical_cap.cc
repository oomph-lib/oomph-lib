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

// The mesh
#include "meshes/two_layer_spine_mesh.h"

//This is required so that the classes are instantiated for the Hijacked elements
#include "fix_vol_axi_int_el.h"

//Use the std namespace
using namespace std;

using namespace oomph;

//============================================================================
/// Specific mesh for the spherical cap problem: Essentially a 
/// standard two-layer mesh with an additional element that applies
/// the volume constraint.
//============================================================================
template <class ELEMENT>
class SphericalCapMesh : public TwoLayerSpineMesh<
 SpineElement<ELEMENT>, 
 FixVolSpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> > >
{

private:

 
 /// \short Pointer to the point element that is used to enforce
 /// conservation of mass
 FiniteElement* Point_element_pt;

public:

 /// Constructor: Pass number of elements in radial direction, number
 /// of elements in the two layers and radius.
 SphericalCapMesh(const unsigned &nr, const unsigned &nh1,
                  const unsigned &nh2, const double & radius);
 
 /// Return pointer to the volumetric constraint element
 FiniteElement* &point_element_pt() {return Point_element_pt;}
 

};


//======================================================================
/// Constructor: Pass number of elements in radial direction, number
/// of elements in the two layers and radius.
//======================================================================
template<class ELEMENT>
SphericalCapMesh<ELEMENT>::SphericalCapMesh(const unsigned &nr,
                                            const unsigned &nh1,
                                            const unsigned &nh2, 
                                            const double & radius) :
 TwoLayerSpineMesh<SpineElement<ELEMENT>, 
                   FixVolSpineAxisymmetricFluidInterfaceElement<
 SpineElement<ELEMENT> > >(nr,nh1,nh2,radius,1.0,1.0)
                                               
{
 //Reorder the elements
this-> element_reorder();

 // Last interface element:
 FixVolSpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >* el_pt=
 dynamic_cast<FixVolSpineAxisymmetricFluidInterfaceElement<
  SpineElement<ELEMENT> >*>(this->Interface_element_pt[this->Nx-1]);

 //Create the volume constraint edge element
 Point_element_pt = el_pt->make_edge_element(1);

 //Push it back onto the stack
 this->Element_pt.push_back(Point_element_pt);
}



//============================================================================
///A Problem class that solves the Navier--Stokes equations
///in an axisymmetric geometry
///N.B. we usually template problems by element-type and a timestepper in
///any problems where different elements or timesteppers MAY be used
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

 /// The normal to the wall of contact
 Vector<double> Wall_normal;

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

 /// Overload access function for the mesh
 SphericalCapMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<SphericalCapMesh<ELEMENT>*>(Problem::mesh_pt());}
 
 /// Update the spine mesh after every Newton step
 void actions_before_newton_convergence_check() {mesh_pt()->node_update();}

 /// No other actions to be taken before or after each solve step
 void actions_after_newton_solve() {}
 void actions_before_newton_solve() {}
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
 Wall_normal.resize(2); Wall_normal[0] = 1.0; Wall_normal[1] = 0.0;

 //Construct our mesh
 Problem::mesh_pt() = 
  new SphericalCapMesh<ELEMENT>(Nr,Nz1,Nz2,width);

 //Set the linear solver (frontal solver) (function defined in Problem class)
 //linear_solver_pt() = new HSL_MA42;

 //We are going to hijack one of the pressure values and make a copy of it
 //in global data
 flush_global_data();


 //Hijack one of the pressure values in the upper fluid
 add_global_data(dynamic_cast<ELEMENT*>(mesh_pt()->upper_layer_element_pt(0))
  ->hijack_internal_value(0,0));
 
 //Loop over the elements on the free surface
 unsigned Ninterface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<Ninterface;e++)
  {
   //Cast to a 1D element
   FixVolSpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> > 
    *el_pt = 
    dynamic_cast<
    FixVolSpineAxisymmetricFluidInterfaceElement<SpineElement<ELEMENT> >*>
    (mesh_pt()->interface_element_pt(e));
   //Set the Capillary number
   el_pt->ca_pt() = &Ca;

   //Add the pressure dof traded for the volume constraint 
   //to the element as global data
   el_pt->set_traded_pressure_data(global_data_pt(0));
  }

 {
 //Finally do the point elementm, adding in the pressure
 SpineAxisymmetricVolumeConstraintPointElement<SpineElement<ELEMENT> >* 
  el_pt = dynamic_cast<
  SpineAxisymmetricVolumeConstraintPointElement<SpineElement<ELEMENT> >*>
  (mesh_pt()->point_element_pt());

   el_pt->volume_pt() = &Volume;
   el_pt->set_traded_pressure_data(global_data_pt(0));
 }    

 //Set the boundary conditions
 //Pin all velocity components on the bottom, wall and top
 for(unsigned b=0;b<3;b++)
  {
   //Find the number of nodes on the boundary
   unsigned n_boundary_node = mesh_pt()->nboundary_node(b);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     for(unsigned i=0;i<3;i++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(i);
      }
    }
  }

 //Axis boundary, pin u (component 0) and v (component 2), but leave w free
 {
  //Find the number of nodes on the boundary
  unsigned Nboundary_node = mesh_pt()->nboundary_node(3);
  //Loop over the nodes on the boundary
  for(unsigned n=0;n<Nboundary_node;n++)
   {
    mesh_pt()->boundary_node_pt(3,n)->pin(0);
    mesh_pt()->boundary_node_pt(3,n)->pin(2);
   }
 }
 
 //Set the contact angle on the RHS
 dynamic_cast<FluidInterfaceEdgeElement*>
  (mesh_pt()->element_pt(mesh_pt()->nelement()-1))->set_contact_angle(&Angle);
 //Set the wall normal on the rhs
 dynamic_cast<FluidInterfaceEdgeElement*>
  (mesh_pt()->element_pt(mesh_pt()->nelement()-1))->wall_normal_pt() = 
  &Wall_normal;;
 

 //Pin a single pressure value
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);

 //Setup all the equation numbering and look-up schemes 
 //(function defined in Problem class)
 cout << assign_eqn_numbers() << std::endl; 
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
 mesh_pt()->output(file,5);
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
   mesh_pt()->output(file,5);
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








