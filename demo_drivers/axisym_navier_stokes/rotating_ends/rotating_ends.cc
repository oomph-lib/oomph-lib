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
//A simple axisymmetric Navier--Stokes problem:
//
//The lid and base of a cylinder of aspect ratio (height/width) 2
//are counter-rotated. The Reynolds number is increased from zero to 50.
//in two steps.

//OOMPH-LIB include files
#include "generic.h"
#include "axisym_navier_stokes.h"
#include "meshes/rectangular_quadmesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;


//======================================================================
/// A Class for the mesh to be used in the rotating Disc equations
/// The mesh is a modification of the standard RectangularQuadMesh. We
/// define x to be the radial (r) coordinate and y to be the axial (z)
/// coordinate in the axisymmetric coordinate system.
/// The mesh consists of upper and lower regions and left
/// and right regions, with uniform element spacing in each region. 
/// The idea is that one might want more closely-spaced elements 
/// near the downstream end of the domain and different numbers of elements
/// in the upper and lower regions, to prevent there being any
/// inherent symmetry in the mesh
//======================================================================
template <class ELEMENT>
class CylinderMesh : public RectangularQuadMesh<ELEMENT>
{
private:

 //Data to hold the number of elements in each section of the mesh
 unsigned Nx1, Nx2, Ny1, Ny2;

 //Data that holds the fractions along the disc in which to introduce
 //the second zones
 double Xfraction, Yfraction;

 public:
 /// Constructor, which "builds" the mesh. The arguments are the number
 /// of elements in each zone, the radius of the disc
 CylinderMesh(const unsigned &nx1, const unsigned &nx2,
              const unsigned &ny1, const unsigned &ny2,
              const double &radius) :
  //Call a protected constructor of the rectangular quad mesh
  //that sets internal variables, but does not actually build the
  //mesh. The arguments are the total number of elements in the x
  //and y directions, xmin, xman, ymin, ymax, periodic_flag and
  //the build flag
  RectangularQuadMesh<ELEMENT>(nx1+nx2,ny1+ny2,0.0,radius,-1.0,1.0,
                               false,false)
  {
   //Set private mesh variables, number of elements in each section
   Nx1 = nx1; Nx2 = nx2; Ny1 = ny1; Ny2 = ny2;

   //Set the fraction, at which the left and right regions are split
   Xfraction=0.8;

   //Set the fraction at which the upper and lower regions are split
   Yfraction=0.5;

   //Call the generic mesh construction routine
   //(Defined in RectangularQuadMesh<ELEMENT>)
   this->build_mesh();
        
   //Now reorder the elements in the best form for the frontal solver
   //This arranges in elements in vertical slices starting from x=0 and moving
   //along the channel to x=radius
   this->element_reorder();
  }

 /// Define the spacing of the nodes in the r (radial) direction
 /// In this mesh, there are two different regions with uniform spacing 
 /// in each. 
 /// Note that this can easily be changed to generate non-uniform meshes
 /// if desired
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode)
  {
   //Find the minimum and maximum values of the X (radial) coordinate
   double x_min = this->Xmin;
   double x_max = this->Xmax;
   //Find the number of nodes in the element in the radial direction
   unsigned np = this->Np;
   //Set up some spacing parameters
   //Region one starts at Xmin 
   //Region two starts at Xmin + Xfraction(Xmax-Xmin)
   double x1init = x_min, x2init = x_min + Xfraction*(x_max-x_min);
   //Calculate the spacing between the nodes in each region
   //Assuming uniform spacing
   //Region one has a length Xfraction*(Xmax-Xmin)
   double x1step = Xfraction*(x_max-x_min)/((np-1)*Nx1);
   //Region two has a length (1.0-Xfraction)*(Xmax-Xmin)
   double x2step = (1.0-Xfraction)*(x_max-x_min)/((np-1)*Nx2);
   
   //Now set up the particular spacing 
   //(it's different in the two different regions)
   if(xelement < Nx1) {return (x1init + x1step*((np-1)*xelement + xnode));}
   else {return (x2init + x2step*((np-1)*(xelement-Nx1) + xnode));}
  }

 /// Defines the spacing of the nodes in the z 
 /// (axial) direction.In this mesh there are two different regions
 /// with uniform spacing in each
 /// Note that this can easily be changed to generate non-uniform meshes
 /// if desired. This may be desirable if sharp boundary, or shear layers,
 /// develop. 
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode)
  {
   //Find the minimum and maximum values of the Y (axial) coordinate
   double y_min = this->Ymin;
   double y_max = this->Ymax;
   //Read out the number of nodes in the y direction
   unsigned np = this->Np;

   //Set up some spacing parameters
   //The lower region starts at Ymin
   //The upper region starts at Ymin + Yfraction*(Ymax-Ymin)
   double y1init = y_min, y2init = y_min + Yfraction*(y_max-y_min);
   //Calculate the space between each node in each region,
   //Assumming uniform spacing
   //Region one has a length Yfraction(Ymax-Ymin)
   double y1step = Yfraction*(y_max-y_min)/((np-1)*Ny1);
   //Region two has a length (1.0-Yfraction)*(Ymax-Ymin)
   double y2step = (1.0-Yfraction)*(y_max-y_min)/((np-1)*Ny2);
 
   //Now return the actual node position, it's different in the two
   //regions, of course
   if(yelement < Ny1) {return (y1init + y1step*((np-1)*yelement + ynode));}
   else {return (y2init + y2step*((np-1)*(yelement-Ny1) + ynode));}
  }

};

//==========================================================================
/// Solve the Axisymmetric Navier Stokes equations for a cylinder with 
/// counter-rotating ends and aspect ratio 2.
//==========================================================================
template<class ELEMENT>
class RotatingProblem : public Problem
{
private:

 /// The Reynolds number will be private member data
 double Re;

public:
 /// Constructor:
 /// Nr: Number of elements in the r (radial) direction
 /// Nz: Number of elements in the z (axial) direction
 RotatingProblem(const unsigned &Nr1, const unsigned &Nr2,
                 const unsigned &Nz1, const unsigned &Nz2);

 /// Set boundary conditions on the walls
 void set_boundary_conditions();

 /// Function that is used to run the parameter study
 void solve_system();

 /// Return a pointer to the specific mesh used
 CylinderMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<CylinderMesh<ELEMENT>*>(Problem::mesh_pt());}
 
 /// No actions to be taken after each solve step 
 void actions_after_newton_solve() {}

 /// No actions to be taken before each solve step
 void actions_before_newton_solve() {}
};

//============================================================================
/// Constructor:
/// Nr1: Number of elements in the left region in r (radial) direction
/// Nr2: Number of elements in the right region in r (radial) direction
/// Nz1: Number of elements in the lower region in z (axial) direction
/// Nz2: Number of elements in the upper region in z (axial) direction
//============================================================================
template<class ELEMENT>
RotatingProblem<ELEMENT>::RotatingProblem
(const unsigned &Nr1, const unsigned &Nr2,
 const unsigned &Nz1, const unsigned &Nz2) :
 Re(0.0) //Initialise value of Re to zero
{
 //Now create the mesh, a generic square mesh, the boundaries of the mesh
 //are labelled 0 to 3 starting from the "bottom" and proceeding in an
 //anti-clockwise direction. The parameters of the constructor are passed 
 //directly to the mesh constructor, and the width of the channel is 
 //always 1.0. The mesh is defined in stdmesh.h
 Problem::mesh_pt() = new CylinderMesh<ELEMENT>(Nr1,Nr2,Nz1,Nz2,1);
  
 //Loop over all the (fluid) elements 
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   //Cast to the particular element type, this is necessary because
   //the base elements don't have the member functions that we're about
   //to call!
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //There is no need for ALE
   el_pt->disable_ALE();

   //Set the Reynolds number for each element 
   //(yes we could have different Reynolds number in each element!!)
   el_pt->re_pt() = &Re;
  }

 //Let this problem be conventional form by setting gamma to zero
 ELEMENT::Gamma[0] = 0.0; //r-momentum
 ELEMENT::Gamma[1] = 0.0; //z-momentum

 //Set the boundary conditions
 //Loop over all four boundaries
 for(unsigned i=0;i<4;i++)
  {
   //Find the number of nodes on the boundary
   unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<Nboundary_node;n++)
    {
     //Pin the u (radial) velocity on all cylinder boundaries
     mesh_pt()->boundary_node_pt(i,n)->pin(0);
     //Pin the v velocity on all cylinder boundaries
     mesh_pt()->boundary_node_pt(i,n)->pin(2);
     //If not on the axis, pin the w-velocity
     if(i!=3) {mesh_pt()->boundary_node_pt(i,n)->pin(1);}
    }
  }
 
 //Pin a single pressure value 
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
 
 //Setup all the equation numbering and look-up schemes 
 cout << assign_eqn_numbers() << std::endl; 
}

//========================================================================
/// Set the boundary conditions
//========================================================================
template<class ELEMENT>
void RotatingProblem<ELEMENT>::set_boundary_conditions()
{
 //NOTE: The default value of all parameters is zero, so we need only 
 //set the values that are non-zero on the boundaries, i.e. the swirl
 
 //Top boundary
 {
  unsigned i=2;
  //Find the number of nodes on the boundary
  unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
  //Loop over the nodes on the boundary
  for(unsigned n=0;n<Nboundary_node;n++)
   {
    //Get the radial values
    double r = mesh_pt()->boundary_node_pt(i,n)->x(0);
    //Set the value of the v-velocity (negative rotation)
    mesh_pt()->boundary_node_pt(i,n)->set_value(2,-1.0*r);
   }
 }
 
 //Bottom boundary
 {
  unsigned i=0;
  //Find the number of nodes on the boundary
  unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
  //Loop over the nodes on the boundary
  for(unsigned n=0;n<Nboundary_node;n++)
   {
    //Get the radial values
    double r = mesh_pt()->boundary_node_pt(i,n)->x(0);
    //Set the value of the v-velocity (positive rotation)
    mesh_pt()->boundary_node_pt(i,n)->set_value(2,1.0*r);
   }
 }
}

//==========================================================================
/// Solve the system for a number of different values of the Reynolds number
//==========================================================================
template<class ELEMENT>
void RotatingProblem<ELEMENT>::solve_system()
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
 ofstream file("Re0.dat");
 mesh_pt()->output(file,5);
 file.close();
 
 //Increase the Reynolds number in steps of 25
 for(unsigned i=1;i<3;i++)
  {
   Re += 25.0;
   //Solve the problem
   steady_newton_solve();
   
   //Output data at each step
   //Create the filename, including the array index
   sprintf(filename,"Re%g.dat",Re);
   //Actually, write the data
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();
  }
}

//Main driver loop
int main()
{
 //Construct and solve the problem
 RotatingProblem<AxisymmetricQTaylorHoodElement> problem(4,4,4,4);
 problem.solve_system();
}








