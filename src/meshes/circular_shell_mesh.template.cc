//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
#ifndef OOMPH_CIRCULAR_SHELL_MESH_TEMPLATE_CC
#define OOMPH_CIRCULAR_SHELL_MESH_TEMPLATE_CC


#include "circular_shell_mesh.template.h"
#include "rectangular_quadmesh.template.cc"


namespace oomph
{



//=======================================================================
/// Mesh constructor
/// Argument list:
/// nx  : number of elements in the axial direction
/// ny : number of elements in the azimuthal direction
/// lx  : length in the axial direction
/// ly  : length in theta direction
//=======================================================================
template <class ELEMENT>
CircularCylindricalShellMesh<ELEMENT>::CircularCylindricalShellMesh(
 const unsigned &nx,
 const unsigned &ny,
 const double &lx,
 const double &ly,
 TimeStepper* time_stepper_pt) :
 RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Now in this case it is the Lagrangian coordinates that we want to set,
 //so we have to loop over all nodes and set them to the Eulerian
 //coordinates that are set by the generic mesh generator
 for(unsigned i=0;i<n_node;i++)
  {
   node_pt(i)->xi(0) = node_pt(i)->x(0);
   node_pt(i)->xi(1) = node_pt(i)->x(1);
  }


 //Assign gradients, etc for the Lagrangian coordinates of
 //Hermite-type elements

 //Read out number of position dofs
 unsigned n_position_type = finite_element_pt(0)->nnodal_position_type();

 //If this is greater than 1 set the slopes, which are the distances between
 //nodes. If the spacing were non-uniform, this part would be more difficult
 if(n_position_type > 1)
  {
   double xstep = (this->Xmax - this->Xmin)/((this->Np-1)*this->Nx);
   double ystep = (this->Ymax - this->Ymin)/((this->Np-1)*this->Ny);
   for(unsigned n=0;n<n_node;n++)
    {
     //The factor 0.5 is because our reference element has length 2.0
     node_pt(n)->xi_gen(1,0) = 0.5*xstep;
     node_pt(n)->xi_gen(2,1) = 0.5*ystep;
    }
  }
}


//=======================================================================
/// Set the undeformed coordinates of the nodes
//=======================================================================
template <class ELEMENT>
void CircularCylindricalShellMesh<ELEMENT>::assign_undeformed_positions(
 GeomObject* const &undeformed_midplane_pt)
{
 //Find out how many nodes there are
 unsigned n_node = nnode();

 //Loop over all the nodes
 for(unsigned n=0;n<n_node;n++)
  {
   //Get the Lagrangian coordinates
   Vector<double> xi(2);
   xi[0] = node_pt(n)->xi(0);
   xi[1] = node_pt(n)->xi(1);

   //Assign memory for values of derivatives, etc
   Vector<double> R(3);
   DenseMatrix<double> a(2,3);
   RankThreeTensor<double>  dadxi(2,2,3);

   //Get the geometrical information from the geometric object
   undeformed_midplane_pt->d2position(xi,R,a,dadxi);

   //Loop over coordinate directions
   for(unsigned i=0;i<3;i++)
    {
     //Set the position
     node_pt(n)->x_gen(0,i) = R[i];

     //Set the derivative wrt Lagrangian coordinates
     //Note that we need to scale by the length of each element here!!
     node_pt(n)->x_gen(1,i) = 0.5*a(0,i)*((this->Xmax - this->Xmin)/this->Nx);
     node_pt(n)->x_gen(2,i) = 0.5*a(1,i)*((this->Ymax - this->Ymin)/this->Ny);

     //Set the mixed derivative
     //(symmetric so doesn't matter which one we use)
     node_pt(n)->x_gen(3,i) = 0.25*dadxi(0,1,i);
    }
  }
}


}
#endif
