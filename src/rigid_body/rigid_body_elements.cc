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
//Non-inline functions for Rigid Body Elements
#include "rigid_body_elements.h"

namespace oomph
{

//=======================================================================
/// \short Constructor: Specify coordinates of a point inside the hole
/// and a vector of pointers to TriangleMeshPolyLines
/// that define the boundary segments of the polygon.
/// Each TriangleMeshPolyLine has its own boundary ID and can contain
/// multiple (straight-line) segments. The optional final argument
/// is a pointer to a Data object whose three values represent 
/// the two displacements of and the rotation angle about the polygon's 
/// centre of mass.
//=======================================================================
RigidBodyTriangleMeshInternalPolygon::RigidBodyTriangleMeshInternalPolygon(
 const Vector<double>& hole_center,
 const Vector<TriangleMeshPolyLine*>& 
 boundary_polyline_pt,
 TimeStepper* const &time_stepper_pt,
 Data* const &centre_displacement_data_pt) :
 TriangleMeshPolygon(boundary_polyline_pt),
 TriangleMeshInternalClosedCurve(hole_center),
 TriangleMeshInternalPolygon(hole_center,boundary_polyline_pt), 
 RigidBodyElement(time_stepper_pt,centre_displacement_data_pt)
{  
 // Original rotation angle is zero
 Initial_Phi=0.0;

 // Compute coordinates of centre of gravity etc
 Vector<double> r_left(2);
 Vector<double> r_right(2);
 Mass=0.0;
 Initial_centre_of_mass[0]=0.0;
 Initial_centre_of_mass[1]=0.0;
 double inertia_x=0.0;
 double inertia_y=0.0;

 // Loop over polylines
 unsigned nboundary=boundary_polyline_pt.size();
 for (unsigned i=0;i<nboundary;i++)
  {
   // Loop over the segments to get the vertex coordinates
   unsigned nseg=boundary_polyline_pt[i]->nsegment();
   for(unsigned j=0;j<nseg;j++)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[i]->vertex_coordinate(j);
     r_right=boundary_polyline_pt[i]->vertex_coordinate(j+1);
   
     // Mass (area)
     Mass+=0.5*(r_left[0]*r_right[1]-r_right[0]*r_left[1]);

     // Centroid
     Initial_centre_of_mass[0]+=(r_left[0]+r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
     Initial_centre_of_mass[1]+=(r_left[1]+r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
    }
   if (nboundary==1)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[0]->vertex_coordinate(nseg);
     r_right=boundary_polyline_pt[0]->vertex_coordinate(0);
   
     // Mass (area)
     Mass+=0.5*(r_left[0]*r_right[1]-r_right[0]*r_left[1]);

     // Centroid
     Initial_centre_of_mass[0]+=(r_left[0]+r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
     Initial_centre_of_mass[1]+=(r_left[1]+r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
    }
  }
   
 // Normalise
 Initial_centre_of_mass[0]/=(6.0*Mass);
 Initial_centre_of_mass[1]/=(6.0*Mass);
   
 // Another loop over polylines for moment of inertia
 for (unsigned i=0;i<nboundary;i++)
  {
   // Loop over the segments to get the vertex coordinates
   unsigned nseg=boundary_polyline_pt[i]->nsegment();
   for(unsigned j=0;j<nseg;j++)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[i]->vertex_coordinate(j);
     r_right=boundary_polyline_pt[i]->vertex_coordinate(j+1);
       
     // Get moment about centroid
     r_left[0]-=Initial_centre_of_mass[0];
     r_left[1]-=Initial_centre_of_mass[1];
     r_right[0]-=Initial_centre_of_mass[0];
     r_right[1]-=Initial_centre_of_mass[1];
       
     // Moment of inertia
     inertia_x+=1.0/12.0*(r_left[1]*r_left[1]+
                          r_left[1]*r_right[1]+
                          r_right[1]*r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
       
     inertia_y+=1.0/12.0*(r_left[0]*r_left[0]+
                          r_left[0]*r_right[0]+
                          r_right[0]*r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);       
    }
     
   if (nboundary==1)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[0]->vertex_coordinate(nseg);
     r_right=boundary_polyline_pt[0]->vertex_coordinate(0);
       
     // Get moment about centroid
     r_left[0]-=Initial_centre_of_mass[0];
     r_left[1]-=Initial_centre_of_mass[1];
     r_right[0]-=Initial_centre_of_mass[0];
     r_right[1]-=Initial_centre_of_mass[1];
       
     // Moment of inertia
     inertia_x+=1.0/12.0*(r_left[1]*r_left[1]+
                          r_left[1]*r_right[1]+
                          r_right[1]*r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
       
     inertia_y+=1.0/12.0*(r_left[0]*r_left[0]+
                          r_left[0]*r_right[0]+
                          r_right[0]*r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);    
    }
  }
   
 // Polar moment of inertia is sum of two orthogonal planar moments
 Moment_of_inertia=inertia_x+inertia_y;
   
//    // Tested for circular and elliptical cross section
//    cout << "Mass             : " << Mass << std::endl;
//    cout << "Moment of inertia: " << Moment_of_inertia << std::endl;
//    cout << "X_c              : " << Initial_centre_of_mass[0] << std::endl;
//    cout << "Y_c              : " << Initial_centre_of_mass[1] << std::endl;
//    pause("done");


 //Assign the intrinsic coordinate
 this->assign_zeta();

//  {
//   unsigned n_poly = this->npolyline();
//   for(unsigned p=0;p<n_poly;++p)
//    {
//     std::cout << "Polyline " << p << "\n";
//     std::cout << "-----------------------\n";
//     unsigned n_vertex = Zeta_vertex[p].size();
//     for(unsigned v=0;v<n_vertex;v++)
//      {
//       std::cout << v << " " << Zeta_vertex[p][v] << "\n";
//      }
//    }
// }

}
 


//===============================================================
/// \short Update the reference configuration by re-setting the original
/// position of the vertices to their current ones, re-set the 
/// original position of the centre of mass, and the displacements 
/// and rotations relative to it
//===============================================================
void RigidBodyTriangleMeshInternalPolygon::reset_reference_configuration()
{
 Vector<double> x_orig(2);
 Vector<double> r(2);
 
 // Loop over the polylines and update their vertex positions
 unsigned npoly=Boundary_polyline_pt.size();
 for (unsigned i=0;i<npoly;i++)
  {
   TriangleMeshPolyLine* poly_line_pt=Boundary_polyline_pt[i];
   unsigned nvertex=poly_line_pt->nvertex();
   for (unsigned j=0;j<nvertex;j++)
    {
     x_orig=poly_line_pt->vertex_coordinate(j);
     this->apply_rigid_body_motion(0,x_orig,r);
     poly_line_pt->vertex_coordinate(j)=r;
    }
  }

 // Update coordinates of hole
 Vector<double> orig_hole_coord(this->internal_point());
 this->apply_rigid_body_motion(0,orig_hole_coord,this->internal_point());

 // Update centre of gravity
 double x_displ=Centre_displacement_data_pt->value(0);
 double y_displ=Centre_displacement_data_pt->value(1);
 double phi_displ=Centre_displacement_data_pt->value(2);
 Initial_centre_of_mass[0]+=x_displ;
 Initial_centre_of_mass[1]+=y_displ;
 Initial_Phi+=phi_displ;

 // Reset displacement and rotation ("non-previous-value" 
 // history values stay)
 TimeStepper* timestepper_pt=Centre_displacement_data_pt->time_stepper_pt();
 unsigned nprev=timestepper_pt->nprev_values();
 for (unsigned t=0;t<nprev;t++)
  {
   Centre_displacement_data_pt->
    set_value(t,0,Centre_displacement_data_pt->value(t,0)-x_displ);

   Centre_displacement_data_pt->
    set_value(t,1,Centre_displacement_data_pt->value(t,1)-y_displ);

   Centre_displacement_data_pt->
    set_value(t,2,Centre_displacement_data_pt->value(t,2)-phi_displ);
  }
}

}
