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
//Demo driver for adaptive bifurcation detection
//based on the flow past a cylinder within a channel
//There is a Hopf bifurcation at Re = 92.
//The study is based on Cliffe & Tavener (2004) JFM 501, 125-133.
#include <iostream>
#include <fstream>
#include <cstdio>

// Generic oomph-lib stuff
#include "generic.h"

// Navier Stokes
#include "navier_stokes.h"
#include "linearised_navier_stokes.h"

#include "assert.h"

using namespace oomph;

using namespace QuadTreeNames;


//===============================================
/// Global parameters
//===============================================
namespace Global_Parameters
{

 // Reynolds number
 double Re=75.0;

 // Blockage ratio
 double B=0.7;

 // Rotation ratio
 double Alpha=0.0;

 // Control Flag that will read in the eigenfunction from disk
 // This is required because we are not allowed to redistribute
 // ARPACK, so cannot assume that it has been installed.
 bool Read_in_eigenfunction_from_disk = true;

 //Set the normalisation for the eigenfunction
 std::complex<double> Eigenfunction_normalisation(1.0,0.0);
}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




//===============================================
/// My own Ellipse class
//===============================================
class GeneralEllipse : public GeomObject
{

public:
 
 ///Constructor
 GeneralEllipse(const double &centre_x, const double &centre_y,
                const double &a, const double &b)
  : GeomObject(1,2), Centre_x(centre_x), Centre_y(centre_y), A(a), B(b)
  {}

 ///Destructor (empty)
 ~GeneralEllipse(){}

 ///Return the position
 void position(const Vector<double> &xi, Vector<double> &r) const
  {
   r[0] = Centre_x + A*cos(xi[0]);
   r[1] = Centre_y + B*sin(xi[0]);
  }


 /// Access to  x-coordinate of centre
 double& centre_x(){return Centre_x;}

 /// Access to  y-coordinate of centre
 double& centre_y(){return Centre_y;}

 /// Access to x-half axis
 double& a(){return A;}

 /// Access to y-half axis
 double& b(){return B;}


private:

 /// x-coordinate of centre
 double Centre_x;

 /// y-coordinate of centre
 double Centre_y;

 /// x-half axis
 double A;

 /// y-half axis
 double B;

};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////



//===========================================================
/// Rectangular domain with circular whole
//===========================================================
class RectangleWithHoleDomain : public Domain
{

public:

 double centre_x() {return Centre_x;}
 double centre_y() {return Centre_y;}

 /// Constructor, Pass pointer to geometric object that 
 /// represents the cylinder, the length and height of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction.
 RectangleWithHoleDomain(GeomObject* cylinder_pt,
                         const double &length, const double &height) :
  Cylinder_pt(cylinder_pt)
  {

   Centre_x = static_cast<GeneralEllipse*>(cylinder_pt)->centre_x();
   Centre_y = static_cast<GeneralEllipse*>(cylinder_pt)->centre_y();

   double lower_wall = -0.5*height;
   double upper_wall = 0.5*height;

   //Vertices of rectangle
   Lower_left.resize(2);
   Lower_left[0] = 0.0;
   Lower_left[1] = lower_wall;
   
   Upper_left.resize(2);
   Upper_left[0] = 0.0;
   Upper_left[1] = upper_wall;

   Lower_right.resize(2);
   Lower_right[0] = length;
   Lower_right[1] = lower_wall;
   
   Upper_right.resize(2);
   Upper_right[0] = length;
   Upper_right[1] = upper_wall;


   double left_cylinder_split = 3.5;
   double right_cylinder_split = 5.5;

   // Coordinates of points where the "radial" lines from central
   // cylinder meet the upper and lower boundaries
   Lower_mid_left.resize(2);
   Lower_mid_left[0] = left_cylinder_split;
   Lower_mid_left[1] = lower_wall;

   Upper_mid_left.resize(2);
   Upper_mid_left[0] = left_cylinder_split;
   Upper_mid_left[1] = upper_wall;

   Lower_mid_right.resize(2);
   Lower_mid_right[0] = right_cylinder_split;
   Lower_mid_right[1] = lower_wall;

   Upper_mid_right.resize(2);
   Upper_mid_right[0] = right_cylinder_split;
   Upper_mid_right[1] = upper_wall;
   

   // Coordinates of internal points for creation of macro elements
   
   double help = (length - right_cylinder_split)/4.;
   
   Upper_1.resize(2);
   Upper_1[0] = right_cylinder_split + help;
   Upper_1[1] = upper_wall;

   Upper_2.resize(2);
   Upper_2[0] = right_cylinder_split + 2.*help;
   Upper_2[1] = upper_wall;

   Upper_3.resize(2);
   Upper_3[0] = right_cylinder_split + 3.*help;
   Upper_3[1] = upper_wall;
   
   Lower_1.resize(2);
   Lower_1[0] = right_cylinder_split + help;
   Lower_1[1] = lower_wall;

   Lower_2.resize(2);
   Lower_2[0] = right_cylinder_split + 2.*help;
   Lower_2[1] = lower_wall;

   Lower_3.resize(2);
   Lower_3[0] = right_cylinder_split + 3.*help;
   Lower_3[1] = lower_wall;
   
   //There are nine macro elements
   Macro_element_pt.resize(9); 

   // Build the 2D macro elements
   for (unsigned i=0;i<9;i++)
    {Macro_element_pt[i]= new QMacroElement<2>(this,i);}
  }



 /// Destructor: Empty; cleanup done in base class
 ~RectangleWithHoleDomain(){}


 ///  Helper function to interpolate linearly between the
 /// "right" and "left" points; \f$ s \in [-1,1] \f$
 void linear_interpolate(Vector<double> left, Vector<double> right,
                         const double &s, Vector<double> &f)
  {
   for(unsigned i=0;i<2;i++)
    {
     f[i] = left[i] + (right[i] - left[i])*0.5*(s+1.0);
    }
  }

   

 ///  Parametrisation of macro element boundaries: f(s) is the position
 /// vector to macro-element m's boundary in the specified direction [N/S/E/W]
 /// at the specfied discrete time level (time=0: present; time>0: previous)
 void macro_element_boundary(const unsigned &time,
                             const unsigned &m,
                             const unsigned &direction,
                             const Vector<double> &s,
                             Vector<double>& f)
 {

#ifdef WARN_ABOUT_SUBTLY_CHANGED_OOMPH_INTERFACES
   // Warn about time argument being moved to the front
   OomphLibWarning(
    "Order of function arguments has changed between versions 0.8 and 0.85",
    "RectangleWithHoleDomain::macro_element_boundary(...)",
    OOMPH_EXCEPTION_LOCATION);
#endif

  // Lagrangian coordinate along surface of cylinder
  Vector<double> xi(1);

  // Point on circle
  Vector<double> point_on_circle(2);

  //Switch on the macro element
  switch(m)
   {

    //Macro element 0, is is immediately left of the cylinder
   case 0:
    
    switch(direction)
     {
     case N:
      xi[0] = 3.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(Upper_mid_left,point_on_circle,s[0],f);
      break;

     case S:  
      xi[0] = -3.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(Lower_mid_left,point_on_circle,s[0],f);
      break;

     case W:
       linear_interpolate(Lower_mid_left,Upper_mid_left,s[0],f);
      break;

     case E:
      xi[0] = 5.0*atan(1.0) - 2.0*atan(1.0)*0.5*(1.0+s[0]);
      Cylinder_pt->position(time,xi,f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    
    break;

   //Macro element 1, is immediately above the cylinder
   case 1:
    
    switch(direction)
     {
     case N:
       linear_interpolate(Upper_mid_left,Upper_mid_right,s[0],f);
      break;
      
     case S:  
      xi[0] = 3.0*atan(1.0) - 2.0*atan(1.0)*0.5*(1.0+s[0]);
      Cylinder_pt->position(time,xi,f);
      break;

     case W:
      xi[0] = 3.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(point_on_circle,Upper_mid_left,s[0],f);
      break;

     case E:
      xi[0] = 1.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(point_on_circle,Upper_mid_right,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    
    break;

    //Macro element 2, is immediately right of the cylinder
   case 2:

    switch(direction)
     {
     case N:
      xi[0] = 1.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(point_on_circle,Upper_mid_right,s[0],f);
      break;
      
     case S:  
      xi[0] = -1.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(point_on_circle,Lower_mid_right,s[0],f);
      break;

     case W:
      xi[0] = -atan(1.0) + 2.0*atan(1.0)*0.5*(1.0+s[0]);
      Cylinder_pt->position(time,xi,f);
      break;

     case E:
       linear_interpolate(Lower_mid_right,Upper_mid_right,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }

    break;
    
    //Macro element 3, is immediately below cylinder
   case 3:
    
    switch(direction)
     {
     case N:
      xi[0] = -3.0*atan(1.0) + 2.0*atan(1.0)*0.5*(1.0+s[0]);
      Cylinder_pt->position(time,xi,f);
      break;
      
     case S:  
       linear_interpolate(Lower_mid_left,Lower_mid_right,s[0],f);
      break;

     case W:
      xi[0] = -3.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(Lower_mid_left,point_on_circle,s[0],f);
      break;

     case E:
      xi[0] = -1.0*atan(1.0);
      Cylinder_pt->position(time,xi,point_on_circle);
      linear_interpolate(Lower_mid_right,point_on_circle,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }  

    break;

    //Macro element 4, is right hand block 1
   case 4:
    
    switch(direction)
     {
     case N:
       linear_interpolate(Upper_mid_right,Upper_1,s[0],f);
      break;
      
     case S:  
       linear_interpolate(Lower_mid_right,Lower_1,s[0],f);
      break;

     case W:
       linear_interpolate(Lower_mid_right,Upper_mid_right,s[0],f);
      break;

     case E:
       linear_interpolate(Lower_1,Upper_1,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }  

    break;

    //Macro element 5, is right hand block 2
   case 5:
    
    switch(direction)
     {
     case N:
       linear_interpolate(Upper_1,Upper_2,s[0],f);
      break;
      
     case S:  
       linear_interpolate(Lower_1,Lower_2,s[0],f);
      break;

     case W:
       linear_interpolate(Lower_1,Upper_1,s[0],f);
      break;

     case E:
       linear_interpolate(Lower_2,Upper_2,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }

    break;

    //Macro element 6, is right hand block 3
   case 6:
    
    switch(direction)
     {
     case N:
       linear_interpolate(Upper_2,Upper_3,s[0],f);
      break;
      
     case S:  
       linear_interpolate(Lower_2,Lower_3,s[0],f);
      break;

     case W:
       linear_interpolate(Lower_2,Upper_2,s[0],f);
      break;

     case E:
       linear_interpolate(Lower_3,Upper_3,s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }

    break;

    //Macro element 7, is right hand block 4
   case 7:
    
    switch(direction)
     {
     case N:
       linear_interpolate(Upper_3,Upper_right,s[0],f);
      break;
      
     case S:  
       linear_interpolate(Lower_3,Lower_right,s[0],f);
      break;

     case W:
       linear_interpolate(Lower_3,Upper_3,s[0],f);
      break;

     case E:
       linear_interpolate(Lower_right,Upper_right,s[0],f);
      break;
      
     default:
      
      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    break;
    
    //Macro element 8, is inlet region
   case 8:
    switch(direction)
     {
     case N:
      linear_interpolate(Upper_left,Upper_mid_left,s[0],f);
      break;
      
     case S:  
      linear_interpolate(Lower_left,Lower_mid_left,s[0],f);
      break;
      
     case W:
      linear_interpolate(Lower_left,Upper_left,s[0],f);
      break;
      
     case E:
      linear_interpolate(Lower_mid_left,Upper_mid_left,s[0],f);
      break;
      
     default:
      
      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
    break;
    
   default:
    
    std::ostringstream error_stream;
    error_stream << "Wrong macro element number" << m << std::endl;
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
 }
 
 GeomObject* const &cylinder_pt() {return Cylinder_pt;}

 
private:
 
 /// Lower left corner of rectangle
 Vector<double> Lower_left;

 /// Lower right corner of rectangle
 Vector<double> Lower_right;

 /// Where the "radial" line from circle meets lower boundary on left
 Vector<double> Lower_mid_left;

 /// Where the "radial" line from circle meets lower boundary on right
 Vector<double> Lower_mid_right;

 /// Upper left corner of rectangle
 Vector<double> Upper_left;

 /// Upper right corner of rectangle
 Vector<double> Upper_right;

 /// Where the "radial" line from circle meets upper boundary on left
 Vector<double> Upper_mid_left;

 /// Where the "radial" line from circle meets upper boundary on right
 Vector<double> Upper_mid_right;
 

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Upper_1;

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Upper_2;

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Upper_3;

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Lower_1;

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Lower_2;

 /// Coordinate of internal point in the boundaries of the domain
 Vector<double> Lower_3;
 

 /// x-coordinate of circle centre
 double Centre_x;

 /// y-coordinate of circle centre
 double Centre_y;

 /// Pointer to geometric object that represents the central cylinder
 GeomObject* Cylinder_pt;

};




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=============================================================
/// Domain-based mesh for rectangular mesh with circular hole
//=============================================================
template<class ELEMENT>
class RectangleWithHoleMesh : public virtual Mesh
{

public:


 /// Constructor: Pass pointer to geometric object that 
 /// represents the cylinder, the length and height of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction. Timestepper defaults to Steady
 /// default timestepper.
 RectangleWithHoleMesh(GeomObject* cylinder_pt, 
                       const double &length, const double &height,
                       TimeStepper* time_stepper_pt=
                       &Mesh::Default_TimeStepper) 
  {
   //Create the domain
   Domain_pt = new RectangleWithHoleDomain(cylinder_pt,length,height);

   Length = length; Height = height;
   
   //Initialise the node counter
   unsigned long node_count=0;

   //Vectors used to get data from domains
   Vector<double> s(2),r(2);
   
   //Setup temporary storage for the Node
   Vector<Node*> Tmp_node_pt;
   
   //Now blindly loop over the macro elements and associate and finite
   //element with each
   unsigned nmacro_element = Domain_pt->nmacro_element();
   for(unsigned e=0;e<nmacro_element;e++)
    {
     //Create the FiniteElement and add to the Element_pt Vector
     Element_pt.push_back(new ELEMENT);
     
     //Read out the number of linear points in the element
     unsigned np = 
      dynamic_cast<ELEMENT*>(finite_element_pt(e))->nnode_1d();
     
     //Loop over nodes in the column
     for(unsigned l1=0;l1<np;l1++)
      {
       //Loop over the nodes in the row
       for(unsigned l2=0;l2<np;l2++)
        {
         //Allocate the memory for the node
         Tmp_node_pt.push_back(finite_element_pt(e)->
                               construct_node(l1*np+l2,time_stepper_pt));
         
         //Read out the position of the node from the macro element
         s[0] = -1.0 + 2.0*(double)l2/(double)(np-1);
         s[1] = -1.0 + 2.0*(double)l1/(double)(np-1);
         Domain_pt->macro_element_pt(e)->macro_map(s,r);
         
         //Set the position of the node
         Tmp_node_pt[node_count]->x(0) = r[0];
         Tmp_node_pt[node_count]->x(1) = r[1];
         
         //Increment the node number
         node_count++;
        }
      }
    } //End of loop over macro elements
 


   //Now the elements have been created, but there will be nodes in 
   //common, need to loop over the common edges and sort it, by reassigning
   //pointers and the deleting excess nodes
   
   //Read out the number of linear points in the element
   unsigned np=dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

   //Edge between Elements 0 and 8
   for(unsigned n=0;n<np;n++)
    {
     finite_element_pt(8)->node_pt(n*np + (np-1)) =
      finite_element_pt(0)->node_pt(n*np);
     //Remover the nodes in element 8 from the temporary node list
     delete Tmp_node_pt[8*np*np + n*np + (np-1)];
     Tmp_node_pt[8*np*np + n*np + (np-1)] = 0;
    }

   //Edge between Elements 0 and 1
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 1 to be the same as in element 0
     finite_element_pt(1)->node_pt(n*np)
      = finite_element_pt(0)->node_pt((np-1)*np+np-1-n);

     //Remove the nodes in element 1 from the temporary node list
     delete Tmp_node_pt[np*np + n*np];
     Tmp_node_pt[np*np + n*np] = 0;
    }
   
   //Edge between Elements 0 and 3
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 3 to be the same as in element 0
     finite_element_pt(3)->node_pt(n*np)
      = finite_element_pt(0)->node_pt(n);

     //Remove the nodes in element 3 from the temporary node list
     delete Tmp_node_pt[3*np*np + n*np];
     Tmp_node_pt[3*np*np + n*np] = 0;
    }

   //Edge between Element 1 and 2
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 2 to be the same as in element 1
     finite_element_pt(2)->node_pt(np*(np-1)+n)
      = finite_element_pt(1)->node_pt(np*n+np-1);

     //Remove the nodes in element 2 from the temporary node list
     delete Tmp_node_pt[2*np*np + np*(np-1)+n];
     Tmp_node_pt[2*np*np + np*(np-1)+n] = 0;
    }


   //Edge between Element 3 and 2
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 2 to be the same as in element 3
     finite_element_pt(2)->node_pt(n)
      = finite_element_pt(3)->node_pt(np*(np-n-1)+np-1);

     //Remove the nodes in element 2 from the temporary node list
     delete Tmp_node_pt[2*np*np + n];
     Tmp_node_pt[2*np*np + n] = 0;
    }


   //Edge between Element 2 and 4
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 4 to be the same as in element 2
     finite_element_pt(4)->node_pt(n*np)
      = finite_element_pt(2)->node_pt(np*n+np-1);

     //Remove the nodes in element 4 from the temporary node list
     delete Tmp_node_pt[4*np*np + n*np];
     Tmp_node_pt[4*np*np + n*np] = 0;
    }

   //Edge between Element 4 and 5
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 5 to be the same as in element 4
     finite_element_pt(5)->node_pt(n*np)
      = finite_element_pt(4)->node_pt(np*n+np-1);

     //Remove the nodes in element 5 from the temporary node list
     delete Tmp_node_pt[5*np*np + n*np];
     Tmp_node_pt[5*np*np + n*np] = 0;
    }

   //Edge between Element 5 and 6
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 6 to be the same as in element 5
     finite_element_pt(6)->node_pt(n*np)
      = finite_element_pt(5)->node_pt(np*n+np-1);

     //Remove the nodes in element 6 from the temporary node list
     delete Tmp_node_pt[6*np*np + n*np];
     Tmp_node_pt[6*np*np + n*np] = 0;
    }

   //Edge between Element 6 and 7
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element 7 to be the same as in element 6
     finite_element_pt(7)->node_pt(n*np)
      = finite_element_pt(6)->node_pt(np*n+np-1);

     //Remove the nodes in element 7 from the temporary node list
     delete Tmp_node_pt[7*np*np + n*np];
     Tmp_node_pt[7*np*np + n*np] = 0;
    }

   //Now set the actual true nodes
   for(unsigned long n=0;n<node_count;n++)
    {
     if(Tmp_node_pt[n]!=0) {Node_pt.push_back(Tmp_node_pt[n]);}
    }

   //Finally set the nodes on the boundaries
   set_nboundary(5);
   
   for(unsigned n=0;n<np;n++)
    {
     //Left hand side
     Node* nod_pt=finite_element_pt(8)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(3,nod_pt);
     
     //Right hand side
     nod_pt=finite_element_pt(7)->node_pt(n*np+np-1);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(1,nod_pt);
     
     //First part of lower boundary
     nod_pt = finite_element_pt(8)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);

     //First part of upper boundary
     nod_pt = finite_element_pt(8)->node_pt(np*(np-1) + n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
     
     //First part of hole boundary
     nod_pt=finite_element_pt(3)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }

   for(unsigned n=1;n<np;n++)
    {
     //Second part of lower boundary
     Node* nod_pt=finite_element_pt(3)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);

     //Third part of lower boundary
     nod_pt=finite_element_pt(4)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);
     
     //Second part of upper boundary
     nod_pt=finite_element_pt(1)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
     
     //Third part of upper boundary                                
     nod_pt=finite_element_pt(4)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
     
     //Next part of hole
     nod_pt=finite_element_pt(2)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }
   
   for(unsigned n=1;n<np;n++)
    {
     //Third part of lower boundary
     Node* nod_pt=finite_element_pt(5)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);
     
     //Third part of upper boundary                                
     nod_pt=finite_element_pt(5)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
    }

   for(unsigned n=1;n<np;n++)
    {
     //Fourth part of lower boundary
     Node* nod_pt=finite_element_pt(6)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);
     
     //Fourth part of upper boundary                                
     nod_pt=finite_element_pt(6)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
    }

   for(unsigned n=1;n<np;n++)
    {
     //Final part of lower boundary
     Node* nod_pt=finite_element_pt(7)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);

     //Final part of upper boundary                                
     nod_pt=finite_element_pt(7)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
     
     //Next part of hole
     nod_pt=finite_element_pt(1)->node_pt(np-n-1);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }

   for(unsigned n=1;n<np-1;n++)
    {
     //Final part of hole
     Node* nod_pt=finite_element_pt(0)->node_pt(np*(np-n-1)+np-1);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }
   

 }

 /// Access function to the domain
 RectangleWithHoleDomain* domain_pt() {return Domain_pt;}

 /// Access function to length
 const double &length() {return Length;}
 /// Access function to height
 const double &height() {return Height;}
 
 
protected:

 /// Pointer to the domain
 RectangleWithHoleDomain* Domain_pt;

 double Length;

 double Height;
 
};


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//===================================================================
/// Refineable version of RectangleWithHoleMesh. Applies one uniform
/// refinement immediately to avoid problems with the automatic
/// applications of boundary conditions in subsequent refinements
//===================================================================
template<class ELEMENT>
class RefineableRectangleWithHoleMesh :
 public RectangleWithHoleMesh<ELEMENT>, public RefineableQuadMesh<ELEMENT>
{
public: 

 /// Constructor. Pass pointer to geometric object that 
 /// represents the cylinder, the length and height of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction. Timestepper defaults to Steady
 /// default timestepper.
 RefineableRectangleWithHoleMesh(GeomObject* cylinder_pt, 
                               const double &length, const double &height,
                               TimeStepper* time_stepper_pt=
                               &Mesh::Default_TimeStepper) :
  RectangleWithHoleMesh<ELEMENT>(cylinder_pt,length,height,time_stepper_pt) 
  {

   // Nodal positions etc. were created in constructor for
   // Cylinder...<...>. Need to setup adaptive information.

   // Loop over all elements and set macro element pointer
   for (unsigned e=0;e<9;e++)
    {
     dynamic_cast<ELEMENT*>(this->element_pt(e))->
      set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }
   
   // Setup quadtree forest for mesh refinement
   this->setup_quadtree_forest();
  }
 
 
 ///  Destructor: Empty
 virtual ~RefineableRectangleWithHoleMesh() {}
 
};


//Make my own NavierStokes class
class MyNavierStokesElement : public RefineableQCrouzeixRaviartElement<2>
{
public:

 //Constructor
 MyNavierStokesElement()
  {}

 //Overload the fill contribution to jacobian
  /// Compute the element's residual Vector and the jacobian matrix
 /// Virtual function can be overloaded by hanging-node version
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_nst(residuals,jacobian,
                                         GeneralisedElement::Dummy_matrix,1);
   //Sort out the Reynolds number as added in
   this->fill_in_jacobian_from_external_by_fd(residuals,jacobian);
  }
 
};




//===================================================================
/// Flow around a cylinder in rectangular domain
//===================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
class FlowAroundCylinderProblem : public Problem
{
 bool Eigenproblem_flag;
 
public:

 /// Constructor: Pass geometric object that represents
 /// central cylinder, and length and height of domain.
 FlowAroundCylinderProblem(GeomObject* cylinder_pt, 
                           const double &length, 
                           const double &height);

 /// Destructor: clean up the memory
 ~FlowAroundCylinderProblem();

 ///  Set the boundary conditions on the cylinder 
 /// and at the inlet
 void set_boundary_conditions()
  {
   using namespace Global_Parameters;
   
   //Set the boundary conditions on the cylinder
   unsigned ibound=4;
   unsigned n_node = Base_flow_mesh_pt->nboundary_node(ibound);
   for(unsigned n=0;n<n_node;n++)
    {
     double x = Base_flow_mesh_pt->boundary_node_pt(ibound,n)->x(0);
     double y = Base_flow_mesh_pt->boundary_node_pt(ibound,n)->x(1);
     
     //Now find the vector distance to the centre
     double len_x = x - base_flow_mesh_pt()->domain_pt()->centre_x();
     double len_y = y - base_flow_mesh_pt()->domain_pt()->centre_y();
     
     //Calculate the angle and radius
     double theta = atan2(len_y,len_x);
     
     double u_x = -Alpha*sin(theta);
     double u_y =  Alpha*cos(theta);

     //Now set the velocities
     Base_flow_mesh_pt->boundary_node_pt(ibound,n)->set_value(0,u_x);
     Base_flow_mesh_pt->boundary_node_pt(ibound,n)->set_value(1,u_y);
    } 
   
   // update the parabolic tangential inflow velocity
   ibound=3;
   double ycoord,uy;
   n_node = Base_flow_mesh_pt->nboundary_node(ibound);
   double um = 1.0/12.0 - 0.25*Domain_height*Domain_height;
   for (unsigned n=0;n<n_node;n++)
    {
     ycoord = Base_flow_mesh_pt->boundary_node_pt(ibound,n)->x(1);
     uy = (ycoord*ycoord - 0.25*Domain_height*Domain_height)/um;

     // Tangential flow
     Base_flow_mesh_pt->boundary_node_pt(ibound,n)->set_value(0,uy);
     //No penetration
     Base_flow_mesh_pt->boundary_node_pt(ibound,n)->set_value(1,0.0);
    }

   if(Eigenproblem_flag)
    {

 unsigned n_boundary = Eigenproblem_mesh_pt->nboundary();
 for(unsigned b=0;b<n_boundary;++b)
  {
   unsigned n_boundary_node = Eigenproblem_mesh_pt->nboundary_node(b);
   for(unsigned n=0;n<n_boundary_node;n++)
    {
     if(b!=1)
      {
       //Pin all four (two real, two imaginary) velocity components on the
       //same boundaries as the base flow problem
       for(unsigned i=0;i<4;i++)
        {
         Eigenproblem_mesh_pt->boundary_node_pt(b,n)->pin(i);
         Eigenproblem_mesh_pt->boundary_node_pt(b,n)->set_value(i,0.0);
        }
      }
    }
  }

 //Pin all normalisation degrees of freedom
 unsigned n_node = Eigenproblem_mesh_pt->nnode();
 for(unsigned n=0;n<n_node;++n)
  {
   for(unsigned i=0;i<4;++i)
    {
     Eigenproblem_mesh_pt->node_pt(n)->pin(4+i);
    }
  }
    }
  }

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
  {
   //Dump out the eigenvalue
   if(Eigenproblem_flag)
    {
        ///Report the eigenvalue
   std::cout << "Eigenvalue is " <<
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(0)
             << " + " << 
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(1)
             << "I\n";

   std::cout << "Critical Reynolds number is " <<
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(2)
             << "\n";
    }
  }

 void actions_before_newton_convergence_check() {set_boundary_conditions();}

 ///  Update the problem specs before solve (empty; all prescribed
 /// velocities are constant along their respective boundares, therefore
 /// their FE interpolation onto the newly created nodes is sufficiently
 /// accurate)
 void actions_before_newton_solve() {}

 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   //Pin both velocities at all boundaries
   //This is required in case the automatic stuff goes wrong
   unsigned num_bound = Base_flow_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned num_nod= Base_flow_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Parallel, axially traction free outflow at downstream end
       if (ibound != 1)
        {
         Base_flow_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
         Base_flow_mesh_pt->boundary_node_pt(ibound,inod)->pin(1);
        }
      }
    }
   
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(Base_flow_mesh_pt->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(Base_flow_mesh_pt->element_pt());

   if(Eigenproblem_flag)
    {
     // Unpin all pressure dofs
     RefineableLinearisedNavierStokesEquations::
      unpin_all_pressure_dofs(Eigenproblem_mesh_pt->element_pt());
     
     // Pin redudant pressure dofs
     RefineableLinearisedNavierStokesEquations::
      pin_redundant_nodal_pressures(Eigenproblem_mesh_pt->element_pt());
     
     // First interaction (base state velocities)
     Multi_domain_functions::
      setup_multi_domain_interaction<BASE_ELEMENT>
      (this,Eigenproblem_mesh_pt,Base_flow_mesh_pt,0);
     
     // Second interaction (base state velocity derivatives w.r.t. space)
     Multi_domain_functions::
      setup_multi_domain_interaction<BASE_ELEMENT>
      (this,Eigenproblem_mesh_pt,Base_flow_mesh_pt,1);
     
     //Reset the boundary conditions to be safe
     this->set_boundary_conditions();
     
     //Loop over all the (fluid) elements 
     unsigned n_element = Eigenproblem_mesh_pt->nelement();
     for(unsigned e=0;e<n_element;e++)
      {
       //Cast to the particular element type, this is necessary because
       //the base elements don't have the member functions that we're about
       //to call!
       PERTURBED_ELEMENT *el_pt = 
        dynamic_cast<PERTURBED_ELEMENT*>(Eigenproblem_mesh_pt->element_pt(e));
       
       //Pin the pressure normalisation dofs (can't be handled automatically
       //because it's a choice we've made)
       el_pt->pin_pressure_normalisation_dofs();
      }

     {
      //Pull out a pointer to the Reynolds number
      double* local_re_pt = Normalisation_mesh_pt->element_pt(0)
       ->internal_data_pt(0)->value_pt(2);
      
      // Pass pointer to Reynolds number to elements
      unsigned nelem=Base_flow_mesh_pt->nelement();
      for (unsigned e=0;e<nelem;e++)
       {
        BASE_ELEMENT* el_pt = 
         dynamic_cast<BASE_ELEMENT*>(Base_flow_mesh_pt->element_pt(e));
        
        //Set the Reynolds number for each element 
        //(yes we could have different Reynolds number in each element!!)
        el_pt->re_pt() = local_re_pt;//&Re;
        
        //Set the product of Reynolds and Strouhal numbers
        el_pt->re_st_pt() = local_re_pt; //&Re;
        
        //and add as external data
        el_pt->add_external_data(
         Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0));
       }
     }
    }

  }

 /// Access function for the specific mesh
 RefineableRectangleWithHoleMesh<BASE_ELEMENT>* base_flow_mesh_pt() 
  {
   return dynamic_cast<RefineableRectangleWithHoleMesh<BASE_ELEMENT>*>
    (Base_flow_mesh_pt);
  }

 /// Return a pointer to the specific mesh used
 RefineableRectangleWithHoleMesh<PERTURBED_ELEMENT>* eigenproblem_mesh_pt() 
  {return dynamic_cast<RefineableRectangleWithHoleMesh<PERTURBED_ELEMENT>*>
    (this->Eigenproblem_mesh_pt);}

 
 void add_eigenproblem()
  {
   Eigenproblem_flag = true;
   using namespace Global_Parameters;
  
   GeomObject* cylinder_pt = base_flow_mesh_pt()->domain_pt()->cylinder_pt();
   double length = base_flow_mesh_pt()->length();
   double height = base_flow_mesh_pt()->height();
   
   // Build mesh
   Eigenproblem_mesh_pt =
    new RefineableRectangleWithHoleMesh<PERTURBED_ELEMENT>(cylinder_pt,
                                                           length,height);
   
   // Set error estimator (not sure this will work with the peturbed problem) 
   Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
   eigenproblem_mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
   // Maximum number of refinements (increase this if you want a finer mesh)
   //eigenproblem_mesh_pt()->max_refinement_level() = max_refinement_level;
   
   // Error targets for adaptive refinement
   //eigenproblem_mesh_pt()->max_permitted_error() = max_error_target; 
   //eigenproblem_mesh_pt()->min_permitted_error() = min_error_target; 

   eigenproblem_mesh_pt()->max_refinement_level() = 10;
   
   //that the natural boundary condition is pseudo-traction-free
   LinearisedNavierStokesEquations::Gamma[0] = 0.0;
   LinearisedNavierStokesEquations::Gamma[1] = 0.0;
   
   //Sort out normalisation
   using namespace Global_Parameters;
   
   LinearisedNavierStokesEigenfunctionNormalisationElement*
    normalisation_element_pt =
    new LinearisedNavierStokesEigenfunctionNormalisationElement(
     &Eigenfunction_normalisation);
   
   Normalisation_mesh_pt = new Mesh;
   Normalisation_mesh_pt->add_element_pt(normalisation_element_pt);
   
   //Pull out a pointer to the Reynolds number
   double* local_re_pt = Normalisation_mesh_pt->element_pt(0)
    ->internal_data_pt(0)->value_pt(2);
   //Set the value
   *local_re_pt = Re;
   
   //Need to set the parameters in the elements in the base mesh
   unsigned n_element = Eigenproblem_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to the particular element type, this is necessary because
     //the base elements don't have the member functions that we're about
     //to call!
     PERTURBED_ELEMENT *el_pt = 
      dynamic_cast<PERTURBED_ELEMENT*>(Eigenproblem_mesh_pt->element_pt(e));
     
     //There is no need for ALE
     el_pt->disable_ALE();
     
     //Set the Reynolds number for each element 
     //(yes we could have different Reynolds number in each element!!)
     el_pt->re_pt() = local_re_pt;//&Re;
     
     //Set the product of Reynolds and Strouhal numbers
     el_pt->re_st_pt() = local_re_pt; //&Re;
     
     //Set the eigenfunction normalisation element
     el_pt->set_eigenfunction_normalisation_element(normalisation_element_pt);
    }
   

   //Also set the Reynolds number in the original elements to be variable
   {
    // Pass pointer to Reynolds number to elements
    unsigned nelem=Base_flow_mesh_pt->nelement();
    for (unsigned e=0;e<nelem;e++)
     {
      BASE_ELEMENT* el_pt = 
       dynamic_cast<BASE_ELEMENT*>(Base_flow_mesh_pt->element_pt(e));
      
      //Set the Reynolds number for each element 
      //(yes we could have different Reynolds number in each element!!)
      el_pt->re_pt() = local_re_pt;//&Re;
      
      //Set the product of Reynolds and Strouhal numbers
      el_pt->re_st_pt() = local_re_pt; //&Re;

      //and add as external data
      el_pt->add_external_data(
       Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0));

     }
   }
   
   //Refine as the other mesh
   dynamic_cast<TreeBasedRefineableMeshBase*>(Eigenproblem_mesh_pt)
    ->refine_base_mesh_as_in_reference_mesh(
     dynamic_cast<TreeBasedRefineableMeshBase*>
     (Base_flow_mesh_pt));
   
   
   //Set the boundary conditions (no slip on the walls)
   //Loop over the nodes on the (only) mesh boundary
   unsigned n_boundary = Eigenproblem_mesh_pt->nboundary();
   for(unsigned b=0;b<n_boundary;++b)
    {
     unsigned n_boundary_node = Eigenproblem_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_boundary_node;n++)
      {
       if(b!=1)
        {
         //Pin all four (two real, two imaginary) velocity components on the
         //same boundaries as the base flow problem
         for(unsigned i=0;i<4;i++)
          {
           Eigenproblem_mesh_pt->boundary_node_pt(b,n)->pin(i);
          }
        }
      }
    }
   
   //Pin all normalisation degrees of freedom
   unsigned n_node = Eigenproblem_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;++n)
    {
     for(unsigned i=0;i<4;++i)
      {
       Eigenproblem_mesh_pt->node_pt(n)->pin(4+i);
      }
    }
   
   //Loop over all the (fluid) elements 
   {unsigned n_element = Eigenproblem_mesh_pt->nelement();
    for(unsigned e=0;e<n_element;e++)
     {
      //Cast to the particular element type, this is necessary because
      //the base elements don't have the member functions that we're about
      //to call!
      PERTURBED_ELEMENT *el_pt = 
       dynamic_cast<PERTURBED_ELEMENT*>(Eigenproblem_mesh_pt->element_pt(e));
      
      //There is no need for ALE
      el_pt->disable_ALE();
      
      //Set the Reynolds number for each element 
      //(yes we could have different Reynolds number in each element!!)
      el_pt->re_pt() = local_re_pt;//&Re;
      
      //Set the product of Reynolds and Strouhal numbers
      el_pt->re_st_pt() = local_re_pt; //&Re;
      
      //Set the eigenfunction normalisation element
      el_pt->set_eigenfunction_normalisation_element(normalisation_element_pt);
      
      //Pin the pressure normalisation dofs
      el_pt->pin_pressure_normalisation_dofs();
     }
   }
   
   
// Pin redudant pressure dofs
   RefineableLinearisedNavierStokesEquations::
    pin_redundant_nodal_pressures(Eigenproblem_mesh_pt->element_pt());
 
 //The setup of the interaction uses a binning approach that
 //works best if there are a reasonably large number of bins
 //in each direction. If there are ever problems, increasing
 //these numbers tends to fix them.

 // Change default number of bins in each dimension
 //Multi_domain_functions::Nx_bin=500;
 //Multi_domain_functions::Ny_bin=500;

 // Don't compute extreme bin coordinates
 //Multi_domain_functions::Compute_extreme_bin_coordinates=false;
 
 // Specify extreme bin coordinates directly
 //Multi_domain_functions::X_min= 1.0/Delta - 1.0;
 //Multi_domain_functions::X_max= 1.0/Delta + 1.0;
 //Multi_domain_functions::Y_min = -1.0;
 //Multi_domain_functions::Y_max= + 1.0;

 // First interaction (base state velocities)
 Multi_domain_functions::
  setup_multi_domain_interaction<BASE_ELEMENT>
  (this,Eigenproblem_mesh_pt,Base_flow_mesh_pt,0);

 // Second interaction (base state velocity derivatives w.r.t. space)
 Multi_domain_functions::
   setup_multi_domain_interaction<BASE_ELEMENT>
  (this,Eigenproblem_mesh_pt,Base_flow_mesh_pt,1);

 //Build the global mesh
 this->add_sub_mesh(Eigenproblem_mesh_pt);
 this->add_sub_mesh(Normalisation_mesh_pt);
 this->rebuild_global_mesh();

 Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->pin(2);
 
 //Setup all the equation numbering and look-up schemes 
 std::cout << assign_eqn_numbers() << std::endl; 
  }


 void set_eigenvalue(const double &real, const double &imag)
  {
   Data* dat_pt =
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0);
   dat_pt->set_value(0,real);
   dat_pt->set_value(1,imag);
  }

 void pin_all_base_flow_dofs()
  {
   const unsigned n_node = Base_flow_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;++n)
    {
     Node* nod_pt = Base_flow_mesh_pt->node_pt(n);
     unsigned n_value = nod_pt->nvalue();
     for(unsigned i=0;i<n_value;++i)
      {
       if(!nod_pt->is_constrained(i))
        {
         nod_pt->pin(i);
        }
      }
    }
   const unsigned n_element = Base_flow_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;++e)
    {
     GeneralisedElement* el_pt = Base_flow_mesh_pt->element_pt(e);
     unsigned n_internal = el_pt->ninternal_data();
     for(unsigned i=0;i<n_internal;++i)
      {
       Data* data_pt = el_pt->internal_data_pt(i);
       unsigned n_value = data_pt->nvalue();
       for(unsigned j=0;j<n_value;++j)
        {
         if(!data_pt->is_constrained(j))
          {
           data_pt->pin(j);
          }
        }
      }
    }
  }

  void unpin_all_base_flow_dofs()
  {
   const unsigned n_node = Base_flow_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;++n)
    {
     Node* nod_pt = Base_flow_mesh_pt->node_pt(n);
     unsigned n_value = nod_pt->nvalue();
     for(unsigned i=0;i<n_value;++i)
      {
       if(!nod_pt->is_constrained(i))
        {
         nod_pt->unpin(i);
        }
      }
    }
   const unsigned n_element = Base_flow_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;++e)
    {
     GeneralisedElement* el_pt = Base_flow_mesh_pt->element_pt(e);
     unsigned n_internal = el_pt->ninternal_data();
     for(unsigned i=0;i<n_internal;++i)
      {
       Data* data_pt = el_pt->internal_data_pt(i);
       unsigned n_value = data_pt->nvalue();
       for(unsigned j=0;j<n_value;++j)
        {
         if(!data_pt->is_constrained(j))
          {
           data_pt->unpin(j);
          }
        }
      }
    }
  }

 
 
 //Transfer eigenfunction onto the mesh
 //N.B. This assumes that the discretisation (and boundary conditions)
 //are exactly the same as in the base mesh
 void transfer_eigenfunction_as_initial_condition(DoubleVector &real,
                                                  DoubleVector &imaginary)
  {
   this->pin_all_base_flow_dofs();

   //Pin the eigenvalues
   Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->pin_all();
   
   //Pin everything apart from the real part of the eigenfunction
   unsigned n_element = Eigenproblem_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;++e)
    {
     PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
      (Eigenproblem_mesh_pt->element_pt(e));

     //Pin the imaginary parts
     el_pt->pin_real_or_imag(1);
    }

   
     //Reassign the equation numbers
     std::cout << "Assigning real part" << this->assign_eqn_numbers()
               << std::endl;

     //Transfer it over
     this->assign_eigenvector_to_dofs(real);
     Eigenproblem_mesh_pt->output("E1.dat",5);
     
     
     //Pin everything apart from the imaginary part of the eigenfunction
     for(unsigned e=0;e<n_element;++e)
      {
       PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
      (Eigenproblem_mesh_pt->element_pt(e));

       el_pt->unpin_real_or_imag(1);
       //Pin the real parts
       el_pt->pin_real_or_imag(0);
      }

     //Reapply the boundary conditions
     this->set_boundary_conditions();
   
     //Reassign the equation numbers
     std::cout << "Assigning imaginary part" << this->assign_eqn_numbers()
               << std::endl;

     //Transfer it over
     this->assign_eigenvector_to_dofs(imaginary);

     //Now we convert things over to the normalisation function
     //Pin everything apart from the imaginary part of the eigenfunction
     for(unsigned e=0;e<n_element;++e)
      {
       PERTURBED_ELEMENT* el_pt = dynamic_cast<PERTURBED_ELEMENT*>
        (Eigenproblem_mesh_pt->element_pt(e));
       //Unpin the real part
       el_pt->unpin_real_or_imag(0);
       //Transfer the normalisation
       el_pt->copy_efunction_to_normalisation();
      }
     
     //Reapply the boundary conditions
     this->set_boundary_conditions();

     //Unpin the eigenvalues
     Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->unpin_all();
     //Pin the Reynolds number
     //Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->pin(2);

     unpin_all_base_flow_dofs();
     this->actions_after_adapt();

     //Reassign the equation numbers
     std::cout << "Re-assigning the equation numbers"
               << this->assign_eqn_numbers()
               << std::endl;

     //Let's output
     Eigenproblem_mesh_pt->output("Initial_mesh.dat",5);
     
     
  }

 void doc_solution()
  {
   Base_flow_mesh_pt->output("final_base_flow.dat",5);
   //Output the eigenproblem
   Eigenproblem_mesh_pt->output("final_eigenfunction.dat",5);
   ///Report the eigenvalue
   std::cout << "Eigenvalue is " <<
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(0)
             << " + " << 
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(1)
             << "I\n";

   std::cout << "Critical Reynolds number is " <<
    Normalisation_mesh_pt->element_pt(0)->internal_data_pt(0)->value(2)
             << "\n";
  }

 
 private:

 Mesh* Base_flow_mesh_pt;

 Mesh* Eigenproblem_mesh_pt;

 Mesh* Normalisation_mesh_pt;
 
 /// Height of the domain
 double Domain_height;
 
 /// Length of the domain
 double Domain_length;
 
 /// The geometric cylinder
 GeomObject *Cylinder_pt;
 
};





//========================================================================
/// Constructor 
//========================================================================
template<class BASE_ELEMENT, class PERTURBED_ELEMENT>
FlowAroundCylinderProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
FlowAroundCylinderProblem(
 GeomObject* cylinder_pt, const double &length, const double &height) :
 Eigenproblem_flag(false), Domain_height(height), Domain_length(length),
 Cylinder_pt(cylinder_pt)
 

{
 //this->linear_solver_pt() = new HSL_MA42;
 
 //Increase the maximum residuals so that we can get convergence on
 //the coarsest mesh
 Max_residuals = 100.0;
 
 if(!Global_Parameters::Read_in_eigenfunction_from_disk)
  {
   this->eigen_solver_pt() = new ARPACK;
   static_cast<ARPACK*>(eigen_solver_pt())->set_shift(50.0);
   static_cast<ARPACK*>(eigen_solver_pt())->narnoldi() = 70;
  }

 // Build mesh
 Base_flow_mesh_pt=
  new RefineableRectangleWithHoleMesh<BASE_ELEMENT>(cylinder_pt,length,height);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 base_flow_mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(base_flow_mesh_pt()->element_pt());
 
 //Make the problem conventional form, so 
 //that the natural boundary condition is pseudo-traction-free
 NavierStokesEquations<2>::Gamma[0] = 0.0;
 NavierStokesEquations<2>::Gamma[1] = 0.0;
 
  using namespace Global_Parameters;

  // Pass pointer to Reynolds number to elements
  unsigned nelem=Base_flow_mesh_pt->nelement();
  for (unsigned e=0;e<nelem;e++)
   {
    dynamic_cast<BASE_ELEMENT*>(Base_flow_mesh_pt->element_pt(e))->re_pt()= &Re;
    dynamic_cast<BASE_ELEMENT*>(Base_flow_mesh_pt->element_pt(e))->re_st_pt()=&Re;
   }

  this->add_sub_mesh(Base_flow_mesh_pt);
  this->build_global_mesh();
  
}

//========================================================================
/// Destructor clean up all memory allocated in the problem
//=========================================================================
template<class BASE_ELEMENT,class PERTURBED_ELEMENT>
FlowAroundCylinderProblem<BASE_ELEMENT,PERTURBED_ELEMENT>::
~FlowAroundCylinderProblem()
{
 //Delete the error estimator
 delete base_flow_mesh_pt()->spatial_error_estimator_pt();
 //Delete the mesh
 delete Base_flow_mesh_pt;
}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=====================================================================
/// Driver
//=====================================================================
int main()
{
 // Length and height of domain
 double length = 14.0;
 double height = 1.0/Global_Parameters::B;

 //Create a new ellipse object with equal semi-minor and semi-major axes
 //to be the central cylinder
 GeneralEllipse* cylinder_pt = 
  new GeneralEllipse(4.5,0.0,0.5,0.5);

 // Create Problem
 FlowAroundCylinderProblem 
  <MyNavierStokesElement/*RefineableQCrouzeixRaviartElement<2>*/,
 RefineableLinearisedQCrouzeixRaviartMultiDomainElement>
  problem(cylinder_pt,length,height);

 //Refine the problem uniformly a couple of times
 problem.refine_uniformly(); problem.refine_uniformly();

 // Solve adaptively with up to two rounds of refinement
 problem.newton_solve(2);

 //Now increase the Reynolds number to 90 in steps of 5
 //without any further adaptation
 for(unsigned i=0;i<3;i++)
  {
   Global_Parameters::Re += 5.0;
   problem.newton_solve();
  }

 //Assign memory for the eigenvalues and eigenvectors
 Vector<DoubleVector> eigenvectors;
 double frequency = 0.0;

 //If we are reading in from the disk
 if(Global_Parameters::Read_in_eigenfunction_from_disk)
  {
   //Read in the  eigenvalue from the data file
   std::ifstream input("eigen.dat");
   input >> frequency;
   
   //Read in the eigenvector from the data file
   const unsigned n_dof = problem.ndof();
   eigenvectors.resize(2);
   LinearAlgebraDistribution dist(problem.communicator_pt(),n_dof,false);
   //Rebuild the vector
   eigenvectors[0].build(&dist,0.0);
   eigenvectors[1].build(&dist,0.0);

   for(unsigned n=0;n<n_dof;n++)
    {
     input >> eigenvectors[0][n];
     input >> eigenvectors[1][n];
    }
   input.close();
  }
 //Otherwise solve the eigenproblem
 else
  {
   Vector<std::complex<double> > eigenvalues;
   //Now solve the eigenproblem
   problem.solve_eigenproblem(6,eigenvalues,eigenvectors);
   frequency = eigenvalues[0].imag();
  }


 //Create the perturbation problem and solve the eigenproblem on it
 //FlowAroundCylinderLinPertProblem 
 // <RefineableQCrouzeixRaviartElement<2>,RefineableLinearisedQCrouzeixRaviartMu//ltiDomainElement
 // > pert_problem(problem.mesh_pt());

 problem.add_eigenproblem();
 
 //Transfer eigenfunction across to perturbation mesh
 problem.transfer_eigenfunction_as_initial_condition(eigenvectors[0],
                                                     eigenvectors[1]);
 
 problem.set_eigenvalue(0.0,frequency);

 problem.newton_solve(2); //Slow!

 problem.doc_solution();
 
 
 exit(1);
 //Try to find the Hopf exactly by using the initial
 //guesses for the eigenvalue and eigenvector from
 //the data file
 problem.activate_hopf_tracking(&Global_Parameters::Re,
                                frequency,
                                eigenvectors[0],      
                                eigenvectors[1]);      
 //Solve the problem 
 problem.newton_solve();
 //Report the value of the bifurcation
 std::cout << "Hopf bifurcation found at  " << Global_Parameters::Re << " " <<
  problem.dof(problem.ndof()-1) << std::endl;

 //Solve it again with three rounds of adaptation
 problem.newton_solve(2);

 std::cout << "Hopf bifurcation found at " << Global_Parameters::Re << " " <<
  problem.dof(problem.ndof()-1) << std::endl;
 problem.mesh_pt()->output("bif_soln.dat");

 //Output the value of the critical Reynolds number into a data file
 std::ofstream trace("trace.dat");
 trace << Global_Parameters::Re << " " <<
  problem.dof(problem.ndof()-1) << "\n";
 trace.close();
}

