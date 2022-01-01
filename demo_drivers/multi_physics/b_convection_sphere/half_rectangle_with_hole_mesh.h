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
#ifndef HALF_RECTANGLE_WITH_MESH_HEADER
#define HALF_RECTANGLE_WITH_MESH_HEADER

#include "generic.h"

namespace oomph
{

using namespace QuadTreeNames;

//===============================================
/// My own Ellipse class
//===============================================
class HalfEllipse : public GeomObject
{

public:
 
 /// Constructor
 HalfEllipse(const double &centre_z,
                const double &a, const double &b)
  : GeomObject(1,2), Centre_z(centre_z), A(a), B(b)
  {}

 /// Destructor (empty)
 ~HalfEllipse(){}

 /// Return the position of the Half Ellipse
 void position(const Vector<double> &xi, Vector<double> &r) const
  {
   r[0] = A*sin(xi[0]);
   r[1] = Centre_z + B*cos(xi[0]);
  }


 /// Access to  z-coordinate of centre
 double& centre_z(){return Centre_z;}

 /// Access to r-half axis
 double& a(){return A;}

 /// Access to z-half axis
 double& b(){return B;}


private:

 /// z-coordinate of centre
 double Centre_z;

 /// r-half axis
 double A;

 /// z-half axis
 double B;

};



/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////



//===========================================================
/// Rectangular domain with Half-elliptical hole
//===========================================================
class HalfRectangleWithHoleDomain : public Domain
{

public:


 /// Constructor, Pass pointer to geometric object that 
 /// represents the cylinder, the length (z) and half-height (r) of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction.
 HalfRectangleWithHoleDomain(GeomObject* cylinder_pt,
                             const double &radius, 
                             const double &length,
                             const double &up_length,
                             const unsigned &nup,
                             const double &down_length,
                             const unsigned &ndown,
                             const double &width_near_cylinder,
                             const unsigned &ncolumn) :
 Cylinder_pt(cylinder_pt), Nup(nup), Ndown(ndown), Ncolumn(ncolumn)
  {   

   //Axial spacing between lines
   double up_spacing = up_length/(double)nup;
   double down_spacing = down_length/(double)ndown;
   
   //Resize the storage for the corners of the squared
   Up.resize(ncolumn+1); Down.resize(ncolumn+1);

   //The first column is special
   {
    unsigned j=0;
    Up[j].resize(nup+1); Up[j].resize(nup+1);
    //Set the coordinates
    for(unsigned i=0;i<(nup+1);i++)
     {
      Up[j][i].resize(2); 
      Up[j][i][0] = 0.0; Up[j][i][1] = i*up_spacing;
     }

     //There are going to be ndown+1 lines in the downstream region
     Down[j].resize(ndown+1); Down[j].resize(ndown+1);
     
     //Set the coordinates
     for(unsigned i=0;i<(ndown+1);i++)
      {
       Down[j][i].resize(2); Down[j][i].resize(2);
       Down[j][i][0] = 0.0; 
       Down[j][i][1] = length - (ndown - i)*down_spacing;
      }
   }

   //What is the column spacing.
   double radial_start = radius;
   double radial_spacing = 0.0;
    
   if(ncolumn > 1)
    {
     radial_start = width_near_cylinder;
     radial_spacing = (radius - width_near_cylinder)/(double)(ncolumn-1);
    }

   //There are going to be nup+1 lines in the upstream region
   for(unsigned j=1;j<(ncolumn+1);j++)
    {
     Up[j].resize(nup+1); Up[j].resize(nup+1);
     //Set the coordinates
     for(unsigned i=0;i<(nup+1);i++)
      {
       Up[j][i].resize(2); 
       Up[j][i][0] = radial_start + (j-1)*radial_spacing; 
       Up[j][i][1] = i*up_spacing;
      }

     //There are going to be ndown+1 lines in the downstream region
     Down[j].resize(ndown+1); Down[j].resize(ndown+1);
     
     //Set the coordinates
     for(unsigned i=0;i<(ndown+1);i++)
      {
       Down[j][i].resize(2); Down[j][i].resize(2);
       Down[j][i][0] = radial_start + (j-1)*radial_spacing; 
       Down[j][i][1] = length - (ndown - i)*down_spacing;
      }
    }

   //There are three + nup + ndown macro elements in the first column
   //Plus the additional columns of 1 + nup + ndowm
   unsigned n_macro = 3 + nup + ndown + (ncolumn-1)*(1 + nup + ndown);
   Macro_element_pt.resize(n_macro); 

   // Build the 2D macro elements
   for (unsigned i=0;i<n_macro;i++)
    {Macro_element_pt[i]= new QMacroElement<2>(this,i);}
  }



 /// Destructor: Empty; cleanup done in base class
 ~HalfRectangleWithHoleDomain() {}

 /// Helper function to interpolate linearly between the
 /// "right" and "left" points; \f$ s \in [-1,1] \f$
 void linear_interpolate(Vector<double> left, Vector<double> right,
                         const double &s, Vector<double> &f)
  {
   for(unsigned i=0;i<2;i++)
    {
     f[i] = left[i] + (right[i] - left[i])*0.5*(s+1.0);
    }
  }

   

 /// Parametrisation of macro element boundaries: f(s) is the position
 /// vector to macro-element m's boundary in the specified direction [N/S/E/W]
 /// at the specfied discrete time level (time=0: present; time>0: previous)
 void macro_element_boundary(const unsigned &time,
                             const unsigned &m,
                             const unsigned &direction,
                             const Vector<double> &s,
                             Vector<double>& f)
 {
  // Lagrangian coordinate along surface of cylinder
  Vector<double> xi(1);

  // Point on circle
  Vector<double> point_on_circle(2);

  //Upstream region all rectangular blocks
  if(m < Nup)
   {
    //Switch on the direction
    
    switch(direction)
     {
     case N:
      linear_interpolate(Up[0][m+1],Up[1][m+1],s[0],f);
      break;
      
     case S:  
       linear_interpolate(Up[0][m],Up[1][m],s[0],f);
      break;

     case W:
      linear_interpolate(Up[0][m],Up[0][m+1],s[0],f);
      break;

     case E:
      linear_interpolate(Up[1][m],Up[1][m+1],s[0],f);
      break;

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(
       error_stream.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }  
   }
  //The special cases around the half-domain
  else if(m < Nup + 3)
   {
    //Scale the macro element number
    unsigned m_mod = m - Nup;
    
    //Switch on the macro element
    switch(m_mod)
     {
      //Macro element 0, is is immediately below the cylinder
     case 0:
      
      switch(direction)
       {
       case N:
        xi[0] = 4.0*atan(1.0) - atan(1.0)*0.5*(1.0+s[0]);
        Cylinder_pt->position(time,xi,f);
        break;
        
       case S:
        linear_interpolate(Up[0][Nup],Up[1][Nup],s[0],f);
        break;
        
       case W:  
        xi[0] = 4.0*atan(1.0);
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(Up[0][Nup],point_on_circle,s[0],f);
        break;
        
       case E:
        xi[0] = 3.0*atan(1.0);
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(Up[1][Nup],point_on_circle,s[0],f);
        break;
        
       default:
        
        std::ostringstream error_stream;
        error_stream << "Direction is incorrect:  " << direction 
                     << std::endl;
        throw OomphLibError(
         error_stream.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
      
      break;
      
      //Macro element 1, is immediately to the right the cylinder
     case 1:
      
      switch(direction)
       {
       case N:
        xi[0] = atan(1.0);
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(point_on_circle,Down[1][0],s[0],f);
        break;
        
       case S:  
        xi[0] = 3.0*atan(1.0);
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(point_on_circle,Up[1][Nup],s[0],f);
        break;
        
       case W:
        xi[0] = 3.0*atan(1.0) - 2.0*atan(1.0)*0.5*(1.0+s[0]);
        Cylinder_pt->position(time,xi,f);
        break;
        
       case E:
        linear_interpolate(Up[1][Nup],Down[1][0],s[0],f);
        break;
        
       default:
        
        std::ostringstream error_stream;
        error_stream << "Direction is incorrect:  " << direction << std::endl;
        throw OomphLibError(
         error_stream.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
      
      break;
      
      //Macro element 2, is immediately above the cylinder
     case 2:
      
      switch(direction)
       {
       case N:
        linear_interpolate(Down[0][0],Down[1][0],s[0],f);
        break;
        
       case S:
        xi[0] = atan(1.0)*0.5*(1.0+s[0]);
        Cylinder_pt->position(time,xi,f);
        break;
        
       case W:  
        xi[0] = 0.0;
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(point_on_circle,Down[0][0],s[0],f);
        break;
        
       case E:
        xi[0] = atan(1.0);
        Cylinder_pt->position(time,xi,point_on_circle);
        linear_interpolate(point_on_circle,Down[1][0],s[0],f);
        break;
        

     default:

      std::ostringstream error_stream;
      error_stream << "Direction is incorrect:  " << direction << std::endl;
      throw OomphLibError(error_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
       }
      break;
       } 
   } //End of macro elements around the half-cylinder
  else
   {
    //Other cases
    if(m < Nup+Ndown+3)
     {
      unsigned m_mod = m - Nup -3;
      
      //Switch on the direction
      
      switch(direction)
       {
       case N:
        linear_interpolate(Down[0][m_mod+1],Down[1][m_mod+1],s[0],f);
        break;
        
       case S:  
        linear_interpolate(Down[0][m_mod],Down[1][m_mod],s[0],f);
        break;
        
       case W:
        linear_interpolate(Down[0][m_mod],Down[0][m_mod+1],s[0],f);
        break;
        
       case E:
        linear_interpolate(Down[1][m_mod],Down[1][m_mod+1],s[0],f);
        break;
        
       default:
        
        std::ostringstream error_stream;
        error_stream << "Direction is incorrect:  " << direction << std::endl;
        throw OomphLibError(
         error_stream.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }  
     }
    else if(m < Nup+Ndown+3 + (Ncolumn-1)*(Nup+Ndown+1))
     {
      //Work out the modified m
      unsigned m_col = m - (Nup+Ndown+3);
      //Work out which column we are in
      unsigned j_col = 1 + m_col/(Nup+Ndown+1);
      //Work out the actual vertical position
      unsigned m_mod = m_col%(Nup+Ndown+1);
     
      //If we're in the upstream region
      if(m_mod < Nup)
       {
        switch(direction)
         {
         case N:
          linear_interpolate(Up[j_col][m_mod+1],Up[j_col+1][m_mod+1],s[0],f);
          break;
          
         case S:  
          linear_interpolate(Up[j_col][m_mod],Up[j_col+1][m_mod],s[0],f);
          break;
          
         case W:
          linear_interpolate(Up[j_col][m_mod],Up[j_col][m_mod+1],s[0],f);
          break;
          
         case E:
          linear_interpolate(Up[j_col+1][m_mod],Up[j_col+1][m_mod+1],s[0],f);
          break;
          
         default:
          
          std::ostringstream error_stream;
          error_stream << "Direction is incorrect:  " << direction << std::endl;
          throw OomphLibError(
           error_stream.str(),
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }  
       }
      //Otherwise central zone
      else if(m_mod==Nup)
       {
        switch(direction)
         {
         case N:
          linear_interpolate(Down[j_col][0],Down[j_col+1][0],s[0],f);
          break;
          
         case S:  
          linear_interpolate(Up[j_col][Nup],Up[j_col+1][Nup],s[0],f);
          break;
          
         case W:
          linear_interpolate(Up[j_col][Nup],Down[j_col][0],s[0],f);
          break;
          
         case E:
          linear_interpolate(Up[j_col+1][Nup],Down[j_col+1][0],s[0],f);
          break;
          
         default:
          
          std::ostringstream error_stream;
          error_stream << "Direction is incorrect:  " << direction << std::endl;
          throw OomphLibError(
           error_stream.str(),
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }  
       }
      else  if(m_mod < Nup+Ndown+1)
       {
        unsigned m_mod2 = m_mod - Nup -1;
      
        //Switch on the direction
        
        switch(direction)
         {
         case N:
          linear_interpolate(Down[j_col][m_mod2+1],Down[j_col+1][m_mod2+1],s[0],f);
          break;
          
         case S:  
          linear_interpolate(Down[j_col][m_mod2],Down[j_col+1][m_mod2],s[0],f);
          break;
          
         case W:
          linear_interpolate(Down[j_col][m_mod2],Down[j_col][m_mod2+1],s[0],f);
          break;
          
         case E:
          linear_interpolate(Down[j_col+1][m_mod2],Down[j_col+1][m_mod2+1],s[0],f);
          break;
          
         default:
          
          std::ostringstream error_stream;
          error_stream << "Direction is incorrect:  " << direction << std::endl;
          throw OomphLibError(
           error_stream.str(),
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }  
       }
     }
    else
     {
      std::ostringstream error_stream;
      error_stream << "Wrong macro element number" << m << std::endl;
      throw OomphLibError(
       error_stream.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 }

private:
 
 /// Left and right side of lines in the upstream region
 Vector<Vector<Vector<double> > > Up;

 /// Left and right side of lines in the downstream region
 Vector<Vector<Vector<double> > > Down;

 /// Pointer to geometric object that represents the central cylinder
 GeomObject* Cylinder_pt;

 /// Number of upstream macro elements
 unsigned Nup;

 /// Number of downstream macro elements
 unsigned Ndown;

 /// Number of columns
 unsigned Ncolumn;

};




/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////


//=============================================================
/// Domain-based mesh for rectangular mesh with circular hole
//=============================================================
template<class ELEMENT>
class HalfRectangleWithHoleMesh : public virtual Mesh
{

public:


 /// Constructor: Pass pointer to geometric object that 
 /// represents the cylinder, the length and height of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction. Timestepper defaults to Steady
 /// default timestepper.
 HalfRectangleWithHoleMesh(GeomObject* cylinder_pt, 
                           const double &radius, const double &length,
                           const double &up_length, const unsigned &nup,
                           const double &down_length, const unsigned &ndown,
                           const double &width_near_cylinder,
                           const unsigned &ncolumn,
                           TimeStepper* time_stepper_pt=
                           &Mesh::Default_TimeStepper) 
  {
   //Create the domain
   Domain_pt = new HalfRectangleWithHoleDomain(cylinder_pt,radius,length,
                                               up_length,nup,down_length,
                                               ndown,width_near_cylinder,
                                               ncolumn);

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

   //Edges between the upstream elements including that immediately below
   //the cylinder
   for(unsigned m=0;m<(nup+1);m++)
    {
     //Edge between Elements m and m+1
     for(unsigned n=0;n<np;n++)
      {
       //Set the nodes in element m+1 to be the same as in element m
       finite_element_pt(m+1)->node_pt(n)
        = finite_element_pt(m)->node_pt(np*(np-1) + n);
       
       //Remove the nodes in element m+1 from the temporary node list
       delete Tmp_node_pt[(m+1)*np*np + n];
       Tmp_node_pt[(m+1)*np*np + n] = 0;
      }
    }

   //Edge between Elements nup and nup+1
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element nup+1 to be the same as in element nup
     finite_element_pt(nup+1)->node_pt(n)
      = finite_element_pt(nup)->node_pt((np-n-1)*np + (np-1));

     //Remove the nodes in element nup+1 from the temporary node list
     delete Tmp_node_pt[(nup+1)*np*np + n];
     Tmp_node_pt[(nup+1)*np*np + n] = 0;
    }

   //Edge between Element nup+1 and nup2
   for(unsigned n=0;n<np;n++)
    {
     //Set the nodes in element nup+2 to be the same as in element nup+1
     finite_element_pt(nup+2)->node_pt(np*n + np-1)
      = finite_element_pt(nup+1)->node_pt(np*(np-1) + n);

     //Remove the nodes in element nup+2 from the temporary node list
     delete Tmp_node_pt[(nup+2)*np*np + np*n + np-1];
     Tmp_node_pt[(nup+2)*np*np + np*n + np-1] = 0;
    }

   //Edges between the downstream elements including that immediately above
   //the cylinder
   for(unsigned m=nup+2;m<(nup+2+ndown);m++)
    {
     //Edge between Elements m and m+1
     for(unsigned n=0;n<np;n++)
      {
       //Set the nodes in element m+1 to be the same as in element m
       finite_element_pt(m+1)->node_pt(n)
        = finite_element_pt(m)->node_pt(np*(np-1) + n);
       
       //Remove the nodes in element m+1 from the temporary node list
       delete Tmp_node_pt[(m+1)*np*np + n];
       Tmp_node_pt[(m+1)*np*np + n] = 0;
      }
    }


   //Now we need to sort out all the edges between the remaining columns 
   //and rows
   bool first_col=true;

   unsigned left_col_offset = 0;
   unsigned right_col_offset = nup+ndown+3;
   //Loop over the columns
   for(unsigned j=1;j<ncolumn;j++)
    {
     for(unsigned i=0;i<(nup+ndown+1);i++)
      {
       //Sort out the left-hand boundary of elements
       for(unsigned n=0;n<np;n++)
        {
         finite_element_pt(right_col_offset)->node_pt(n*np)
          = finite_element_pt(left_col_offset)->node_pt(n*np + np-1);
         
         //Remove the nodes in the element from the temporary node list
         delete Tmp_node_pt[right_col_offset*np*np + n*np];
         Tmp_node_pt[right_col_offset*np*np + n*np] = 0;
        }
       
       
       //Sort out the bottom of the element
       if(i!=(nup+ndown))
        {
         for(unsigned n=1;n<np;n++)
          {
           finite_element_pt(right_col_offset+1)->node_pt(n) 
            = finite_element_pt(right_col_offset)->node_pt((np-1)*np + n);
           
           //Remove the nodes in the element from the temporary node list
           delete Tmp_node_pt[(right_col_offset+1)*np*np + n];
           Tmp_node_pt[(right_col_offset+1)*np*np + n] = 0;
          }
        }

       ++right_col_offset; ++left_col_offset;
       //Add another offset if the first column
       if(first_col==true)
        {
         if((i==nup-1) || (i==nup)) {++left_col_offset;}
        }
      }
     if(j==1) {first_col=false;}
    }
        

   //Now set the actual true nodes
   for(unsigned long n=0;n<node_count;n++)
    {
     if(Tmp_node_pt[n]!=0) {Node_pt.push_back(Tmp_node_pt[n]);}
    }

   //Finally set the nodes on the boundaries
   set_nboundary(5);

   //Find the offset for the right hand side
   unsigned rhs_offset=0;
   if(ncolumn > 1) {rhs_offset = nup+ndown+3 + (ncolumn-2)*(nup+ndown+1);}
   
   for(unsigned n=0;n<np;n++)
    {
     //First part of left hand side
     Node* nod_pt=finite_element_pt(0)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(3,nod_pt);

     //Part of LHS immediately after hole
     nod_pt=finite_element_pt(nup+2)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(3,nod_pt);
     
     //First part of right hand side
     nod_pt=finite_element_pt(rhs_offset)->node_pt(n*np + np-1);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(1,nod_pt);
     
     //First part of lower boundary
     nod_pt=finite_element_pt(0)->node_pt(n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(0,nod_pt);
     
     //First part of upper boundary
     nod_pt=finite_element_pt(nup+ndown+2)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(2,nod_pt);
     
     //First part of hole boundary
     nod_pt=finite_element_pt(nup)->node_pt(np*(np-1)+n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }

   //Upstream section
   for(unsigned n=1;n<np;n++)
    {
     for(unsigned m=1;m<nup;m++)
      {
       //Next part of left hand side
       Node* nod_pt=finite_element_pt(m)->node_pt(n*np);
       convert_to_boundary_node(nod_pt);
       add_boundary_node(3,nod_pt);
       
       //Next part of right hand side
       nod_pt=finite_element_pt(rhs_offset + m)->node_pt(n*np + np-1);
       convert_to_boundary_node(nod_pt);
       add_boundary_node(1,nod_pt);
      }
     
     //Next part of hole
     Node* nod_pt=finite_element_pt(nup+1)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }
   
   //Need to take care with the region around the hole in the special
   //case of a single column
   unsigned rhs_extra_offset=1;
   if(ncolumn>1) {rhs_extra_offset=0;}
   for(unsigned n=1;n<np;n++)
    {
     //Next two parts of left boundary
     Node* nod_pt=finite_element_pt(nup)->node_pt(n*np);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(3,nod_pt);

     // Next part of right boundary
     nod_pt=finite_element_pt(rhs_offset + rhs_extra_offset + nup)
      ->node_pt(n*np + np-1);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(1,nod_pt);
     
     //Final part of hole boundary                                
     nod_pt=finite_element_pt(nup+2)->node_pt(np-1-n);
     convert_to_boundary_node(nod_pt);
     add_boundary_node(4,nod_pt);
    }

   //Downstream section
   for(unsigned n=1;n<np;n++)
    {
     for(unsigned m=nup+3;m<(nup+ndown+3);m++)
      {
       //Next part of left hand side
       Node* nod_pt=finite_element_pt(m)->node_pt(n*np);
       convert_to_boundary_node(nod_pt);
       add_boundary_node(3,nod_pt);
       
       //Next part of right hand side
       if(rhs_offset==0)
        {
         nod_pt=finite_element_pt(m)->node_pt(n*np + np-1);
         convert_to_boundary_node(nod_pt);
         add_boundary_node(1,nod_pt);
        }
      }
     if(rhs_offset!=0)
      {
       for(unsigned m=nup+1;m<(nup+ndown+1);m++)
        {
         Node* nod_pt=finite_element_pt(rhs_offset + m)->node_pt(n*np + np-1);
         convert_to_boundary_node(nod_pt);
         add_boundary_node(1,nod_pt);
        }
      }
    }

   //Upper and Lower boundaries
   unsigned lower_index = nup + ndown + 3;
   for(unsigned j=1;j<ncolumn;j++)
    {
     for(unsigned n=1;n<np;n++)
      {
       Node* nod_pt = finite_element_pt(lower_index)->node_pt(n);
       convert_to_boundary_node(nod_pt);
       add_boundary_node(0,nod_pt);

       nod_pt = 
        finite_element_pt(lower_index+nup+ndown)->node_pt(np*(np-1) + n);
       convert_to_boundary_node(nod_pt);
       add_boundary_node(2,nod_pt);
      }
     lower_index += nup+ndown+1;
    }


   //Setup the info
   this->setup_boundary_element_info();

  }

 /// Access function to the domain
 HalfRectangleWithHoleDomain* domain_pt() {return Domain_pt;}

protected:

 /// Pointer to the domain
 HalfRectangleWithHoleDomain* Domain_pt;

};


/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//===================================================================
/// Refineable version of RectangleWithHoleMesh. Applies one uniform
/// refinement immediately to avoid problems with the automatic
/// applications of boundary conditions in subsequent refinements
//===================================================================
template<class ELEMENT>
class RefineableHalfRectangleWithHoleMesh :
 public HalfRectangleWithHoleMesh<ELEMENT>, public RefineableQuadMesh<ELEMENT>
{
public: 

 /// Constructor. Pass pointer to geometric object that 
 /// represents the cylinder, the length and height of the domain.
 /// The GeomObject must be parametrised such that 
 /// \f$\zeta \in [0,2\pi]\f$ sweeps around the circumference
 /// in anticlockwise direction. Timestepper defaults to Steady
 /// default timestepper.
 RefineableHalfRectangleWithHoleMesh(GeomObject* cylinder_pt, 
                                     const double &radius, 
                                     const double &length,
                                     const double &up_length, 
                                     const unsigned &nup,
                                     const double &down_length, 
                                     const unsigned &ndown,
                                     const double &width_near_cylinder,
                                     const unsigned &ncolumn,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  HalfRectangleWithHoleMesh<ELEMENT>(cylinder_pt,radius,length,up_length,nup,
                                     down_length,ndown,
                                     width_near_cylinder,
                                     ncolumn,time_stepper_pt) 
  {

   // Nodal positions etc. were created in constructor for
   // Cylinder...<...>. Need to setup adaptive information.
  
   // Loop over all elements and set macro element pointer
   for (unsigned e=0;e<(ndown+nup+3);e++)
    {
     dynamic_cast<ELEMENT*>(this->element_pt(e))->
      set_macro_elem_pt(this->Domain_pt->macro_element_pt(e));
    }
   
   // Setup quadtree forest for mesh refinement
   this->setup_quadtree_forest();
   
   // Do one round of uniform refinement to avoid problems
   // with automatic application of boundary conditions
   // in subsequent refinements
   this->refine_uniformly();

   // Automatically setup the boundary element lookup scheme
   this->setup_boundary_element_info();
  }
 
 
 /// Destructor: Empty
 virtual ~RefineableHalfRectangleWithHoleMesh() {}
 
};

}
#endif
