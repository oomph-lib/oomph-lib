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

// Generic oomph-lib headers
#include "generic.h"

// Circle header
#include "circle.h"

using namespace std;

using namespace oomph;

//========================================================================
/// Driver
//========================================================================
int main()
{

 // X-coordinate of the circle's centre
 double x_c=0.5;

 // Y-coordinate of the circle's centre
 double y_c=1.5;

 // Radius
 double R=0.9;

 // Build circle object with specified (const) parameters 
 //-------------------------------------------------------
 GeneralCircle circle0(x_c,y_c,R);


 // Build circle object with Data -- the Data values might be determine
 //--------------------------------------------------------------------
 // "somewhere else", e.g. as part of the solution of another problem 
 //------------------------------------------------------------------

 // The circle's shape is determine by a single Data object whose
 // three values specify x_c, y_c and R:
 Data* circle_data_pt=new Data(3);

 // Set the values
 circle_data_pt->set_value(0,x_c);
 circle_data_pt->set_value(1,y_c);
 circle_data_pt->set_value(2,R);

 // Build the object
 GeneralCircle circle1(circle_data_pt);

 // Number of plot points
 unsigned npts=100;

 // Lagrangian coordinate and position vector (both as vectors)
 Vector<double> xi(1);
 Vector<double> r(2);
 
 // Output circles
 ofstream some_file0, some_file1;
 some_file0.open("circle0.dat");
 some_file1.open("circle1.dat");

 for (unsigned i=0;i<npts;i++)
  {
   xi[0]=2.0*MathematicalConstants::Pi*double(i)/double(npts-1);
   circle0.position(xi,r);  
   some_file0 << r[0] << " " << r[1] << std::endl;
   circle1.position(xi,r);  
   some_file1 << r[0] << " " << r[1] << std::endl;
  }
 some_file0.close();
 some_file1.close();


}


