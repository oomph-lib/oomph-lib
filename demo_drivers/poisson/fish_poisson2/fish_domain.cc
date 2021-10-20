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
/// Driver code for mucking around with GeomObjects, MacroElements
/// and Domains

// Generic oomph-lib headers
#include "generic.h"

// The fish domain
#include "meshes/fish_domain.h"

using namespace std;

using namespace oomph;

//=======================start_of_unit_circle==============================
/// Unit circle in 2D, centred at the origin, parametrised by a single
/// Lagrangian coordinate, the polar angle.
//=========================================================================
class UnitCircle : public GeomObject
{

 public:

  ///  Constructor: Pass the number of Lagrangian
  /// and Eulerian coordinates to the constructor of the
  ///  GeomObject base class.
  UnitCircle() : GeomObject(1,2) {}

  /// Destructor -- emtpy
  virtual ~UnitCircle(){}

  ///  Position vector, r, to the point on the circle identified by  
  /// its 1D Lagrangian coordinate, xi (passed as a 1D Vector):
  void position(const Vector<double>& xi, Vector<double>& r) const
  {
   // Eulerian position vector
   r[0] = cos(xi[0]);
   r[1] = sin(xi[0]);
  }


  ///  Position vector, r, to the point on the circle identified by  
  /// its 1D Lagrangian coordinate, xi (passed as a 1D Vector) at discrete time
  /// level t (t=0: present; t>0: previous). The shape of the object 
  /// is not time-dependent, therefore we forward this call to the 
  /// steady version. 
  void position(const unsigned& t, const Vector<double>& xi, 
                Vector<double>& r) const
  {
   position(xi,r);
  }

}; // end of unit circle class






//===================================================================
/// Driver code for mucking around with GeomObjects, MacroElements
/// and domains
//===================================================================
int main()
{


 // Play around with a GeomObject
 //------------------------------

 // Create a unit circle
 UnitCircle unit_circle;

 // Plot:
 ofstream circle_file("unit_circle.dat");

 // Number of plot points
 unsigned nplot=50;

 // 1D vector for the Lagrangian coordinate
 Vector<double> s(1);

 // 2D vector for the Eulerian position
 Vector<double> r(2);

 for (unsigned i=0;i<nplot;i++)
  {
   // Lagrangian coordinate at plot point
   s[0]=2.0*MathematicalConstants::Pi*double(i)/double(nplot-1);

   // Get Eulerian position vector from GeomObject:
   unit_circle.position(s,r);

   // Plot
   circle_file << r[0] << " " << r[1] << std::endl;
  }
   
   
 // Close output file
 circle_file.close();





 // Build a FishDomain and plot it
 //-------------------------------


 // Create the fish's back as a circle of radius 1, centred at (0.5,0.0)
 double x_c=0.5;
 double y_c=0.0;
 double r_back=1.0;
 GeomObject* back_pt=new Circle(x_c,y_c,r_back);


 // Start and end coordinates of the fish back
 double s_nose=2.6; 
 double s_tail=0.4; 

 // Create the domain
 Domain* domain_pt=new FishDomain(back_pt,s_nose,s_tail);
 

 // Plot the domain

 // Number of plot points in each coordinate direction.
 unsigned npts=10;
 
 // Output domain (= plot its macro elements) and the macro element boundaries
 ofstream domain_file("fish_domain.dat");
 domain_pt->output(domain_file,npts);
 domain_pt->output_macro_element_boundaries(domain_file,npts);
 domain_file.close();


}
