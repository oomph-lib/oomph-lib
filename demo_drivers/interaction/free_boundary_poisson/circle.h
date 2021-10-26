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
#ifndef CIRCLE_HEADER
#define CIRCLE_HEADER


// Generic oomph-lib headers
#include "generic.h"

namespace oomph
{

//=====start_of_circle=====================================================
/// Simple circle in 2D space.
/// \f[ x = X_c + R \cos(\zeta)  \f]
/// \f[ y = Y_c + R \sin(\zeta)  \f]
//=========================================================================
class SimpleCircle : public GeomObject
{

public:

 /// Constructor:  Pass x and y-coords of centre and radius
 SimpleCircle(const double& x_c, const double& y_c, 
               const double& r) : GeomObject(1,2), X_c(x_c), Y_c(y_c), R(r)
  {}

 /// Position Vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position vector
   r[0] = X_c+R*cos(zeta[0]);
   r[1] = Y_c+R*sin(zeta[0]);
  }

 /// Position Vector at Lagrangian coordinate zeta  at time level t
 /// (t=0: present; t>0: previous level). Steady object, so we 
 /// simply forward the call to the steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {position(zeta,r);}

protected:

 /// X-coordinate of centre
 double X_c;

 /// Y-coordinate of centre
 double Y_c;

 /// Radius
 double R;

}; // end of SimpleCircle class




/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=========================================================================
/// GeneralCircle in 2D space.
/// \f[ x = X_c + R \cos(\zeta)  \f]
/// \f[ y = Y_c + R \sin(\zeta)  \f]
/// The three parameters \f$ X_c, Y_c \f$ and \f$ R \f$ are represented
/// by Data and can therefore be unknowns in the problem.
//=========================================================================
class GeneralCircle : public GeomObject
{

public:

 /// Constructor:  Pass x and y-coords of centre and radius (all pinned)
 GeneralCircle(const double& x_c, const double& y_c, 
               const double& r) : GeomObject(1,2)
  {
   // Create Data: We have one Data item with 3 values. The Data object
   // stores the x position of the circle's centre as value 0,
   // the y position of the circle's centre as value 1, and the 
   // circle's radius as value 2.
   Geom_data_pt.resize(1);
   Geom_data_pt[0] = new Data(3);
   
   // Assign data: X_c -- the value is free by default: Need to pin it.

   // Pin the data
   Geom_data_pt[0]->pin(0); 
   // Give it a value: 
   Geom_data_pt[0]->set_value(0,x_c);

   // Assign data: Y_c -- the value is free by default: Need to pin it.

   // Pin the data
   Geom_data_pt[0]->pin(1); 
   // Give it a value: 
   Geom_data_pt[0]->set_value(1,y_c);

   // Assign data: R -- the value is free by default: Need to pin it.

   // Pin the data
   Geom_data_pt[0]->pin(2); 
   // Give it a value: 
   Geom_data_pt[0]->set_value(2,r);

   // I've created the data, I need to clean up
   Must_clean_up=true; 
  }

 /// Alternative constructor:  Pass x and y-coords of centre and radius
 /// (all as part of Data)
 /// \code 
 /// Geom_data_pt[0]->value(0) = X_c;
 /// Geom_data_pt[0]->value(1) = Y_c;
 /// Geom_data_pt[0]->value(2) = R;
 /// \endcode
 GeneralCircle(Data* geom_data_pt) : GeomObject(1,2)
  {
#ifdef PARANOID
   if (geom_data_pt->nvalue()!=3)
    {
     std::ostringstream error_stream;
     error_stream << "Geometric Data must have 3 values, not "
                  << geom_data_pt->nvalue() << std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   Geom_data_pt.resize(1);
   Geom_data_pt[0]=geom_data_pt;

   // Data has been created externally: Must not clean up
   Must_clean_up=false;

  } //end of alternative constructor


 /// Destructor:  Clean up if necessary
 virtual ~GeneralCircle()
  {
   // Do I need to clean up?
   if (Must_clean_up)
    {
     unsigned ngeom_data=Geom_data_pt.size();
     for (unsigned i=0;i<ngeom_data;i++)
      {
       delete Geom_data_pt[i];
       Geom_data_pt[i]=0;
      }
    }
  } // end of destructor


 /// Position Vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Extract data
   double X_c= Geom_data_pt[0]->value(0);
   double Y_c= Geom_data_pt[0]->value(1);
   double R= Geom_data_pt[0]->value(2);

   // Position vector
   r[0] = X_c+R*cos(zeta[0]);
   r[1] = Y_c+R*sin(zeta[0]);

  } // end of position(...)


 /// Position Vector at Lagrangian coordinate zeta  at time level t
 /// (t=0: present; t>0: previous level). Steady object, so we 
 /// simply forward the call to the steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {position(zeta,r);}

 /// Access function to x-coordinate of centre of circle
 double& x_c(){return *Geom_data_pt[0]->value_pt(0);}

 /// Access function to y-coordinate of centre of circle
 double& y_c(){return *Geom_data_pt[0]->value_pt(1);}

 /// Access function to radius of circle
 double& R(){return *Geom_data_pt[0]->value_pt(2);}

 /// How many items of Data does the shape of the object depend on?
 unsigned ngeom_data() const {return Geom_data_pt.size();}
 
 /// Return pointer to the j-th Data item that the object's 
 /// shape depends on 
 Data* geom_data_pt(const unsigned& j) {return Geom_data_pt[j];}
 

protected:

 /// Vector of pointers to Data items that affects the object's shape
 Vector<Data*> Geom_data_pt;

 /// Do I need to clean up?
 bool Must_clean_up;

};

}

#endif
