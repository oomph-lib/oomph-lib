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

#include<fenv.h>

//Generic routines
#include "generic.h"

using namespace std;

using namespace oomph;


//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{


 // Warped disk with specified amplitude and wavenumber for warping
 double epsilon=0.1;
 unsigned n=5;
 WarpedCircularDisk* disk_pt=new WarpedCircularDisk(epsilon,n);


 Vector<double> r(3);
 Vector<double> zeta(2);
 ofstream disk_file;
 disk_file.open("disk.dat");
 unsigned nr=50;
 unsigned nphi=300;
 disk_file << "ZONE I=" << nphi << ", J=" << nr << std::endl;
 for (unsigned i=0;i<nr;i++)
  {
   double radius=double(i)/double(nr-1);
   for (unsigned j=0;j<nphi;j++)
    {
     double phi=double(j)/double(nphi-1)*2.0*MathematicalConstants::Pi;
     zeta[0]=radius*cos(phi);
     zeta[1]=radius*sin(phi);
     disk_pt->position(zeta,r);
     disk_file << r[0] << " " 
               << r[1] << " " 
               << r[2] << " " 
               << zeta[0] <<" " 
               << zeta[1] <<" " 
               << std::endl;
    }
  }
 disk_file.close();


 ofstream two_d_boundaries_file("two_d_boundaries.dat");
 ofstream three_d_boundaries_file("three_d_boundaries.dat");
 ofstream boundaries_tangent_file("boundaries_tangent.dat");
 ofstream boundaries_normal_file("boundaries_normal.dat");
 ofstream boundaries_binormal_file("boundaries_binormal.dat");
 unsigned nplot=100;
 disk_pt->output_boundaries_and_triads(nplot,
                                       two_d_boundaries_file,
                                       three_d_boundaries_file,
                                       boundaries_tangent_file,
                                       boundaries_normal_file,
                                       boundaries_binormal_file);
 two_d_boundaries_file.close();
 three_d_boundaries_file.close();
 boundaries_tangent_file.close();
 boundaries_normal_file.close();
 boundaries_binormal_file.close();

}
