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
// Doc sparse MacroElement-based node update

 
// Generic oomph-lib headers
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The fish mesh 
#include "meshes/fish_mesh.h"

// Circle as generalised element:
#include "circle_as_generalised_element.h"

using namespace std;

using namespace oomph;

//==================start_of_main=========================================
/// Driver to document sparse MacroElement-based node update.
//========================================================================
int main()
{

 // Shorthand for element type
 typedef MacroElementNodeUpdateElement<RefineableQPoissonElement<2,3> > 
  ELEMENT;

 // Set coordinates and radius for the circle that will become the fish back
 double x_c=0.5;
 double y_c=-0.2;
 double r_back=1.0;

 // Build geometric object that will become the fish back
  ElasticallySupportedRingElement* Fish_back_pt=
   new ElasticallySupportedRingElement(x_c,y_c,r_back);

 // Build fish mesh with geometric object that specifies the fish back 
 MacroElementNodeUpdateRefineableFishMesh<ELEMENT>* Fish_mesh_pt=new 
  MacroElementNodeUpdateRefineableFishMesh<ELEMENT>(Fish_back_pt);



 // Number of plot points in each coordinate direction.
 unsigned npts=11; 

 ofstream some_file;
 char filename[100];

 // Output initial mesh
 unsigned count=0;
 sprintf(filename,"RESLT/soln%i.dat",count);
 some_file.open(filename);
 Fish_mesh_pt->output(some_file,npts);
 some_file.close();
 count++; 



 // Increment y_c
 Fish_back_pt->y_c()+=0.2;



 // Adjust each node in turn and doc
 unsigned nnod=Fish_mesh_pt->nnode();
 for (unsigned i=0;i<nnod;i++)
  {
   // Update individual nodal position
   Fish_mesh_pt->node_pt(i)->node_update();

   // Doc mesh
   sprintf(filename,"RESLT/soln%i.dat",count);
   some_file.open(filename);
   Fish_mesh_pt->output(some_file,npts);
   some_file.close();
   count++; 
  }

} // end of main


