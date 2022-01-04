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

//Generic routines
#include "generic.h"

// Navier--Stokes equations
#include "navier_stokes.h"

// Get the mesh
#include "oomph_tetgen.h"
#include "meshes/tetgen_mesh.h"
#include "meshes/refineable_tetgen_mesh.h"

using namespace std;

using namespace oomph;

namespace Global_Parameters
{
 double Re = 0.0;
 double Visc_Ratio = 2.0;
 double Box_width = 3.0; 
 double Box_length = 10.0;
}


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////


//=============================================================
/// TetMeshFacetedSurface that defines inner boundary
//=============================================================
class SphericalTetMeshFacetedSurface : public virtual 
TetMeshFacetedClosedSurface
{


public: 

 /// Constructor
 SphericalTetMeshFacetedSurface()
  {
   
   //Golden ratio
   const double phi = 0.5*(1.0 + sqrt(5.0));
   
   // Make vertices
   unsigned n_vertex=12;
   Vertex_pt.resize(n_vertex);

   //Set basic icosahedron points
   Vector<double> icosa_point(3);

   icosa_point[0] = 0.0;
   icosa_point[1] = 1.0;
   icosa_point[2] = phi;   
   Vertex_pt[0]=new TetMeshVertex(icosa_point);

   icosa_point[0] = 0.0;
   icosa_point[1] = -1.0;
   icosa_point[2] = phi;
   Vertex_pt[1]=new TetMeshVertex(icosa_point);

   icosa_point[0] = 0.0;
   icosa_point[1] = 1.0;
   icosa_point[2] = -phi;
   Vertex_pt[2]=new TetMeshVertex(icosa_point);

   icosa_point[0] = 0.0;
   icosa_point[1] = -1.0;
   icosa_point[2] = -phi;
   Vertex_pt[3]=new TetMeshVertex(icosa_point);

   icosa_point[0] = 1.0;
   icosa_point[1] = phi;
   icosa_point[2] = 0.0;
   Vertex_pt[4]=new TetMeshVertex(icosa_point);

   icosa_point[0] = -1.0;
   icosa_point[1] = phi;
   icosa_point[2] = 0.0;
   Vertex_pt[5]=new TetMeshVertex(icosa_point);

   icosa_point[0] = 1.0;
   icosa_point[1] = -phi;
   icosa_point[2] = 0.0;
   Vertex_pt[6]=new TetMeshVertex(icosa_point);

   icosa_point[0] = -1.0;
   icosa_point[1] = -phi;
   icosa_point[2] = 0.0;
   Vertex_pt[7]=new TetMeshVertex(icosa_point);

   icosa_point[0] = phi;
   icosa_point[1] = 0.0;
   icosa_point[2] = 1.0;
   Vertex_pt[8]=new TetMeshVertex(icosa_point);

   icosa_point[0] = phi;
   icosa_point[1] = 0.0;
   icosa_point[2] = -1.0;
   Vertex_pt[9]=new TetMeshVertex(icosa_point);

   icosa_point[0] = -phi;
   icosa_point[1] = 0.0;
   icosa_point[2] = 1.0;
   Vertex_pt[10]=new TetMeshVertex(icosa_point);

   icosa_point[0] = -phi;
   icosa_point[1] = 0.0;
   icosa_point[2] = -1.0;
   Vertex_pt[11]=new TetMeshVertex(icosa_point);

   // Make facets
   unsigned n_facet=20;
   Facet_pt.resize(n_facet);
   
   unsigned n_vertex_on_facet=3;
   Facet_pt[0]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[0]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[0]->set_vertex_pt(1,Vertex_pt[1]);
   Facet_pt[0]->set_vertex_pt(2,Vertex_pt[8]);
   
   // icosa_facet[0][0] = 0;
   // icosa_facet[0][1] = 1;
   // icosa_facet[0][2] = 8;
   
   Facet_pt[1]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[1]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[1]->set_vertex_pt(1,Vertex_pt[10]);
   Facet_pt[1]->set_vertex_pt(2,Vertex_pt[1]);
   
   // icosa_facet[1].resize(3);
   // icosa_facet[1][0] = 0;
   // icosa_facet[1][1] = 10;
   // icosa_facet[1][2] = 1;
   
   Facet_pt[2]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[2]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[2]->set_vertex_pt(1,Vertex_pt[5]);
   Facet_pt[2]->set_vertex_pt(2,Vertex_pt[10]);
   
   // icosa_facet[2].resize(3);
   // icosa_facet[2][0] = 0;
   // icosa_facet[2][1] = 5;
   // icosa_facet[2][2] = 10;
   
   Facet_pt[3]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[3]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[3]->set_vertex_pt(1,Vertex_pt[4]);
   Facet_pt[3]->set_vertex_pt(2,Vertex_pt[5]);
   
   // icosa_facet[3].resize(3);
   // icosa_facet[3][0] = 0;
   // icosa_facet[3][1] = 4;
   // icosa_facet[3][2] = 5;
   
   Facet_pt[4]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[4]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[4]->set_vertex_pt(1,Vertex_pt[8]);
   Facet_pt[4]->set_vertex_pt(2,Vertex_pt[4]);
   
   // icosa_facet[4].resize(3);
   // icosa_facet[4][0] = 0;
   // icosa_facet[4][1] = 8;
   // icosa_facet[4][2] = 4;
   
   Facet_pt[5]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[5]->set_vertex_pt(0,Vertex_pt[5]);
   Facet_pt[5]->set_vertex_pt(1,Vertex_pt[11]);
   Facet_pt[5]->set_vertex_pt(2,Vertex_pt[10]);

   // icosa_facet[5].resize(3);
   // icosa_facet[5][0] = 5;
   // icosa_facet[5][1] = 11;
   // icosa_facet[5][2] = 10;
   
   Facet_pt[6]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[6]->set_vertex_pt(0,Vertex_pt[5]);
   Facet_pt[6]->set_vertex_pt(1,Vertex_pt[2]);
   Facet_pt[6]->set_vertex_pt(2,Vertex_pt[11]);
   
   // icosa_facet[6].resize(3);
   // icosa_facet[6][0] = 5;
   // icosa_facet[6][1] = 2;
   // icosa_facet[6][2] = 11;
   
   Facet_pt[7]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[7]->set_vertex_pt(0,Vertex_pt[4]);
   Facet_pt[7]->set_vertex_pt(1,Vertex_pt[2]);
   Facet_pt[7]->set_vertex_pt(2,Vertex_pt[5]);
   
   // icosa_facet[7].resize(3);
   // icosa_facet[7][0] = 4;
   // icosa_facet[7][1] = 2;
   // icosa_facet[7][2] = 5;
   
   Facet_pt[8]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[8]->set_vertex_pt(0,Vertex_pt[4]);
   Facet_pt[8]->set_vertex_pt(1,Vertex_pt[9]);
   Facet_pt[8]->set_vertex_pt(2,Vertex_pt[2]);

   // icosa_facet[8].resize(3);
   // icosa_facet[8][0] = 4;
   // icosa_facet[8][1] = 9;
   // icosa_facet[8][2] = 2;
   
   Facet_pt[9]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[9]->set_vertex_pt(0,Vertex_pt[8]);
   Facet_pt[9]->set_vertex_pt(1,Vertex_pt[9]);
   Facet_pt[9]->set_vertex_pt(2,Vertex_pt[4]);
   
   // icosa_facet[9].resize(3);
   // icosa_facet[9][0] = 8;
   // icosa_facet[9][1] = 9;
   // icosa_facet[9][2] = 4;
   
   Facet_pt[10]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[10]->set_vertex_pt(0,Vertex_pt[6]);
   Facet_pt[10]->set_vertex_pt(1,Vertex_pt[9]);
   Facet_pt[10]->set_vertex_pt(2,Vertex_pt[8]);
   
   // icosa_facet[10].resize(3);
   // icosa_facet[10][0] = 6;
   // icosa_facet[10][1] = 9;
   // icosa_facet[10][2] = 8;
   
   Facet_pt[11]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[11]->set_vertex_pt(0,Vertex_pt[1]);
   Facet_pt[11]->set_vertex_pt(1,Vertex_pt[6]);
   Facet_pt[11]->set_vertex_pt(2,Vertex_pt[8]);
   
   // icosa_facet[11].resize(3);
   // icosa_facet[11][0] = 1;
   // icosa_facet[11][1] = 6;
   // icosa_facet[11][2] = 8;
   
   Facet_pt[12]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[12]->set_vertex_pt(0,Vertex_pt[1]);
   Facet_pt[12]->set_vertex_pt(1,Vertex_pt[7]);
   Facet_pt[12]->set_vertex_pt(2,Vertex_pt[6]);
   
   // icosa_facet[12].resize(3);
   // icosa_facet[12][0] = 1;
   // icosa_facet[12][1] = 7;
   // icosa_facet[12][2] = 6;
   
   Facet_pt[13]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[13]->set_vertex_pt(0,Vertex_pt[10]);
   Facet_pt[13]->set_vertex_pt(1,Vertex_pt[7]);
   Facet_pt[13]->set_vertex_pt(2,Vertex_pt[1]);
   
   // icosa_facet[13].resize(3);
   // icosa_facet[13][0] = 10;
   // icosa_facet[13][1] = 7;
   // icosa_facet[13][2] = 1;
   
   Facet_pt[14]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[14]->set_vertex_pt(0,Vertex_pt[10]);
   Facet_pt[14]->set_vertex_pt(1,Vertex_pt[11]);
   Facet_pt[14]->set_vertex_pt(2,Vertex_pt[7]);
   
   // icosa_facet[14].resize(3);
   // icosa_facet[14][0] = 10;
   // icosa_facet[14][1] = 11;
   // icosa_facet[14][2] = 7;
   
   Facet_pt[15]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[15]->set_vertex_pt(0,Vertex_pt[11]);
   Facet_pt[15]->set_vertex_pt(1,Vertex_pt[3]);
   Facet_pt[15]->set_vertex_pt(2,Vertex_pt[7]);
   
   // icosa_facet[15].resize(3);
   // icosa_facet[15][0] = 11;
   // icosa_facet[15][1] = 3;
   // icosa_facet[15][2] = 7;
   
   Facet_pt[16]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[16]->set_vertex_pt(0,Vertex_pt[7]);
   Facet_pt[16]->set_vertex_pt(1,Vertex_pt[3]);
   Facet_pt[16]->set_vertex_pt(2,Vertex_pt[6]);
   
   // icosa_facet[16].resize(3);
   // icosa_facet[16][0] = 7;
   // icosa_facet[16][1] = 3;
   // icosa_facet[16][2] = 6;
   
   Facet_pt[17]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[17]->set_vertex_pt(0,Vertex_pt[6]);
   Facet_pt[17]->set_vertex_pt(1,Vertex_pt[3]);
   Facet_pt[17]->set_vertex_pt(2,Vertex_pt[9]);
   
   // icosa_facet[17].resize(3);
   // icosa_facet[17][0] = 6;
   // icosa_facet[17][1] = 3;
   // icosa_facet[17][2] = 9;
   
   Facet_pt[18]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[18]->set_vertex_pt(0,Vertex_pt[9]);
   Facet_pt[18]->set_vertex_pt(1,Vertex_pt[3]);
   Facet_pt[18]->set_vertex_pt(2,Vertex_pt[2]);
   
   // icosa_facet[18].resize(3);
   // icosa_facet[18][0] = 9;
   // icosa_facet[18][1] = 3;
   // icosa_facet[18][2] = 2;
   
   Facet_pt[19]=new TetMeshFacet(n_vertex_on_facet);
   Facet_pt[19]->set_vertex_pt(0,Vertex_pt[2]);
   Facet_pt[19]->set_vertex_pt(1,Vertex_pt[3]);
   Facet_pt[19]->set_vertex_pt(2,Vertex_pt[11]);
   
   // icosa_facet[19].resize(3);
   // icosa_facet[19][0] = 2;
   // icosa_facet[19][1] = 3;
   // icosa_facet[19][2] = 11;
   
   // Set one-based boundary IDs
   unsigned one_based_boundary_id=1;
   for (unsigned f=0;f<n_facet;f++)
    {
     Facet_pt[f]->set_one_based_boundary_id(one_based_boundary_id);
    }
   
   // Identify point in hole
   Vector<double> inner_point(3,0.0);
   set_hole_for_tetgen(inner_point);
   
  }

};


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////


//=============================================================
/// TetMeshFacetedSurface that defines outer boundary
//=============================================================
class CubicTetMeshFacetedSurface : public virtual TetMeshFacetedClosedSurface
{


public: 

 /// Constructor
 CubicTetMeshFacetedSurface()
  {
   const double box_width = Global_Parameters::Box_width;
   const double box_length = Global_Parameters::Box_length;

   // Make vertices
   unsigned n_vertex=8;
   Vertex_pt.resize(n_vertex);
   Vector<double> box_point(3);

   box_point[0] = -box_width;
   box_point[1] = -box_width;
   box_point[2] = -box_length;
   Vertex_pt[0]=new TetMeshVertex(box_point);
  
   box_point[0] = -box_width;
   box_point[1] =  box_width;
   box_point[2] = -box_length;
   Vertex_pt[1]=new TetMeshVertex(box_point);

   box_point[0] = -box_width;
   box_point[1] =  box_width;
   box_point[2] =  box_length;
   Vertex_pt[2]=new TetMeshVertex(box_point);
   
   box_point[0] = -box_width;
   box_point[1] = -box_width;
   box_point[2] =  box_length;
   Vertex_pt[3]=new TetMeshVertex(box_point);

   box_point[0] =  box_width;
   box_point[1] = -box_width;
   box_point[2] = -box_length;
   Vertex_pt[4]=new TetMeshVertex(box_point);
   
   box_point[0] =  box_width;
   box_point[1] =  box_width;
   box_point[2] = -box_length;
   Vertex_pt[5]=new TetMeshVertex(box_point);
   
   box_point[0] =  box_width;
   box_point[1] =  box_width;
   box_point[2] =  box_length;
   Vertex_pt[6]=new TetMeshVertex(box_point);
   
   box_point[0] =  box_width;
   box_point[1] =  -box_width;
   box_point[2] =  box_length;
   Vertex_pt[7]=new TetMeshVertex(box_point);
   

   // Make facets
   unsigned n_facet=6;
   Facet_pt.resize(n_facet);
   
   unsigned n_vertex_on_facet=4;
   Facet_pt[0]=new TetMeshFacet(n_vertex_on_facet);
   unsigned one_based_boundary_id=2;
   Facet_pt[0]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[0]->set_vertex_pt(0,Vertex_pt[0]);
   Facet_pt[0]->set_vertex_pt(1,Vertex_pt[4]);
   Facet_pt[0]->set_vertex_pt(2,Vertex_pt[7]);
   Facet_pt[0]->set_vertex_pt(3,Vertex_pt[3]);

   Facet_pt[1]=new TetMeshFacet(n_vertex_on_facet);
   one_based_boundary_id=3;
   Facet_pt[1]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[1]->set_vertex_pt(0,Vertex_pt[4]);
   Facet_pt[1]->set_vertex_pt(1,Vertex_pt[5]);
   Facet_pt[1]->set_vertex_pt(2,Vertex_pt[6]);
   Facet_pt[1]->set_vertex_pt(3,Vertex_pt[7]);

   // top
   Facet_pt[2]=new TetMeshFacet(n_vertex_on_facet);
   one_based_boundary_id=4;
   Facet_pt[2]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[2]->set_vertex_pt(0,Vertex_pt[3]);
   Facet_pt[2]->set_vertex_pt(1,Vertex_pt[7]);
   Facet_pt[2]->set_vertex_pt(2,Vertex_pt[6]);
   Facet_pt[2]->set_vertex_pt(3,Vertex_pt[2]);

   Facet_pt[3]=new TetMeshFacet(n_vertex_on_facet);
   one_based_boundary_id=5;
   Facet_pt[3]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[3]->set_vertex_pt(0,Vertex_pt[6]);
   Facet_pt[3]->set_vertex_pt(1,Vertex_pt[5]);
   Facet_pt[3]->set_vertex_pt(2,Vertex_pt[1]);
   Facet_pt[3]->set_vertex_pt(3,Vertex_pt[2]);

   Facet_pt[4]=new TetMeshFacet(n_vertex_on_facet);
   one_based_boundary_id=6;
   Facet_pt[4]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[4]->set_vertex_pt(0,Vertex_pt[3]);
   Facet_pt[4]->set_vertex_pt(1,Vertex_pt[2]);
   Facet_pt[4]->set_vertex_pt(2,Vertex_pt[1]);
   Facet_pt[4]->set_vertex_pt(3,Vertex_pt[0]);

   // bottom
   Facet_pt[5]=new TetMeshFacet(n_vertex_on_facet);
   one_based_boundary_id=7;
   Facet_pt[5]->set_one_based_boundary_id(one_based_boundary_id);
   Facet_pt[5]->set_vertex_pt(0,Vertex_pt[5]);
   Facet_pt[5]->set_vertex_pt(1,Vertex_pt[1]);
   Facet_pt[5]->set_vertex_pt(2,Vertex_pt[0]);
   Facet_pt[5]->set_vertex_pt(3,Vertex_pt[4]);
  }

};


/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//====================================================================
/// Micky mouse  problem.
//====================================================================
template<class ELEMENT> 
class FallingBlockProblem : public Problem
{

public:


 /// Constructor
 FallingBlockProblem();
  
 /// Destructor (empty)
 ~FallingBlockProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

 /// Snap the boundary nodes onto the sphere
 void snap_onto_sphere()
  {
   this->mesh_pt()->output("pre_moved.dat",5);
   std::ofstream nodes("moved_nodes.dat");

   unsigned n_bound = mesh_pt()->nboundary_node(0);
   for(unsigned n=0;n<n_bound;++n)
    {
     Node* nod_pt = mesh_pt()->boundary_node_pt(0,n);
     double x = nod_pt->x(0);
     double y = nod_pt->x(1);
     double z = nod_pt->x(2);

     nodes << x << " " << y << " " << z << "  ";

     //Now let's snap by calculating the angle
     double r = sqrt(x*x + y*y + z*z);
     double theta = acos(z/r);
     double phi = atan2(y,x);
     
     nodes << r << " " << theta << " " << phi << " ";

     //Do the snapping
     double R_new = sqrt(1.0 + 0.5*0.5*(1.0 + sqrt(5.0))*(1.0 + sqrt(5.0)));
     nod_pt->x(0) = R_new*sin(theta)*cos(phi);
     nod_pt->x(1) = R_new*sin(theta)*sin(phi);
     nod_pt->x(2) = R_new*cos(theta);

     nodes << nod_pt->x(0) << " " << nod_pt->x(1) << " " << nod_pt->x(2) << "\n";
    }
   nodes.close();
   this->mesh_pt()->output("post_moved.dat",5);
  }
      
     
 /// Totally new mesh, need to fix it
 void actions_after_adapt()
  {
   // Set the boundary conditions for this problem 
   // Only on the "outer boundaries"
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned final_index = 3;
     //Do no pin the outlet z-velocity
     if(ibound==3) {final_index = 2;}
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<final_index;++i)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
        }
      }
    }

   // Complete the build of all elements so they are fully functional
   
   //Find number of elements in mesh
   unsigned n_element = mesh_pt()->nelement();
   
   // Loop over the elements to set up element-specific 
   // things that cannot be handled by constructor
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from GeneralElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
     
     //Set the source function pointer
     el_pt->re_pt() = &Global_Parameters::Re;
    }
  }
 


 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   //Loop over the boundaries 
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     //Don't boundary 3's z-coordinates
     unsigned final_index = 3;
     if(ibound==3) {final_index = 2;}

     // Loop over the nodes on boundary
     unsigned num_nod=mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
       for(unsigned i=0;i<final_index;++i) {nod_pt->set_value(i,0.0);}
      }
    }

   //Now set a Poiseuille like inflow
   {
    using namespace Global_Parameters;
    
    // Loop over the nodes on boundary
    unsigned num_nod=mesh_pt()->nboundary_node(6);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
       Node* nod_pt=mesh_pt()->boundary_node_pt(6,inod);
       double x = nod_pt->x(0);
       double y = nod_pt->x(1);
       double u = (Box_width-x)*(Box_width+x)*(Box_width-y)*(Box_width+y);

       nod_pt->set_value(2,u);
     }
   }
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}


 //Access function for the specific mesh
RefineableTetgenMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableTetgenMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

 
 /// Calculate the fluid dissipation
 double get_dissipation()
  {
   double dissipation=0.0;
   const unsigned n_element = this->mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *el_pt = 
      dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(e));
     //Add to the dissipation
     dissipation += el_pt->dissipation(); 
    }
   return dissipation;
  }

 /// Storage for the outer boundary object
 TetMeshFacetedClosedSurface* Outer_boundary_pt;

 Vector<TetMeshFacetedSurface*> Inner_boundary_pt;

};



//========================================================================
/// Constructor for FallingBlock problem
//========================================================================
template<class ELEMENT>
FallingBlockProblem<ELEMENT>::FallingBlockProblem()
{ 
 //Let's have stupidly high tolerance
 //Newton_solver_tolerance = 1000;

 //Add a steady time stepper
 this->add_time_stepper_pt(new Steady<0>);

 //Make the outer boundary object
 Outer_boundary_pt=new CubicTetMeshFacetedSurface;
 Outer_boundary_pt->output("outer.dat");

 //Create the inner boundary object
 Inner_boundary_pt.resize(1);
 Inner_boundary_pt[0] = new SphericalTetMeshFacetedSurface;
 Inner_boundary_pt[0]->output("inner.dat");

 Problem::mesh_pt() = 
  new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                    Inner_boundary_pt,2.0,
                                    this->time_stepper_pt(),
                                    true);

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 mesh_pt()->max_permitted_error()=0.005;
 mesh_pt()->min_permitted_error()=0.001; 
 mesh_pt()->max_element_size()=1.0;
 mesh_pt()->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   mesh_pt()->max_element_size()=2.0;
   mesh_pt()->min_element_size()=0.1; 
  }

 // Set the boundary conditions for this problem 
 // Only on the "outer boundaries"
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
  unsigned final_index = 3;
  //Do no pin the outlet z-velocity
  if(ibound==3) {final_index = 2;}
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    for(unsigned i=0;i<final_index;++i)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
     }
   }
 }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->re_pt() = &Global_Parameters::Re;
  }

 //Loop over the elements in region 1
 /*unsigned n_inner = mesh_pt()->nregion_element(1);
 for(unsigned e=0;e<n_inner;++e)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->region_element_pt(1,e));

   //Set the source function pointer
   el_pt->viscosity_ratio_pt() = &Global_Parameters::Visc_Ratio;
   }*/

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FallingBlockProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                           DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Doc local node numbering
 //-------------------------
 sprintf(filename,"%s/node_numbering%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 FiniteElement* el_pt=mesh_pt()->finite_element_pt(0);
 unsigned nnode=el_pt->nnode();
 unsigned ndim=el_pt->node_pt(0)->ndim();
 for (unsigned j=0;j<nnode;j++)
  {
   for (unsigned i=0;i<ndim;i++)
    {
     some_file << el_pt->node_pt(j)->x(i) << " " ;
    }
   some_file << j << std::endl;
  }
 some_file.close();

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,nplot);
 some_file.close();



} // end of doc




//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=2;

 //Output trace file
 std::ofstream dissipation("RESLT/diss.dat");

 // Do the problem with quadratic elements
 //---------------------------------------
 {
  FallingBlockProblem<ProjectableCrouzeixRaviartElement<
  TCrouzeixRaviartElement<3> > > problem;
  //problem.snap_onto_sphere();

  for(unsigned n=0;n<1;++n)
   {
    // Solve the problem
    problem.steady_newton_solve(1);
    
    //Output solution with 5 points per edge
    nplot=5;
    problem.doc_solution(nplot,doc_info);
    
    //Increment counter for solutions 
    doc_info.number()++;
    //Output the dissipation
    dissipation << " " << Global_Parameters::Re << "  "
                << problem.get_dissipation() << std::endl;

    //Increase the Reynolds number
    Global_Parameters::Re += 0.1;
   }

  dissipation.close();
 }


}




