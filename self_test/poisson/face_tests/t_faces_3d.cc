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
/// Driver code for a progam that exhaustively tests the faces of
/// Tetrahedral elements.
/// Note that the deformation of the cubic mesh into parabolic sides
/// will not be exactly representable by the TElement's because
/// the cross terms x^{2}y^{2} are not in the approximation space.
/// In order to construct an exactly representable surface, a single
/// TElement is created and rotated so that each face is "sloping"
/// in turn.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// Triangular mesh
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for functions that deform the cubic mesh into
/// a shape with parabolic sides.
//====================================================================
namespace MeshDeformation
{
 /// Deform the cubic mesh so that its six sides are all parabolic
 void deform_mesh(Mesh* const &mesh_pt)
 {

 //Deform each side of the mesh into a parabola, that
 //can be exactly represented by quadratic elements
 unsigned n_node = mesh_pt->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   Node* nod_pt = mesh_pt->node_pt(n);
   double x = nod_pt->x(0) - 1.0;
   double y = nod_pt->x(1) - 1.0;
   double z = nod_pt->x(2) - 1.0;
   
   double X_new = 1 + 0.2*(z-1)*(z+1)*(y-1)*(y+1);
   double Y_new = 1 - 0.1*(x-1)*(x+1)*(z-1)*(z+1);
   double Z_new = 1 + 0.2*(x-1)*(x+1)*(y-1)*(y+1);

   nod_pt->x(0) = X_new*x  + 1;
   nod_pt->x(1) = Y_new*y  + 1;
   nod_pt->x(2) = Z_new*z  + 1;
  }
 }

 /// Return the exact normal on the given face of the cubic mesh
 /// N.B. This is slightly different than the case in q_faces_3d.cc
 /// because the faces are labelled differently (DOH).
 void exact_normal(const unsigned &face,
                   const Vector<double> &x, Vector<double> &n)
 {
  double N[3] = {0.0,0.0,0.0};
  switch(face)
   {
   case 4:
    N[0] = -0.4*x[1]*(x[1]-2)*(x[0]-1);
    N[1] = -0.4*x[0]*(x[0]-2)*(x[1]-1);
    N[2] = -1.0;
    break;

   case 2:
    N[0] = 0.2*x[2]*(x[2]-2)*(x[0]-1); 
    N[1] = -1.0;
    N[2] = 0.2*x[0]*(x[0]-2)*(x[2]-1);  
    break;

   case 1:
    N[0] = 1.0;
    N[1] = -0.4*x[2]*(x[2]-2)*(x[1]-1);
    N[2] = -0.4*x[1]*(x[1]-2)*(x[2]-1);
    break;

   case 3:
    N[0] = 0.2*x[2]*(x[2]-2)*(x[0]-1); 
    N[1] = 1.0;
    N[2] = 0.2*x[0]*(x[0]-2)*(x[2]-1);  
    break;

   case 0:
    N[0] = -1.0;
    N[1] = -0.4*x[2]*(x[2]-2)*(x[1]-1);
    N[2] = -0.4*x[1]*(x[1]-2)*(x[2]-1);
    break;

   case 5:
    N[0] = -0.4*x[1]*(x[1]-2)*(x[0]-1);
    N[1] = -0.4*x[0]*(x[0]-2)*(x[1]-1);
    N[2] = 1.0;
    break;
   }
  double length = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
  n[0] = N[0]/length;
  n[1] = N[1]/length;
  n[2] = N[2]/length;
    
 }
}

//============================================================
/// Helper namespace for the single element test case
//===========================================================
namespace SingleElement
{
 //Return the normal to the face based on the global x-coordinate
 void face_normal(const unsigned &face,
                  const Vector<double> &x, Vector<double> &n)
 {
  double N[3] = {0.0,0.0,0.0};
  switch(face)
   {
   case 0:
    N[0] = -1.0;
    N[1] = 0.0;
    N[2] = 0.0;
    break;

   case 1:
    N[0] = 0.0;
    N[1] = -1.0;
    N[2] = 0.0;
    break;

   case 2:
    N[0] = 0.0;
    N[1] = 0.0;
    N[2] = -1.0;
    break;

   case 3:
    Vector<double> s(3);
    //Calculate the local "bulk" (unrotated) coordinate on the face
    for(unsigned i=0;i<3;i++)
     {
      s[0] = -0.5 + 0.5*sqrt(4.0*x[0] - 3);
      s[1] = -0.5 + 0.5*sqrt(4.0*x[1] - 3);
      s[2] = -0.5 + 0.5*sqrt(4.0*x[2] - 3);
     }
      
    //Hence compute the outer normal
    N[0] = (1.0 + 2.0*s[1])*(1.0 + 2.0*s[2]);
    N[1] = (1.0 + 2.0*s[0])*(1.0 + 2.0*s[2]);
    N[2] = (1.0 + 2.0*s[0])*(1.0 + 2.0*s[1]);
    break;
   }

  double length = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
  n[0] = N[0]/length;
  n[1] = N[1]/length;
  n[2] = N[2]/length;
 }

}


/// Little helper class that is used to convert each element to 
/// a unique filestem
template<class ELEMENT>
class FileStem
{
 static string Stem;
public:

 explicit FileStem() {}

 inline const string &operator()() const {return Stem;}

};

/// The TPoissonElement<3,2> is linear
template<>
string FileStem<TPoissonElement<3,2> >::Stem = "linear"; 
/// The TPoissonElement<3,3> is quadratic
template<>
string FileStem<TPoissonElement<3,3> >::Stem = "quadratic"; 
/// The TPoissonElement<3,4> is cubic
template<>
string FileStem<TPoissonElement<3,4> >::Stem = "cubic"; 



//====================================================================
/// TFaceTest problem.
//====================================================================
template<class ELEMENT> 
class TFaceTestProblem : public Problem
{

public:

 /// Constructor
 TFaceTestProblem(const unsigned& h_power, const unsigned &rotate,
                  const bool &single_element=true);
 
 /// Destructor (empty)
 ~TFaceTestProblem(){};

 /// Empty actions before solve (we never solve a problem)
 void actions_before_newton_solve() {}

 /// Empty actions after solve (we never solve a problem)
 void actions_after_newton_solve()
  {}

private:

 /// Exponent for h scaling
 unsigned H_power;

};


//========================================================================
/// Constructor for TFaceTest problem
///
//========================================================================
template<class ELEMENT>
TFaceTestProblem<ELEMENT>::
TFaceTestProblem(const unsigned& h_power, const unsigned &rotate,
                 const bool &single_element) : 
 H_power(h_power)
{ 
 //Create a unique stem for the files, depending on the element type
 string filestem = FileStem<ELEMENT>()();

 //Permutation (by rotation) of the faces
 //The "old" face represented by the "new" face i
 //is given by perm[rotate][i]
 unsigned perm[4][4] ={{0,1,2,3},{3,0,2,1},{1,3,2,0},{2,1,3,0}};

 //If we are just computing the single element (the default)
 if(single_element)
  {
   //Make a single element
   ELEMENT* local_elem_pt = new ELEMENT;
   //Find the number of nodess
   unsigned n_node = local_elem_pt->nnode();
   //Local storage for the nodes
   Vector<Node*> local_node_pt(n_node);
   //The local coordinate of each node
   Vector<double> s(3);
   //Create the local nodes
   for(unsigned n=0;n<n_node;n++)
    {
     //Make the node
     local_node_pt[n] = local_elem_pt->construct_node(n);
     //Get the local coordinate of the node
     local_elem_pt->local_coordinate_of_node(n,s);
     //Now set the nodal values to vary quadratically
     for(unsigned i=0;i<3;i++)
      {
       local_node_pt[n]->x(i) = 1.0 + s[i] + s[i]*s[i];
      }
    }

   //Permute the nodes so that we can test all the faces
   //These correspond to rotations
   if(dynamic_cast<TPoissonElement<3,2>*>(local_elem_pt))
    {
     switch(rotate)
      {
      case 0:
       break;

      case 1:
       local_elem_pt->node_pt(0) = local_node_pt[3];
       local_elem_pt->node_pt(1) = local_node_pt[0];
       local_elem_pt->node_pt(3) = local_node_pt[1];
       break;
      
      case 2:
       local_elem_pt->node_pt(0) = local_node_pt[1];
       local_elem_pt->node_pt(1) = local_node_pt[3];
       local_elem_pt->node_pt(3) = local_node_pt[0];
       break;

      case 3:
       local_elem_pt->node_pt(0) = local_node_pt[2];
       local_elem_pt->node_pt(2) = local_node_pt[3];
       local_elem_pt->node_pt(3) = local_node_pt[0];
       break;
      }
    }
   else if(dynamic_cast<TPoissonElement<3,3>*>(local_elem_pt))
    {
     switch(rotate)
      {
      case 0:
       break;

      case 1:
       local_elem_pt->node_pt(0) = local_node_pt[3];
       local_elem_pt->node_pt(1) = local_node_pt[0];
       local_elem_pt->node_pt(3) = local_node_pt[1];

       local_elem_pt->node_pt(4) = local_node_pt[6];
       local_elem_pt->node_pt(5) = local_node_pt[8];
       local_elem_pt->node_pt(6) = local_node_pt[9];
       local_elem_pt->node_pt(7) = local_node_pt[5];
       local_elem_pt->node_pt(8) = local_node_pt[7];
       local_elem_pt->node_pt(9) = local_node_pt[4];
       break;
      
      case 2:
       local_elem_pt->node_pt(0) = local_node_pt[1];
       local_elem_pt->node_pt(1) = local_node_pt[3];
       local_elem_pt->node_pt(3) = local_node_pt[0];

       local_elem_pt->node_pt(4) = local_node_pt[9];
       local_elem_pt->node_pt(5) = local_node_pt[7];
       local_elem_pt->node_pt(6) = local_node_pt[4];
       local_elem_pt->node_pt(7) = local_node_pt[8];
       local_elem_pt->node_pt(8) = local_node_pt[5];
       local_elem_pt->node_pt(9) = local_node_pt[6];
       break;

      case 3:
       local_elem_pt->node_pt(0) = local_node_pt[2];
       local_elem_pt->node_pt(2) = local_node_pt[3];
       local_elem_pt->node_pt(3) = local_node_pt[0];

       local_elem_pt->node_pt(4) = local_node_pt[7];
       local_elem_pt->node_pt(5) = local_node_pt[8];
       local_elem_pt->node_pt(6) = local_node_pt[5];
       local_elem_pt->node_pt(7) = local_node_pt[9];
       local_elem_pt->node_pt(8) = local_node_pt[6];
       local_elem_pt->node_pt(9) = local_node_pt[4];
       break;
      }
    }       

   //Construct face elements on all the boundaries of the single element
   for(unsigned b=0;b<4;b++)
    {
     double error=0.0;
     const unsigned n_pts = 5;
     
     // Build the corresponding prescribed-flux element
     PoissonFluxElement<ELEMENT>* flux_element_pt = new 
      PoissonFluxElement<ELEMENT>(local_elem_pt,b);
     
     //Local storage for the positions (from local and bulk representations)
     Vector<double> x(3);
     Vector<double> X(3);
     //Local storage for the exact and computed normals
     Vector<double> n(3);
     Vector<double> N(3);
     //Local storage for the local and "bulk local" coordinates
     Vector<double> s(2);
     Vector<double> s_bulk(3);
    
     //Create an output file
     char filename[100];
     sprintf(filename,"%s_normals%i_%i.dat",filestem.c_str(),b,rotate);
     ofstream output(filename);
     output << flux_element_pt->tecplot_zone_string(n_pts);
     double local_error = 0.0;
     unsigned num_plot_pts = flux_element_pt->nplot_points(n_pts);
     //Loop over the plot points
     for(unsigned i=0;i<num_plot_pts;i++)
      {
       //Get the local coordinate
       flux_element_pt->get_s_plot(i,n_pts,s);
       //Get the Eulerian position from the face element
       flux_element_pt->interpolated_x(s,x);
       //Get the outer unit normal from the face
       flux_element_pt->outer_unit_normal(s,n);
       //Get the local coordinate in the bulk
       flux_element_pt->get_local_coordinate_in_bulk(s,s_bulk);
       //Get the excat solution for the face normal
       SingleElement::face_normal(perm[rotate][b],x,N);
       //Get the Eulerian position from the bulk element
       local_elem_pt->interpolated_x(s_bulk,X);
       //Print the results
       output << " " << x[0] << " " << x[1] << " " << x[2] << " "
              << X[0] << " " << X[1] << " " << X[2] << " "  << 
        n[0] << " " << n[1] << "  " << n[2] << " " 
              << N[0] << " " << N[1] << " " << N[2] << "\n";
       
       //Calculate the error in the normals
       error += sqrt((N[0]-n[0])*(N[0]-n[0]) + (N[1]-n[1])*(N[1]-n[1])
                     + (N[2] - n[2])*(N[2]-n[2]));
      }
     local_error /= num_plot_pts;
     error += local_error;
     flux_element_pt->write_tecplot_zone_footer(output,n_pts);
     //Delete the storage
     delete flux_element_pt;
     output.close();
     std::cout << "Boundary: " << b << " error " << error << "\n";
    }
   
   //Delete the element and all the nodes
   delete local_elem_pt;
   for(unsigned n=0;n<n_node;n++) {delete local_node_pt[n];}
  }
 //Otherwise we are doing the full mesh NOT included in the tests
 else
  {
   //Create mesh and assign element lengthscale h
   unsigned nx=unsigned(pow(2.0,int(h_power)));
   unsigned ny=unsigned(pow(2.0,int(h_power)));
   unsigned nz=unsigned(pow(2.0,int(h_power)));
   double lx=2.0;
   double ly=2.0;
   double lz=2.0;
   
   // Build and assign mesh
   Problem::mesh_pt() = new SimpleCubicTetMesh<ELEMENT>(nx,ny,nz,lx,ly,lz);
   //Setup the boundary information
   mesh_pt()->setup_boundary_element_info();

   //Deform the mesh
   MeshDeformation::deform_mesh(mesh_pt());

   //Storage for the faces on each boundary
 Vector<std::list<FaceElement*> > faces_on_boundary(6);
 
 //Construct face elements on all the boundaries
 for(unsigned b=0;b<6;b++)
  {
   double error=0.0;
   unsigned n_el = mesh_pt()->nboundary_element(b);
   const unsigned n_pts = 5;
   for(unsigned e=0;e<n_el;e++)
    {
     //Get Pointer to element adjacent to the boundary
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      mesh_pt()->boundary_element_pt(b,e));
     
     //What is the index of face of the element e along boundary b
     int face_index = mesh_pt()->face_index_at_boundary(b,e);
          
     // Build the corresponding prescribed-flux element
     PoissonFluxElement<ELEMENT>* flux_element_pt = new 
      PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

     faces_on_boundary[b].push_back(flux_element_pt);
    } //end of loop over bulk elements adjacent to boundary b
   
   //Local storage for position, local coordinate and the normals
   Vector<double> x(3);
   Vector<double> n(3);
   Vector<double> N(3);
   Vector<double> s(2);
   
   char filename[100];
   sprintf(filename,"%s_normals%i.dat",filestem.c_str(),b);
   ofstream output(filename);
   for(std::list<FaceElement*>::iterator it=faces_on_boundary[b].begin();
       it!=faces_on_boundary[b].end();++it)
    {
     output << (*it)->tecplot_zone_string(n_pts);
     double local_error = 0.0;
     unsigned num_plot_pts = (*it)->nplot_points(n_pts);
     for(unsigned i=0;i<num_plot_pts;i++)
      {
       (*it)->get_s_plot(i,n_pts,s);
       (*it)->interpolated_x(s,x);
       (*it)->outer_unit_normal(s,n);
       MeshDeformation::exact_normal(b,x,N);
       output << b << " " << x[0] << " " << x[1] << " " << x[2] << " " <<
        n[0] << " " << n[1] << "  " << n[2] << " " 
              << N[0] << " " << N[1] << " " << N[2] << "\n";
       
       error += sqrt((N[0]-n[0])*(N[0]-n[0]) + (N[1]-n[1])*(N[1]-n[1])
                     + (N[2] - n[2])*(N[2]-n[2]));
      }
     local_error /= num_plot_pts;
     error += local_error;
     (*it)->write_tecplot_zone_footer(output,n_pts);
     //Delete the storage
     delete *it;
    }
 output.close();
 error/= n_el;
 std::cout << "Boundary: " << b << " error " << error << "\n";
  }
  }
}

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//========================================================================
/// Run the test for linear and quadratic TPoisson elements
//========================================================================
int main(int argc, char* argv[])
{
 //Set a reasonable (not too large) power
 unsigned h_power = 2;

 cout << "Linear elements " << std::endl;
 cout << "================" << std::endl << std::endl;
 
 //Loop over all four possible rotations
 for(unsigned rotate=0;rotate<4;rotate++)
 {
 //Set up the problem
  cout << "Rotation " << rotate << std::endl;
  TFaceTestProblem<TPoissonElement<3,2> >  problem(h_power,rotate);
 }

 
 cout << "Quadratic elements " << std::endl;
 cout << "==================" << std::endl << std::endl;
 
 //Loop over all four possible rotations
 for(unsigned rotate=0;rotate<4;rotate++)
 {
  cout << "Rotation " << rotate << std::endl;
 //Set up the problem
 TFaceTestProblem<TPoissonElement<3,3> > problem(h_power,rotate);

 }

}


