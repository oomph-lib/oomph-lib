//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Driver function for an exhaustive test of the face element
//construcion of three-dimensional quad-type elements.

//The idea is to take the cubic mesh and deform each side
//into a parabola, which can be exactly represented by quadratic (or higher)
//interpolation. The normals computed by the FaceElements should then be
//exact for QElement<3,3>.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// Cubic mesh
#include "meshes/simple_cubic_mesh.h"


using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for the functions associated with the mesh deformation
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

 /// The exact normal vector for each face at global coordinate x
 void exact_normal(const unsigned &face,
                   const Vector<double> &x, Vector<double> &n)
 {
  double N[3] = {0.0,0.0,0.0};
  switch(face)
   {
   case 0:
    N[0] = -0.4*x[1]*(x[1]-2)*(x[0]-1);
    N[1] = -0.4*x[0]*(x[0]-2)*(x[1]-1);
    N[2] = -1.0;
    break;

   case 1:
    N[0] = 0.2*x[2]*(x[2]-2)*(x[0]-1); 
    N[1] = -1.0;
    N[2] = 0.2*x[0]*(x[0]-2)*(x[2]-1);  
    break;

   case 2:
    N[0] = 1.0;
    N[1] = -0.4*x[2]*(x[2]-2)*(x[1]-1);
    N[2] = -0.4*x[1]*(x[1]-2)*(x[2]-1);
    break;

   case 3:
    N[0] = 0.2*x[2]*(x[2]-2)*(x[0]-1); 
    N[1] = 1.0;
    N[2] = 0.2*x[0]*(x[0]-2)*(x[2]-1);  
    break;

   case 4:
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

/// The QPoissonElement<3,2> is linear
template<>
string FileStem<QPoissonElement<3,2> >::Stem = "linear"; 
/// The QPoissonElement<3,3> is quadratic
template<>
string FileStem<QPoissonElement<3,3> >::Stem = "quadratic"; 
/// The QPoissonElement<3,4> is cubic
template<>
string FileStem<QPoissonElement<3,4> >::Stem = "cubic"; 


//====================================================================
/// Face Test problem.
//====================================================================
template<class ELEMENT> 
class QFaceTestProblem : public Problem
{

public:

 /// Constructor
  QFaceTestProblem(const unsigned& h_power);

 /// Destructor (empty)
 ~QFaceTestProblem(){};

 /// Empty actions before solve (we don't solve a problem)
 void actions_before_newton_solve() {}

 /// Empty actions after solve (we don't solve a problem)
 void actions_after_newton_solve() {}

private:

 /// Exponent for h scaling
 unsigned H_power;

};


//========================================================================
/// Constructor for QFaceTest problem
//========================================================================
template<class ELEMENT>
QFaceTestProblem<ELEMENT>::
QFaceTestProblem(const unsigned& h_power) : H_power(h_power)
{ 
 //Create mesh and assign element lengthscale h
 unsigned nx=unsigned(pow(2.0,int(h_power)));
 unsigned ny=unsigned(pow(2.0,int(h_power)));
 unsigned nz=unsigned(pow(2.0,int(h_power)));
 double lx=2.0;
 double ly=2.0;
 double lz=2.0;

 // Build and assign mesh
 Problem::mesh_pt() = new SimpleCubicMesh<ELEMENT>(nx,ny,nz,lx,ly,lz);
 //Deform the mesh
 MeshDeformation::deform_mesh(mesh_pt());

 //Storage for the faces on each boundary of the mesh
 Vector<std::list<FaceElement*> > faces_on_boundary(6);

  //Create a unique stem for the files, depending on the element type
 string filestem = FileStem<ELEMENT>()();
 //Output for the errors on each face
 char error_file[100];
 sprintf(error_file,"%s_errors.dat",filestem.c_str());
 ofstream output_error(error_file);

 //Construct face elements on all the boundaries
 for(unsigned b=0;b<6;b++)
  {
   //Number of output points
   const unsigned n_pts = 5;
   //The RMS error between computed and exact normals on the boundary
   double error=0.0;
   //Loop over the elements adjacent to each face
   unsigned n_el = mesh_pt()->nboundary_element(b);
   for(unsigned e=0;e<n_el;e++)
    {
     //Get Pointer to element adjacent to the boundary
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      mesh_pt()->boundary_element_pt(b,e));
     
     //What is the index of the face of element e along boundary b
     int face_index = mesh_pt()->face_index_at_boundary(b,e);
          
     // Build the corresponding prescribed-flux element
     PoissonFluxElement<ELEMENT>* flux_element_pt = new 
      PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

     faces_on_boundary[b].push_back(flux_element_pt);
    } //end of loop over bulk elements adjacent to boundary b
   
   //Storage for local and global position and the computed and exact normals
   Vector<double> x(3), n(3), N(3),s(2);
   //Open a filename with a unique name for each different element type
   char filename[100];
   sprintf(filename,"%s_normals%i.dat",filestem.c_str(),b);
   ofstream output(filename);
   
   //Now loop over the face elements
   for(std::list<FaceElement*>::iterator it=faces_on_boundary[b].begin();
       it!=faces_on_boundary[b].end();++it)
    {
     //Loop over n_pts plot points in each direction
     for(unsigned i=0;i<n_pts;i++)
      {
       s[0] = -1.0 + i*(2.0/(double)(n_pts-1));
       for(unsigned j=0;j<n_pts;j++)
        {
         s[1] = -1.0 + j*(2.0/(double)(n_pts-1));

         //Calculate and output position and the normals
         (*it)->interpolated_x(s,x);
         (*it)->outer_unit_normal(s,n);
         MeshDeformation::exact_normal(b,x,N);
         output << b << " " << x[0] << " " << x[1] << " " << x[2] << " " <<
          n[0] << " " << n[1] << "  " << n[2] << " " 
                << N[0] << " " << N[1] << " " << N[2] << "\n";
         
         error += sqrt((N[0]-n[0])*(N[0]-n[0]) + (N[1]-n[1])*(N[1]-n[1])
                       + (N[2] - n[2])*(N[2]-n[2]));
        }
      }
     //Delete the face element
     delete *it;
    }
   output.close();
   error/= n_el*n_pts*n_pts;

   //Output the error
   output_error << b << " " << error << "\n";
   std::cout << "Boundary: " << b << " error " << error << "\n";
  }

 //Close the output file
 output_error.close();
}

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//========================================================================
/// Demonstrate how to solve QFaceTest problem
//========================================================================
int main(int argc, char* argv[])
{
 //Choose a resonable resolution
 unsigned h_power = 2;

 {
  cout << "Linear elements " << std::endl;
  cout << "================" << std::endl << std::endl;
  //Set up the problem
  QFaceTestProblem<QPoissonElement<3,2> >  problem(h_power);
 }

 {
  cout << "Quadratic elements " << std::endl;
  cout << "==================" << std::endl << std::endl;
  //Set up the problem
  QFaceTestProblem<QPoissonElement<3,3> > problem(h_power);
 }
 
 {
  cout << "Cubic elements " << std::endl;
  cout << "================" << std::endl << std::endl;
  //Set up the problem
  QFaceTestProblem<QPoissonElement<3,4> > problem(h_power);
 }

}


