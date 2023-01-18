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
//Driver function for an exhaustive test of the face element
//construcion of two-dimensional quad-type elements.

//The idea is to take the rectangular quad mesh and deform each side
//into a parabola, which can be exactly represented by quadratic (or higher)
//interpolation. The normals computed by the FaceElements should then be
//exact for QElement<2,3> and QElement<2,4>.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// Triangular mesh
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for the functions associated with the mesh deformation
//====================================================================
namespace MeshDeformation
{
 
 //Deform each side of the mesh into a parabola, that
 //can be exactly represented by quadratic elements
 void deform_mesh(Mesh* const &mesh_pt)
 {
  unsigned n_node = mesh_pt->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    Node* nod_pt = mesh_pt->node_pt(n);
    double x = nod_pt->x(0) - 4;
    double y = nod_pt->x(1) - 4;
    
    double X_new = 5 - y*y/16;
    double Y_new = 3 + x*x/16;
    
    nod_pt->x(0) = X_new*x/4.0 + 4;
    nod_pt->x(1) = Y_new*y/4.0 + 4;
   }
 }


 /// The exact normal vector for each face at global coordinate x.
 void exact_normal(const unsigned &face,
                   const Vector<double> &x, Vector<double> &n)
 {
  double N[2] = {0.0,0.0};
  switch(face)
   {
   case 0:
    N[0] = -(x[0]-4.0)/8.0;
    N[1] = -1.0;
    break;

   case 1:
    N[0] = 1.0;
    N[1] = (x[1]-4)/8.0;
    break;

   case 2:
    N[0] = -(x[0]-4)/8.0;
    N[1] = 1.0;
    break;

   case 3:
    N[0] = -1.0;
    N[1] = (x[1]-4)/8.0; 
    break;
   }
  double length = sqrt(N[0]*N[0] + N[1]*N[1]);
  n[0] = N[0]/length;
  n[1] = N[1]/length;
    
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

/// The QPoissonElement<2,2> is linear
template<>
string FileStem<QPoissonElement<2,2> >::Stem = "linear"; 
/// The QPoissonElement<2,3> is quadratic
template<>
string FileStem<QPoissonElement<2,3> >::Stem = "quadratic"; 
/// The QPoissonElement<2,4> is cubic
template<>
string FileStem<QPoissonElement<2,4> >::Stem = "cubic"; 



//====================================================================
/// QFaceTest problem.
//====================================================================
template<class ELEMENT> 
class QFaceTestProblem : public Problem
{

public:

 /// Constructor
  QFaceTestProblem(const unsigned& h_power);

 /// Destructor (empty)
 ~QFaceTestProblem(){};

 /// Empty actions before solve (we never solve the problem in this test)
 void actions_before_newton_solve() {}

 /// Empty actions after solve
 void actions_after_newton_solve() {}

private:

 /// Exponent for h scaling
 unsigned H_power;

};


//========================================================================
/// Constructor for QFaceTest problem
///
//========================================================================
template<class ELEMENT>
QFaceTestProblem<ELEMENT>::
QFaceTestProblem(const unsigned& h_power) : H_power(h_power)
{ 
 //Create mesh and assign element lengthscale h
 unsigned nx=unsigned(pow(2.0,int(h_power)));
 unsigned ny=unsigned(pow(2.0,int(h_power)));
 double lx=8.0;
 double ly=8.0;
 
 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly);

 //Deform the mesh
 MeshDeformation::deform_mesh(mesh_pt());

 //Storage for the faces on each boundary of the mesh
 Vector<std::list<FaceElement*> > faces_on_boundary(4);

 //Create a unique stem for the files, depending on the element type
 string filestem = FileStem<ELEMENT>()();
 //Output for the errors on each face
 char error_file[100];
 sprintf(error_file,"%s_errors.dat",filestem.c_str());
 ofstream output_error(error_file);
 
 //Construct face elements on all the boundaries
 for(unsigned b=0;b<4;b++)
  {
   //Number of output points
   const unsigned n_pts = 5;
      
   //The error between the computed and exact normals on each boundary
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
   

   //Storage for global position, computed normal and exact normal
   Vector<double> x(2), n(2), N(2);
   //Open a filename with a unique name for each different element type
   char filename[100];
   sprintf(filename,"%s_normals%i.dat",filestem.c_str(),b);
   ofstream output(filename);

   //Loop over the face elements
   for(std::list<FaceElement*>::iterator it=faces_on_boundary[b].begin();
       it!=faces_on_boundary[b].end();++it)
    {
     //Loop over n_pts plot points and calculate the
     //global position, normal and exact normal
     for(unsigned i=0;i<n_pts;i++)
      {
       Vector<double> s(1,-1.0 + i*(2.0/(double)(n_pts-1)));
       (*it)->interpolated_x(s,x);
       (*it)->outer_unit_normal(s,n);
       MeshDeformation::exact_normal(b,x,N);
       //Send the results to the output
       output << b << " " << x[0] << " " << x[1] << " " <<
        n[0] << " " << n[1] << "  "<< N[0] << " " << N[1] << "\n";
       
       //Compute the mean square error
       error += sqrt((N[0]-n[0])*(N[0]-n[0]) + (N[1]-n[1])*(N[1]-n[1]));
      }
     //Delete the face element
     delete *it;
    }
   output.close();
   error/= n_el*n_pts;
   
   //Output the error
   output_error << b << " " << error << "\n";
   std::cout << "Boundary: " << b << " error " << error << "\n";
  }
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
 //Let's choose a reasonable level of refinement
 unsigned h_power = 3;

 {
  cout << "Linear elements " << std::endl;
  cout << "================" << std::endl << std::endl;
  //Set up the problem
  QFaceTestProblem<QPoissonElement<2,2> >  problem(h_power);
 }

 {
 cout << "Quadratic elements " << std::endl;
 cout << "==================" << std::endl << std::endl;
 //Set up the problem
 QFaceTestProblem<QPoissonElement<2,3> > problem(h_power);
 }

 {
 cout << "Cubic elements " << std::endl;
 cout << "================" << std::endl << std::endl;
 //Set up the problem
 QFaceTestProblem<QPoissonElement<2,4> > problem(h_power);
 }

}


