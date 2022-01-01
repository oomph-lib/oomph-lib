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
//Driver function to test the calculation of outer unit normals
//for face elements of triangular elements TElement<2,.>
//We actually construct TPoisson elements because we already
//have all the machinary in place to calculate the Flux (face) elements.

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// Triangular mesh
#include "meshes/simple_rectangular_tri_mesh.h"

using namespace std;

using namespace oomph;

//====================================================================
/// Namespace for function that represents the exact 
/// (quadratically-varying) outer unit normal.
/// and the function the deforms the mesh
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

 /// The outer unit normal of the mesh boundaries
 /// that have been deformed into parabolas. The normal is different
 /// for each boundary, so the input arguments are the boundary and
 /// the global Cartesian coordinate on that boundary, x. The outer
 /// unit normal is returned in the vector n.
 void exact_normal(const unsigned &boundary,
                   const Vector<double> &x, Vector<double> &n)
 {
  //Temporary local storage for the non-unit normals
  double N[2] = {0.0,0.0};

  //Determine the boundary
  switch(boundary)
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

   default:
    //Throw an error
    std::ostringstream error_stream;
    error_stream << "There are only four boundaries in this mesh.\n";
    error_stream << "Called with boundary number " << boundary << "\n";
    
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    break;
   }

  //Normalise the vector
  double length = sqrt(N[0]*N[0] + N[1]*N[1]);
  n[0] = N[0]/length;
  n[1] = N[1]/length;
 }

}

//====================================================================
/// TriangleFaceTest problem.
//====================================================================
template<class ELEMENT> 
class TriangleFaceTestProblem : public Problem
{

public:

 /// Constructor
 TriangleFaceTestProblem(const unsigned& h_power,const bool &permute=false);
 
 /// Destructor (empty)
 ~TriangleFaceTestProblem(){};

 /// Empty actions before solve (we never solve the problem)
 void actions_before_newton_solve() {}

 /// Empty actions after solve (we never solve the problem)
 void actions_after_newton_solve() {}

private:

 /// Exponent for h scaling
 unsigned H_power;

};



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
string FileStem<TPoissonElement<2,2> >::Stem = "linear"; 
/// The QPoissonElement<2,3> is quadratic
template<>
string FileStem<TPoissonElement<2,3> >::Stem = "quadratic"; 
/// The QPoissonElement<2,4> is cubic
template<>
string FileStem<TPoissonElement<2,4> >::Stem = "cubic"; 


//========================================================================
/// Constructor for TriangleFaceTest problem
///
//========================================================================
template<class ELEMENT>
TriangleFaceTestProblem<ELEMENT>::
TriangleFaceTestProblem(const unsigned& h_power,const bool &permute) : 
 H_power(h_power)
{ 
 //Create mesh and assign element lengthscale h
 unsigned nx=unsigned(pow(2.0,int(h_power)));
 unsigned ny=unsigned(pow(2.0,int(h_power)));
 double lx=8.0;
 double ly=8.0;
 
 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularTriMesh<ELEMENT>(nx,ny,lx,ly);

 //Permute the nodes of all triangle
 //If we do not do such a permutation then not all boundaries of the
 //triangles are tested.
 if(permute)
 {
  unsigned n_element = mesh_pt()->nelement();
  for(unsigned e=0;e<n_element;e++)
   {
    FiniteElement* el_pt = mesh_pt()->finite_element_pt(e);
    unsigned n_node = el_pt->nnode();
    Node* local_node_pt[n_node];
    for(unsigned n=0;n<n_node;n++)
     {
      local_node_pt[n] = el_pt->node_pt(n);
     }
    
    //Now do the permutation
    if(dynamic_cast<TPoissonElement<2,2>* >(el_pt))
     {
      el_pt->node_pt(0) = local_node_pt[1];
      el_pt->node_pt(1) = local_node_pt[2];
      el_pt->node_pt(2) = local_node_pt[0];
     }

    if(dynamic_cast<TPoissonElement<2,3>* >(el_pt))
     {
      el_pt->node_pt(0) = local_node_pt[1];
      el_pt->node_pt(1) = local_node_pt[2];
      el_pt->node_pt(2) = local_node_pt[0];
      el_pt->node_pt(3) = local_node_pt[4];
      el_pt->node_pt(4) = local_node_pt[5];
      el_pt->node_pt(5) = local_node_pt[3];
     }

    if(dynamic_cast<TPoissonElement<2,4>* >(el_pt))
     {
      el_pt->node_pt(0) = local_node_pt[1];
      el_pt->node_pt(1) = local_node_pt[2];
      el_pt->node_pt(2) = local_node_pt[0];
      el_pt->node_pt(3) = local_node_pt[5];
      el_pt->node_pt(4) = local_node_pt[6];
      el_pt->node_pt(5) = local_node_pt[7];
      el_pt->node_pt(6) = local_node_pt[8];
      el_pt->node_pt(7) = local_node_pt[3];
      el_pt->node_pt(8) = local_node_pt[4];
     }
   }
 }

 // Now Setup the boundary info (after any potential permutation)
 this->mesh_pt()->setup_boundary_element_info();

 // Now deform the mesh
 MeshDeformation::deform_mesh(mesh_pt());
 
 //Storage for the face elements formed from a boundary of the mesh
 Vector<std::list<FaceElement*> > faces_on_boundary(4);


 //Create a unique stem for the files, depending on the element type
 string filestem = FileStem<ELEMENT>()();
 //Output for the errors on each face
 char error_file[100];
 sprintf(error_file,"%s_errors_%i.dat",filestem.c_str(),permute);
 ofstream output_error(error_file);

 if(permute) {std::cout << "Permuted Nodes" << std::endl;}
  
 //Construct face elements on all the boundaries
 for(unsigned b=0;b<4;b++)
  {
   double error=0.0;
   unsigned n_el = mesh_pt()->nboundary_element(b);
   const unsigned n_pts = 5;
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

   //Now print the outer unit normals halfway along each edge
   Vector<double> x(2);
   Vector<double> n(2);
   Vector<double> N(2);
   char filename[100];
   sprintf(filename,"%s_normals%i_%i.dat",
           filestem.c_str(),b,permute);
   ofstream output(filename);
     for(std::list<FaceElement*>::iterator it=faces_on_boundary[b].begin();
         it!=faces_on_boundary[b].end();++it)
      {
       for(unsigned i=0;i<n_pts;i++)
        {
         Vector<double> s(1,i/(double)(n_pts-1));
         (*it)->interpolated_x(s,x);
         (*it)->outer_unit_normal(s,n);
         MeshDeformation::exact_normal(b,x,N);
         output << b << " " << x[0] << " " << x[1] << " " <<
          n[0] << " " << n[1] << "  "<< N[0] << " " << N[1] << "\n";
         
         error += sqrt((N[0]-n[0])*(N[0]-n[0]) + (N[1]-n[1])*(N[1]-n[1]));
        }
       //Delete the storage
       delete *it;
      }
     output.close();
     error/= n_el*n_pts;
     output_error << b << " " << error << "\n";
   std::cout << "Boundary: " << b << " error " << error << "\n";
   }
 output_error.close();
}

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//========================================================================
/// Demonstrate how to solve TriangleFaceTest problem
//========================================================================
int main(int argc, char* argv[])
{
 unsigned h_power = 3;

 {
  cout << "Linear elements " << std::endl;
  cout << "================" << std::endl << std::endl;
  //Set up the problem
  TriangleFaceTestProblem<TPoissonElement<2,2> >  problem(h_power);
  //Permuted version
  TriangleFaceTestProblem<TPoissonElement<2,2> >  problem2(h_power,true);
 }

 {
 cout << "Quadratic elements " << std::endl;
 cout << "==================" << std::endl << std::endl;
 //Set up the problem
 TriangleFaceTestProblem<TPoissonElement<2,3> > problem(h_power);
  //Permuted version
 TriangleFaceTestProblem<TPoissonElement<2,3> >  problem2(h_power,true);
 }
 
 {
  cout << "Cubic elements " << std::endl;
  cout << "================" << std::endl << std::endl;
  //Set up the problem
  TriangleFaceTestProblem<TPoissonElement<2,4> > problem(h_power);
  //Permuted version
  TriangleFaceTestProblem<TPoissonElement<2,4> >  problem2(h_power,true);
 }
 
}


