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
//Driver for a simple 2D poisson problem with adaptive mesh refinement
#include <typeinfo>

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

using namespace oomph;

using namespace std;

//================================================================
/// Class for sorting the elements into lexicographic order based
/// on the coordinate at the centre of the element. This is the only
/// "rational" way to preserve the ordering given rotations of the elements.
//================================================================
template<class ELEMENT>
class CompareElementCoordinate 
{
public:
 
 ///The actual comparison operator
 int operator() (GeneralisedElement* const &element1_pt,
                 GeneralisedElement* const &element2_pt)
  {
   
   //Dynamic cast the elements
   ELEMENT* cast_element1_pt = dynamic_cast<ELEMENT*>(element1_pt);
   ELEMENT* cast_element2_pt = dynamic_cast<ELEMENT*>(element2_pt);
   
   //Check that we managed to successfully cast the elements
#ifdef PARANOID
     if (cast_element1_pt==0)
      {
       std::ostringstream error_message;
       error_message 
        << "Failed to cast element1_pt to an ELEMENT" 
        << std::endl;
       throw OomphLibError(error_message.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }

     if (cast_element2_pt==0)
      {
       std::ostringstream error_message;
       error_message 
        << "Failed to cast element2_pt to an ELEMENT"
        << std::endl;
       throw OomphLibError(error_message.str(),
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
#endif

     //Find the dimension of the element
     unsigned dim = cast_element1_pt->dim();
     //Check that the dimensions of bothelements are the same
     {
      unsigned dim2 = cast_element2_pt->dim();
      
      if(dim != dim2)
       {
        std::ostringstream warn_message;
        warn_message
        << "Warning: Two elements do not have the same dimension"
        << dim << " and " << dim2 << ". They cannot be compared\n";
        OomphLibWarning(warn_message.str(),
                        "CompareElementCoordinate::()",
                        OOMPH_EXCEPTION_LOCATION);
       }
     }

     //Find the positions of the centres of the elements
     Vector<double> x1(dim), x2(dim);
     //Not general, only works for quads or octs (centre 0,0,0)
     Vector<double> s(dim,0.0);
     cast_element1_pt->interpolated_x(s,x1);
     cast_element2_pt->interpolated_x(s,x2);

     //This is the Stroustrup-approved way to do lexicographical ordering
     //Loop over the components until they are not equal
     //to within a given tolerance
     {
      unsigned i=0; double tol = 1.0e-14;
      while(i!=dim && (std::abs(x1[i]-x2[i]) < tol)){ ++i;}
      //If we've reached the end, the coordinates are equal, return false
      if(i==dim) {return 0;}
      //Otherwise, return the ordering on the final component
      return x1[i] < x2[i];
     }
    }
  };
 

//=============================================================
/// Class the overloaded RefineablePoissonElements so that the
/// elements that can be rotated about their centre
//=============================================================
template<class ELEMENT>
class Rotateable : 
 public ELEMENT
{
private:

 /// Integer to store the rotation angle
 unsigned Rotated;

public:

 ///Constructor, initialise rotation to NULL (default)
 Rotateable() : RefineableElement(), ELEMENT(), Rotated(0) { }
 
 /// Overload the further build function to pass the rotate information
 ///to the sons
 void further_build()
  {
   ELEMENT::further_build();
   this->Rotated = dynamic_cast<Rotateable<ELEMENT>*>
    (this->father_element_pt())->Rotated;
  }

  
 /// Rotate the element by a given angle:
 /// 0 (0), 1(90 degrees), 2(180 degrees), 3 (270 degrees)
 void rotate(const unsigned &angle)
  {
   Rotated = angle;
   //Get the nodes and make a copy
   unsigned n_node = this->nnode();
   Vector<Node*> elemental_nodes_pt(n_node);
   for(unsigned n=0;n<n_node;n++)
   {
    elemental_nodes_pt[n] = this->node_pt(n);
   }

   //Face rotation four possibilities
   unsigned face_rot = angle%4;
   
   //We now merely permute the nodes
   unsigned n_p = this->nnode_1d();
   //The permutation depends on the angle
   switch(face_rot)
    {
     //FIRST WE JUST ROTATE ABOUT THE Z-AXIS
     //No permutation
    case 0:
     break;
     
     //Rotation by 90 degrees (i and j are swapped and i is reversed)
    case 1:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[j + n_p*(n_p-1-i) + k*n_p*n_p];
          }
        }
      }
     break;

     //Rotation by 180 degrees (i and j are reversed)
    case 2:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[(n_p-1-i) + n_p*(n_p-1-j) + k*n_p*n_p];
          }
        }
      }
     break;

     //Rotation by 270 degrees (i and j are swapped and  j is reversed)
    case 3:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[(n_p-1-j) + n_p*i + k*n_p*n_p];
          }
        }
      }
     break;
    }
 
   //Now get the nodes again, after the first rotation
   for(unsigned n=0;n<n_node;n++)
    {
     elemental_nodes_pt[n] = this->node_pt(n);
    }

   //We now perform the second rotation to move the face into position
   unsigned body_rot = angle/4;
   switch(body_rot)
    {
     //No permutation
    case 0:
     break;

     //ROTATE ABOUT X-AXIS
     //Rotation by 90 degrees (k and j are swapped and k is reversed)
    case 1:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[i + n_p*(n_p-1-k) + j*n_p*n_p];
          }
        }
      }
     break;

     //Rotation by 180 degrees (k and j are reversed)
    case 2:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[i + n_p*(n_p-1-j) + (n_p-1-k)*n_p*n_p];
          }
        }
      }
     break;

     //Rotation by 270 degrees (k and j are swapped and j is reversed)
    case 3:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[i + n_p*k + (n_p-1-j)*n_p*n_p];
          }
        }
      }
     break;
 
     //ROTATE ABOUT Y-AXIS
     //Rotation by 90 degrees (i and k are swapped and i is reversed)
    case 4:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[k + n_p*j + (n_p-1-i)*n_p*n_p];
          }
        }
      }
     break;
     
     //Rotation by 270 degrees (i and k are swapped and k is reversed)
    case 5:
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         for(unsigned k=0;k<n_p;k++)
          {
           this->node_pt(i + j*n_p + k*n_p*n_p)
            = elemental_nodes_pt[(n_p-1-k) + n_p*j + i*n_p*n_p];
          }
        }
      }
     break;
    }  
  }

 ///  Rotate the local coordinates so that that output is
 /// will be consistent, irrespective of the rotation of the element
 void rotate_local_coordinates(Vector<double> &s)
  {
   //Now rotate about the face
   Vector<double> s_rot(3);
   unsigned face_rot = Rotated%4;
   switch(face_rot)
    {
    case 0:
     s_rot[0] = s[0]; s_rot[1] = s[1]; s_rot[2] = s[2];
     break;
     
    case 1:
     s_rot[0] = -s[1]; s_rot[1] = s[0]; s_rot[2] = s[2];
     break;

    case 2:
     s_rot[0] = -s[0]; s_rot[1] = -s[1]; s_rot[2] = s[2];
     break;

    case 3:
     s_rot[0] = s[1]; s_rot[1] = -s[0]; s_rot[2] = s[2];
     break;
    }

   //Rotate the face into position
   unsigned body_rot = Rotated/4;
   Vector<double> s_rot2(3);
   switch(body_rot)
    {
    case 0:
     s_rot2[0] = s_rot[0]; s_rot2[1] = s_rot[1]; s_rot2[2] = s_rot[2];
     break;
     
    case 1:
     s_rot2[0] = s_rot[0]; s_rot2[1] = s_rot[2]; s_rot2[2] = -s_rot[1];
     break;

    case 2:
     s_rot2[0] = s_rot[0]; s_rot2[1] = -s_rot[1]; s_rot2[2] = -s_rot[2];
     break;
     
    case 3:
     s_rot2[0] = s_rot[0]; s_rot2[1] = -s_rot[2]; s_rot2[2] = s_rot[1];
     break;

    case 4:
     s_rot2[0] = -s_rot[2]; s_rot2[1] = s_rot[1]; s_rot2[2] = s_rot[0];
     break;

    case 5:
     s_rot2[0] = s_rot[2]; s_rot2[1] = s_rot[1]; s_rot2[2] = -s_rot[0];
     break;
    }

   //Set the input to the rotated coordinate
   s[0] = s_rot2[0]; s[1] = s_rot2[1]; s[2] = s_rot2[2];
  }

   

///  Output function overloaded to produce identical output
/// under rotation
void output(std::ostream &outfile, 
            const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(3);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points = this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   //Rotate the local coordinates
   this->rotate_local_coordinates(s);
   //Print the output
   for(unsigned i=0;i<3;i++) 
    {
     outfile << this->interpolated_x(s,i) << " ";
    }
   outfile << this->interpolated_u_poisson(s) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);

}

 ///  Output exact solution: Overloaded to produce identical
 /// output under rotation
 void output_fct(std::ostream &outfile, 
                 const unsigned &nplot, 
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
 //Vector of local coordinates
 Vector<double> s(3);
 
 // Vector for coordintes
 Vector<double> x(3);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 // Loop over plot points
 unsigned num_plot_points = this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   //Rotate the local coordinates
   this->rotate_local_coordinates(s);
   
   // Get x position as Vector
   this->interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<3;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << std::endl;  
  }
 
 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}


/// Error computation stuff
/// Overloaded to be identical under rotation
void compute_error(std::ostream &outfile, 
                   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                   double& error, double& norm)
{ 
 // Initialise
 error=0.0;
 norm=0.0;
 
 //Vector of local coordinates
 Vector<double> s(3);
 
 // Vector for coordintes
 Vector<double> x(3);
 
 //Find out how many nodes there are in the element
 unsigned n_node = this->nnode();
 
 Shape psi(n_node);
 
 //Set the value of n_intpt
 unsigned n_intpt = this->integral_pt()->nweight();
  
 // Tecplot 
 outfile << "ZONE" << std::endl;
 
 // Exact solution Vector (here a scalar)
 Vector<double> exact_soln(1);
 
 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   
   //Assign values of s
   for(unsigned i=0;i<3;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Rotate the local coordinates
   this->rotate_local_coordinates(s);
   
   //Get the integral weight
   double w = this->integral_pt()->weight(ipt);
   
   // Get jacobian of mapping
   double J = this->J_eulerian(s);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // Get x position as Vector
   this->interpolated_x(s,x);
   
   // Get FE function value
   double u_fe = this->interpolated_u_poisson(s);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,error
   for(unsigned i=0;i<3;i++)
    {
     outfile << x[i] << " ";
    }
   outfile << exact_soln[0] << " " << exact_soln[0]-u_fe << std::endl;  
   
   // Add to error and norm
   norm+=exact_soln[0]*exact_soln[0]*W;
   error+=(exact_soln[0]-u_fe)*(exact_soln[0]-u_fe)*W;
  }
  
}

};


//==============================start_of_mesh======================
/// Refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
/// The base mesh is 3x3 and of size 1x1. We can then rotate
/// the middle element to produce all possible variations of rotated
/// QuadTrees.
//=================================================================
template<class ELEMENT>
class TestRefineableCubicMesh : 
 public virtual SimpleCubicMesh<ELEMENT>,  
 public RefineableBrickMesh<ELEMENT>
{ 

public: 

 ///  Pass the angle of rotation and the timestepper to
 /// the Mesh.
 TestRefineableCubicMesh(const unsigned &angle,
                         TimeStepper* time_stepper_pt=
                         &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(3,3,3,1.0,1.0,1.0,time_stepper_pt)
  {

   //To Check rotations, we shall rotate the central element
   dynamic_cast<ELEMENT*>(this->finite_element_pt(13))->rotate(angle);
   
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_octree_forest();

  } // end of constructor
 

 /// Sort the elements into lexographical order
 void sort_elements()
  {
   std::sort(this->element_pt().begin(),this->element_pt().end(),
             CompareElementCoordinate<ELEMENT>());
  }

 ///  Destructor: Empty -- all cleanup gets handled in the base
 /// classes
 virtual ~TestRefineableCubicMesh() {}

}; // end of mesh



//=============start_of_namespace=====================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of step
 double Alpha=50.0;

 /// Orientation (non-normalised x-component of unit vector in direction
 /// of step plane)
 double N_x=-1.0;

 /// Orientation (non-normalised y-component of unit vector in direction
 /// of step plane)
 double N_y=-1.0;

 /// Orientation (non-normalised z-component of unit vector in direction
 /// of step plane)
 double N_z=1.0;


 /// Orientation (x-coordinate of step plane) 
 double X_0=0.0;

 /// Orientation (y-coordinate of step plane) 
 double Y_0=0.0;

 /// Orientation (z-coordinate of step plane) 
 double Z_0=0.0;


 // Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
                     N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
                     N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  u = tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-Y_0)*
                     N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*
                     N_z/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)));
 }


 /// Source function to make it an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {

  double s1,s2,s3,s4;

  s1 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_x*N_x/(N_x*N_x+N_y*N_y+N_z*N_z);
      s3 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_y*N_y/(N_x*N_x+N_y*N_y+N_z*N_z);
      s4 = -2.0*tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z)))*(1.0-pow(tanh(Alpha*((x[0]-X_0)*N_x/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[1]-
Y_0)*N_y/sqrt(N_x*N_x+N_y*N_y+N_z*N_z)+(x[2]-Z_0)*N_z/sqrt(N_x*N_x+N_y*N_y+N_z*
N_z))),2.0))*Alpha*Alpha*N_z*N_z/(N_x*N_x+N_y*N_y+N_z*N_z);
      s2 = s3+s4;
      source = s1+s2;
 }


} // end of namespace




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_class_definition====================================
/// Poisson problem in refineable eighth of a sphere mesh.
//====================================================================
template<class ELEMENT>
class TestPoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 TestPoissonProblem(const unsigned &angle,
  PoissonEquations<3>::PoissonSourceFctPt source_fct_pt);

 /// Destructor: Empty
 ~TestPoissonProblem(){}

 ///  Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 TestRefineableCubicMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<TestRefineableCubicMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()  {}

 ///  Update the problem specs before solve: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_newton_solve()
 {
  //Loop over the boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
   {
    // Loop over the nodes on boundary
    unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
     double u;
     Vector<double> x(3);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     x[2]=nod_pt->x(2);
     TanhSolnForPoisson::get_exact_u(x,u);
     nod_pt->set_value(0,u);
    }
   }
 }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:

 /// Pointer to source function
 PoissonEquations<3>::PoissonSourceFctPt Source_fct_pt;

}; // end of class definition





//====================start_of_constructor================================
/// Constructor for Poisson problem on eighth of a sphere mesh
//========================================================================
template<class ELEMENT>
TestPoissonProblem<ELEMENT>::TestPoissonProblem(const unsigned &angle,
 PoissonEquations<2>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 //Create the mesh
 Problem::mesh_pt() = new TestRefineableCubicMesh<ELEMENT>(angle);

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 //Default parmaters
 mesh_pt()->max_permitted_error()=0.01;
 mesh_pt()->min_permitted_error()=0.001;

 //Doc the mesh boundaries
 //ofstream some_file;
 //some_file.open("boundaries.dat");
 //mesh_pt()->output_boundaries(some_file);
 //some_file.close();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here (all the nodes on the boundary)
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  } // end of pinning


 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }


 //Manually adapt the mesh to our solution
 Vector<unsigned> to_be_refined(10);
 to_be_refined[0] = 0;
 to_be_refined[1] = 1;
 to_be_refined[2] = 3;
 to_be_refined[3] = 9;
 to_be_refined[4] = 10;
 to_be_refined[5] = 12;
 to_be_refined[6] = 13;
 to_be_refined[7] = 18;
 to_be_refined[8] = 22;
 to_be_refined[9] = 26; 

 

 mesh_pt()->refine_selected_elements(to_be_refined);

 // Setup equation numbering 
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TestPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 


 // Output solution 
 //-----------------
 mesh_pt()->sort_elements();
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
//  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
//  some_file.close();


//  // Doc error
//  //----------
//  double error,norm;
//  sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
//          doc_info.number());
//  some_file.open(filename);
//  mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
//                           error,norm); 
//  some_file.close();
//  cout << "error: " << sqrt(error) << std::endl; 
//  cout << "norm : " << sqrt(norm) << std::endl << std::endl;

} // end of doc

///Global function that is used to run the rotation tests for
///different elements
template<class ELEMENT>
void run_test(const unsigned &i)
{
 char dirname[100];
 sprintf(dirname,"RESLT%i",i);
 
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory(dirname);
 
 // Step number
 doc_info.number()=0;
 
 //Set up the problem
 //------------------
 for(unsigned i=0;i<24;i++)
  {
   //Create the problem with rotation i 
   TestPoissonProblem<ELEMENT> 
    problem(i,&TanhSolnForPoisson::get_source);
   
   // Solve the problem, performing up to 2 adaptive refinements
   problem.newton_solve();
   
   //Output the solution
   problem.doc_solution(doc_info);
   
   doc_info.number()++;
  }
}

//===== start_of_main=====================================================
/// Driver code for 2D Poisson problem
//========================================================================
int main()
{
 run_test<Rotateable<RefineableQPoissonElement<3,2> > >(1);
 run_test<Rotateable<RefineableQPoissonElement<3,3> > >(2);
 run_test<Rotateable<RefineableQSpectralPoissonElement<3,4> > >(3);
} //end of main







