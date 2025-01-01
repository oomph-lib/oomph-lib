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

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

using namespace oomph;

using namespace std;

namespace oomph
{


// //===========================================================================
// /// Class for sorting the elements into lexicographic order based
// /// on the coordinate at the centre of the element. This is the only
// /// "rational" way to preserve the ordering given rotations of the elements.
// //===========================================================================
// template<class ELEMENT>
// class CompareElementCoordinate 
// {
// public:
 
//  ///The actual comparison operator
//  int operator() (GeneralisedElement* const &element1_pt,
//                  GeneralisedElement* const &element2_pt)
//   {
   
//    //Dynamic cast the elements
//    ELEMENT* cast_element1_pt = dynamic_cast<ELEMENT*>(element1_pt);
//    ELEMENT* cast_element2_pt = dynamic_cast<ELEMENT*>(element2_pt);
   
//    //Check that we managed to successfully cast the elements
// #ifdef PARANOID
//      if (cast_element1_pt==0)
//       {
//        std::ostringstream error_message;
//        error_message 
//         << "Failed to cast element1_pt to an ELEMENT" 
//         << std::endl;
//        throw OomphLibError(error_message.str(),
//                            OOMPH_CURRENT_FUNCTION,
//                            OOMPH_EXCEPTION_LOCATION);
//       }

//      if (cast_element2_pt==0)
//       {
//        std::ostringstream error_message;
//        error_message 
//         << "Failed to cast element2_pt to an ELEMENT"
//         << std::endl;
//        throw OomphLibError(error_message.str(),
//                            OOMPH_CURRENT_FUNCTION,
//                            OOMPH_EXCEPTION_LOCATION);
//       }
// #endif

//      //Find the dimension of the element
//      unsigned dim = cast_element1_pt->dim();
//      //Check that the dimensions of bothelements are the same
//      {
//       unsigned dim2 = cast_element2_pt->dim();
      
//       if(dim != dim2)
//        {
//         std::ostringstream warn_message;
//         warn_message
//         << "Warning: Two elements do not have the same dimension"
//         << dim << " and " << dim2 << ". They cannot be compared\n";
//         OomphLibWarning(warn_message.str(),
//                         "CompareElementCoordinate::()",
//                         OOMPH_EXCEPTION_LOCATION);
//        }
//      }

//      //Find the positions of the centres of the elements
//      Vector<double> x1(dim), x2(dim);
//      //Not general, only works for quads or octs (centre 0,0,0)
//      Vector<double> s(dim,0.0);
//      cast_element1_pt->interpolated_x(s,x1);
//      cast_element2_pt->interpolated_x(s,x2);

//      //This is the Stroustrup-approved way to do lexicographical ordering
//      //Loop over the components until they are not equal
//      //to within a given tolerance
//      {
//       unsigned i=0; double tol = 1.0e-14;
//       while(i!=dim && (std::abs(x1[i]-x2[i]) < tol)){ ++i;}
//       //If we've reached the end, the coordinates are equal, return false
//       if(i==dim) {return 0;}
//       //Otherwise, return the ordering on the final component
//       return x1[i] < x2[i];
//      }
//     }
//   };
 

//=============================================================
/// Class to overload an refineable QElement<3,...> so that it
/// can be rotated about its centre.
//=============================================================
template<class ELEMENT>
class Rotateable : 
 public ELEMENT
{
private:

 /// Integer to store the rotation angle
 unsigned Rotated;

public:

 /// Constructor, initialise rotation to zero (default)
 Rotateable() : RefineableElement(), ELEMENT(), Rotated(0) { }
 

 /// Overload the further build function to pass the rotate information
 /// to the sons
 void further_build()
  {
   ELEMENT::further_build();
   this->Rotated = dynamic_cast<Rotateable<ELEMENT>*>
    (this->father_element_pt())->Rotated;
  }

  
 /// Rotate the element
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

 /// Rotate the local coordinates so that that output is
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

};


//==============================start_of_mesh======================
/// 3x3x3 SimpleCubicMesh upgraded to become refineable
/// and with the ability to rotate one of its elements.
//=================================================================
template<class ELEMENT>
class RotatableRefineableCubicMesh : 
 public virtual SimpleCubicMesh<ELEMENT>,  
 public RefineableBrickMesh<ELEMENT>
{ 

public: 

 /// Constructor: Pass the rotation index, the number of the 
 /// rotated element and the timestepper to the Mesh.
 RotatableRefineableCubicMesh(const unsigned &angle, 
                              const unsigned& rotated_element,
                              TimeStepper* time_stepper_pt=
                              &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(3,3,3,1.0,1.0,1.0,time_stepper_pt)
  {

   // Sanity check
   if(rotated_element>this->nelement())
    {
     std::ostringstream error_stream;
     error_stream << "rotated_element>this->nelement(): " 
                  << rotated_element << " " << this->nelement()
                  << std::endl;
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }

   //To check rotations, we rotate the specified element
   dynamic_cast<ELEMENT*>(this->finite_element_pt(rotated_element))
    ->rotate(angle);
   
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_octree_forest();

  } // end of constructor
 

 /// Destructor: Empty
 virtual ~RotatableRefineableCubicMesh(){}

}; // end of mesh


}


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




/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_class_definition====================================
/// Poisson problem on rotatable, refineable cubic brick mesh
//====================================================================
template<class ELEMENT>
class TestPoissonProblem : public Problem
{

public:

 /// Constructor: Pass index of rotation, number of rotated
 /// element, and pointer to source function. 
 TestPoissonProblem(const unsigned &angle, const unsigned& rotated_element,
  PoissonEquations<3>::PoissonSourceFctPt source_fct_pt);

 /// Destructor: close trace file
 ~TestPoissonProblem()
  {
   Trace_file.close();
  }

 /// Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 RotatableRefineableCubicMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RotatableRefineableCubicMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Impose specific refinement on mesh (up to nrefine refinements)
 void impose_specific_refinement(const unsigned& nrefine);

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve: 
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

 /// Trace file
 ofstream Trace_file;

}; // end of class definition





//====================start_of_constructor================================
/// Constructor for for test Poisson problem.
//========================================================================
template<class ELEMENT>
TestPoissonProblem<ELEMENT>::TestPoissonProblem(
 const unsigned &angle,
 const unsigned& rotated_element,
 PoissonEquations<2>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 


 // Open trace file (append)
 Trace_file.open("trace.dat",std::ios_base::app);
 Trace_file << angle << " " << rotated_element << " ";


 //Create the mesh
 Problem::mesh_pt() = 
  new RotatableRefineableCubicMesh<ELEMENT>(angle,rotated_element);

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 //Default parmaters
 mesh_pt()->max_permitted_error()=0.01;
 mesh_pt()->min_permitted_error()=0.001;

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

 // Document the mesh adaptation
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 mesh_pt()->doc_info_pt()=&doc_info;

 // Impose the specific refinement
 impose_specific_refinement(2);

 // Setup equation numbering 
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor








//========================================================================
/// Impose a specific refinement pattern with up to 3 non-uniform
/// refinements that produces all sorts of neighbours. Triply
/// refined next to unrefined, double refined next to triply refined
/// etc.
//========================================================================
template<class ELEMENT>
void TestPoissonProblem<ELEMENT>::impose_specific_refinement(
 const unsigned& nrefine)
{

 // First four block elements get triply refined
 set<RefineableElement*> triply_refined_el_pt;
 if (nrefine>2)
  {
   triply_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(7)));
   triply_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(8)));
   triply_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(16)));
   triply_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(17)));
  }


 // First four block elements get doubly refined
 set<RefineableElement*> doubly_refined_el_pt;
 if (nrefine>1)
  {
   doubly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(0)));
   doubly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(1)));
   doubly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(3)));
   doubly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                                mesh_pt()->element_pt(4)));
  }


 // Next four block elements get singly refined
 set<RefineableElement*> singly_refined_el_pt;
 singly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                              mesh_pt()->element_pt(2)));
 singly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                              mesh_pt()->element_pt(5)));
 singly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                              mesh_pt()->element_pt(11)));
 singly_refined_el_pt.insert(dynamic_cast<RefineableElement*>(
                              mesh_pt()->element_pt(14)));

 // Loop over refinement levels
 Vector<RefineableElement*> ref_el_pt;
 for (unsigned i=0;i<nrefine;i++)
  {
   
   // Form Vector of elements that have to be refined
   ref_el_pt.resize(0);
   typedef set<RefineableElement*>::iterator IT;
   
   // Triple refinement
   for (IT it=triply_refined_el_pt.begin();
        it!=triply_refined_el_pt.end();
        ++it)
    {
     ref_el_pt.push_back(*it);
    }
   
   
   // Double refinement
   if (i<2)
    {
     for (IT it=doubly_refined_el_pt.begin();
          it!=doubly_refined_el_pt.end();
          ++it)
      {
       ref_el_pt.push_back(*it);
      }
    }
   
   // Only add singly refined ones during first round
   if (i==0)
    {
     for (IT it=singly_refined_el_pt.begin();
          it!=singly_refined_el_pt.end();
          ++it)
      {
       ref_el_pt.push_back(*it);
      }
     
     // Add central element
     ref_el_pt.push_back(dynamic_cast<RefineableElement*>(
                          mesh_pt()->element_pt(13)));
    }
   
   // Total number of elements to be refined at least once
   unsigned n_total=ref_el_pt.size();
   
   
   // Refine these elements
   refine_selected_elements(ref_el_pt);
   
   
   // Now that they've been refined, remove them from the set
   // of elements that have to be refined next time but add
   // their sons instead
   Vector<RefineableQElement<3>*> aux_el_pt;
   for (IT it=triply_refined_el_pt.begin();
          it!=triply_refined_el_pt.end();
          ++it)
    {
     aux_el_pt.push_back(dynamic_cast<RefineableQElement<3>*>(*it));
    }

   unsigned nr=aux_el_pt.size();
   for (unsigned ii=0;ii<nr;ii++)
    {
     // Add sons
     for (unsigned s=0;s<8;s++)
      {
       triply_refined_el_pt.insert(aux_el_pt[ii]->octree_pt()->
                                   son_pt(s)->object_pt());
      }
     // Remove element itself
     triply_refined_el_pt.erase(aux_el_pt[ii]);
    }
   
 
   
   if (i<2)
    {

     // Now that they've been refined, remove them from the set
     // of elements that have to be refined next time but add
     // their sons instead
     aux_el_pt.resize(0);
     for (IT it=doubly_refined_el_pt.begin();
          it!=doubly_refined_el_pt.end();
          ++it)
      {
       aux_el_pt.push_back(dynamic_cast<RefineableQElement<3>*>(*it));
      }
     
     nr=aux_el_pt.size();
     for (unsigned ii=0;ii<nr;ii++)
      {
       // Add sons
       for (unsigned s=0;s<8;s++)
        {
         doubly_refined_el_pt.insert(aux_el_pt[ii]->octree_pt()->
                                     son_pt(s)->object_pt());
        }
       // Remove element itself
       doubly_refined_el_pt.erase(aux_el_pt[ii]);
      }
     // Finally: Last son in central element
     doubly_refined_el_pt.insert(dynamic_cast<RefineableQElement<3>*>(
                                  ref_el_pt[n_total-1])->octree_pt()->
                                 son_pt(7)->object_pt());
    }   
  }
}




//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TestPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 ofstream some_file,some_file2;
 char filename[100];

 

 unsigned success_flag=0;

 // Check for repeated nodes in mesh
 if (mesh_pt()->check_for_repeated_nodes()!=0) success_flag=1;
  
 // Do octree self test
 if (dynamic_cast<OcTreeForest*>(mesh_pt()->forest_pt())->self_test()!=0)
  {
   success_flag=1;
  }

 // Doc success flag to trace file
 Trace_file  << success_flag << std::endl; 
 
 
 // Limited output if run as self test
 if ((CommandLineArgs::Argc==1)||(doc_info.number()==0))
  {
   // Number of plot points
   unsigned npts;
   npts=5; 
   
   // Doc orientation of the root elements
   sprintf(filename,"%s/orientation%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
  
   sprintf(filename,"%s/root_elements%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file2.open(filename);
  
   // Loop over octree roots to display orientation
   Vector<double> s(3);
   Vector<double> x_centre(3),x_end(3);
   unsigned ntree=mesh_pt()->forest_pt()->ntree();
   for (unsigned i=0;i<ntree;i++)
    {
     RefineableQElement<3>* el_pt=
      dynamic_cast<RefineableQElement<3>*>(
       mesh_pt()->forest_pt()->tree_pt(i)->object_pt());  
    
     // Centre of element
     s[0]=0.0;
     s[1]=0.0;
     s[2]=0.0;
     el_pt->interpolated_x(s,x_centre);
    
    
     // End of "up" in  element
     s[0]=0.0;
     s[1]=1.0;
     s[2]=0.0;
     el_pt->interpolated_x(s,x_end);
    
     some_file << "ZONE I=1, C=\"RED\"" << std::endl;
     some_file << x_centre[0] << " " 
               << x_centre[1] << " " 
               << x_centre[2] << " " 
               << x_end[0]-x_centre[0] << " " 
               << x_end[1]-x_centre[1] << " " 
               << x_end[2]-x_centre[2] << " " << std::endl;
    
     // End of "right" in  element
     s[0]=1.0;
     s[1]=0.0;
     s[2]=0.0;
     el_pt->interpolated_x(s,x_end);
    
     some_file << "ZONE I=1, C=\"GREEN\"" << std::endl;
     some_file << x_centre[0] << " " 
               << x_centre[1] << " " 
               << x_centre[2] << " " 
               << x_end[0]-x_centre[0] << " " 
               << x_end[1]-x_centre[1] << " " 
               << x_end[2]-x_centre[2] << " " << std::endl;
    
     el_pt->output_corners(some_file2,"BLUE");
    }
   some_file.close();
   some_file2.close();
  
   // Doc neighbours
   doc_info.enable_doc();
   mesh_pt()->forest_pt()->check_all_neighbours(doc_info);
  
  
   // Output solution 
   //-----------------
   //mesh_pt()->sort_elements();
   sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   mesh_pt()->output(some_file,npts);
   some_file.close();
  
 }

} // end of doc







//===== start_of_main=====================================================
/// Driver code 
//========================================================================
int main(int argc, char* argv[])
{

 // Save command line arguments
 CommandLineArgs::setup(argc,argv);

 // Write header to trace file and close (also wipes any existing ones
 ofstream trace_file("trace.dat");
 trace_file << "rotation_index rotated_element success_flag" << std::endl;
 trace_file.close();

 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;
 
 //Set up the problems and doc
 //---------------------------

 // Loop over rotations
 for(unsigned r=1;r<24;r++)
  {
   // Loop over elements that are rotated
   for(unsigned rotated_element=0;rotated_element<27;rotated_element++)
    {
     
     std::cout << "Case: r, rotated_element" 
               << r << " " << rotated_element << std::endl;
     
    //Create the problem with rotation r for element rotated_element
    TestPoissonProblem<Rotateable<RefineableQPoissonElement<3,2> > > 
     problem(r,rotated_element,&TanhSolnForPoisson::get_source);
    
    //Output the solution
    problem.doc_solution(doc_info);
    
    doc_info.number()++;
    
    // We're not interested in the solution...
//    problem.newton_solve();

   }
  }
}







