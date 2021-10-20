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
#include "meshes/simple_rectangular_quadmesh.h"

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
   //Get the nodes and make a copy
   unsigned n_node = this->nnode();
   Vector<Node*> elemental_nodes_pt(n_node);
   for(unsigned n=0;n<n_node;n++)
   {
    elemental_nodes_pt[n] = this->node_pt(n);
   }
   
   //We now merely permute the nodes
   unsigned n_p = this->nnode_1d();
   //The permutation depends on the angle
   switch(angle)
    {
     //No permutation
    case 0:
     Rotated = 0;
     break;
     
     //Rotation by 90 degrees (i and j are swapped and i is reversed)
    case 1:
     Rotated = 1;
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         this->node_pt(i + j*n_p)
          = elemental_nodes_pt[j + n_p*(n_p-1-i)];
        }
      }
     break;

     //Rotation by 180 degrees (i and j are reversed)
    case 2:
     Rotated = 2;
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         this->node_pt(i + j*n_p)
          = elemental_nodes_pt[(n_p-1-i) + n_p*(n_p-1-j)];
        }
      }
     break;

     //Rotation by 270 degrees (i and j are swapped and new j is reversed)
    case 3:
     Rotated = 3;
     for(unsigned i=0;i<n_p;i++)
      {
       for(unsigned j=0;j<n_p;j++)
        {
         this->node_pt(i + j*n_p)
          = elemental_nodes_pt[(n_p-1-j) + n_p*i];
        }
      }
     break;

    default:
     Rotated = 0;
    }
  }
   
 
 ///  Rotate the local coordinates so that that output is
 /// will be consistent, irrespective of the rotation of the element
 void rotate_local_coordinates(Vector<double> &s)
  {
   //Rotated coordinates
   Vector<double> s_rot(2);
   //Do the four different cases
   switch(Rotated)
    {
    case 0:
     s_rot[0] = s[0]; s_rot[1] = s[1];
     break;
     
    case 1:
     s_rot[0] = -s[1]; s_rot[1] = s[0];
     break;

    case 2:
     s_rot[0] = -s[0]; s_rot[1] = -s[1];
     break;

    case 3:
     s_rot[0] = s[1]; s_rot[1] = -s[0];
     break;
    }
   
   //Now do the rotation
   s[0] = s_rot[0]; s[1] = s_rot[1];
  }

 

///  Output function overloaded to produce identical output
/// under rotation
void output(std::ostream &outfile, 
            const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points = this->nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   //Rotate the coordinates
   this->rotate_local_coordinates(s);

   for(unsigned i=0;i<2;i++) 
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
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
 
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
   //Rotate the coordinates
   this->rotate_local_coordinates(s);
 
   // Get x position as Vector
   this->interpolated_x(s,x);
   
   // Get exact solution at this point
   (*exact_soln_pt)(x,exact_soln);
   
   //Output x,y,...,u_exact
   for(unsigned i=0;i<2;i++)
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
 Vector<double> s(2);
 
 // Vector for coordintes
 Vector<double> x(2);
 
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
   for(unsigned i=0;i<2;i++)
    {
     s[i] = this->integral_pt()->knot(ipt,i);
    }

   //Rotate the coordinates
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
   for(unsigned i=0;i<2;i++)
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
class TestRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 ///  Pass the angle of rotation and the timestepper to
 /// the Mesh.
 TestRefineableRectangularQuadMesh(const unsigned &angle,
                                   TimeStepper* time_stepper_pt=
                                   &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(3,3,1.0,1.0,time_stepper_pt)
  {

   //To Check rotations, we shall rotate the central element
   dynamic_cast<ELEMENT*>(this->finite_element_pt(4))->rotate(angle);
   
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();

  } // end of constructor
 

 /// Sort the elements into lexographical order
 void sort_elements()
  {
   std::sort(this->element_pt().begin(),this->element_pt().end(),
             CompareElementCoordinate<ELEMENT>());
  }

 ///  Destructor: Empty -- all cleanup gets handled in the base
 /// classes
 virtual ~TestRefineableRectangularQuadMesh() {}

}; // end of mesh



//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=50.0;

 /// Parameter for angle Phi of "step" (45 degrees)
 double TanPhi=1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }
 
} // end of namespace


//====== start_of_problem_class=======================================
/// 2D Poisson problem on rectangular domain, discretised with
/// refineable 2D QPoisson elements. The specific type of element is
/// specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class TestRefineablePoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 TestRefineablePoissonProblem(const unsigned &angle,
PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty -- all the cleanup is done in base class)
 ~TestRefineablePoissonProblem(){};

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 ///  Overloaded version of the Problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 TestRefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<TestRefineableRectangularQuadMesh<ELEMENT>*>
    (Problem::mesh_pt());
  }

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class




//=====start_of_constructor===============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
TestRefineablePoissonProblem<ELEMENT>::
TestRefineablePoissonProblem(
 const unsigned &angle,
 PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
 :  Source_fct_pt(source_fct_pt)
{ 

 // Build and assign mesh
 Problem::mesh_pt() = 
  new TestRefineableRectangularQuadMesh<ELEMENT>(angle);
 
 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet conditions
 // here. 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//========================================start_of_actions_before_newton_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void TestRefineablePoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned num_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // How many nodes are there on this boundary?
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);

   // Loop over the nodes on boundary
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     // Compute the value of the exact solution at the nodal point
     Vector<double> u(1);
     TanhSolnForPoisson::get_exact_u(x,u);

     // Assign the value to the one (and only) nodal value at this node
     nod_pt->set_value(0,u[0]);
    }
  } 
}  // end of actions before solve



//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void TestRefineablePoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output the nodes
 //sprintf(filename,"%s/nodes%i.dat",doc_info.directory().c_str(),
 //        doc_info.number());
 //some_file.open(filename);
 // {
 // unsigned n_node = mesh_pt()->nnode();
 // for(unsigned n=0;n<n_node;n++)
 //  {
 //   some_file << mesh_pt()->node_pt(n)->x(0) << " "
 //             << mesh_pt()->node_pt(n)->x(1) << " "
 //             << mesh_pt()->node_pt(n)->value(0) << std::endl;
 //  }
 //}
 //some_file.close();

 // Output solution 
 //-----------------
 this->mesh_pt()->sort_elements();
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();


 // Output exact solution 
 //----------------------
 //sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
 //        doc_info.number());
 //some_file.open(filename);
 //mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 //some_file.close();

 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 //double error,norm;
 //sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
 //        doc_info.number());
 //some_file.open(filename);
 //mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
 //                         error,norm); 
 //some_file.close();

 // Doc L2 error and norm of solution
 //cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 //cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

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
 for(unsigned i=0;i<4;i++)
  {
   //Create the problem with rotation i 
   TestRefineablePoissonProblem<ELEMENT> 
    problem(i,&TanhSolnForPoisson::get_source);
   
   //Attach the info to the mesh
   //problem.mesh_pt()->doc_info_pt() = &doc_info;

   
   // Check if we're ready to go:
   //----------------------------
   cout << "\n\n\nProblem self-test ";
   if (problem.self_test()==0) 
    {
     cout << "passed: Problem can be solved." << std::endl;
    }
     else 
      {
       throw OomphLibError("Self test failed",
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
   
   // Solve the problem, performing up to 3 adaptive refinements
   problem.newton_solve(3);
   
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
 //Run the test three times for three different types of element
 run_test<Rotateable<RefineableQPoissonElement<2,2> > >(1);
 run_test<Rotateable<RefineableQPoissonElement<2,3> > >(2);
 run_test<Rotateable<RefineableQSpectralPoissonElement<2,4> > >(3);
} 









