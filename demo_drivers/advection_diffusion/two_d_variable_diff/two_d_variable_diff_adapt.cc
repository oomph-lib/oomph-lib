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
//Driver for a simple 2D adv diff problem with adaptive mesh refinement

//Generic routines
#include "generic.h"

// The Poisson equations
#include "advection_diffusion.h"

// The mesh 
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//======start_of_namespace============================================
/// Namespace for exact solution for AdvectionDiffusion equation 
/// with "sharp" step 
//====================================================================
namespace RachelsAdvectionDiffusion
{

 /// Peclet number
 double Peclet=0.01;

 /// Swimming speed
 double U = 0.7;

 //Diagonal components of the diffusivity tensor 
 //(used in the asympotic solution)
 double D11 = 0.1;
 double D22 = 0.2;

 /// The value of pi
 const double pi = MathematicalConstants::Pi;

 /// The asymptotic solution
 void get_exact_u(const Vector<double> &x, Vector<double> &u)
 {
  double gamma = U/D22;
  //Leading order term
  u[0] = std::exp(gamma*x[1]);

  //Get the surd in the solution of the ODE
  double surd = sqrt(0.25*gamma*gamma + pi*pi*D11/D22);
  //Powers of the exponentials
  double alpha = 0.5*gamma + surd;
  double beta  = 0.5*gamma - surd;

  //Coefficients of the exact solution
  double Gamma = 2.0*D22*pi*gamma - U*pi;
  double Kappa = D22*(gamma*gamma - pi*pi) - U*gamma - D11*pi*pi;
  
  double temp = pi*gamma/(Gamma*Gamma + Kappa*Kappa);

  double A = - Kappa*temp;
  double B =  Gamma*temp;

  //Now the final messy bit
  double C = (gamma - alpha)*exp(alpha);
  double D = (gamma -  beta)*exp(beta);

  double b = (pi*A*exp(gamma) - C*B)/(C - D);
  double a = -B -b;

  //Done, I hope
  double f= a*exp(alpha*x[1]) + b*exp(beta*x[1]) + 
   (A*sin(pi*x[1]) + B*cos(pi*x[1]))*exp(gamma*x[1]);

  /*double fprime = a*alpha*exp(alpha*x[1]) + b*beta*exp(beta*x[1])
   + ((pi*A + B*gamma)*cos(pi*x[1]) + (gamma*A- pi*B)*sin(pi*x[1]))
   *exp(gamma*x[1]);

  double fpp = a*alpha*alpha*exp(alpha*x[1]) + b*beta*beta*exp(beta*x[1])
   + (-pi*(pi*A + B*gamma)*sin(pi*x[1])
      + pi*(gamma*A - pi*B)*cos(pi*x[1]))*exp(gamma*x[1])
   + ((pi*A + B*gamma)*cos(pi*x[1]) + (gamma*A- pi*B)*sin(pi*x[1]))
   *gamma*exp(gamma*x[1]);*/

  //Add the first-order term
  u[0] += Peclet*cos(pi*x[0])*f;
 }

 /// Wind (the underlying velocity field)
 void wind_function(const Vector<double>& x, Vector<double>& wind)
 {
  wind[0]=  pi*cos(pi*x[1])*sin(pi*x[0]);
  wind[1]= -pi*sin(pi*x[1])*cos(pi*x[0]);
 }

 /// Conserved bit (the swiming velocity)
 void swimming(const Vector<double> &x, Vector<double> &swim)
 {
  swim[0] = 0.0;
  swim[1] = U;
 }
 
 /// Diffusivity tensor
 void diff_function(const Vector<double> &x, DenseMatrix<double> &D)
 {
  D(0,0) = D11;
  D(0,1) = 0.0;
  D(1,0) = 0.0;
  D(1,1) = D22;
 }  

} // end of namespace

/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////


//====== start_of_problem_class=======================================
/// 2D GeneralisedAdvectionDiffusion problem on a 
/// rectangular domain, discretised 
/// with refineable 2D QGeneralisedAdvectionDiffusion elements. 
/// The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class RefineableAdvectionDiffusionProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source, wind and diffusivity functions
 RefineableAdvectionDiffusionProblem(
  GeneralisedAdvectionDiffusionEquations<2>::
  GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt,
  GeneralisedAdvectionDiffusionEquations<2>::
  GeneralisedAdvectionDiffusionWindFctPt swim_fct_pt,
  GeneralisedAdvectionDiffusionEquations<2>::
  GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt);

 /// Destructor. Delete the mesh and error estimator
 ~RefineableAdvectionDiffusionProblem()
  {
   delete this->mesh_pt();
  }

 /// Update the problem specs before solve
 void actions_before_newton_solve() {} 

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt: Document the solution
 void actions_before_adapt()
  {
   // Doc the solution
   doc_solution();
   
   // Increment label for output files
   Doc_info.number()++;
  }

 /// Doc the solution.
 void doc_solution();

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

private:

 /// DocInfo object
 DocInfo Doc_info;

 /// Pointer to diffusivity function
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionDiffFctPt Diff_fct_pt;

 /// Pointer to wind function
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionWindFctPt Wind_fct_pt;

 /// Pointer to a swimming speed
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionWindFctPt Swim_fct_pt;
 

}; // end of problem class



//=====start_of_constructor===============================================
/// Constructor for AdvectionDiffusion problem: Pass pointer to 
/// source function.
//========================================================================
template<class ELEMENT>
RefineableAdvectionDiffusionProblem<ELEMENT>::RefineableAdvectionDiffusionProblem(
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt,
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionWindFctPt swim_fct_pt,
 GeneralisedAdvectionDiffusionEquations<2>::
 GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt)
 :  Diff_fct_pt(diff_fct_pt), Wind_fct_pt(wind_fct_pt),
    Swim_fct_pt(swim_fct_pt)
{ 

 // Set output directory
 Doc_info.set_directory("RESLT");

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet 
 // conditions here. That is those on the lower wall
 for(unsigned ibound=0;ibound<1;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
     mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,1.0); 
    }
  } // end loop over boundaries
   
 // Complete the build of all elements so they are fully functional 

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by the (argument-free!) ELEMENT 
 // constructor: Pass pointer to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the diffusivity function pointer
   el_pt->diff_fct_pt() = Diff_fct_pt;

   //Set the wind function pointer
   el_pt->wind_fct_pt() = Wind_fct_pt;

   //Set the swimming speed
   el_pt->conserved_wind_fct_pt() = Swim_fct_pt;

   // Set the Peclet number
   el_pt->pe_pt() = &RachelsAdvectionDiffusion::Peclet;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//===============start_of_doc=============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableAdvectionDiffusionProblem<ELEMENT>::doc_solution()
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 
 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,RachelsAdvectionDiffusion::get_exact_u); 
 some_file.close();


  // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,RachelsAdvectionDiffusion::get_exact_u,
                          error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

} // end of doc


//===== start_of_main=====================================================
/// Driver code for 2D AdvectionDiffusion problem
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 // Create the problem with 2D nine-node refineable elements from the
 // RefineableQuadAdvectionDiffusionElement family. Pass pointer to 
 // source and wind function. 
 RefineableAdvectionDiffusionProblem<
  RefineableQGeneralisedAdvectionDiffusionElement<2,
  3> > problem(&RachelsAdvectionDiffusion::wind_function,
               &RachelsAdvectionDiffusion::swimming,
               &RachelsAdvectionDiffusion::diff_function);
 
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

 // Solve the problem, performing up to 4 adptive refinements
 problem.newton_solve(4);

 //Output the solution
 problem.doc_solution();
 
} // end of main









