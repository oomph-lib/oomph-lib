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
//Driver for a 3D poisson problem

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/eighth_sphere_mesh.h"

using namespace std;

using namespace oomph;

//=============start_of_namespace=====================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of step
 double Alpha=1;

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
class EighthSpherePoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 EighthSpherePoissonProblem(
  PoissonEquations<3>::PoissonSourceFctPt source_fct_pt);

 /// Destructor: Empty
 ~EighthSpherePoissonProblem(){}

 ///  Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 RefineableEighthSphereMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableEighthSphereMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 void actions_after_adapt()
  {
   //Doc the mesh boundaries
   ofstream some_file;
   some_file.open("boundaries_end.dat");
   mesh_pt()->output_boundaries(some_file);
   some_file.close();
  }

 void actions_before_adapt()  {}

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
     // Pin
     nod_pt->pin(0);
     // Set value
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
EighthSpherePoissonProblem<ELEMENT>::EighthSpherePoissonProblem(
   PoissonEquations<3>::PoissonSourceFctPt source_fct_pt) : 
         Source_fct_pt(source_fct_pt)
{

 // Change solver to CG
 IterativeLinearSolver* solver_pt = new CG<CRDoubleMatrix>;
 linear_solver_pt()=solver_pt;
 
 //// Specify preconditioner
 //solver_pt->preconditioner_pt() = new ILUZeroPreconditioner<CRDoubleMatrix>;

 // Setup parameters for exact tanh solution
 // Steepness of step
 TanhSolnForPoisson::Alpha=10.0;

 /// Create mesh for sphere of radius 5
 double radius=5.0; 
 Problem::mesh_pt() = new RefineableEighthSphereMesh<ELEMENT>(radius);

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Adjust error targets for adaptive refinement
 if (CommandLineArgs::Argc>1)
  {
   // Validation: Relax tolerance to get nonuniform refinement during
   // first step
   mesh_pt()->max_permitted_error()=0.1;
   mesh_pt()->min_permitted_error()=0.01;
  }
 else
  {
   mesh_pt()->max_permitted_error()=0.01;
   mesh_pt()->min_permitted_error()=0.001;
  } // end adjustment

 //Doc the mesh boundaries
 ofstream some_file;
 some_file.open("boundaries.dat");
 mesh_pt()->output_boundaries(some_file);
 some_file.close();

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

 // Setup equation numbering 
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void EighthSpherePoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 


 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 std::cout << "Output file: " << filename << std::endl;


 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
 some_file.close();


 // Doc error
 //----------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->compute_error(some_file,TanhSolnForPoisson::get_exact_u,
                          error,norm); 
 some_file.close();
 cout << "error: " << sqrt(error) << std::endl; 
 cout << "norm : " << sqrt(norm) << std::endl << std::endl;

} // end of doc


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=========start_of_main=============================================
/// Driver for 3D Poisson problem in eighth of a sphere. Solution 
/// has a sharp step. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only a single adaptation.
//===================================================================
int main(int argc, char *argv[]) 
{ 

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Set up the problem with 27-node brick elements, pass pointer to 
 // source function
 EighthSpherePoissonProblem<PRefineableQPoissonElement<3> >
  problem(&TanhSolnForPoisson::get_source);

 // Setup labels for output
  DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;

 // Check if we're ready to go
 std::cout << "Self test: " << problem.self_test() << std::endl;
 
 // Solve the problem
 problem.newton_solve();
 
 //Output solution
 problem.doc_solution(doc_info);

 //Increment counter for solutions 
 doc_info.number()++;

 // Perform initial mesh refinement and solve
 problem.refine_uniformly();
 problem.refine_uniformly();
 problem.p_refine_uniformly();
 
 //std::cout << "Self test: " << problem.self_test() << std::endl;
 problem.newton_solve();
 problem.doc_solution(doc_info);
 doc_info.number()++;
 
 // Adapt and solve
 problem.adapt();
 //std::cout << "Self test: " << problem.self_test() << std::endl;
 problem.newton_solve();
 problem.doc_solution(doc_info);
 doc_info.number()++;
 
 problem.p_adapt();
 //std::cout << "Self test: " << problem.self_test() << std::endl;
 problem.newton_solve();
 problem.doc_solution(doc_info);
 doc_info.number()++;
 
} // end of main









