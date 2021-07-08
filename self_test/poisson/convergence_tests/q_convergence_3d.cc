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

 ///\short Orientation (non-normalised x-component of unit vector in direction
 /// of step plane)
 double N_x=-1.0;

 ///\short Orientation (non-normalised y-component of unit vector in direction
 /// of step plane)
 double N_y=-1.0;

 ///\short Orientation (non-normalised z-component of unit vector in direction
 /// of step plane)
 double N_z=1.0;


 ///\short Orientation (x-coordinate of step plane) 
 double X_0=0.0;

 ///\short Orientation (y-coordinate of step plane) 
 double Y_0=0.0;

 ///\short Orientation (z-coordinate of step plane) 
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
  PoissonEquations<3>::PoissonSourceFctPt source_fct_pt,
  const unsigned& h_power);

 /// Destructor: Empty
 ~EighthSpherePoissonProblem(){}

 /// \short Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 RefineableEighthSphereMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableEighthSphereMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()  {}

 /// \short Update the problem specs before solve: 
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
 void doc_solution(DocInfo& doc_info,ostream& convergence_file);

private:

 /// Pointer to source function
 PoissonEquations<3>::PoissonSourceFctPt Source_fct_pt;

 /// Exponent for h scaling
 unsigned H_power;


}; // end of class definition





//====================start_of_constructor================================
/// Constructor for Poisson problem on eighth of a sphere mesh
//========================================================================
template<class ELEMENT>
EighthSpherePoissonProblem<ELEMENT>::EighthSpherePoissonProblem(
   PoissonEquations<2>::PoissonSourceFctPt source_fct_pt,
   const unsigned& h_power) : 
         Source_fct_pt(source_fct_pt), H_power(h_power)

{ 

 Problem::Problem_is_nonlinear=false;

 // Setup parameters for exact tanh solution
 // Steepness of step
 TanhSolnForPoisson::Alpha=5.0;

 /// Create mesh for sphere of radius 5
 double radius=5.0; 
 Problem::mesh_pt() = new RefineableEighthSphereMesh<ELEMENT>(radius);

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

 // Do refinement
 for (unsigned i=1;i<h_power;i++) refine_uniformly();

} // end of constructor



//========================start_of_doc====================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void EighthSpherePoissonProblem<ELEMENT>::doc_solution(
DocInfo& doc_info,
ostream& convergence_file)
{
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 



 // Switch on full validation during self test
 if (CommandLineArgs::Argc>1)
  {

   // Output boundaries
   //------------------
   sprintf(filename,"%s/boundaries.dat",doc_info.directory().c_str());
   some_file.open(filename);
   mesh_pt()->output_boundaries(some_file);
   some_file.close();

   // Output solution 
   //-----------------
   sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   mesh_pt()->output(some_file,npts);
   some_file.close();
   
   
   // Output exact solution 
   //----------------------
   sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   mesh_pt()->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
   some_file.close();
  }
 

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

 double log_h=log(1.0/double(pow(2.0,int(H_power))));
 convergence_file << log_h << " "  
                  << log(sqrt(error)) << " ";

 unsigned nnode_1d=mesh_pt()->finite_element_pt(0)->nnode_1d();
 switch (nnode_1d)
  {
  case 2:
   convergence_file << 2.0*log_h+4.0 << " ";
   break;

  case 3:
   convergence_file << 3.0*log_h+4.0 << " ";
   break;

  case 4:
   convergence_file << 4.0*log_h+4.0 << " ";
   break;

  default:
   throw OomphLibError("Never get here",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

  convergence_file << sqrt(norm) << std::endl;


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

 /// Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");

 

 // Output for convergence data
 ofstream convergence_file("RESLT/convergence.dat");

 convergence_file << "VARIABLES=\"log h\",\"log e\",\"theory\",\"log u\" " 
                  << std::endl;

 unsigned nstep;



 cout << "Linear element" << std::endl;
 cout << "==============" << std::endl << std::endl;



 // Number of steps (reduced for validation)
 nstep = 5;  
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }

 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Linear element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   EighthSpherePoissonProblem<RefineableQPoissonElement<3,2> >
    problem(&TanhSolnForPoisson::get_source,h_power);

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }



 cout << "Quadratic element" << std::endl;
 cout << "=================" << std::endl << std::endl;

 // Number of steps (reduced for validation)
 nstep = 4;  
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }
 
 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Quadratic element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   EighthSpherePoissonProblem<RefineableQPoissonElement<3,3> >
    problem(&TanhSolnForPoisson::get_source,h_power);

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }




 cout << "Cubic element" << std::endl;
 cout << "=============" << std::endl << std::endl;

 // Number of steps (reduced for validation)
 nstep = 4;  
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }

 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Cubic element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   EighthSpherePoissonProblem<RefineableQPoissonElement<3,4> >
    problem(&TanhSolnForPoisson::get_source,h_power);

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }









} // end of main









