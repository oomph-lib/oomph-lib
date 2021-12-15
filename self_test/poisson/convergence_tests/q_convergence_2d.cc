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
//Driver function for a simple test poisson problem using a triangle grid

//Generic routines
#include "generic.h"

// Poisson equations
#include "poisson.h"

// Triangular mesh
#include "meshes/simple_rectangular_quadmesh.h"


using namespace std;

using namespace oomph;

using namespace MathematicalConstants;

//====================================================================
/// Namespace for exact solution for Poisson equation with sharp step 
//====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of step
 double Alpha;

 /// Parameter for angle of step
 double Beta;
 
 // Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  u=tanh(1.0-Alpha*(Beta*x[0]-x[1]));
 }

 /// Source function to make it an exact solution 
 void get_source(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
             (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*
             Alpha*Alpha*Beta*Beta+2.0*tanh(-1.0+Alpha*(Beta*x[0]-x[1]))*
            (1.0-pow(tanh(-1.0+Alpha*(Beta*x[0]-x[1])),2.0))*Alpha*Alpha;
 }

}


//====================================================================
/// Poisson problem.
//====================================================================
template<class ELEMENT> 
class PoissonProblem : public Problem
{

public:

 /// Constructor
  PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt, 
                 const unsigned& h_power);

 /// Destructor (empty)
 ~PoissonProblem(){};

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   //Loop over the boundaries
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned long ibound=0;ibound<num_bound;ibound++)
    {
     // Loop over the nodes on boundary
     unsigned num_nod=mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
       double u;
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       TanhSolnForPoisson::get_exact_u(x,u);
       nod_pt->set_value(0,u);
      }
    }
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve()
  {}

  /// Doc the solution
 void doc_solution(DocInfo& doc_info, ostream& convergence_file);

private:

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

 /// Exponent for h scaling
 unsigned H_power;

};


//========================================================================
/// Constructor for Poisson problem
///
//========================================================================
template<class ELEMENT>
PoissonProblem<ELEMENT>::
PoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt, 
               const unsigned& h_power) : Source_fct_pt(source_fct_pt),
                                        H_power(h_power)
                                        
{ 
 
 Problem::Problem_is_nonlinear=false;

 // Setup parameters for exact tanh solution

 // Steepness of step
 TanhSolnForPoisson::Alpha=1.0; 

 // Orientation of step
 TanhSolnForPoisson::Beta=1.4;

 //Create mesh and assign element lengthscale h
 unsigned nx=unsigned(pow(2.0,int(h_power)));
 unsigned ny=unsigned(pow(2.0,int(h_power)));
 double lx=8.0;
 double ly=8.0;
// Problem::mesh_pt() = new SimpleRectangularTriMesh<ELEMENT>(nx,ny,lx,ly);
 
 // Build and assign mesh
 Problem::mesh_pt() = new SimpleRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly);

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here. 
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned long ibound=0;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
     }
   }
  
  // Complete the build of all elements so they are fully functional
  
  
  //Find number of elements in mesh
  unsigned n_element = mesh_pt()->nelement();
  
  // Loop over the elements to set up element-specific 
  // things that cannot be handled by constructor
  for(unsigned i=0;i<n_element;i++)
   {
    // Upcast from GeneralElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
    
    //Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;

   }
  
  // Setup equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
}


//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info, 
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
   //----------------
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
}

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 /// Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;

 //Set step orientation
 TanhSolnForPoisson::Beta=0.4;

 
 // Number of steps (reduced for validation)
 unsigned nstep = 7;  
 if (CommandLineArgs::Argc>1)
  {
   nstep=3;
  }

 // Output for convergence data
 ofstream convergence_file("RESLT/convergence.dat");

 convergence_file << "VARIABLES=\"log h\",\"log e\",\"theory\",\"log u\" " 
                  << std::endl;


 cout << "Linear elements " << std::endl;
 cout << "================" << std::endl << std::endl;


 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Linear element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   PoissonProblem<QPoissonElement<2,2> > 
    problem(&TanhSolnForPoisson::get_source,h_power);
   

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }



 cout << "Quadratic elements " << std::endl;
 cout << "==================" << std::endl << std::endl;


 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Quadratic element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   PoissonProblem<QPoissonElement<2,3> > 
    problem(&TanhSolnForPoisson::get_source,h_power);


//    // Create an instance of a lower-order integration scheme
//    Gauss<2,2>* int_pt=new Gauss<2,2>;

//   //Find number of elements in mesh
//   unsigned n_element = problem.mesh_pt()->nelement();
  
//   // Loop over the elements to set up element-specific 
//   // things that cannot be handled by constructor
//   for(unsigned i=0;i<n_element;i++)
//    {
//     //Set the a different integration scheme
//     problem.mesh_pt()->finite_element_pt(i)->set_integration_scheme(int_pt); 
//    }

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }



 cout << "Cubicc elements " << std::endl;
 cout << "================" << std::endl << std::endl;


 //Tecplot output separator for convergence results
 convergence_file << "ZONE T=\"Cubic element\"" << std::endl;
 for (unsigned h_power=1;h_power<nstep+1;h_power++)
  {
   //Set up the problem
   PoissonProblem<QPoissonElement<2,4> > 
    problem(&TanhSolnForPoisson::get_source,h_power);
   

//    // Create an instance of a lower-order integration scheme
//    Gauss<2,3>* int_pt=new Gauss<2,3>;

//   //Find number of elements in mesh
//   unsigned n_element = problem.mesh_pt()->nelement();
  
//   // Loop over the elements to set up element-specific 
//   // things that cannot be handled by constructor
//   for(unsigned i=0;i<n_element;i++)
//    {
//     //Set the a different integration scheme
//     problem.mesh_pt()->finite_element_pt(i)->set_integration_scheme(int_pt); 
//    }

   // Solve the problem
   problem.newton_solve();
   
   //Output solution & convergence data
   problem.doc_solution(doc_info,convergence_file);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }

 convergence_file.close();

}


