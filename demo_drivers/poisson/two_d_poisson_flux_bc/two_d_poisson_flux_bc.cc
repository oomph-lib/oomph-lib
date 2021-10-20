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
//Driver for a simple 2D Poisson problem with flux boundary conditions

//Generic routines
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The mesh
#include "meshes/simple_rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//===== start_of_namespace=============================================
/// Namespace for exact solution for Poisson equation with "sharp step" 
//=====================================================================
namespace TanhSolnForPoisson
{

 /// Parameter for steepness of "step"
 double Alpha=1.0;

 /// Parameter for angle Phi of "step"
 double TanPhi=0.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = 2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*
   Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-1.0+Alpha*(TanPhi*x[0]-x[1]))*
   (1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*Alpha;
 }

 /// Flux required by the exact solution on a boundary on which x is fixed
 void prescribed_flux_on_fixed_x_boundary(const Vector<double>& x, 
                                          double& flux)
 {
  //The outer unit normal to the boundary is (1,0)
  double N[2] = {1.0, 0.0};
  //The flux in terms of the normal is
  flux = 
   -(1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*TanPhi*N[0]+(
   1.0-pow(tanh(-1.0+Alpha*(TanPhi*x[0]-x[1])),2.0))*Alpha*N[1];
 }
 
} // end of namespace





//========= start_of_problem_class=====================================
/// 2D Poisson problem on rectangular domain, discretised with
/// 2D QPoisson elements. Flux boundary conditions are applied
/// along boundary 1 (the boundary where x=L). The specific type of 
/// element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class FluxPoissonProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source function
 FluxPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~FluxPoissonProblem(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

private:

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 ///  Create Poisson flux elements on the b-th boundary of the 
 /// problem's mesh
 void create_flux_elements(const unsigned &b);

 ///  Number of Poisson "bulk" elements (We're attaching the flux 
 /// elements to the bulk mesh --> only the first Npoisson_elements elements
 /// in the mesh are bulk elements!)
 unsigned Npoisson_elements;

 /// Pointer to source function
 PoissonEquations<2>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Poisson problem: Pass pointer to source function.
//========================================================================
template<class ELEMENT>
FluxPoissonProblem<ELEMENT>::
      FluxPoissonProblem(PoissonEquations<2>::PoissonSourceFctPt source_fct_pt)
       :  Source_fct_pt(source_fct_pt)
{ 

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build and assign mesh
 Problem::mesh_pt()=new SimpleRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Store number of Poisson bulk elements (= number of elements so far).
 Npoisson_elements=mesh_pt()->nelement();
 
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1 and add them to the (single) global mesh
 create_flux_elements(1);

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   //Leave nodes on boundary 1 free
   if(b!=1)
    {
     unsigned n_node= mesh_pt()->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       mesh_pt()->boundary_node_pt(b,n)->pin(0); 
      }
    }
  }

 // Loop over the Poisson bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to source
 // function
 for(unsigned e=0;e<Npoisson_elements;e++)
  {
   // Upcast from GeneralisedElement to Poisson bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;
  }

 // Total number of elements:
 unsigned n_element=mesh_pt()->nelement();
 
 // Loop over the flux elements (located at the "end" of the
 // mesh) to pass function pointer to prescribed-flux function.
 for (unsigned e=Npoisson_elements;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Poisson flux element
   PoissonFluxElement<ELEMENT> *el_pt = 
    dynamic_cast<PoissonFluxElement<ELEMENT>*>(mesh_pt()->element_pt(e));

   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForPoisson::prescribed_flux_on_fixed_x_boundary;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor


//============start_of_create_flux_elements==============================
/// Create Poisson Flux Elements on the b-th boundary of the Mesh
//=======================================================================
template<class ELEMENT>
void FluxPoissonProblem<ELEMENT>::create_flux_elements(const unsigned &b)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = mesh_pt()->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    mesh_pt()->boundary_element_pt(b,e));
   
   // What is the index of the face of the bulk element at the boundary
   int face_index = mesh_pt()->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   PoissonFluxElement<ELEMENT>* flux_element_pt = new 
   PoissonFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the mesh
   mesh_pt()->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements



//====================start_of_actions_before_newton_solve================
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void FluxPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the boundaries
 for(unsigned i=0;i<n_bound;i++)
  {
   // Only update Dirichlet nodes
   if (i!=1)
    {
     // How many nodes are there on this boundary?
     unsigned n_node = mesh_pt()->nboundary_node(i);
     
     // Loop over the nodes on boundary
     for (unsigned n=0;n<n_node;n++)
      {
       // Get pointer to node
       Node* nod_pt = mesh_pt()->boundary_node_pt(i,n);
       
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
  }
} // end of actions before solve



//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void FluxPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for(unsigned e=0;e<Npoisson_elements;e++)
  {
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   el_pt->output_fct(some_file,npts,TanhSolnForPoisson::get_exact_u); 
  }
 some_file.close();

 // Can't use black-box error computatation routines because
 // the mesh contains two different types of elements. 
 // error function hasn't been implemented for the prescribed
 // flux elements...

} // end of doc

 

//==========start_of_main=================================================
/// Demonstrate how to solve 2D Poisson problem with flux boundary 
/// conditions
//========================================================================
int main()
{

 //Set up the problem
 //------------------

 //Set up the problem with 2D nine-node elements from the
 //QPoissonElement family. Pass pointer to source function. 
 FluxPoissonProblem<QPoissonElement<2,3> > 
  problem(&TanhSolnForPoisson::source_function);
 

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Step number
 doc_info.number()=0;



 // Check that we're ready to go:
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

 
 // Set the orientation of the "step" to 45 degrees
 TanhSolnForPoisson::TanPhi=1.0;
 
 // Initial value for the steepness of the "step"
 TanhSolnForPoisson::Alpha=1.0; 

 // Do a couple of solutions for different forcing functions
 //---------------------------------------------------------
 unsigned nstep=4;
 for (unsigned istep=0;istep<nstep;istep++)
  {
   // Increase the steepness of the step:
   TanhSolnForPoisson::Alpha+=2.0;

   cout << "\n\nSolving for TanhSolnForPoisson::Alpha="
        << TanhSolnForPoisson::Alpha << std::endl << std::endl;

   // Solve the problem
   problem.newton_solve();

   //Output solution
   problem.doc_solution(doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++; 
  }

} //end of main









