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
// Driver for an adaptive 2D advection diffusion problem with flux boundary 
// conditions using two separate meshes for the bulk and surface meshes

//Generic routines
#include "generic.h"

// The Advection Diffusion equations
#include "advection_diffusion.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//======start_of_namespace============================================
/// Namespace for exact solution for AdvectionDiffusion equation 
/// with "sharp" step 
//====================================================================
namespace TanhSolnForAdvectionDiffusion
{

 /// Peclet number
 double Peclet=200.0;

 /// Parameter for steepness of step
 double Alpha=10.0;

 /// Parameter for angle of step
 double TanPhi=1.0;

 /// Exact solution as a Vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }

 /// Exact solution as a scalar
 void get_exact_u(const Vector<double>& x, double& u)
 {
  u=tanh(1.0-Alpha*(TanPhi*x[0]-x[1]));
 }


 /// Source function required to make the solution above an exact solution 
 void source_function(const Vector<double>& x_vect, double& source)
 {
  double x=x_vect[0];
  double y=x_vect[1];
  source =
2.0*tanh(-0.1E1+Alpha*(TanPhi*x-y))*(1.0-pow(tanh(-0.1E1+Alpha*(
TanPhi*x-y)),2.0))*Alpha*Alpha*TanPhi*TanPhi+2.0*tanh(-0.1E1+Alpha*(TanPhi*x-y)
)*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha*Alpha-Peclet*(-sin(6.0*y
)*(1.0-pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha*TanPhi+cos(6.0*x)*(1.0-
pow(tanh(-0.1E1+Alpha*(TanPhi*x-y)),2.0))*Alpha);
//  Use this one for the funky unforced solution
//  source=0.0;
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
//  Use this one for the funky unforced solution
//  flux=-1.0; 
 }
 
 /// Wind
 void wind_function(const Vector<double>& x, Vector<double>& wind)
 {
  wind[0]=sin(6.0*x[1]);
  wind[1]=cos(6.0*x[0]);
 }

} //end of namespace



//========= start_of_problem_class=====================================
/// 2D AdvectionDiffusion problem on rectangular domain, discretised with
/// 2D QAdvectionDiffusion elements. Flux boundary conditions are applied
/// along boundary 1 (the boundary where x=L). The specific type of 
/// element is specified via the template parameter.
//====================================================================
template<class ELEMENT> 
class TwoMeshFluxAdvectionDiffusionProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source and wind functions
 TwoMeshFluxAdvectionDiffusionProblem(
  AdvectionDiffusionEquations<2>::AdvectionDiffusionSourceFctPt source_fct_pt,
  AdvectionDiffusionEquations<2>::AdvectionDiffusionWindFctPt wind_fct_pt);

 /// Destructor (empty)
 ~TwoMeshFluxAdvectionDiffusionProblem(){}

 ///  Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);


private:

 ///  Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();

 ///  Create Advection Diffusion flux elements on boundary b of 
 /// the Mesh pointed to by bulk_mesh_pt and add them to the Mesh 
 /// object pointed to by surface_mesh_pt
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const &surface_mesh_pt);

 ///  Delete Advection Diffusion flux elements and wipe the surface mesh
 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 /// Pointer to the "bulk" mesh
 RefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// Pointer to source function
 AdvectionDiffusionEquations<2>::AdvectionDiffusionSourceFctPt Source_fct_pt;

 /// Pointer to wind function
 AdvectionDiffusionEquations<2>::AdvectionDiffusionWindFctPt Wind_fct_pt;

}; // end of problem class




//=======start_of_constructor=============================================
/// Constructor for AdvectionDiffusion problem: Pass pointer to source and
/// wind functions.
//========================================================================
template<class ELEMENT>
TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::
TwoMeshFluxAdvectionDiffusionProblem(
 AdvectionDiffusionEquations<2>::AdvectionDiffusionSourceFctPt source_fct_pt,
 AdvectionDiffusionEquations<2>::AdvectionDiffusionWindFctPt wind_fct_pt)
 :  Source_fct_pt(source_fct_pt), Wind_fct_pt(wind_fct_pt)
{ 
 // Setup "bulk" mesh

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;
 
 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=2.0;

 // Build "bulk" mesh
 Bulk_mesh_pt=new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Create "surface mesh" that will contain only the prescribed-flux 
 // elements. The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_pt = new Mesh;

 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1, but add them to a separate mesh.
 // Note that this is exactly the same function as used in the 
 // single mesh version of the problem, we merely pass different Mesh pointers.
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   //Leave nodes on boundary 1 free
   if (b!=1)
    {
     unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
     for (unsigned n=0;n<n_node;n++)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 
      }
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the AdvectionDiffusion bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // source function
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to AdvectionDiffusion bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = Source_fct_pt;

   //Set the wind function pointer
   el_pt->wind_fct_pt() = Wind_fct_pt;

   // Set the Peclet number
   el_pt->pe_pt() = &TanhSolnForAdvectionDiffusion::Peclet;
  }

 // Loop over the flux elements to pass pointer to prescribed flux function
 n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to AdvectionDiffusion flux element
   AdvectionDiffusionFluxElement<ELEMENT> *el_pt = 
    dynamic_cast< AdvectionDiffusionFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));

   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForAdvectionDiffusion::prescribed_flux_on_fixed_x_boundary;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//====================start_of_actions_before_newton_solve=======================
/// Update the problem specs before solve: Reset boundary conditions
/// to the values from the exact solution.
//========================================================================
template<class ELEMENT>
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::actions_before_newton_solve()
{
 // How many boundaries are in the bulk mesh?
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 
 //Loop over the boundaries in the bulk mesh
 for(unsigned i=0;i<n_bound;i++)
  {
   // Only update the Dirichlet conditions!
   if (i!=1)
    {
     // How many nodes are there on this boundary?
     unsigned n_node = Bulk_mesh_pt->nboundary_node(i);
     
     // Loop over the nodes on boundary
     for (unsigned n=0;n<n_node;n++)
      {
       // Get pointer to node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(i,n);
       
       // Extract nodal coordinates from node:
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       
       // Compute the value of the exact solution at the nodal point
       Vector<double> u(1);
       TanhSolnForAdvectionDiffusion::get_exact_u(x,u);
       
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
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::
doc_solution(DocInfo& doc_info)
{ 

 // Doc refinement levels in bulk mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 Bulk_mesh_pt->get_refinement_levels(min_refinement_level,
                                     max_refinement_level); 
 cout << "Ultimate min/max. refinement levels in bulk mesh : " 
      << min_refinement_level << " " 
      << max_refinement_level << std::endl;

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
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,TanhSolnForAdvectionDiffusion::get_exact_u); 
 some_file.close();


 // Doc error and return of the square of the L2 error
 //---------------------------------------------------
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,TanhSolnForAdvectionDiffusion::get_exact_u,
                               error,norm); 
 some_file.close();

 // Doc L2 error and norm of solution
 cout << "\nNorm of error   : " << sqrt(error) << std::endl; 
 cout << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;


} // end of doc
 
//============start_of_create_flux_elements==============================
/// Create AdvectionDiffusion Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::
create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed-flux element
   AdvectionDiffusionFluxElement<ELEMENT>* flux_element_pt = new 
   AdvectionDiffusionFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements

//============start_of_delete_flux_elements==============================
/// Delete Advection Diffusion Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();

 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();

} // end of delete_flux_elements

//============start_of_actions_before_adapt==============================
// Remove the flux elements from the mesh.
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the flux elements and wipe surface mesh
 delete_flux_elements(Surface_mesh_pt);
 
 // Rebuild the global mesh
 rebuild_global_mesh();
} // end of actions_before_adapt

//============start_of_actions_after_adapt==============================
// Attach flux elements to the mesh.
//=======================================================================
template<class ELEMENT>
void TwoMeshFluxAdvectionDiffusionProblem<ELEMENT>::actions_after_adapt()
{
 // Create prescribed-flux elements from all elements that are 
 // adjacent to boundary 1 and add them to surfac mesh
 create_flux_elements(1,Bulk_mesh_pt,Surface_mesh_pt);
 
 // Rebuild the global mesh
 rebuild_global_mesh();
 
 // Loop over the flux elements to pass pointer to prescribed flux function
 unsigned n_element=Surface_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Advection Diffusion flux element
   AdvectionDiffusionFluxElement<ELEMENT> *el_pt = 
    dynamic_cast< AdvectionDiffusionFluxElement<ELEMENT>*>(
     Surface_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() = 
    &TanhSolnForAdvectionDiffusion::prescribed_flux_on_fixed_x_boundary;
  }
 
 // Doc refinement levels in bulk mesh
 unsigned min_refinement_level;
 unsigned max_refinement_level;
 Bulk_mesh_pt->get_refinement_levels(min_refinement_level,
                                     max_refinement_level); 
 cout << "Min/max. refinement levels in bulk mesh: " 
      << min_refinement_level << " " 
      << max_refinement_level << std::endl;

} // end of actions_after_adapt



//==========start_of_run_it===============================================
/// Run 2D AdvectionDiffusion problem with flux boundary 
/// conditions, using the iterative linear solver and
/// preconditioner specified via the pointers.
//========================================================================
void run_it(LinearSolver* linear_solver_pt)
{

 //Set up the problem
 //------------------

 //Set up the problem with 2D nine-node elements from the
 //RefineableQuadAdvectionDiffusionElement family. 
 //Pass pointer to source function. 
 TwoMeshFluxAdvectionDiffusionProblem<RefineableQAdvectionDiffusionElement<2,
  3> > problem(&TanhSolnForAdvectionDiffusion::source_function,
               &TanhSolnForAdvectionDiffusion::wind_function);
 
 // Set linear solver
 problem.linear_solver_pt()=linear_solver_pt;

 // Change max. number of iterations (some of the solvers used here
 // are rubbish!
 IterativeLinearSolver* iter_solver_pt=
  dynamic_cast<IterativeLinearSolver*>(linear_solver_pt);
 if (iter_solver_pt!=0)
  {
   iter_solver_pt->max_iter()=1000;
  }

 // Get Jacobian matrix and create dummy rhs for re-solve
 // (solution from re-solve should be a vector of ones)
 unsigned n_dof=problem.ndof();
 CRDoubleMatrix matrix;
 DoubleVector residual;
 DoubleVector other_rhs;
 DoubleVector one;
 DoubleVector solution;
 problem.get_jacobian(residual,matrix);
 one.build(residual.distribution_pt(),0.0);
 other_rhs.build(residual.distribution_pt(),0.0);
 one.initialise(1.0);
 matrix.multiply(one,other_rhs);
 linear_solver_pt->enable_resolve();


 // Normal solve
 //-------------

 // Create label for output
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");

 // Solve the problem
 problem.newton_solve();
 
 //Output solution
 problem.doc_solution(doc_info);


 // Re-solve
 //---------

 // Open convergence history for resolve
 IterativeLinearSolver* it_solver_pt=
  dynamic_cast<IterativeLinearSolver*>(linear_solver_pt);
 if (it_solver_pt!=0)
  {
   it_solver_pt->
    open_convergence_history_file_stream("RESLT/resolve_convergence.dat");
  }
 
 // Resolve
 linear_solver_pt->resolve(other_rhs,solution);

 // Get error
 double error=0.0;
 for (unsigned i=0;i<n_dof;i++)
  {
   error+=pow((solution[i]-1.0),2);
  }
 
 // Doc
 ofstream error_file("RESLT/resolve_error.dat");
 error_file << "Error for resolve: " 
            << sqrt(error/double(n_dof)) 
            << std::endl;
 std::cout  << "Error for resolve: " 
            << sqrt(error/double(n_dof)) 
            << std::endl;
 error_file.close();

 // Switch of ability to resolve
 linear_solver_pt->disable_resolve();



 // Linear-algebra solve
 //---------------------

 if (it_solver_pt!=0)
  {
   it_solver_pt->
    open_convergence_history_file_stream("RESLT/la_solve_convergence.dat");
  }
 
 DoubleVector solution2;
 linear_solver_pt->solve(&matrix,other_rhs,solution2);

 // Get error
 error=0.0;
 for (unsigned i=0;i<n_dof;i++)
  {
   error+=pow((solution2[i]-1.0),2);
  }
 
 // Doc
 ofstream error_file2("RESLT/la_solve_error.dat");
 error_file2 << "Error for LA solve: " 
             << sqrt(error/double(n_dof)) 
             << std::endl;
 std::cout << "Error for LA solve: " 
           << sqrt(error/double(n_dof)) 
           << std::endl;
 error_file2.close();

} //end of run_it




//==========start_of_main=================================================
/// Check out various iterative linear solvers on an advection
/// diffusion problem. Command line arg specifies the solver.
//========================================================================
int main(int argc, char *argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 //Check
 if (CommandLineArgs::Argc!=2)
  {
   throw OomphLibError(
    "Code needs one argument to specify which solver to use",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 // Default value for Peclet number -- can be reset to zero for
 // solvers that require symmetric matrices
 TanhSolnForAdvectionDiffusion::Peclet=200.0;

 // Tolerance for iterative linear solvers
 double tol=1.0e-12;

 // Pointer to linear solver
 LinearSolver* linear_solver_pt=new SuperLUSolver;

 // Pointer to iterative linear solver
 IterativeLinearSolver* it_linear_solver_pt=0;

 // Preconditioner
 Preconditioner* prec_pt=new ILUZeroPreconditioner<CRDoubleMatrix>;


 /// Possible endless loop for memory leak testing
 bool doit=true;
 while (doit)
  {
   
   // Switch based on one and only command line argument
   int solver_flag=atoi(CommandLineArgs::Argv[1]);
   switch (solver_flag)
    {
     
     
     // SuperLU
     //--------
    case 0:
     
     // This solver can handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=200.0;
     
     std::cout << "Solving with SuperLU direct solver" 
               << std::endl;  
     break;

     
     // GMRES
     //------
    case 1:
     
     // This solver can handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=200.0;
     
     // Build solver
     it_linear_solver_pt = new GMRES<CRDoubleMatrix>;
     
     // Build/pass preconditioner
     it_linear_solver_pt->preconditioner_pt() = prec_pt;
     
     // Switch on doc of convergence history and adjust tolerance
     it_linear_solver_pt->tolerance()=tol;
     it_linear_solver_pt->
      open_convergence_history_file_stream("RESLT/convergence.dat");
     
     // Set solver flag
     linear_solver_pt=it_linear_solver_pt;
     
     // doc
     std::cout << "Solving with GMRES" << std::endl;  
     
     break;
     
     
     
     
     // BiCGStab
     //---------
    case 2:
   
     // This solver can handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=200.0;

     // Build solver
     it_linear_solver_pt = new BiCGStab<CRDoubleMatrix>;
   
     // Build/pass preconditioner
     it_linear_solver_pt->preconditioner_pt() = prec_pt;
   
     // Switch on doc of convergence history and adjust tolerance
     it_linear_solver_pt->tolerance()=tol;
     it_linear_solver_pt->
      open_convergence_history_file_stream("RESLT/convergence.dat");
   
     // Set solver flag
     linear_solver_pt=it_linear_solver_pt;
   
     // doc
     std::cout << "Solving with BiCGStab" << std::endl;  
  
     break;



     // CG
     //----
    case 3:
   
     // This solver can't handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=0.0;

     // Build solver
     it_linear_solver_pt = new CG<CRDoubleMatrix>;
   
     // Build/pass preconditioner
     it_linear_solver_pt->preconditioner_pt() = prec_pt;
   
     // Switch on doc of convergence history and adjust tolerance
     it_linear_solver_pt->tolerance()=tol;
     it_linear_solver_pt->
      open_convergence_history_file_stream("RESLT/convergence.dat");
   
     // Set solver flag
     linear_solver_pt=it_linear_solver_pt;
   
     // doc
     std::cout << "Solving with CG" << std::endl;  
  
     break;


     // GS
     //----
    case 4:
      
     // This solver can't handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=0.0;

     // Build solver 
     it_linear_solver_pt = new GS<CRDoubleMatrix>;
   
     // Note: Doesn't have a preconditioner
   
     // Switch on doc of convergence history and adjust tolerance
     it_linear_solver_pt->tolerance()=tol;
     it_linear_solver_pt->
      open_convergence_history_file_stream("RESLT/convergence.dat");
   
     // Set solver flag
     linear_solver_pt=it_linear_solver_pt;
   
     // doc
     std::cout << "Solving with GS" << std::endl;  
  
     break;
     
     // Damped Jacobi
     //--------------
    case 5:
      
     // This solver can't handle non-symmetric matrices
     TanhSolnForAdvectionDiffusion::Peclet=0.0;

     // Build solver 
     it_linear_solver_pt = new DampedJacobi<CRDoubleMatrix>;
   
     // Note: Doesn't have a preconditioner
   
     // Switch on doc of convergence history and adjust tolerance
     it_linear_solver_pt->tolerance()=tol;
     it_linear_solver_pt->
      open_convergence_history_file_stream("RESLT/convergence.dat");
   
     // Set solver flag
     linear_solver_pt=it_linear_solver_pt;
   
     // doc
     std::cout << "Solving with damped Jacobi" << std::endl;  
  
     break;


     // Wrong
     //------
    default:

     throw OomphLibError("Invalid command line argument: solver flag: ",
			 OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
     break;
   
    }

   // Build problem and solve it with specified linear solver
   run_it(linear_solver_pt);
 
   delete it_linear_solver_pt;

   // Comment this out if you want an endless loop to test
   // for memory leaks.
   doit=false;

  }

} //end of main
