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
// Helmoholtz equation with one PML and forcing on the opposite side and
// periodic BC on all other sides
#include <fenv.h>
#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The Helmholtz equations and complex-valued multigrid machinery
#include "pml_helmholtz.h"

// The mesh
#include "meshes/simple_cubic_mesh.h"

using namespace std;
using namespace oomph;

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//================================================start_of_namespace======
/// Namespace for the Helmholtz problem parameters
//========================================================================
namespace GlobalParameters
{
 /// Solver specific parameters:
 /// ----------------------------------- 
 /// The number of nodes in one direction (default=2)
 unsigned Nnode_1d=2;
 
 /// The minimum level of uniform refinement
 unsigned Min_refinement_level=1;
 
 /// The additional levels of uniform refinement 
 unsigned Add_refinement_level=0;
 
 /// The number of adaptations allowed by the Newton solver
 unsigned N_adaptations=1;

 /// The choice of whether or not to use adaptation
 ///    0 = Uniform refinement
 ///    1 = Adaptive refinement
 unsigned Use_adaptation_flag=0;

 /// The choice of pre-smoother:
 ///    0 = Automatic (GMRES as a smoother on levels where kh>0.5)
 ///    1 = Damped Jacobi on all levels with a constant omega value
 unsigned Pre_smoother_flag=0;
 
 /// The choice of post-smoother:
 ///    0 = Automatic (GMRES as a smoother on levels where kh>0.5)
 ///    1 = Damped Jacobi on all levels with a constant omega value
 unsigned Post_smoother_flag=0;

 /// The choice of linear solver
 ///    0 = SuperLU
 ///    1 = Multigrid
 unsigned Linear_solver_flag=1;
 
 /// The MG solver allows for five different levels of output:
 ///    0 = Outputs everything
 ///    1 = Outputs everything except the smoother timings 
 ///    2 = Outputs setup information but no V-cycle timings
 ///    3 = Suppresses all output
 unsigned Output_management_flag=0;
  
 /// Variable used to decide whether or not convergence information
 /// is displayed:
 ///    0 = Don't display convergence information
 ///    1 = Display convergence information
 unsigned Doc_convergence_flag=0;
 
 /// DocInfo object used for documentation of the solution
 DocInfo Doc_info;
 
 // Pointer to the output stream -- defaults to oomph_info
 std::ostream* Stream_pt;
  
 /// Problem specific parameters:
 /// ------------------------------------
 /// Length of the cube in each direction
 double Lx=1.0;
 double Ly=1.0;
 double Lz=1.0;

 /// Number of elements in each direction (used by SimpleCubicMesh)
 unsigned Nx=7;
 unsigned Ny=7;
 unsigned Nz=7;
 
 /// The element width
 double Element_width=Lx/double(Nx);
 
 /// Length of cube in each direction
 double Pml_thickness=Element_width;
 
 /// Store the value of Pi
 double Pi=MathematicalConstants::Pi;

 /// Choose the value of the shift to create the complex-shifted
 /// Laplacian preconditioner (CSLP)
 double Alpha_shift=0.0;
 
 /// Square of the wavenumber (also known as k^2)
 double K_squared=20.0;

 /// Wavenumber (also known as k),k=omega/c
 double Wavenumber=sqrt(K_squared);

 /// Update the parameters passed in at the command line
 void update_parameters()
 {
  /// Wavenumber (also known as k), k=omega/c
  Wavenumber=sqrt(K_squared);
 }
 
 /// The x and y coordinate of the centre of the cube 
 double Centre=Lx/2.0;

 /// Get the exact solution, u, at the spatial position, x
 void get_simple_exact_u(const Vector<double>& x,Vector<double>& u)
 {
  // Initialise a variable to store the radial distance
  double r=std::sqrt((x[0]-Centre)*(x[0]-Centre)
		     +(x[1]-Centre)*(x[1]-Centre)
		     +(x[2]-Centre)*(x[2]-Centre));

  // Scale the radial distance by the wavenumber
  double kr=Wavenumber*r;

  // The solution is singular at the centre so set it to zero 
  if (r==0.0)
  {
   // Set the real part of the solution value
   u[0]=0.0;
   
   // Set the imaginary part of the solution value
   u[1]=0.0;
  }
  // Otherwise set the correct solution value
  else
  {
   // Set the real part of the solution value
   u[0]=cos(kr)/kr;
   
   // Set the imaginary part of the solution value
   u[1]=sin(kr)/kr;
  }
 } // End of get_simple_exact_u

 // Set the exact solution pointer to the get_simple_exact_u function above
 FiniteElement::SteadyExactSolutionFctPt simple_exact_u_pt=&get_simple_exact_u;

 /// New mapping function that makes the mapping independent of the
 /// PML thickness
 class TestPMLMapping : public virtual PMLMapping
 {

 public:

  /// Default constructor (empty)
  TestPMLMapping(){};

  /// Overwrite the pure PML mapping coefficient function to return the
  /// coeffcients proposed by Bermudez et al
  std::complex<double> gamma(const double& nu_i,
			     const double& pml_width_i,
			     const double& k_squared_local,
			     const double& alpha_shift)
   {
    // The "effective k^2" is shifted, so we shift the k used in the
    // transformation too
    std::complex<double> k_shifted=
     sqrt(k_squared_local*std::complex<double>(1.0,alpha_shift));

    // Return the gamma in J++, with the shifted k
    return (1.0/k_shifted)*std::complex<double>
     (0.0,1.0/(std::fabs(pml_width_i-nu_i)));    
   } // End of gamma
 }; // End of TestPMLMapping

 /// Set the new PML mapping
 TestPMLMapping* Test_pml_mapping_pt=new TestPMLMapping;
 
 /// The choice of whether or not to enable the new test mapping
 ///    1 = Enable test mapping
 ///    0 = Disable test mapping
 unsigned Enable_test_pml_mapping_flag=0;
 
 /// The tolerance for a point relative to the bounding inner square
 double Eps=1.0e-12;
   
 /// Function to determine whether or not a point lies in the centre
 /// of the mesh (in the pinned region)
 bool is_in_pinned_region(const Vector<double>& x)
 {
  // Check if the element lies in the central cube region
  return (abs(x[0]-GlobalParameters::Centre)<
	  (0.5*GlobalParameters::Element_width+Eps)&&
	  abs(x[1]-GlobalParameters::Centre)<
	  (0.5*GlobalParameters::Element_width+Eps)&&
	  abs(x[2]-GlobalParameters::Centre)<
	  (0.5*GlobalParameters::Element_width+Eps));
 } // End of is_in_pinned_region 
} // End of namespace

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Returns a pointer to a smoother of the appropriate type
//========================================================================
namespace Smoother_Factory_Function_Helper
{
 /// The value of the damping factor for the damped Jacobi smoother
 double Omega=0.4;
 
 /// Returns a pointer to a Smoother object which is to be used as
 /// the pre-smoother
 HelmholtzSmoother* set_pre_smoother()
 {
  // Create a new DampedJacobi object
  return new ComplexDampedJacobi<CRDoubleMatrix>(Omega);
 } 
 
 /// Returns a pointer to a Smoother object which is to be used as
 /// the post-smoother
 HelmholtzSmoother* set_post_smoother()
 {
  // Create a new DampedJacobi object
  return new ComplexDampedJacobi<CRDoubleMatrix>(Omega);
 }
} // End of Smoother_Factory_Function_Helper

/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////

//============================================start_of_problem_class======
/// Problem class
//========================================================================
template<class ELEMENT>
class PMLStructuredCubicHelmholtz : public HelmholtzMGProblem
{

public:

 /// Constructor
 PMLStructuredCubicHelmholtz();

 /// Destructor (empty)
 ~PMLStructuredCubicHelmholtz();

 /// Doc the solution
 void doc_solution();

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
 {
  // Document the solution
  doc_solution();
 }

 /// Actions before adapt: (empty)
 void actions_before_adapt(){}

 /// Actions after adapt:(empty)
 void actions_after_adapt();

 /// Set GMRES preconditioner by multigrid as the linear solver
 void set_gmres_multigrid_solver();
 
 /// Enable the PML mapping function for all nodes in the PML region
 void enable_pmls();

 // Apply boundary conditions
 void apply_boundary_conditions();

private:

 /// Pointer to the "bulk" mesh
 RefineableSimpleCubicMesh<ELEMENT>* Bulk_mesh_pt;

 /// Overload the make_new_problem function to return an object of this class
 HelmholtzMGProblem* make_new_problem()
 {
  // Return a new problem pointer
  return new PMLStructuredCubicHelmholtz<ELEMENT>;
 }

 /// Overload the mg_bulk_mesh_pt function to return a pointer to the
 /// "refineable" portion of the mesh
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
 {
  // Return the pointer to the bulk mesh
  return Bulk_mesh_pt;
 }
 
 /// Trace file
 ofstream Trace_file;
}; // End of PMLStructuredCubicHelmholtz class

//==============================================start_of_constructor======
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLStructuredCubicHelmholtz<ELEMENT>::PMLStructuredCubicHelmholtz()
{
 // Indicate that the problem is nonlinear to ensure the residual is
 // calculated at the end of the iteration
 problem_is_nonlinear(true);

 // Set the number of Newton iterations to one
 max_newton_iterations()=1;

 // Set up solver specific information:
 //------------------------------------
 // If we're choosing to use GMRES & MG as our linear solver
 if (GlobalParameters::Linear_solver_flag==1)
 {
  // Set the solver
  set_gmres_multigrid_solver();
 }
 
 // Open trace file
 Trace_file.open("RESLT/trace.dat");

 // Build the mesh using the specified parameters:
 //-----------------------------------------------
 // Build the "bulk" mesh
 Bulk_mesh_pt=new RefineableSimpleCubicMesh<ELEMENT>(
  GlobalParameters::Nx,GlobalParameters::Ny,GlobalParameters::Nz,
  GlobalParameters::Lx,GlobalParameters::Ly,GlobalParameters::Lz);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Choose error tolerances
 Bulk_mesh_pt->max_permitted_error()=1.0e-03;
 Bulk_mesh_pt->min_permitted_error()=1.0e-06;
 
 // Create the main mesh
 add_sub_mesh(Bulk_mesh_pt);

 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // How many elements in the mesh?
 unsigned n_element=mesh_pt()->nelement();

 // Loop over the elements and pass a pointer to the value of k^2
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Set the wavenumber function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;
   
   // If we're using Jonathon's new test mapping
   if (GlobalParameters::Enable_test_pml_mapping_flag)
   {
    // Set the PML mapping function
    el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
   }
  } // if (el_pt!=0)
 } // for (unsigned e=0;e<n_element;e++)

 // Apply the boundary conditions, both in the central region and on the
 // outer boundary (since these nodes are PML nodes)
 apply_boundary_conditions();
 
 // Enable the PML mapping in elements in the PML region
 enable_pmls();
 
 // Setup equation numbering scheme
 assign_eqn_numbers(); 
} // End of constructor

//===============================================start_of_destructor======
/// Destructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLStructuredCubicHelmholtz<ELEMENT>::~PMLStructuredCubicHelmholtz()
{
 // If we're using GMRES & MG as the linear solver
 if (GlobalParameters::Linear_solver_flag==1)
 {
  // Delete the MG solver pointers
  delete dynamic_cast<HelmholtzFGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt();

  // Set the pointer to null
  dynamic_cast<HelmholtzFGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt()=0;
    
  // Delete the MG solver pointers
  delete linear_solver_pt();

  // Set the pointer to null
  linear_solver_pt()=0;    
 }
   
 // Delete the error estimator
 delete Bulk_mesh_pt->spatial_error_estimator_pt();
 
 // Set the pointer to null
 Bulk_mesh_pt->spatial_error_estimator_pt()=0;

 // Delete the "bulk" mesh
 delete Bulk_mesh_pt;

 // Set the pointer to null
 Bulk_mesh_pt=0;
 
} // End of ~PMLStructuredCubicHelmholtz

//=======set_gmres_multigrid_solver=======================================
/// Build and set GMRES preconditioner by multigrid as the linear solver
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::set_gmres_multigrid_solver()
{
 // Create linear solver
 HelmholtzFGMRESMG<CRDoubleMatrix>* solver_pt=
  new HelmholtzFGMRESMG<CRDoubleMatrix>;

 // Set the number of iterations
 solver_pt->max_iter()=200;

 // Set the tolerance (to ensure the Newton solver converges in one step)
 solver_pt->tolerance()=1.0e-10;
   
 // If the user wishes to document the convergence information
 if (GlobalParameters::Doc_convergence_flag)
 {
  // Create a file to record the convergence history
  solver_pt->open_convergence_history_file_stream("RESLT/conv.dat");
 }
 
 // Create linear solver 
 linear_solver_pt()=solver_pt;
 
 // This preconditioner uses multigrid on the block version of the full
 // matrix. 2 V-cycles will be used here per preconditioning step
 HelmholtzMGPreconditioner<3>* prec_pt=new HelmholtzMGPreconditioner<3>(this);

 // Set preconditioner
 solver_pt->preconditioner_pt()=prec_pt;
  
 // Set the shift
 prec_pt->alpha_shift()=GlobalParameters::Alpha_shift;
   
 // If the user wants to use damped Jacobi on every level as a smoother
 if (GlobalParameters::Pre_smoother_flag==1)
 {
  // Set the pre-smoother factory function
  prec_pt->set_pre_smoother_factory_function
   (Smoother_Factory_Function_Helper::set_pre_smoother);
 }

 // If the user wants to use damped Jacobi on every level as a smoother
 if (GlobalParameters::Post_smoother_flag==1)
 {
  // Set the post-smoother factory function
  prec_pt->set_post_smoother_factory_function
   (Smoother_Factory_Function_Helper::set_post_smoother);
 }
 
 // Suppress certain timings
 if (GlobalParameters::Output_management_flag==1)
 {
  prec_pt->disable_doc_time();
 }
 else if (GlobalParameters::Output_management_flag==2)
 {
  prec_pt->disable_v_cycle_output();
 }
 else if (GlobalParameters::Output_management_flag==3)
 {
  prec_pt->disable_output();
 }
} // End of set_gmres_multigrid_solver

//================================start_of_apply_boundary_conditions======
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::apply_boundary_conditions()
{
 // Find the number of elements in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Vector to hold the local coordinates of a point in any given element 
 Vector<double> s(3,0.0);
 
 // Vector to hold the (Eulerian) coordinates of a point
 Vector<double> x(3,0.0);

 // Vector to hold the real and imaginary part of the solution
 Vector<double> u(2,0.0);

 // If the user wishes to silence everything
 if (GlobalParameters::Output_management_flag==3)
 {
  // Store the output stream pointer
  GlobalParameters::Stream_pt=oomph_info.stream_pt();

  // Now set the oomph_info stream pointer to the null stream to
  // disable all possible output
  oomph_info.stream_pt()=&oomph_nullstream;
 }
    
 // Loop over the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);

   // Check if the element lies in the central cube region
   if (GlobalParameters::is_in_pinned_region(x))
   {
    // Calculate the number of nodes in the element
    unsigned nnode=el_pt->nnode();

    // Loop over all of the nodes in the element
    for (unsigned i=0;i<nnode;i++)
    {
     // Create a node pointer to store the i-th node in the element
     Node* node_pt=el_pt->node_pt(i);

     // Get the spatial position of this node
     for (unsigned k=0;k<3;k++)
     {
      // Store the k-th coordinate value in the vector, x
      x[k]=node_pt->x(k);
     }

     // Get the exact solution at this (Eulerian) position
     GlobalParameters::get_simple_exact_u(x,u);

     // Make sure each dof at this point is pinned (real and imaginary)
     node_pt->pin(0);
     node_pt->pin(1);

     // Set the solution value at this point
     node_pt->set_value(0,u[0]);
     node_pt->set_value(1,u[1]);
    }
   } // if(abs(x[0]-GlobalParameters::Centre) < 0.51 ... 
  } // if (el_pt!=0)
 } // for (unsigned e=0;e<n_element;e++)
        
 // If the user wishes to suppress all output restore the the stream pointer
 if (GlobalParameters::Output_management_flag==3)
 {
  // Now set the oomph_info stream pointer to the null stream to
  // disable all possible output
  oomph_info.stream_pt()=GlobalParameters::Stream_pt;
 }
  
 // Find the number of boundaries in the mesh
 unsigned n_bound=Bulk_mesh_pt->nboundary();
 
 // Loop over all boundaries
 for (unsigned b=0;b<n_bound;b++)
 {
  // Find the number of nodes on the b-th boundary
  unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

  // Loop over the nodes on the b-th boundary
  for(unsigned n=0;n<n_node;n++)
  {
   // All of these nodes sides are PMLs so pin to 0
   Node* boundary_node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

   // Pin the (real) dof at this node
   boundary_node_pt->pin(0);
   
   // Pin the (imaginary) dof at this node
   boundary_node_pt->pin(1);

   // Set the solution value at this point (real part)
   boundary_node_pt->set_value(0,0.0);
   
   // Set the solution value at this point (imaginary part)
   boundary_node_pt->set_value(1,0.0);
  }
 } // for(unsigned b=0;b<n_bound;b++)
} // End of apply_boundary_conditions

//==============================================start_of_enable_pmls======
/// Enable the PML mapping function for each node in the PML region
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::enable_pmls()
{
 // Find the number of elements in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Vector to hold the local coordinates of a point in any given element 
 Vector<double> s(3,0.0);
 
 // Vector to hold the (Eulerian) coordinates of a point
 Vector<double> x(3,0.0);

 // Vector to hold the real and imaginary part of the solution
 Vector<double> u(2,0.0);

 // Store the required coordinate of the inner boundary of the left PML; 
 // in any given direction this will be the value of Pml_thickness
 double left_boundary=GlobalParameters::Pml_thickness;
   
 // Store the required coordinate of the inner boundary of the right
 // PML; in the x-direction this will be the value of Lx-Pml_thickness
 // (or Ly-Pml_thickness in the y-direction and Lz-Pml_thickness in
 // the z-direction) but we assume the PML has the same thickness in
 // all directions
 double right_boundary=GlobalParameters::Lx-GlobalParameters::Pml_thickness;

 // Loop over the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // If we're using Jonathon's new test mapping
   if (GlobalParameters::Enable_test_pml_mapping_flag)
   {
    // Set the PML mapping function
    el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
   }
  
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);
   
   // If it's in the left (x-direction) PML region
   if (x[0]<=left_boundary)
    el_pt->enable_pml(0,left_boundary,0.0);

   // If it's in the right (x-direction) PML region
   if (x[0]>=right_boundary)
    el_pt->enable_pml(0,right_boundary,GlobalParameters::Lx);

   // If it's in the left (y-direction) PML region
   if (x[1]<=left_boundary)
    el_pt->enable_pml(1,left_boundary,0.0);

   // If it's in the right (y-direction) PML region
   if (x[1]>=right_boundary)
    el_pt->enable_pml(1,right_boundary,GlobalParameters::Ly);

   // If it's in the left (z-direction) PML region
   if (x[2]<=left_boundary)
    el_pt->enable_pml(2,left_boundary,0.0);

   // If it's in the right (z-direction) PML region
   if (x[2]>=right_boundary)
    el_pt->enable_pml(2,right_boundary,GlobalParameters::Lz);
  }
 } // for (unsigned e=0;e<n_element;e++)
} // End of enable_pmls

//======================================start_of_actions_after_adapt======
/// Actions after adapt: Re-apply the boundary conditions
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::actions_after_adapt()
{
 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // How many elements in the mesh?
 unsigned n_element=mesh_pt()->nelement();

 // Loop over the elements and pass a pointer to the value of k^2
 for (unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Set the wavenumber function pointer
   el_pt->k_squared_pt()=&GlobalParameters::K_squared;
  }
 } // for (unsigned e=0;e<n_element;e++)
 
 // Re-apply boundary conditions
 apply_boundary_conditions();

 // Re-enable the PML mapping in elements in the PML region
 enable_pmls();

 // Rebuild the mesh
 rebuild_global_mesh();
} // End of actions_after_adapt

//======================================================start_of_doc======
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLStructuredCubicHelmholtz<ELEMENT>::doc_solution()
{
 // Tell the user
 oomph_info << "\nDocumentation step: "
	    << GlobalParameters::Doc_info.number() << std::endl;
 
 // Create an output stream
 ofstream some_file;

 // Create space for the file name
 char filename[100];

 // Number of plot points
 unsigned npts=5;
 
 // Number of plot points in the coarse solution
 unsigned npts_coarse=2;

 // Output solution
 //-----------------
 sprintf(filename,"%s/soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Ouput exact solution
 //---------------------
 sprintf(filename,"%s/exact_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,npts,GlobalParameters::simple_exact_u_pt);
 some_file.close();

 // Output coarse solution
 //-----------------------
 sprintf(filename,"%s/coarse_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Compute error
 //--------------
 sprintf(filename,"%s/error%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
	 GlobalParameters::Doc_info.number());
 some_file.open(filename);
 
 //---------------------------------------------------------------------
 // To compute the norm of the error norm we first need to loop over all
 // of the elements in the mesh. Before we compute the norm of the error
 // in any element we need to make sure it doesn't lie in the PML region
 // or in the pinned region
 //--------------------------------------------------------------------- 
 // Variables to hold the L2 norm of the error in the solution
 double error=0.0;
 
 // Variable to hold the L2 norm of the solution
 double norm=0.0;
 
 // Vector to hold the local coordinates of a point in an element
 Vector<double> s(3,0.0);

 // Vector to hold the spatial position of a point in an element
 Vector<double> x(3,0.0);

 // Store the required coordinate of the inner boundary of the left PML; 
 // in any given direction this will be the value of Pml_thickness
 double left_boundary=GlobalParameters::Pml_thickness;
   
 // Store the required coordinate of the inner boundary of the right
 // PML; in the x-direction this will be the value of Lx-Pml_thickness
 // (or Ly-Pml_thickness in the y-direction and Lz-Pml_thickness in
 // the z-direction) but we assume the PML has the same thickness in
 // all directions
 double right_boundary=GlobalParameters::Lx-GlobalParameters::Pml_thickness;

 // Find out how many elements there are in the mesh
 unsigned n_element=Bulk_mesh_pt->nelement();

 // Loop over all of the elements in the mesh
 for (unsigned e=0;e<n_element;e++)
 { 
  // Variables to hold the L2 norm of the error in the elemental solution
  double el_error=0.0;
 
  // Variable to hold the L2 norm of the elemental solution
  double el_norm=0.0;
 
  // Upcast from GeneralisedElement to Helmholtz bulk element
  ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

  // If the upcast was successful
  if (el_pt!=0)
  {
   // Get the (Eulerian) coordinates of the centre of the element
   el_pt->get_x(s,x);

   // We only take the contribution from this element if it does
   // not lie in the PML region 
   if(x[0]<=left_boundary) continue;
   if(x[0]>=right_boundary) continue;
   if(x[1]<=left_boundary) continue;
   if(x[1]>=right_boundary) continue;
   if(x[2]<=left_boundary) continue;
   if(x[2]>=right_boundary) continue;

   // If it's in the (pinned) central region, ignore it 

   // Check if the element lies in the central cube region
   if (GlobalParameters::is_in_pinned_region(x))
   {
    // Skip to the next element
    continue;
   }

   // Otherwise, compute the L2 norm of the error over this element
   el_pt->compute_error(some_file,
			GlobalParameters::get_simple_exact_u,
			el_error,
			el_norm);

   // Update the global error norm value
   error+=el_error;
   
   // Update the global norm value
   norm+=el_norm;
  }
 } // for(unsigned e=0;e<n_element;e++)

 // Now close the file
 some_file.close();

 // Output the L2 norm of the error and the solution and then output
 // the relative error of the solution
 oomph_info << "\nSolution norm : " << norm
	    << "\nAbsolute error: " << error	   
	    << "\nRelative error: " << error/norm
	    << std::endl;
 
 // Write the L2 norm of the solution to the trace file
 Trace_file << norm << std::endl;

 // Increment the documentation number
 GlobalParameters::Doc_info.number()++;
} // End of doc_solution


//=====================================================start_of_main======
/// Solve 3D Helmholtz problem for a point source in a unit cube
//========================================================================
int main(int argc,char **argv)
{
 //------------------------
 // Command line arguments
 //------------------------
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Choose the number of nodes in one direction of an element;
 // Values of nnode_1d:
 //        2: Bilinear interpolation 
 //        3: Biquadratic interpolation 
 //        4: Bicubic interpolation 
 CommandLineArgs::specify_command_line_flag(
  "--nnode_1d",&GlobalParameters::Nnode_1d);

 // Choose the minimum level of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--min_ref",&GlobalParameters::Min_refinement_level);
 
 // Choose the additional levels of uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--add_ref",&GlobalParameters::Add_refinement_level);
 
 // Choose the maximum number of adaptive refinements
 CommandLineArgs::specify_command_line_flag(
  "--n_adapt",&GlobalParameters::N_adaptations);
 
 // Choose how many additional levels of uniform refinement to use
 CommandLineArgs::specify_command_line_flag(
  "--use_adapt",&GlobalParameters::Use_adaptation_flag);
  
 // Choose the value of k^2
 CommandLineArgs::specify_command_line_flag(
  "--k_sq",&GlobalParameters::K_squared);
 
 // Choose the value of the shift in the CSLP
 CommandLineArgs::specify_command_line_flag(
  "--alpha",&GlobalParameters::Alpha_shift);
 
 // Choose the value of the damping factor in the damped Jacobi solver
 CommandLineArgs::specify_command_line_flag(
  "--omega",&Smoother_Factory_Function_Helper::Omega);
  
 // Choose the pre-smoother
 CommandLineArgs::specify_command_line_flag(
  "--presmoother",&GlobalParameters::Pre_smoother_flag);
  
 // Choose the post-smoother
 CommandLineArgs::specify_command_line_flag(
  "--postsmoother",&GlobalParameters::Post_smoother_flag);
  
 // Choose the linear solver
 CommandLineArgs::specify_command_line_flag(
  "--linear_solver",&GlobalParameters::Linear_solver_flag);
 
 // Decide whether or not to suppress all or some of the MG solver output
 CommandLineArgs::specify_command_line_flag(
  "--output_flag",&GlobalParameters::Output_management_flag);
     
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--conv_flag",&GlobalParameters::Doc_convergence_flag);
 
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--test_pml_mapping",&GlobalParameters::Enable_test_pml_mapping_flag);
 
 // Parse command line
 CommandLineArgs::parse_and_assign();
  
 // Document what has been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Update any parameters that need to be updated
 GlobalParameters::update_parameters();
 
 //--------------------------------
 // Set the documentation directory
 //--------------------------------
 // Set output directory
 GlobalParameters::Doc_info.set_directory("RESLT");

 //-------------------
 // Set up the problem
 //-------------------
 // Initialise a null pointer to the class Problem 
 Problem* problem_pt=0;
 
 // Set the problem pointer depending on the input (defaulted to nnode_1d=2)
 if (GlobalParameters::Nnode_1d==2)
 {  
  // Set up the problem with refineable 3D eight-noded elements from the
  // QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<3,2> ELEMENT;
    
  // Set the problem pointer
  problem_pt=new PMLStructuredCubicHelmholtz<ELEMENT>;
 }
 else if (GlobalParameters::Nnode_1d==3)
 {
  // Set up the problem with refineable 3D twenty-seven-noded elements from
  // the QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<3,3> ELEMENT;
  
  // Set the problem pointer
  problem_pt=new PMLStructuredCubicHelmholtz<ELEMENT>;
 }  
 else if (GlobalParameters::Nnode_1d==4)
 {
  // Set up the problem with refineable 3D sixty-four-noded elements from
  // the QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<3,4> ELEMENT;
  
  // Set the problem pointer
  problem_pt=new PMLStructuredCubicHelmholtz<ELEMENT>;
 }
 else
 {
  // Throw an error otherwise
  throw OomphLibError("nnode_1d can only be 2,3 or 4.",
		      OOMPH_CURRENT_FUNCTION,
		      OOMPH_EXCEPTION_LOCATION);
 }

 //------------------ 
 // Solve the problem
 //------------------ 
 // If the user wishes to use adaptive refinement then we use the Newton
 // solver with a given argument to indicate how many adaptations to use
 if (GlobalParameters::Use_adaptation_flag)
 {
  // If the user wishes to silence everything
  if (GlobalParameters::Output_management_flag==3)
  {
   // Store the output stream pointer
   GlobalParameters::Stream_pt=oomph_info.stream_pt();

   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=&oomph_nullstream;
  }
    
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::Min_refinement_level;i++)
  { 
   oomph_info << "\n===================="
	      << "Initial Refinement"
	      << "====================\n"
	      << std::endl;

   // Add additional refinement
   problem_pt->refine_uniformly();
  }

  // If we silenced the adaptation, allow output again
  if (GlobalParameters::Output_management_flag==3)
  {
   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=GlobalParameters::Stream_pt;
  }
  
  // Solve the problem
  problem_pt->newton_solve();
   
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::N_adaptations;i++)
  { 
   // If the user wishes to silence everything
   if (GlobalParameters::Output_management_flag==3)
   {
    // Store the output stream pointer
    GlobalParameters::Stream_pt=oomph_info.stream_pt();

    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=&oomph_nullstream;
   }
   
   // Adapt the problem
   problem_pt->adapt();
   
   // If we silenced the adaptation, allow output again
   if (GlobalParameters::Output_management_flag==3)
   {
    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=GlobalParameters::Stream_pt;
   }
  
   // Solve the problem
   problem_pt->newton_solve();
  }
 }
 // If the user instead wishes to use uniform refinement
 else
 {
  // If the user wishes to silence everything
  if (GlobalParameters::Output_management_flag==3)
  {
   // Store the output stream pointer
   GlobalParameters::Stream_pt=oomph_info.stream_pt();

   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=&oomph_nullstream;
  }
  
  // Keep refining until the minimum refinement level is reached
  for (unsigned i=0;i<GlobalParameters::Min_refinement_level;i++)
  { 
   oomph_info << "\n===================="
	      << "Initial Refinement"
	      << "====================\n"
	      << std::endl;

   // Add additional refinement
   problem_pt->refine_uniformly();
  }
 
  // If we silenced the adaptation, allow output again
  if (GlobalParameters::Output_management_flag==3)
  {
   // Now set the oomph_info stream pointer to the null stream to
   // disable all possible output
   oomph_info.stream_pt()=GlobalParameters::Stream_pt;
  }
  
  // Solve the problem
  problem_pt->newton_solve();
 
  // Refine and solve until the additional refinements have been completed
  for (unsigned i=0;i<GlobalParameters::Add_refinement_level;i++)
  {
   // If the user wishes to silence everything
   if (GlobalParameters::Output_management_flag==3)
   {
    // Store the output stream pointer
    GlobalParameters::Stream_pt=oomph_info.stream_pt();

    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=&oomph_nullstream;
   }
  
   oomph_info << "==================="
	      << "Additional Refinement"
	      << "==================\n"
	      << std::endl;
 
   // Add additional refinement
   problem_pt->refine_uniformly();
  
   // If we silenced the adaptation, allow output again
   if (GlobalParameters::Output_management_flag==3)
   {
    // Now set the oomph_info stream pointer to the null stream to
    // disable all possible output
    oomph_info.stream_pt()=GlobalParameters::Stream_pt;
   }
  
   // Solve the problem
   problem_pt->newton_solve();
  }
 } // if (GlobalParameters::Use_adaptation_flag)
 
 // Delete the problem pointer
 delete problem_pt;

 // Make it a null pointer
 problem_pt=0;
} // End of main
