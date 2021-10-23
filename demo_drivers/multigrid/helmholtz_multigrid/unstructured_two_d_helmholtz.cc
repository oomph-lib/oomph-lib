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
// Driver for an adaptive 2D Helmholtz problem with perfectly matched layer
// treatment for the exterior boundaries. The default solution is a Bessel
// solution. The alternative is a plane wave incident at angle 'alpha'.

#include <fenv.h>
#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The Helmholtz equations (purely for the Hankel functions namespace)
#include "helmholtz.h"

// The pml Helmholtz equations and complex-valued multigrid machinery
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/quad_from_triangle_mesh.h"

using namespace oomph;
using namespace std;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//======start_of_namespace================================================
/// Namespace for the Helmholtz problem parameters
//========================================================================
namespace GlobalParameters
{
 /// Mesh specific parameters:
 ///---------------------------------
 /// Length of the square in each direction
 double Lx=1.0;
 double Ly=1.0;

 /// The x and y coordinate of the centre of the square 
 double Centre=Lx/2.0;
 
 /// The number of elements in the PML layer
 unsigned N_pml_element=1;

 /// The thickness of the PML layer defaults to 0.2, so 10% of the
 /// size of the physical domain
 double Pml_thickness=1.0e-01;
 
 /// The choice of whether or not to disable PMLs (the default is to
 /// enable them)
 ///    0 = Enable PMLs
 ///    1 = Disable PMLs
 unsigned Disable_pml_flag=0;
 
 /// The choice of whether or not to enable the new test mapping
 ///    1 = Enable test mapping
 ///    0 = Disable test mapping
 unsigned Enable_test_pml_mapping_flag=0;
 
 /// Problem specific parameters:
 ///------------------------------------
 /// Square of the wavenumber (also known as k^2)
 double K_squared=1600.0;

 /// Wavenumber (also known as k), k=omega/c
 double Wavenumber=sqrt(K_squared);

 /// Choose the value of the shift to create the complex-shifted
 /// Laplacian preconditioner (CSLP)
 double Alpha_shift=0.5;
 
 /// Update the parameters passed in at the command line
 void update_parameters()
 {
  /// Wavenumber (also known as k), k=omega/c
  Wavenumber=sqrt(K_squared);
 }
 
 /// Value of the solution on the boundary of the obstacle (here
 /// we assume the solution is a plane wave incident at angle alpha)
 void get_exact_u(const Vector<double>& x,
		  Vector<double>& u,
		  const double& alpha=0.25*MathematicalConstants::Pi)
 {
  // Set the first entry
  u[0]=-cos(Wavenumber*(x[0]*cos(alpha)+x[1]*sin(alpha)));

  // Set the second entry
  u[1]=-sin(Wavenumber*(x[0]*cos(alpha)+x[1]*sin(alpha)));
 } // End of get_exact_u
 
 // The exact solution to the Green’s function for the two-
 // dimensional Helmholtz equation [see C.M.Linton - The Green’s
 // function for the two-dimensional Helmholtz equation in
 // periodic domains]. The solution is G=-(i/4)*H_0^{(1)}(kr)
 // or equivalently, G=(1/4)*Y_0(kr)-(i/4)*J_0(kr), where
 // J_0 is the Bessel function of the first kind and Y_0 is the
 // Bessel function of the second kind
 void get_exact_u_bessel(const Vector<double>& x,Vector<double>& u)
 {
  // The radius in polar coordinates
  double r;
   
  // Switch to polar coordinates
  r=sqrt((x[0]-Centre)*(x[0]-Centre)+(x[1]-Centre)*(x[1]-Centre));

  // Tolerance for the point being in the centre
  double tol=1.0e-06;

  // Check if the point lies practically in the centre
  if (r<tol)
  {
   // Set the real and imaginary parts of the solution. The real
   // part is set arbitrarily since the Bessel function of the
   // second kind is unbounded at the centre of the mesh. The
   // real part of the solution is (1/4)*Y_0(0) and the imaginary
   // part is -(i/4)*J_0(0). We know Y_0(0)=-inf and J_0(0)=1
   // so the imaginary part of the solution is precisely -0.25i
   // but the real part we can set arbitrarily. Note, setting it
   // to -DBL_MAX might cause some problems so we choose some
   // finite value
   u[0]=-0.5;
   u[1]=-0.25;
  }
  else
  {
   // Scale r by the wavenumber
   double rr=Wavenumber*r;

   // Number of first-order Hankel function terms (starts from n=0
   // in H^(1)_n(...))
   unsigned n_terms=1;

   // Hankel function and the its derivative
   Vector<std::complex<double> > h(1);
   Vector<std::complex<double> > hp(1);
   
   // Calculate the solution (i/4)*H^(1)_0(kr) at kr
   Hankel_functions_for_helmholtz_problem::Hankel_first(n_terms,rr,h,hp);
   
   // Calculate the real and imaginary parts of the solution
   u[0]=0.25*imag(h[0]);
   u[1]=-0.25*real(h[0]);
  }
 } // End of get_exact_u_bessel
 
 /// Default value of the solution on the boundary of the obstacle.
 /// This represents a Bessel solution
 void default_get_exact_u(const Vector<double>& x,Vector<double>& u)
 {
  // Set the first entry
  u[0]=0.1;

  // Set the second entry
  u[1]=0.0;
 } // End of default_get_exact_u
 
 // Set the exact solution pointer to the get_simple_exact_u function above
 FiniteElement::SteadyExactSolutionFctPt simple_exact_u_pt=&get_exact_u_bessel;
 
 /// The number of nodes in one direction (default=2)
 unsigned Nnode_1d=2;
 
 /// The minimum level of uniform refinement
 unsigned Min_refinement_level=2;
 
 /// The additional levels of uniform refinement 
 unsigned Add_refinement_level=0;
 
 /// The number of segments used to define the circular boundary
 unsigned N_boundary_segment=6;
 
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
 
 // Pointer to the output stream -- defaults to std::cout
 std::ostream* Stream_pt;
 
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
    // If we're disabling PMLs
    if (!Disable_pml_flag)
    {   
     // The "effective k^2" is shifted, so we shift the k used in the
     // transformation too
     std::complex<double> k_shifted=
      sqrt(k_squared_local*std::complex<double>(1.0,alpha_shift));
         
     // Return the gamma in J++, with the shifted k
     return (1.0/k_shifted)*std::complex<double>
      (0.0,1.0/(std::fabs(pml_width_i-nu_i)));
    }
    else
    {
     // Otherwise just return the value 1.0
     return 1.0;
    }
   } // End of gamma
 }; // End of TestPMLMapping

 /// Set the new PML mapping
 TestPMLMapping* Test_pml_mapping_pt=new TestPMLMapping;
} // End of namespace


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//========= start_of_problem_class========================================
/// HelmholtzMGProblem class to demonstrate use of perfectly matched
/// layers for Helmholtz problems.
//========================================================================
template<class ELEMENT>
class PMLHelmholtzMGProblem : public HelmholtzMGProblem
{

public:

 /// Constructor
 PMLHelmholtzMGProblem();

 /// Destructor (empty)
 ~PMLHelmholtzMGProblem();

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
  {
   // Document the solution
   doc_solution();
  }

 /// Doc the solution. DocInfo object stores flags/labels for where
 /// the output gets written to
 void doc_solution();

 /// Create PML meshes
 void create_pml_meshes();

 // Apply boundary conditions
 void apply_boundary_conditions();

 /// Actions before adapt: Wipe the PML meshes
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the PML meshes
 void actions_after_adapt();

 /// Set GMRES preconditioner by multigrid as the linear solver
 void set_gmres_multigrid_solver();
 
private:

 /// Pointer to the refineable "bulk" mesh
 RefineableQuadFromTriangleMesh<ELEMENT>* Bulk_mesh_pt;

 /// Overload the make_new_problem function to return an object of this class
 HelmholtzMGProblem* make_new_problem()
  {
   // Return a new problem pointer
   return new PMLHelmholtzMGProblem<ELEMENT>;
  }

 /// Overload the mg_bulk_mesh_pt function to return a pointer to the
 /// "refineable" portion of the mesh
 TreeBasedRefineableMeshBase* mg_bulk_mesh_pt()
  {
   // Return the pointer to the bulk mesh
   return Bulk_mesh_pt;
  }
 
 /// Pointer to the right PML mesh
 Mesh* PML_right_mesh_pt;

 /// Pointer to the top PML mesh
 Mesh* PML_top_mesh_pt;

 /// Pointer to the left PML mesh
 Mesh* PML_left_mesh_pt;

 /// Pointer to the bottom PML mesh
 Mesh* PML_bottom_mesh_pt;

 /// Pointer to the top right corner PML mesh
 Mesh* PML_top_right_mesh_pt;

 /// Pointer to the top left corner PML mesh
 Mesh* PML_top_left_mesh_pt;

 /// Pointer to the bottom right corner PML mesh
 Mesh* PML_bottom_right_mesh_pt;

 /// Pointer to the bottom left corner PML mesh
 Mesh* PML_bottom_left_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // End of problem class

//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLHelmholtzMGProblem<ELEMENT>::PMLHelmholtzMGProblem()
{
 // Create boundary information: Outer boundary
 //--------------------------------------------
 TriangleMeshClosedCurve* outer_boundary_pt=0;

 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(4);

 // Each polyline only has three vertices, provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord(2);
 for(unsigned i=0;i<2;i++)
 {
  vertex_coord[i].resize(2);
 }

 // First polyline
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=GlobalParameters::Lx;
 vertex_coord[1][1]=0.0;

 // Build the 1st boundary polyline
 unsigned boundary_id=2;
 outer_boundary_line_pt[0]=new TriangleMeshPolyLine(vertex_coord,boundary_id);

 // Second boundary polyline
 vertex_coord[0][0]=GlobalParameters::Lx;
 vertex_coord[0][1]=0.0;
 vertex_coord[1][0]=GlobalParameters::Lx;
 vertex_coord[1][1]=GlobalParameters::Ly;

 // Build the 2nd boundary polyline
 boundary_id=3;
 outer_boundary_line_pt[1]=new TriangleMeshPolyLine(vertex_coord,boundary_id);

 // Third boundary polyline
 vertex_coord[0][0]=GlobalParameters::Lx;
 vertex_coord[0][1]=GlobalParameters::Ly;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=GlobalParameters::Ly;

 // Build the 3rd boundary polyline
 boundary_id=4;
 outer_boundary_line_pt[2]=new TriangleMeshPolyLine(vertex_coord,boundary_id);

 // Fourth boundary polyline
 vertex_coord[0][0]=0.0;
 vertex_coord[0][1]=GlobalParameters::Ly;
 vertex_coord[1][0]=0.0;
 vertex_coord[1][1]=0.0;

 // Build the 4th boundary polyline
 boundary_id=5;
 outer_boundary_line_pt[3]=new TriangleMeshPolyLine(vertex_coord,boundary_id);

 // Create the triangle mesh polygon for outer boundary
 outer_boundary_pt=new TriangleMeshPolygon(outer_boundary_line_pt);

 // Create boundary information: Inner boundary
 //--------------------------------------------
 // Create circle representing inner boundary. We are comparing this problem
 // to a problem with a pinned square. Instead of choosing the value of r
 // such that the diameter of the circle is the same as the side-length of
 // the pinned square we find the value of r such that the area of both
 // regions are equivalent so r=sqrt(0.2^{2}/pi)=0.1128...
 double a=0.1128;
 double x_c=GlobalParameters::Centre;
 double y_c=GlobalParameters::Centre;
 Circle* inner_circle_pt=new Circle(x_c,y_c,a);

 // Create storage for curves which will represent the obstacle
 Vector<TriangleMeshCurveSection*> inner_boundary_line_pt(2);

 // # of segments for each curve (defining the obstacle in the coarsest mesh)
 unsigned n_segments=GlobalParameters::N_boundary_segment;
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start=0.0;
 double s_end=MathematicalConstants::Pi;
 boundary_id=0;
 inner_boundary_line_pt[0]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);

 // The intrinsic coordinates for the beginning and end of the curve
 s_start=MathematicalConstants::Pi;
 s_end=2.0*MathematicalConstants::Pi;
 boundary_id=1;
 inner_boundary_line_pt[1]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);

 // Combine to hole
 //----------------
 Vector<TriangleMeshClosedCurve*> hole_pt(1);
 Vector<double> hole_coords(2);
 hole_coords[0]=GlobalParameters::Centre;
 hole_coords[1]=GlobalParameters::Centre;
 hole_pt[0]=new TriangleMeshClosedCurve(inner_boundary_line_pt,hole_coords);

 // Use the TriangleMeshParameters object for helping on the manage
 // of the TriangleMesh parameters. The only parameter that needs to take
 // is the outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt()=hole_pt;

 // Target element size in bulk mesh
 triangle_mesh_parameters.element_area()=0.10;

 // Build the mesh using the specified parameters
 //----------------------------------------------
 // Build adaptive "bulk" mesh
 Bulk_mesh_pt=new RefineableQuadFromTriangleMesh<ELEMENT>(
  triangle_mesh_parameters);

 // Output
 Bulk_mesh_pt->output("RESLT/mesh.dat",2);
 
 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=0.00004;
 Bulk_mesh_pt->max_permitted_error()=0.0001;
 
 // Create the main triangular mesh
 add_sub_mesh(Bulk_mesh_pt);
 
 // If we're using at least one PML element
 if (0!=GlobalParameters::N_pml_element)
 {
  // Create PML meshes and add them to the global mesh
  create_pml_meshes();
 }

 // Build the entire mesh from its submeshes
 build_global_mesh();
 
 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // Find out how many elements there are in the mesh
 unsigned n_element=this->mesh_pt()->nelement();

 // Loop over the elements in the mesh
 for(unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to Helmholtz bulk element
  PMLHelmholtzEquations<2>* el_pt =
   dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

  // Set the k_squared double pointer
  el_pt->k_squared_pt()=&GlobalParameters::K_squared;

  // If we're using Jonathon's new test mapping
  if (GlobalParameters::Enable_test_pml_mapping_flag)
  {
   // Set the PML mapping function
   el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
  }
  
  // If we're using GMRES & MG as the linear solver
  if (GlobalParameters::Disable_pml_flag)
  {
   // Disable the PML-ification in these layers
   dynamic_cast<PMLHelmholtzEquations<2>*>
    (mesh_pt()->element_pt(e))->disable_pml();
  }
 } // for(unsigned e=0;e<n_element;e++)

 // Apply boundary conditions
 apply_boundary_conditions();
  
 // Setup equation numbering scheme
 assign_eqn_numbers();
 
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
} // End of the constructor, PMLHelmholtzMGProblem

//=======start_of_destructor==============================================
/// Destructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLHelmholtzMGProblem<ELEMENT>::~PMLHelmholtzMGProblem()
{
 // If we're using GMRES & MG as the linear solver
 if (GlobalParameters::Linear_solver_flag==1)
 {
  // Delete the MG solver pointers
  delete dynamic_cast<HelmholtzGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt();

  // Set the pointer to null
  dynamic_cast<HelmholtzGMRESMG<CRDoubleMatrix>* >
   (linear_solver_pt())->preconditioner_pt()=0;
  
  // Delete the MG solver pointers
  delete linear_solver_pt();

  // Set the pointer to null
  linear_solver_pt()=0;    
 }
   
 // Delete the error estimator (allocated as a 'new' object)
 delete Bulk_mesh_pt->spatial_error_estimator_pt();
 
 // Set the pointer to null
 Bulk_mesh_pt->spatial_error_estimator_pt()=0;

 // Delete the "bulk" mesh
 delete Bulk_mesh_pt;

 // Set the pointer to null
 Bulk_mesh_pt=0;
 
} // End of ~PMLHelmholtzMGProblem


//=======set_gmres_multigrid_solver=======================================
/// Build and set GMRES preconditioner by multigrid as the linear solver
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::set_gmres_multigrid_solver()
{
 // Create linear solver
 HelmholtzGMRESMG<CRDoubleMatrix>* solver_pt=
  new HelmholtzGMRESMG<CRDoubleMatrix>;

 // Use RHS preconditioning
 solver_pt->set_preconditioner_RHS();
 
 // Set the maximum number of iterations
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
 HelmholtzMGPreconditioner<2>* prec_pt=new HelmholtzMGPreconditioner<2>(this);

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
  // Disable smoother timings
  prec_pt->disable_doc_time();
 }
 else if (GlobalParameters::Output_management_flag==2)
 {
  // Disable all output from the V-cycle
  prec_pt->disable_v_cycle_output();
 }
 else if (GlobalParameters::Output_management_flag==3)
 {
  // Suppress all output
  prec_pt->disable_output();
 }
} // End of set_gmres_multigrid_solver

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::actions_before_adapt()
{
 // If we're using at least one PML element
 if (0!=GlobalParameters::N_pml_element)
 {
  // Before adapting the added PML meshes must be removed
  // as they are not refineable and are to be rebuilt from the
  // newly refined triangular mesh
  delete PML_right_mesh_pt;
  PML_right_mesh_pt=0;
  delete PML_top_mesh_pt;
  PML_top_mesh_pt=0;
  delete PML_left_mesh_pt;
  PML_left_mesh_pt=0;
  delete PML_bottom_mesh_pt;
  PML_bottom_mesh_pt=0;
  delete PML_top_right_mesh_pt;
  PML_top_right_mesh_pt=0;
  delete PML_top_left_mesh_pt;
  PML_top_left_mesh_pt=0;
  delete PML_bottom_right_mesh_pt;
  PML_bottom_right_mesh_pt=0;
  delete PML_bottom_left_mesh_pt;
  PML_bottom_left_mesh_pt=0;

  // Rebuild the HelmholtzMGProblem's global mesh from its various sub-meshes
  // but first flush all its submeshes
  flush_sub_meshes();

  // Then add the triangular mesh back
  add_sub_mesh(Bulk_mesh_pt);

  //  Rebuild the global mesh such that it now stores
  // the triangular mesh only
  rebuild_global_mesh();
 }
} // End of actions_before_adapt

//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::actions_after_adapt()
{
 // Re-build the full mesh with PMLs:
 //----------------------------------
 // If we're using at least one PML element
 if (0!=GlobalParameters::N_pml_element)
 {
  // Create PML meshes and add them to the global mesh
  create_pml_meshes(); 
 }
 
 // Build the entire mesh from its submeshes
 rebuild_global_mesh();
  
 // Complete the build of all elements so they are fully functional: 
 //-----------------------------------------------------------------
 // Calculate the number of elements in the mesh
 unsigned n_element=this->mesh_pt()->nelement();

 // Loop over the elements
 for(unsigned e=0;e<n_element;e++)
 {
  // Upcast from GeneralisedElement to PMLHelmholtz bulk element
  PMLHelmholtzEquations<2> *el_pt =
   dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

  // Set the frequency function pointer
  el_pt->k_squared_pt()=&GlobalParameters::K_squared;
  
  // If we're using Jonathon's new test mapping
  if (GlobalParameters::Enable_test_pml_mapping_flag)
  {
   // Set the PML mapping function
   el_pt->pml_mapping_pt()=GlobalParameters::Test_pml_mapping_pt;
  }
  
  // If we're using GMRES & MG as the linear solver
  if (GlobalParameters::Disable_pml_flag)
  {
   // Disable the PML-ification in these layers
   dynamic_cast<PMLHelmholtzEquations<2>*>
    (mesh_pt()->element_pt(e))->disable_pml();
  }
 } // for(unsigned e=0;e<n_element;e++)
 
 // Re-apply boundary conditions
 apply_boundary_conditions();
} // End of actions_after_adapt

//==================start_of_apply_boundary_conditions====================
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::apply_boundary_conditions()
{
 // Boundary conditions are set on the surface of the circle
 // as a constant nonzero Dirichlet boundary condition
 unsigned n_bound=Bulk_mesh_pt->nboundary();

 // Vector to hold the coordinates of the node
 Vector<double> x(2,0.0);
 
 // Vector to hold the solution at the given node (first entry holds
 // the real part and the second entry holds the imaginary part)
 Vector<double> u(2,0.0);

 // Loop over the boundaries
 for(unsigned b=0;b<n_bound;b++)
 {
  // Find the number of nodes on the b-th boundary
  unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

  // Loop over the nodes on boundary b
  for (unsigned n=0;n<n_node;n++)
  {
   // If we're on the boundary of the obstacle (remembering that
   // the circle has been broken up into two semi-circle boundaries)
   if ((0==b) || (1==b))
   {
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Get the coordinates of the node
    for (unsigned k=0;k<2;k++)
    {
     // Store the k-th coordinate value
     x[k]=nod_pt->x(k);
    }

    // Get the solution at this point (default setting -- Bessel solution)
    GlobalParameters::simple_exact_u_pt(x,u);
    
    // Set the values at each dof
    nod_pt->set_value(0,u[0]);
    nod_pt->set_value(1,u[1]);
   }
  } // for (unsigned n=0;n<n_node;n++)
 } // for(unsigned b=0;b<n_bound;b++)
 
 // If we're not using any PML elements
 if (0==GlobalParameters::N_pml_element)
 {
  // Loop over the outer boundaries (not the circle boundary)
  for(unsigned b=2;b<n_bound;b++)
  {
   // Find the number of nodes on the b-th boundary
   unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

   // Loop over the nodes on boundary b
   for (unsigned n=0;n<n_node;n++)
   {
    // Grab the n-th node on the b-th boundary
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

    // Pin the first dof at this node
    nod_pt->pin(0);

    // Pin the second dof at this node
    nod_pt->pin(1);

    // Set the value of the first dof
    nod_pt->set_value(0,0.0);
     
    // Set the value of the second dof
    nod_pt->set_value(1,0.0);
   } // for (unsigned n=0;n<n_node;n++)
  } // for(unsigned b=0;b<n_bound;b++)
 } // if (GlobalParameters::Disable_pml_flag)
} // End of apply_boundary_conditions

//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::doc_solution()
{
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
 
 // If we're using at least one PML element
 if (0!=GlobalParameters::N_pml_element)
 {  
  // Output solution within pml domains
  //-----------------------------------
  sprintf(filename,"%s/pml_soln%i.dat",
	  GlobalParameters::Doc_info.directory().c_str(),
	  GlobalParameters::Doc_info.number());
  some_file.open(filename);
  PML_top_mesh_pt->output(some_file,npts);
  PML_right_mesh_pt->output(some_file,npts);
  PML_bottom_mesh_pt->output(some_file,npts);
  PML_left_mesh_pt->output(some_file,npts);
  PML_top_right_mesh_pt->output(some_file,npts);
  PML_bottom_right_mesh_pt->output(some_file,npts);
  PML_top_left_mesh_pt->output(some_file,npts);
  PML_bottom_left_mesh_pt->output(some_file,npts);
  some_file.close();
 }

 // Output coarse solution
 //-----------------------
 sprintf(filename,"%s/coarse_soln%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
         GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // If we're using at least one PML element
 if (0!=GlobalParameters::N_pml_element)
 {  
  // Output coarse solution within pml domains
  //------------------------------------------
  sprintf(filename,"%s/coarse_pml_soln%i.dat",
	  GlobalParameters::Doc_info.directory().c_str(),
	  GlobalParameters::Doc_info.number());
  some_file.open(filename);
  PML_top_mesh_pt->output(some_file,npts_coarse);
  PML_right_mesh_pt->output(some_file,npts_coarse);
  PML_bottom_mesh_pt->output(some_file,npts_coarse);
  PML_left_mesh_pt->output(some_file,npts_coarse);
  PML_top_right_mesh_pt->output(some_file,npts_coarse);
  PML_bottom_right_mesh_pt->output(some_file,npts_coarse);
  PML_top_left_mesh_pt->output(some_file,npts_coarse);
  PML_bottom_left_mesh_pt->output(some_file,npts_coarse);
  some_file.close();
 }
 
 // Increment the documentation counter
 GlobalParameters::Doc_info.number()++;

 // Doc error and solution norm
 //---------------------------- 
 double error,norm;
 sprintf(filename,"%s/error%i.dat",
	 GlobalParameters::Doc_info.directory().c_str(),
         GlobalParameters::Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
			     GlobalParameters::simple_exact_u_pt,
			     error,norm); 
 some_file.close();
 
 // Document the L2 norm of the error and the L2 norm of the solution
 cout << "\nNorm of error         : " << error
      << "\nNorm of solution      : " << norm
      << "\nNormalised error norm : " << error/norm << std::endl;

 // Write the L2 norm of the solution to the trace file
 Trace_file << norm << std::endl;
} // End of doc

//============start_of_create_pml_meshes==================================
/// Create PML meshes and add them to the problem's sub-meshes
//========================================================================
template<class ELEMENT>
void PMLHelmholtzMGProblem<ELEMENT>::create_pml_meshes()
{
 // Bulk mesh bottom boundary id
 unsigned int bottom_boundary_id=2;

 // Bulk mesh right boundary id
 unsigned int right_boundary_id=3;

 // Bulk mesh top boundary id
 unsigned int top_boundary_id=4;

 // Bulk mesh left boundary id
 unsigned int left_boundary_id=5;

 // PML width in elements for the right layer
 unsigned n_x_right_pml=GlobalParameters::N_pml_element;

 // PML width in elements for the top layer
 unsigned n_y_top_pml=GlobalParameters::N_pml_element;

 // PML width in elements for the left layer
 unsigned n_x_left_pml=GlobalParameters::N_pml_element;

 // PML width in elements for the left layer
 unsigned n_y_bottom_pml=GlobalParameters::N_pml_element;

 // Outer physical length of the PML layers
 double width_x_right_pml=GlobalParameters::Pml_thickness;
 double width_y_top_pml=GlobalParameters::Pml_thickness;
 double width_x_left_pml=GlobalParameters::Pml_thickness;
 double width_y_bottom_pml=GlobalParameters::Pml_thickness;

 // Build the PML meshes based on the new adapted triangular mesh
 PML_right_mesh_pt=
  TwoDimensionalPMLHelper::create_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt,right_boundary_id,
   n_x_right_pml, width_x_right_pml);
 PML_top_mesh_pt=
  TwoDimensionalPMLHelper::create_top_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt, top_boundary_id,
   n_y_top_pml, width_y_top_pml);
 PML_left_mesh_pt=
  TwoDimensionalPMLHelper::create_left_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt, left_boundary_id,
   n_x_left_pml, width_x_left_pml);
 PML_bottom_mesh_pt=
  TwoDimensionalPMLHelper::create_bottom_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt, bottom_boundary_id,
   n_y_bottom_pml, width_y_bottom_pml);

 // Add submeshes to the global mesh
 add_sub_mesh(PML_right_mesh_pt);
 add_sub_mesh(PML_top_mesh_pt);
 add_sub_mesh(PML_left_mesh_pt);
 add_sub_mesh(PML_bottom_mesh_pt);

 // Rebuild corner PML meshes
 PML_top_right_mesh_pt=
  TwoDimensionalPMLHelper::create_top_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_right_mesh_pt, PML_top_mesh_pt,
   Bulk_mesh_pt, right_boundary_id);

 PML_bottom_right_mesh_pt=
  TwoDimensionalPMLHelper::create_bottom_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_right_mesh_pt, PML_bottom_mesh_pt,
   Bulk_mesh_pt, right_boundary_id);

 PML_top_left_mesh_pt=
  TwoDimensionalPMLHelper::create_top_left_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_left_mesh_pt, PML_top_mesh_pt,
   Bulk_mesh_pt, left_boundary_id);

 PML_bottom_left_mesh_pt=
  TwoDimensionalPMLHelper::create_bottom_left_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_left_mesh_pt, PML_bottom_mesh_pt,
   Bulk_mesh_pt, left_boundary_id);

 // Add submeshes to the global mesh
 add_sub_mesh(PML_top_right_mesh_pt);
 add_sub_mesh(PML_bottom_right_mesh_pt);
 add_sub_mesh(PML_top_left_mesh_pt);
 add_sub_mesh(PML_bottom_left_mesh_pt);

} // End of create_pml_meshes

//==========start_of_main=================================================
/// Solve 2D Helmholtz problem
//========================================================================
int main(int argc, char **argv)
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

 // Set the number of elements in the PML layer
 CommandLineArgs::specify_command_line_flag(
  "--npml_element",&GlobalParameters::N_pml_element);
 
 // Set the thickness of the pml
 CommandLineArgs::specify_command_line_flag(
  "--pml_thickness",&GlobalParameters::Pml_thickness);
 
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--nboundary_segment",&GlobalParameters::N_boundary_segment);
 
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
 
 // Decide whether or not to display convergence information
 CommandLineArgs::specify_command_line_flag(
  "--disable_pml",&GlobalParameters::Disable_pml_flag); 
 
 // Parse command line
 CommandLineArgs::parse_and_assign();
  
 // Document what has been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Update parameters
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
  // Set up the problem with refineable 2D four-node elements from the
  // QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<2,2> ELEMENT;
    
  // Set the problem pointer
  problem_pt=new PMLHelmholtzMGProblem<ELEMENT>;
 }
 else if (GlobalParameters::Nnode_1d==3)
 {
  // Set up the problem with refineable 2D nine-node elements from the
  // QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<2,3> ELEMENT;
  
  // Set the problem pointer
  problem_pt=new PMLHelmholtzMGProblem<ELEMENT>;
 }  
 else if (GlobalParameters::Nnode_1d==4)
 {
  // Set up the problem with refineable 2D sixteen-node elements from the
  // QPMLHelmholtzElement family
  typedef RefineableQPMLHelmholtzElement<2,4> ELEMENT;
  
  // Set the problem pointer
  problem_pt=new PMLHelmholtzMGProblem<ELEMENT>;
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
	      << "===================="
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
	      << "===================="
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
	      << "=================="
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
