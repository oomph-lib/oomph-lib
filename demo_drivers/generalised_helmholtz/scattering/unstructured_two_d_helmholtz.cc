//LIC//======================================================================
///LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
//LIC//
//LIC//======================================================================
/// \short Driver for a specific 2D Helmholtz problem with 
/// perfectly matched layer treatment for the exterior boundaries 

#include<fenv.h>

#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The generalised Helmholtz equations
#include "generalised_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;

/// \short Adaptivity works but for scattering problems in general
/// this is a questionable upgrade
//#define ADAPTIVE

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{
 /// \short Number of terms used in the computation 
 /// of the exact solution
 unsigned N_fourier=10;

 /// Radius of outer boundary (must be a circle!)
 double Outer_radius=1.5;

 /// Imaginary unit 
 std::complex<double> I(0.0,1.0);
 
 /// \short Exact solution for scattered field 
 /// (vector returns real and impaginary parts).
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Switch to polar coordinates
  double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta;
  theta=atan2(x[1],x[0]);
  
  // Argument for Bessel/Hankel functions
//  double rr=sqrt(K_squared)*r;  
 
//   // Evaluate Bessel/Hankel functions
//   complex <double > u_ex(0.0,0.0);
//   Vector<double> jn(N_fourier+1), yn(N_fourier+1),
//    jnp(N_fourier+1), ynp(N_fourier+1);
//   Vector<double> jn_a(N_fourier+1),yn_a(N_fourier+1),
//    jnp_a(N_fourier+1), ynp_a(N_fourier+1);
//   Vector<complex<double> > h(N_fourier+1),h_a(N_fourier+1),
//    hp(N_fourier+1), hp_a(N_fourier+1);

  // We want to compute N_fourier terms but the function
  // may return fewer than that.
//int n_actual=0;
//   CRBond_Bessel::bessjyna(N_fourier,sqrt(K_squared),n_actual,
//                           &jn_a[0],&yn_a[0],
//                           &jnp_a[0],&ynp_a[0]); 

  // Shout if things went wrong  
//!! Cleanup #ifdef PARANOID
//   if (n_actual!=int(N_fourier))
//    {
//     std::ostringstream error_stream; 
//     error_stream << "CRBond_Bessel::bessjyna() only computed "
//                  << n_actual << " rather than " << N_fourier 
//                  << " Bessel functions.\n";    
//     throw OomphLibError(error_stream.str(),
//                         "GlobalParameters::prescribed_incoming_flux()",
//                         OOMPH_EXCEPTION_LOCATION);
//    }
// #endif

  // Evaluate Hankel at actual radius
  //Hankel_functions_for_helmholtz_problem::Hankel_first(N_fourier,rr,h,hp);

  // Evaluate Hankel at inner (unit) radius
//Hankel_functions_for_helmholtz_problem::Hankel_first(N_fourier
//                                                        ,sqrt(K_squared),
//                                                        h_a,hp_a);
  
  // Compute the sum: Separate the computation of the negative 
  // and positive terms
//!! Cleanup   for (unsigned i=0;i<N_fourier;i++)
//    {
//     u_ex-=pow(I,i)*h[i]*((jnp_a[i])/hp_a[i])*pow(exp(I*theta),i);
//    }
//   for (unsigned i=1;i<N_fourier;i++)
//    {
//     u_ex-=pow(I,i)*h[i]*((jnp_a[i])/hp_a[i])*pow(exp(-I*theta),i);
//    }
  
//   // Get the real & imaginary part of the result
//   u[0]=real(u_ex);
//   u[1]=imag(u_ex);
  
 }// end of get_exact_u
 


 /// \short Flux (normal derivative) on the unit disk
 /// for a planar incoming wave
 void prescribed_incoming_flux(const Vector<double>& x, 
                               complex<double>& flux)
 {
//   // Switch to polar coordinates
// double r;
//   r=sqrt(x[0]*x[0]+x[1]*x[1]);
//   double theta;
//   theta=atan2(x[1],x[0]);
  
//   // Argument of the Bessel/Hankel fcts
//   double rr=sqrt(K_squared)*r;  
  
//   // Compute Bessel/Hankel functions
//   Vector<double> jn(N_fourier+1), yn(N_fourier+1),
//    jnp(N_fourier+1), ynp(N_fourier+1);

//   // We want to compute N_fourier terms but the function
//   // may return fewer than that.
//   int n_actual=0;
//   CRBond_Bessel::bessjyna(N_fourier,rr,n_actual,&jn[0],&yn[0],
//                           &jnp[0],&ynp[0]);
  
//   // Shout if things went wrong...
// #ifdef PARANOID
//   if (n_actual!=int(N_fourier))
//    {
//     std::ostringstream error_stream; 
//     error_stream << "CRBond_Bessel::bessjyna() only computed "
//                  << n_actual << " rather than " << N_fourier 
//                  << " Bessel functions.\n";    
//     throw OomphLibError(error_stream.str(),
//                         "GlobalParameters::prescribed_incoming_flux()",
//                         OOMPH_EXCEPTION_LOCATION);
//    }
// #endif
  
//   // Compute the sum: Separate the computation of the negative and 
//   // positive terms
//   flux=std::complex<double>(0.0,0.0);
//   for (unsigned i=0;i<N_fourier;i++)
//    {
//     flux+=pow(I,i)*(sqrt(K_squared))*pow(exp(I*theta),i)*jnp[i];
//    }
//   for (unsigned i=1;i<N_fourier;i++)
//    {
//     flux+=pow(I,i)*(sqrt(K_squared))*pow(exp(-I*theta),i)*jnp[i];
//    }


 }// end of prescribed_incoming_flux 

 /// Setting of frequency across the entire domain
 double Omega = sqrt(500.0);

 /// \short Variable soundspeed across the domain
 /// defaults set to constant across the domain 
 void c_function(const Vector<double>& x, double& c)
 {
  c = sqrt(10.0);
 }

 /// \short Absorption field function, 
 /// defaults set to constant across the domain 
 void alpha_function(const Vector<double>& x, double& alpha)
 {
  alpha = 0.0;
 }

} // end of namespace


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class to compute scattering of planar wave from unit disk
//=====================================================================
template<class ELEMENT> 
class ScatteringProblem : public Problem
{

public:
 
 /// Constructor
 ScatteringProblem();
 
 /// Destructor (empty)
 ~ScatteringProblem(){}

 /// \short Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve(){} 

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve()
  {
   DoubleVector residuals;
   get_residuals(residuals);
   residuals.output("residuals.dat");
  }

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();
 
 /// \short Create Helmholtz flux elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified surface Mesh 
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const & helmholtz_inner_boundary_mesh_pt);
 
 /// \short Delete boundary face elements and wipe the surface mesh
 void delete_face_elements( Mesh* const & boundary_mesh_pt);
 
 /// \short Set pointer to prescribed-flux function for all
 /// elements in the surface mesh on the surface of the unit disk
 void set_prescribed_incoming_flux_pt();

#ifdef ADAPTIVE

 /// Pointer to the "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* bulk_mesh_pt;

#else

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* bulk_mesh_pt;

#endif

 /// Pointer to the right PML mesh
 Mesh* pml_right_mesh_pt;

 /// Pointer to the top PML mesh
 Mesh* pml_top_mesh_pt;

 /// Pointer to the left PML mesh
 Mesh* pml_left_mesh_pt;

 /// Pointer to the bottom PML mesh
 Mesh* pml_bottom_mesh_pt;

 /// Pointer to the top right corner PML mesh
 Mesh* pml_top_right_mesh_pt;

 /// Pointer to the top left corner PML mesh
 Mesh* pml_top_left_mesh_pt;

 /// Pointer to the bottom right corner PML mesh
 Mesh* pml_bottom_right_mesh_pt;

 /// Pointer to the bottom left corner PML mesh
 Mesh* pml_bottom_left_mesh_pt;
 
 /// \short Pointer to the mesh containing 
 /// the Helmholtz inner boundary condition elements 
 // Mesh* Helmholtz_inner_boundary_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
ScatteringProblem<ELEMENT>::
ScatteringProblem()
{ 

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Setup "bulk" mesh
  
 // Inner radius
 double a=0.2;
 
 // Thickness of annular computational domain
 double h=0.5; 

 // Set outer radius
 GlobalParameters::Outer_radius=a+h;
 
 // Create circles representing inner and outer boundary
 double x_c=0.0;
 double y_c=0.0;
 Circle* inner_circle_pt=new Circle(x_c,y_c,a);

 // Outer boundary
 //---------------
 TriangleMeshClosedCurve* outer_boundary_pt=0;
 
 unsigned n_segments = 20;
 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(4);
 
 /// \short Each polyline only has three vertices, provide storage for their
 /// coordinates
 Vector<Vector<double> > vertex_coord(2);
 // for(unsigned i=0;i<3;i++)
 for(unsigned i=0;i<2;i++)
  {
   vertex_coord[i].resize(2);
  }

 // First polyline
 vertex_coord[0][0]=-2.0;
 vertex_coord[0][1]=-2.0;
 vertex_coord[1][0]=-2.0;
 vertex_coord[1][1]=2.0;
 
 // Build the 1st boundary polyline
 unsigned boundary_id=2;
 outer_boundary_line_pt[0] = new TriangleMeshPolyLine(vertex_coord,
                                                      boundary_id);
 
 // Second boundary polyline
 vertex_coord[0][0]=-2.0;
 vertex_coord[0][1]=2.0;
 vertex_coord[1][0]=2.0;
 vertex_coord[1][1]=2.0;
 
 // Build the 2nd boundary polyline
 boundary_id=3;
 outer_boundary_line_pt[1] = new TriangleMeshPolyLine(vertex_coord,
                                                      boundary_id);
 
 // Third boundary polyline
 vertex_coord[0][0]=2.0;
 vertex_coord[0][1]=2.0;
 vertex_coord[1][0]=2.0;
 vertex_coord[1][1]=-2.0;
 
 // Build the 3rd boundary polyline
 boundary_id=4;
 outer_boundary_line_pt[2] = new TriangleMeshPolyLine(vertex_coord,
                                                      boundary_id);
 
 // Fourth boundary polyline
 vertex_coord[0][0]=2.0;
 vertex_coord[0][1]=-2.0;
 vertex_coord[1][0]=-2.0;
 vertex_coord[1][1]=-2.0;
 
 // Build the 4th boundary polyline
 boundary_id=5;
 outer_boundary_line_pt[3] = new TriangleMeshPolyLine(vertex_coord,
                                                      boundary_id);
 
 // Create the triangle mesh polygon for outer boundary
 outer_boundary_pt = new TriangleMeshPolygon(outer_boundary_line_pt);

 outer_boundary_pt->output(cout);
 
 // Inner circular boundary
 //------------------------
 Vector<TriangleMeshCurveSection*> inner_boundary_line_pt(2);
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start = 0.0;
 double s_end   = MathematicalConstants::Pi;
 boundary_id = 0;
 inner_boundary_line_pt[0]=
  new TriangleMeshCurviLine(inner_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = MathematicalConstants::Pi;
 s_end   = 2.0*MathematicalConstants::Pi;
 boundary_id = 1;
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
 hole_coords[0]=0.0;
 hole_coords[1]=0.0;
 hole_pt[0]=new TriangleMeshClosedCurve(inner_boundary_line_pt,
                                        hole_coords);
 
 

 /// \short Use the TriangleMeshParameters object for helping on the manage 
 /// of the TriangleMesh parameters. The only parameter that needs to take 
 /// is the outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;

#ifdef ADAPTIVE
 
 // Build "bulk" mesh
 bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Create/set error estimator
 bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Choose error tolerances to force some uniform refinement
 bulk_mesh_pt->min_permitted_error()=0.00004;
 bulk_mesh_pt->max_permitted_error()=0.0001;

#else

 // Refine mesh here, 0.05 is a good default for the given application
 triangle_mesh_parameters.element_area() = 0.05;
 
 // Build "bulk" mesh
 bulk_mesh_pt=new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

#endif
 
 // Create the main triangular mesh
 add_sub_mesh(bulk_mesh_pt);

 // Bulk mesh left boundary id
 unsigned int left_boundary_id = 2;

 // Bulk mesh top boundary id
 unsigned int top_boundary_id = 3;

 // Bulk mesh right boundary id
 unsigned int right_boundary_id = 4;
 
 // Bulk mesh bottom boundary id
 unsigned int bottom_boundary_id = 5;

 // PML width in elements for the right layer
 unsigned n_x_right_pml = 3;

 // PML width in elements for the top layer
 unsigned n_y_top_pml = 3;

 // PML width in elements for the left layer
 unsigned n_x_left_pml = 3; 

 // PML width in elements for the left layer
 unsigned n_y_bottom_pml = 3; 

 /// \short Outer physical length of the PML layers
 /// defaults to 0.2, so 10% of the size of the 
 /// physical domain
 double width_x_right_pml  = 0.2;
 double width_y_top_pml    = 0.2;
 double width_x_left_pml   = 0.2;
 double width_y_bottom_pml = 0.2;

 // Start building PML meshes
 pml_right_mesh_pt = PMLHelper::create_right_pml_mesh(bulk_mesh_pt,
                                                      right_boundary_id, 
                                                      n_x_right_pml, 
                                                      width_x_right_pml);
 pml_top_mesh_pt   = PMLHelper::create_top_pml_mesh(bulk_mesh_pt,
                                                    top_boundary_id, 
                                                    n_y_top_pml, 
                                                    width_y_top_pml);
 pml_left_mesh_pt  = PMLHelper::create_left_pml_mesh(bulk_mesh_pt,
                                                     left_boundary_id, 
                                                     n_x_left_pml, 
                                                     width_x_left_pml);
 pml_bottom_mesh_pt= PMLHelper::create_bottom_pml_mesh(bulk_mesh_pt,
                                                       bottom_boundary_id, 
                                                       n_y_bottom_pml, 
                                                       width_y_bottom_pml);
 
 // Add new submeshes to the global mesh
 add_sub_mesh(pml_right_mesh_pt);
 add_sub_mesh(pml_top_mesh_pt);
 add_sub_mesh(pml_left_mesh_pt);
 add_sub_mesh(pml_bottom_mesh_pt);

 
 // Build corner PML meshes
 pml_top_right_mesh_pt    = 
  PMLHelper::create_top_right_pml_mesh(pml_right_mesh_pt, 
                                       pml_top_mesh_pt, 
                                       bulk_mesh_pt,
                                       right_boundary_id);
 
 pml_bottom_right_mesh_pt = 
  PMLHelper::create_bottom_right_pml_mesh(pml_right_mesh_pt, 
                                          pml_bottom_mesh_pt, 
                                          bulk_mesh_pt,
                                          right_boundary_id);
 
 pml_top_left_mesh_pt     = 
  PMLHelper::create_top_left_pml_mesh(pml_left_mesh_pt, 
                                      pml_top_mesh_pt, 
                                      bulk_mesh_pt,
                                      left_boundary_id);
 
 pml_bottom_left_mesh_pt  = 
  PMLHelper::create_bottom_left_pml_mesh(pml_left_mesh_pt, 
                                         pml_bottom_mesh_pt, 
                                         bulk_mesh_pt,
                                         left_boundary_id);
 
 // Add new submeshes to the global mesh
 add_sub_mesh(pml_top_right_mesh_pt);
 add_sub_mesh(pml_bottom_right_mesh_pt);
 add_sub_mesh(pml_top_left_mesh_pt);
 add_sub_mesh(pml_bottom_left_mesh_pt);
 
 // End building corner PML meshes
 
 // Build the entire mesh from its submeshes
 build_global_mesh();
 
 /// \short Output mesh to visualize the correct construction 
 /// and alignment of the submeshes 
 mesh_pt()->output("mesh.dat");
 
 // Complete the build of all elements so they are fully functional
 
 /// \short Loop over the entire mesh elements to set up element-specific 
 /// things that cannot be handled by constructor: 
 
 unsigned n_element = this->mesh_pt()->nelement();
 
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   
   GeneralisedHelmholtzEquations<2> *el_pt = 
    dynamic_cast<GeneralisedHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));
   
   //Set the frequency function pointer
   el_pt->omega_pt() = &GlobalParameters::Omega;
    
   //Set the soundspeed function pointer
   el_pt->c_fct_pt() =  &GlobalParameters::c_function;
   
   //Set the absorption function pointer
   el_pt->alpha_fct_pt() = &GlobalParameters::alpha_function;
  }
 
 // Set pointer to prescribed flux function for flux elements
 set_prescribed_incoming_flux_pt();
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end of constructor

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::actions_before_adapt()
{
 /// \short Before adapting the added PML meshes must be removed
 /// as they are not refineable and are to be rebuilt from the
 /// newly refined triangular mesh

 // Pointer to the right PML mesh
 delete pml_right_mesh_pt;

 // Pointer to the top PML mesh
 delete pml_top_mesh_pt;

 // Pointer to the left PML mesh
 delete pml_left_mesh_pt;

 // Pointer to the bottom PML mesh
 delete pml_bottom_mesh_pt;

 // Pointer to the top right corner PML mesh
 delete pml_top_right_mesh_pt;

 // Pointer to the top left corner PML mesh
 delete pml_top_left_mesh_pt;

 // Pointer to the bottom right corner PML mesh
 delete pml_bottom_right_mesh_pt;

 // Pointer to the bottom left corner PML mesh
 delete pml_bottom_left_mesh_pt;
 
 // Kill the flux elements and wipe the boundary meshes
// if (!GlobalParameters::Rectangular_outer_boundary) 
//   {
//    delete_face_elements(Helmholtz_outer_boundary_mesh_pt);
//   }
//  delete_face_elements(Helmholtz_inner_boundary_mesh_pt);


 /// \short Rebuild the Problem's global mesh from its various sub-meshes
 /// but first flush all its submeshes
 flush_sub_meshes();
 
 // Then add the triangular mesh back
 add_sub_mesh(bulk_mesh_pt);

 /// \short Rebuild the global mesh such that it now stores
 /// the triangular mesh only
 rebuild_global_mesh();

}// end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::actions_after_adapt()
{

 // Bulk mesh left boundary id
 unsigned int left_boundary_id = 2;

 // Bulk mesh top boundary id
 unsigned int top_boundary_id = 3;

 // Bulk mesh right boundary id
 unsigned int right_boundary_id = 4;
 
 // Bulk mesh bottom boundary id
 unsigned int bottom_boundary_id = 5;

 // PML width in elements for the right layer
 unsigned n_x_right_pml = 3;

 // PML width in elements for the top layer
 unsigned n_y_top_pml = 3;

 // PML width in elements for the left layer
 unsigned n_x_left_pml = 3; 

 // PML width in elements for the left layer
 unsigned n_y_bottom_pml = 3; 

 /// \short Outer physical length of the PML layers
 /// defaults to 0.2, so 10% of the size of the 
 /// physical domain
 double width_x_right_pml  = 0.2;
 double width_y_top_pml    = 0.2;
 double width_x_left_pml   = 0.2;
 double width_y_bottom_pml = 0.2;
 
 // Rebuild the PML meshes based on the new adapted triangular mesh
 pml_right_mesh_pt = PMLHelper::create_right_pml_mesh(bulk_mesh_pt,
                                                      right_boundary_id, 
                                                      n_x_right_pml, 
                                                      width_x_right_pml);
 pml_top_mesh_pt   = PMLHelper::create_top_pml_mesh(bulk_mesh_pt,
                                                    top_boundary_id, 
                                                    n_y_top_pml, 
                                                    width_y_top_pml);
 pml_left_mesh_pt  = PMLHelper::create_left_pml_mesh(bulk_mesh_pt,
                                                     left_boundary_id, 
                                                     n_x_left_pml, 
                                                     width_x_left_pml);
 pml_bottom_mesh_pt= PMLHelper::create_bottom_pml_mesh(bulk_mesh_pt,
                                                       bottom_boundary_id, 
                                                       n_y_bottom_pml, 
                                                       width_y_bottom_pml);
 
 // Add submeshes to the global mesh
 add_sub_mesh(pml_right_mesh_pt);
 add_sub_mesh(pml_top_mesh_pt);
 add_sub_mesh(pml_left_mesh_pt);
 add_sub_mesh(pml_bottom_mesh_pt);
 
 // Rebuild corner PML meshes
 pml_top_right_mesh_pt    = 
  PMLHelper::create_top_right_pml_mesh(pml_right_mesh_pt, 
                                       pml_top_mesh_pt, 
                                       bulk_mesh_pt,
                                       right_boundary_id);
 
 pml_bottom_right_mesh_pt = 
  PMLHelper::create_bottom_right_pml_mesh(pml_right_mesh_pt, 
                                          pml_bottom_mesh_pt, 
                                          bulk_mesh_pt,
                                          right_boundary_id);
 
 pml_top_left_mesh_pt     = 
  PMLHelper::create_top_left_pml_mesh(pml_left_mesh_pt, 
                                      pml_top_mesh_pt, 
                                      bulk_mesh_pt,
                                      left_boundary_id);
 
 pml_bottom_left_mesh_pt  = 
  PMLHelper::create_bottom_left_pml_mesh(pml_left_mesh_pt, 
                                         pml_bottom_mesh_pt, 
                                         bulk_mesh_pt,
                                         left_boundary_id);
 
 // Add submeshes to the global mesh
 add_sub_mesh(pml_top_right_mesh_pt);
 add_sub_mesh(pml_bottom_right_mesh_pt);
 add_sub_mesh(pml_top_left_mesh_pt);
 add_sub_mesh(pml_bottom_left_mesh_pt);
 
 // End rebuilding corner PML meshes
 
 // Build the entire mesh from its submeshes
 rebuild_global_mesh();
 
 // Complete the build of all elements so they are fully functional
 
 /// \short Loop over the entire mesh elements to set up element-specific 
 /// things that cannot be handled by constructor
 unsigned n_element = this->mesh_pt()->nelement();
 
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to GeneralisedHelmholtz bulk element
   GeneralisedHelmholtzEquations<2> *el_pt = 
    dynamic_cast<GeneralisedHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

   //Set the frequency function pointer
   el_pt->omega_pt() = &GlobalParameters::Omega;

   //Set the soundspeed function pointer
   el_pt->c_fct_pt() =  &GlobalParameters::c_function;

   //Set the absorption function pointer
   el_pt->alpha_fct_pt() = &GlobalParameters::alpha_function;
  }
 
 // Reset boundary conditions
 set_prescribed_incoming_flux_pt(); 
  
 
}// end of actions_after_adapt


//==================start_of_set_prescribed_incoming_flux_pt==============
/// Set pointer to prescribed incoming-flux function for all
/// elements in the inner boundary
//========================================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::set_prescribed_incoming_flux_pt()
{
 // Loop over the flux elements to pass pointer to prescribed 
 // flux function for the inner boundary 
//  unsigned n_element=Helmholtz_inner_boundary_mesh_pt->nelement();
//  for(unsigned e=0;e<n_element;e++)
//   {
//    // Upcast from GeneralisedElement to Helmholtz flux element
//    HelmholtzFluxElement<ELEMENT> *el_pt = 
//     dynamic_cast<  HelmholtzFluxElement<ELEMENT>*>(
//      Helmholtz_inner_boundary_mesh_pt->element_pt(e));
   
//    // Set the pointer to the prescribed flux function
//    el_pt->flux_fct_pt() =
//     &GlobalParameters::prescribed_incoming_flux;

//   }

 /// \short Boundary conditions are set on the surface of the circle
 /// as a constant nonzero Dirichlet boundary condition 
 unsigned n_bound = bulk_mesh_pt->nboundary();
 
 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = bulk_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     if ((b==0) || (b==1)) {
      bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
      bulk_mesh_pt->boundary_node_pt(b,n)->pin(1); 
     }
    }
  }
 
 for(unsigned b=0;b<n_bound;b++)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = bulk_mesh_pt->nboundary_node(b);
   
   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,n);
     
     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);
     
     // inner circle boundaries
     if ((b==0) || (b==1))
      { 
       nod_pt->set_value(0,0.1);
       nod_pt->set_value(1,0.0);
      }  
    }
  }
 
}// end of set prescribed_incoming_flux pt


//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::doc_solution(DocInfo& 
                                              doc_info) 
{ 

 ofstream some_file,some_file2;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=5; 

 // Total radiated power
 double power=0.0;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 this->mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Output exact solution 
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt->output_fct(some_file,npts,GlobalParameters::get_exact_u); 
 some_file.close();
 
 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt->compute_error(some_file,GlobalParameters::get_exact_u,
                             error,norm); 
 some_file.close();
 
 // Doc L2 error and norm of solution
 oomph_info << "\nNorm of error   : " << sqrt(error) << std::endl; 
 oomph_info << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;
 

 // Write power to trace file
 Trace_file  << power << std::endl;

 // Do animation of Helmholtz solution
 //-----------------------------------
//  unsigned nstep=40;
//  for (unsigned i=0;i<nstep;i++)
//   {
//    sprintf(filename,"%s/helmholtz_animation%i_frame%i.dat",
//            doc_info.directory().c_str(),
//            doc_info.number(),i);
//    some_file.open(filename);
//    sprintf(filename,"%s/exact_helmholtz_animation%i_frame%i.dat",
//            doc_info.directory().c_str(),
//            doc_info.number(),i);
//    some_file2.open(filename);
//    double phi=2.0*MathematicalConstants::Pi*double(i)/double(nstep-1);
//    unsigned nelem=bulk_mesh_pt->nelement();
//    for (unsigned e=0;e<nelem;e++)
//     {
//      ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
//       bulk_mesh_pt->element_pt(e));
//      el_pt->output_real(some_file,phi,npts);    
//      el_pt->output_real_fct(some_file2,phi,npts,
//                             GlobalParameters::get_exact_u); 
//     }
//    some_file.close();
//    some_file2.close();
//   }

} // end of doc

//============start_of_create_flux_elements==================================
/// Create Helmholtz inner Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object  pointed to by helmholtz_inner_boundary_mesh_pt
//============================================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::
create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const & helmholtz_inner_boundary_mesh_pt)
{

 assert(false);

//  // Loop over the bulk elements adjacent to boundary b
//  unsigned n_element = bulk_mesh_pt->nboundary_element(b);
//  for(unsigned e=0;e<n_element;e++)
//   {
//    // Get pointer to the bulk element that is adjacent to boundary b
//    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
//     bulk_mesh_pt->boundary_element_pt(b,e));
   
//    //Find the index of the face of element e along boundary b 
//    int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

//    // Build the corresponding prescribed incoming-flux element
//    GeneralisedHelmholtzFluxElement<ELEMENT>* flux_element_pt = new 
//     GeneralisedHelmholtzFluxElement<ELEMENT>(bulk_elem_pt,face_index);
   
//    //Add the prescribed incoming-flux element to the surface mesh
//    helmholtz_inner_boundary_mesh_pt->add_element_pt(flux_element_pt);
   
//   } //end of loop over bulk elements adjacent to boundary b
 
} // end of create_flux_elements

//============start_of_delete_face_elements================
/// Delete face elements and wipe the boundary mesh
//==========================================================
template<class ELEMENT>
void ScatteringProblem<ELEMENT>::
delete_face_elements(Mesh* const & boundary_mesh_pt)
{
 // Loop over the surface elements
 unsigned n_element = boundary_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete  boundary_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 boundary_mesh_pt->flush_element_and_node_storage();
 
} // end of delete_outer_face_elements



//==========start_of_main=================================================
/// Solve 2D Helmholtz problem for scattering of a planar wave from a 
/// unit disk 
//========================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags(); 
 
 //Set up the problem
 //------------------
 
#ifdef ADAPTIVE
 
 /// \shortSet up the problem with 2D six-node elements from the
 /// TGeneralisedHelmholtzElement family.
 ScatteringProblem<ProjectableGeneralisedHelmholtzElement
  <TGeneralisedHelmholtzElement<2,3> > > problem;

 /// \short Set up the problem with 2D ten-node elements from the
 /// TGeneralisedHelmholtzElement family. 
//  ScatteringProblem<ProjectableGeneralisedHelmholtzElement
//   <TGeneralisedHelmholtzElement<2,4> > > problem;

#else
 
 /// \short Set up the problem with 2D six-node elements from the
 /// TGeneralisedHelmholtzElement family. 
 ScatteringProblem<TGeneralisedHelmholtzElement<2,3> >  problem;

 /// \short Set up the problem with 2D ten-node elements from the
 /// TGeneralisedHelmholtzElement family. 
//   ScatteringProblem<TGeneralisedHelmholtzElement<2,4> >  problem;

#endif
 
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");


#ifdef ADAPTIVE
 
 // Max. number of adaptations
 unsigned max_adapt=3;
 
 /// \short Solve the problem with Newton's method, allowing
 /// up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);
 
#else
 
 // Solve the problem with Newton's method
 problem.newton_solve();
 
#endif
 
 //Output solution
 problem.doc_solution(doc_info);
    
} //end of main

