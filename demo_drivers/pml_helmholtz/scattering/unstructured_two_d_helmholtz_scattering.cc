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
// Driver for a specific 2D Helmholtz problem with
// perfectly matched layer treatment for the exterior boundaries

#include<fenv.h>


#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The Helmholtz equations
#include "helmholtz.h"

// The pml Helmholtz equations
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{

 /// Wavenumber (also known as k), k=omega/c
 double Wavenumber = sqrt(10.0);

 /// Square of the wavenumber (also known as k^2)
 double K_squared = Wavenumber * Wavenumber;

 /// Number of terms used in the computation
 /// of the exact solution
 unsigned N_fourier=100;

 /// Imaginary unit
 std::complex<double> I(0.0,1.0);

 /// Exact solution for scattered field
 /// (vector returns real and impaginary parts).
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  // Switch to polar coordinates
  double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta;
  theta=atan2(x[1],x[0]);

  // Argument for Bessel/Hankel functions
  double rr=Wavenumber*r;

  // Evaluate Bessel/Hankel functions
  complex <double > u_ex(0.0,0.0);
  Vector<double> jn(N_fourier+1), yn(N_fourier+1),
   jnp(N_fourier+1), ynp(N_fourier+1);
  Vector<double> jn_a(N_fourier+1),yn_a(N_fourier+1),
   jnp_a(N_fourier+1), ynp_a(N_fourier+1);
  Vector<complex<double> > h(N_fourier+1),h_a(N_fourier+1),
   hp(N_fourier+1), hp_a(N_fourier+1);

  // We want to compute N_fourier terms but the function
  // may return fewer than that.
  int n_actual=0;
  CRBond_Bessel::bessjyna(N_fourier,Wavenumber,n_actual,
                          &jn_a[0],&yn_a[0],
                          &jnp_a[0],&ynp_a[0]);

  // Shout if things went wrong
#ifdef PARANOID
  if (n_actual!=int(N_fourier))
   {
    std::ostringstream error_stream;
    error_stream << "CRBond_Bessel::bessjyna() only computed "
                 << n_actual << " rather than " << N_fourier
                 << " Bessel functions.\n";
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Evaluate Hankel at actual radius
  Hankel_functions_for_helmholtz_problem::Hankel_first(N_fourier,rr,h,hp);

  // Evaluate Hankel at inner (unit) radius
  Hankel_functions_for_helmholtz_problem::Hankel_first(N_fourier
                                                       ,Wavenumber,
                                                       h_a,hp_a);

  // Compute the sum: Separate the computation of the negative
  // and positive terms
  // ALH: The construction with the static cast to a double is to avoid
  // a floating point exception when running unoptimised on my machine
  for (unsigned i=0;i<N_fourier;i++)
   {
    u_ex-=pow(I,i)*h[i]*jnp_a[i]*pow(exp(I*theta),static_cast<double>(i))/hp_a[i];
   }
  for (unsigned i=1;i<N_fourier;i++)
   {
    u_ex-=pow(I,i)*h[i]*jnp_a[i]*pow(exp(-I*theta),static_cast<double>(i))/hp_a[i];
   }

  // Get the real & imaginary part of the result
  u[0]=real(u_ex);
  u[1]=imag(u_ex);

 }// end of get_exact_u



 /// Flux (normal derivative) on the unit disk
 /// for a planar incoming wave
 void prescribed_incoming_flux(const Vector<double>& x,
                               complex<double>& flux)
 {
  // Switch to polar coordinates
  double r;
  r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double theta;
  theta=atan2(x[1],x[0]);

  // Argument of the Bessel/Hankel fcts
  double rr=Wavenumber*r;

  // Compute Bessel/Hankel functions
  Vector<double> jn(N_fourier+1), yn(N_fourier+1),
   jnp(N_fourier+1), ynp(N_fourier+1);

  // We want to compute N_fourier terms but the function
  // may return fewer than that.
  int n_actual=0;
  CRBond_Bessel::bessjyna(N_fourier,rr,n_actual,&jn[0],&yn[0],
                          &jnp[0],&ynp[0]);

  // Shout if things went wrong...
#ifdef PARANOID
  if (n_actual!=int(N_fourier))
   {
    std::ostringstream error_stream;
    error_stream << "CRBond_Bessel::bessjyna() only computed "
                 << n_actual << " rather than " << N_fourier
                 << " Bessel functions.\n";
    throw OomphLibError(error_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif

  // Compute the sum: Separate the computation of the negative and
  // positive terms
  flux=std::complex<double>(0.0,0.0);
  for (unsigned i=0;i<N_fourier;i++)
   {
    flux+=pow(I,i)*(Wavenumber)*pow(exp(I*theta),i)*jnp[i];
   }
  for (unsigned i=1;i<N_fourier;i++)
   {
    flux+=pow(I,i)*(Wavenumber)*pow(exp(-I*theta),i)*jnp[i];
   }

 }// end of prescribed_incoming_flux

 class TestPMLMapping : public PMLMapping
 {

   public:

   /// Default constructor (empty)
   TestPMLMapping(){};

   /// Overwrite the pure PML mapping coefficient function to return the
   /// coeffcients proposed by Bermudez et al
   std::complex<double> gamma(const double& nu_i,
                              const double& pml_width_i,
                              const double& k_squared_local,
                              const double& alpha_shift=0.0)
   {
      // (return) gamma = 1 + (1/k) * (i/|outer_boundary - x|)
      /*return 1.0 + (1.0 / std::complex<double> (sqrt(k_squared_local), 0))
       * std::complex<double>
        (0.0, 1.0/(std::fabs(pml_width_i - nu_i)));*/
     return 1.0 + (1.0 / std::complex<double> (sqrt(k_squared_local), 0))
          * std::complex<double>
           (0.0, 2.0/(std::fabs(pml_width_i - nu_i)));
   }

 };

TestPMLMapping* Test_pml_mapping_pt = new TestPMLMapping;

} // end of namespace


/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//========= start_of_problem_class=====================================
/// Problem class to demonstrate use of perfectly matched layers
/// for Helmholtz problems.
//=====================================================================
template<class ELEMENT>
class PMLProblem : public Problem
{

public:

 /// Constructor
 PMLProblem();

 /// Destructor (empty)
 ~PMLProblem(){};

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){} ;

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){};

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// Create PML meshes
 void create_pml_meshes();

 /// Create Helmholtz flux elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified surface Mesh
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const & helmholtz_inner_boundary_mesh_pt);

 /// Create Helmholtz power elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified surface Mesh
 void create_power_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                            Mesh* const & helmholtz_power_boundary_mesh_pt);

#ifdef ADAPTIVE

 /// Actions before adapt: Wipe the PML meshes
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the PML meshes
 void actions_after_adapt();



#endif


#ifdef ADAPTIVE

private:

 /// Pointer to the refineable "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#else

private:

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif


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

 /// Pointer to the mesh containing
 /// the Helmholtz inner boundary condition elements
 Mesh* Helmholtz_inner_boundary_mesh_pt;

 /// Pointer to mesh of elements that compute the radiated power
 Mesh* Helmholtz_power_boundary_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT>
PMLProblem<ELEMENT>::PMLProblem()
{

 // Open trace file
 Trace_file.open("RESLT/trace.dat");

 // Create circle representing inner boundary
 double a=1.0;
 double x_c=0.0;
 double y_c=0.0;
 Circle* inner_circle_pt=new Circle(x_c,y_c,a);

 // Outer boundary
 //---------------
 TriangleMeshClosedCurve* outer_boundary_pt=0;

 unsigned n_segments = 20;
 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(4);

 // Each polyline only has three vertices, provide storage for their
 // coordinates
 Vector<Vector<double> > vertex_coord(2);
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


 // Use the TriangleMeshParameters object for helping on the manage
 // of the TriangleMesh parameters. The only parameter that needs to take
 // is the outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;

 // Target element size in bulk mesh
 triangle_mesh_parameters.element_area() = 0.05;

#ifdef ADAPTIVE

 // Build adaptive "bulk" mesh
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=0.00004;
 Bulk_mesh_pt->max_permitted_error()=0.0001;

#else

 // Build "bulk" mesh
 Bulk_mesh_pt=new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

#endif

 // Pointer to mesh containing the Helmholtz inner boundary condition
 // elements. Specify outer radius
 Helmholtz_inner_boundary_mesh_pt = new Mesh;

 // Create prescribed-flux elements from all elements that are
 // adjacent to the inner boundary , but add them to a separate mesh.
 create_flux_elements(0,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 create_flux_elements(1,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);

 // Pointer to mesh containing the Helmholtz power condition
 // elements.
 Helmholtz_power_boundary_mesh_pt = new Mesh;

 // Create power elements from all elements that are
 // adjacent to the inner boundary , but add them to a separate mesh.
 create_power_elements(0,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);
 create_power_elements(1,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);

 // Create the main triangular mesh
 add_sub_mesh(Bulk_mesh_pt);

 // Create PML meshes and add them to the global mesh
 create_pml_meshes();

 // Add the flux mesh
 add_sub_mesh(Helmholtz_inner_boundary_mesh_pt);

 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional
 unsigned n_element = this->mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Pml Helmholtz bulk element
   PMLHelmholtzEquations<2> *el_pt =
    dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

   if (el_pt!=0)
    {
     //Set the k_squared function pointer
     el_pt->k_squared_pt() = &GlobalParameters::K_squared;
     // Pass in a pointer to the class containing the PML mapping function
     //el_pt->pml_mapping_pt() = GlobalParameters::Test_pml_mapping_pt;
    }

  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;

} // end of constructor


#ifdef ADAPTIVE

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::actions_before_adapt()
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

 // Rebuild the Problem's global mesh from its various sub-meshes
 // but first flush all its submeshes
 flush_sub_meshes();

 // Then add the triangular mesh back
 add_sub_mesh(Bulk_mesh_pt);

 //  Rebuild the global mesh such that it now stores
 // the triangular mesh only
 rebuild_global_mesh();

}// end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::actions_after_adapt()
{
 create_inner_bc_elements(0,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);
 create_inner_bc_elements(1,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);
 // Build PML meshes  and add them to the global mesh
 create_pml_meshes();

 // Complete the build of all elements so they are fully functional

 // Loop over the entire mesh elements to set up element-specific
 // things that cannot be handled by constructor
 unsigned n_element = this->mesh_pt()->nelement();

 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to PMLHelmholtz bulk element
   PMLHelmholtzEquations<2> *el_pt =
    dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

   if (el_pt!=0)
    {
     //Set the k_squared double pointer
     el_pt->k_squared_pt() = &GlobalParameters::K_squared;
    }
  }

 // Create prescribed-flux elements and BC elements
 // from all elements that are adjacent to the boundaries and add them to
 // Helmholtz_boundary_meshes
 create_flux_elements(0,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 create_flux_elements(1,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
er_eleme
 // Create power elements from all elements that are
 // adjacent to the inner boundary , but add them to a separate mesh.
 create_power_elements(0,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);
 create_power_elements(1,Bulk_mesh_pt,Helmholtz_power_boundary_mesh_pt);

 // Add the flux mesh
 add_sub_mesh(Helmholtz_inner_boundary_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_after_adapt

#endif

//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

 ofstream some_file,some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 // Accumulate contribution from elements
 double power=0.0;
 unsigned nn_element=Helmholtz_power_boundary_mesh_pt->nelement();
 for(unsigned e=0;e<nn_element;e++)
  {
   PMLHelmholtzPowerElement<ELEMENT> *el_pt =
    dynamic_cast<PMLHelmholtzPowerElement<ELEMENT>*>(
     Helmholtz_power_boundary_mesh_pt->element_pt(e));
   power += el_pt->global_power_contribution(some_file);
  }
 some_file.close();
 oomph_info << "Total radiated power: " << power << std::endl;

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
 Bulk_mesh_pt->output_fct(some_file,npts,GlobalParameters::get_exact_u);
 some_file.close();

 double error,norm;
 sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,GlobalParameters::get_exact_u,
                             error,norm);
 some_file.close();

 // Doc L2 error and norm of solution
 oomph_info << "\nNorm of error   : " << sqrt(error) << std::endl;
 oomph_info << "Norm of solution: " << sqrt(norm) << std::endl << std::endl;

 // Output PML layers
 //-----------------
 sprintf(filename,"%s/pml_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 PML_right_mesh_pt->output(some_file,npts);
 PML_left_mesh_pt->output(some_file,npts);
 PML_top_mesh_pt->output(some_file,npts);
 PML_bottom_mesh_pt->output(some_file,npts);
 PML_bottom_right_mesh_pt->output(some_file,npts);
 PML_top_right_mesh_pt->output(some_file,npts);
 PML_bottom_left_mesh_pt->output(some_file,npts);
 PML_top_left_mesh_pt->output(some_file,npts);
 some_file.close();



 // Write norm of solution to trace file
 double norm2=0.0;
 Bulk_mesh_pt->compute_norm(norm2);
 Trace_file  << norm2 << std::endl;

 // //Do animation of Helmholtz solution
 // //-----------------------------------
 // unsigned nstep=40;
 // for (unsigned i=0;i<nstep;i++)
 //  {
 //   sprintf(filename,"%s/helmholtz_animation%i_frame%i.dat",
 //           doc_info.directory().c_str(),
 //           doc_info.number(),i);
 //   some_file.open(filename);
 //   double phi=2.0*MathematicalConstants::Pi*double(i)/double(nstep-1);
 //   unsigned nelem=Bulk_mesh_pt->nelement();
 //   for (unsigned e=0;e<nelem;e++)
 //    {
 //     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
 //      Bulk_mesh_pt->element_pt(e));
 //     el_pt->output_real(some_file,phi,npts);
 //    }
 //   some_file.close();
 //  }

} // end of doc

//============start_of_create_flux_elements==================================
/// Create Helmholtz inner Flux Elements on the b-th boundary of
/// the Mesh object pointed to by bulk_mesh_pt and add the elements
/// to the Mesh object  pointed to by helmholtz_inner_boundary_mesh_pt
//============================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::
create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const & helmholtz_inner_boundary_mesh_pt)
{
 // Loop over the bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //Find the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed incoming-flux element
   PMLHelmholtzFluxElement<ELEMENT>* flux_element_pt = new
    PMLHelmholtzFluxElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed incoming-flux element to the surface mesh
   helmholtz_inner_boundary_mesh_pt->add_element_pt(flux_element_pt);

   // Set the pointer to the prescribed flux function
   flux_element_pt->flux_fct_pt() =
    &GlobalParameters::prescribed_incoming_flux;

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_flux_elements

//============start_of_create_power_elements==================================
/// Create Helmholtz inner Flux Elements on the b-th boundary of
/// the Mesh object pointed to by bulk_mesh_pt and add the elements
/// to the Mesh object  pointed to by helmholtz_power_boundary_mesh_pt
//============================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::
create_power_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                     Mesh* const & helmholtz_power_boundary_mesh_pt)
{
 // Loop over the bulk elements adjacent to boundary b
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //Find the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding prescribed power element
   PMLHelmholtzPowerElement<ELEMENT>* power_element_pt = new
    PMLHelmholtzPowerElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed power element to the surface mesh
   helmholtz_power_boundary_mesh_pt->add_element_pt(power_element_pt);

  } //end of loop over bulk elements adjacent to boundary b

} // end of create_power_elements


//============start_of_create_pml_meshes======================================
/// Create PML meshes and add them to the problem's sub-meshes
//============================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::create_pml_meshes()
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

 // PML width in elements for the bottom layer
 unsigned n_y_bottom_pml = 3;

 // Outer physical length of the PML layers
 // defaults to 0.2, so 10% of the size of the
 // physical domain
 double width_x_right_pml  = 0.2;
 double width_y_top_pml    = 0.2;
 double width_x_left_pml   = 0.2;
 double width_y_bottom_pml = 0.2;

 // Build the PML meshes based on the new adapted triangular mesh
 PML_right_mesh_pt =
  TwoDimensionalPMLHelper::create_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt,right_boundary_id,
   n_x_right_pml, width_x_right_pml);
 PML_top_mesh_pt   =
  TwoDimensionalPMLHelper::create_top_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (Bulk_mesh_pt, top_boundary_id,
   n_y_top_pml, width_y_top_pml);
 PML_left_mesh_pt  =
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
 PML_top_right_mesh_pt    =
  TwoDimensionalPMLHelper::create_top_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_right_mesh_pt, PML_top_mesh_pt,
   Bulk_mesh_pt, right_boundary_id);

 PML_bottom_right_mesh_pt =
  TwoDimensionalPMLHelper::create_bottom_right_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_right_mesh_pt, PML_bottom_mesh_pt,
   Bulk_mesh_pt, right_boundary_id);

 PML_top_left_mesh_pt     =
  TwoDimensionalPMLHelper::create_top_left_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_left_mesh_pt, PML_top_mesh_pt,
   Bulk_mesh_pt, left_boundary_id);

 PML_bottom_left_mesh_pt  =
  TwoDimensionalPMLHelper::create_bottom_left_pml_mesh
  <PMLLayerElement<ELEMENT> >
  (PML_left_mesh_pt, PML_bottom_mesh_pt,
   Bulk_mesh_pt, left_boundary_id);

 // Add submeshes to the global mesh
 add_sub_mesh(PML_top_right_mesh_pt);
 add_sub_mesh(PML_bottom_right_mesh_pt);
 add_sub_mesh(PML_top_left_mesh_pt);
 add_sub_mesh(PML_bottom_left_mesh_pt);

} // end of create_pml_meshes



//==========start_of_main=================================================
/// Solve 2D Helmholtz problem
//========================================================================
int main(int argc, char **argv)
{
 //Set up the problem
 //------------------

#ifdef ADAPTIVE

 // Set up the problem with projectable 2D six-node elements from the
 // TPMLHelmholtzElement family.
 PMLProblem<ProjectablePMLHelmholtzElement
  <TPMLHelmholtzElement<2,3> > > problem;

#else

 // Set up the problem with 2D six-node elements from the
 // TPMLHelmholtzElement family.
 PMLProblem<TPMLHelmholtzElement<2,3> >  problem;
#endif

 // Create label for output
 //------------------------
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");


#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=1;

 // Solve the problem with the adaptive Newton's method, allowing
 // up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);

#else

 // Solve the problem with Newton's method
 problem.newton_solve();

#endif

 //Output solution
 problem.doc_solution(doc_info);

} //end of main
