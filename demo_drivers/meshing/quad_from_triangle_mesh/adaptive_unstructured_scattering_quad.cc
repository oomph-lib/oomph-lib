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
//Driver for a specific  2D Helmholtz problem with flux boundary conditions
//uses two separate meshes for bulk and surface mesh

#include "math.h"
#include <complex>

//Generic routines
#include "generic.h"

// The Helmholtz equations
#include "helmholtz.h"

// The mesh
#include "meshes/quad_from_triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

// Get the Bessel functions
#include "oomph_crbond_bessel.h"

using namespace oomph;
using namespace std;

#ifndef ADAPTIVE
#define ADAPTIVE

/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////

//===== start_of_namespace=============================================
/// Namespace for the Helmholtz problem parameters
//=====================================================================
namespace GlobalParameters
{

 /// Square of the wavenumber
 double K_squared=10.0; 
 
 /// Number of terms used in the computation 
 /// of the exact solution
 unsigned N_fourier=10;
 
 /// Flag to choose the Dirichlet to Neumann BC
 /// or ABC BC
 bool DtN_BC=false;

 /// Flag to choose wich order to use
 // in the ABCs BC: 1 for ABC 1st order...
 unsigned ABC_order=3;

 /// Radius of outer boundary (must be a circle!)
 double Outer_radius=1.5;

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
  double rr=sqrt(K_squared)*r;  
 
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
  CRBond_Bessel::bessjyna(N_fourier,sqrt(K_squared),n_actual,
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
                                                       ,sqrt(K_squared),
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
  double rr=sqrt(K_squared)*r;  
  
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
    flux+=pow(I,i)*(sqrt(K_squared))*pow(exp(I*theta),i)*jnp[i];
   }
  for (unsigned i=1;i<N_fourier;i++)
   {
    flux+=pow(I,i)*(sqrt(K_squared))*pow(exp(-I*theta),i)*jnp[i];
   }


 }// end of prescribed_incoming_flux 

} // end of namespace




/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////



//========= start_of_problem_class=====================================
/// Problem class to compute scattering of planar wave from unit disk
//=====================================================================
template<class ELEMENT,class MESH> 
class ScatteringProblem : public Problem
{

public:
 
 /// Constructor
 ScatteringProblem();
 
 /// Destructor (empty)
 ~ScatteringProblem(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){} 

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
 
 /// Recompute gamma integral before checking Newton residuals
 void actions_before_newton_convergence_check()
  {
   if (GlobalParameters::DtN_BC)
    {
     Helmholtz_outer_boundary_mesh_pt->setup_gamma();
    }
  }

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();
 
 /// Actions after adapt: Rebuild the mesh of prescribed flux elements
 void actions_after_adapt();
 
 /// Create BC elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified survace Mesh 
 void create_outer_bc_elements(
  const unsigned &b, Mesh* const &bulk_mesh_pt,
  Mesh* const & helmholtz_outer_boundary_mesh_pt);
 
 /// Create Helmholtz flux elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the specified surface Mesh 
 void create_flux_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                           Mesh* const & helmholtz_inner_boundary_mesh_pt);
 
 /// Delete boundary face elements and wipe the surface mesh
 void delete_face_elements( Mesh* const & boundary_mesh_pt);
 
 /// Set pointer to prescribed-flux function for all
 /// elements in the surface mesh on the surface of the unit disk
 void set_prescribed_incoming_flux_pt();

 /// Set up boundary condition elements on outer boundary
 void setup_outer_boundary();

 /// Pointer to the "bulk" mesh
 MESH* Bulk_mesh_pt;

 /// Pointer to mesh containing the DtN (or ABC) boundary
 /// condition elements
 HelmholtzDtNMesh<ELEMENT>* Helmholtz_outer_boundary_mesh_pt;
 
 /// Pointer to the mesh containing 
 /// the Helmholtz inner boundary condition elements 
 Mesh* Helmholtz_inner_boundary_mesh_pt;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//=======start_of_constructor=============================================
/// Constructor for Helmholtz problem
//========================================================================
template<class ELEMENT,class MESH>
ScatteringProblem<ELEMENT,MESH>::
ScatteringProblem()
{ 

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Setup "bulk" mesh
  
 // Inner radius
 double a=1.0;
 
 // Thickness of annular computational domain
 double h=0.5; 

 // Set outer radius
 GlobalParameters::Outer_radius=a+h;
 
 // Create circles representing inner and outer boundary
 double x_c=0.0;
 double y_c=0.0;
 Circle* inner_circle_pt=new Circle(x_c,y_c,a);
 Circle* outer_circle_pt=new Circle(x_c,y_c,GlobalParameters::Outer_radius);

 // Outer boundary
 //---------------
 TriangleMeshClosedCurve* outer_boundary_pt=0;

 unsigned n_segments=40;
 Vector<TriangleMeshCurveSection*> outer_boundary_line_pt(2);
 
 // The intrinsic coordinates for the beginning and end of the curve
 double s_start = 0.0;
 double s_end   = MathematicalConstants::Pi;
 unsigned boundary_id = 2;
 outer_boundary_line_pt[0]=
  new TriangleMeshCurviLine(outer_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);
 
 // The intrinsic coordinates for the beginning and end of the curve
 s_start = MathematicalConstants::Pi;
 s_end   = 2.0*MathematicalConstants::Pi;
 boundary_id = 3;
 outer_boundary_line_pt[1]=
  new TriangleMeshCurviLine(outer_circle_pt,
                            s_start,
                            s_end,
                            n_segments,
                            boundary_id);
 
 // Create closed curve for outer boundary
 outer_boundary_pt=new TriangleMeshClosedCurve(outer_boundary_line_pt);

 // Inner circular boundary
 //------------------------
 Vector<TriangleMeshCurveSection*> inner_boundary_line_pt(2);

 // The intrinsic coordinates for the beginning and end of the curve
 s_start = 0.0;
 s_end   = MathematicalConstants::Pi;
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
 
 

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters. The only parameter that needs to take is the
 // outer boundary.
 TriangleMeshParameters triangle_mesh_parameters(outer_boundary_pt);

 // Specify the closed curve using the TriangleMeshParameters object
 triangle_mesh_parameters.internal_closed_curve_pt() = hole_pt;

 // Build "bulk" mesh
 Bulk_mesh_pt=new MESH(triangle_mesh_parameters);
 
#ifdef ADAPTIVE
 
 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Choose error tolerances to force some uniform refinement
 Bulk_mesh_pt->min_permitted_error()=0.00004;
 Bulk_mesh_pt->max_permitted_error()=0.0001;

#endif
 
 // Pointer to mesh containing the Helmholtz outer boundary condition
 // elements. Specify outer radius and number of Fourier terms to be
 // used in gamma integral
 Helmholtz_outer_boundary_mesh_pt = 
  new HelmholtzDtNMesh<ELEMENT>(a+h,GlobalParameters::N_fourier);
 
 // Create outer boundary elements from all elements that are 
 // adjacent to the outer boundary , but add them to a separate mesh.
 create_outer_bc_elements(2,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);
 create_outer_bc_elements(3,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);

 // Pointer to mesh containing the Helmholtz inner boundary condition
 // elements. Specify outer radius
 Helmholtz_inner_boundary_mesh_pt = new Mesh;
   
 // Create prescribed-flux elements from all elements that are 
 // adjacent to the inner boundary , but add them to a separate mesh.
 create_flux_elements(0,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 create_flux_elements(1,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 
 // Add the several  sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Helmholtz_outer_boundary_mesh_pt); 
 add_sub_mesh(Helmholtz_inner_boundary_mesh_pt);   
  
 // Build the Problem's global mesh from its various sub-meshes
 build_global_mesh();

 // Complete the build of all elements so they are fully functional
 
 // Loop over the Helmholtz bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wave number squared
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the k_squared  pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Set up elements on outer boundary
 setup_outer_boundary();

 // Set pointer to prescribed flux function for flux elements
 set_prescribed_incoming_flux_pt();

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of face elements
//========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::actions_before_adapt()
{ 
 // Kill the flux elements and wipe the boundary meshs
 delete_face_elements(Helmholtz_outer_boundary_mesh_pt);
 delete_face_elements(Helmholtz_inner_boundary_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt


//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the face element meshes
//========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::actions_after_adapt()
{


 // Complete the build of all elements so they are fully functional
 
 // Loop over the Helmholtz bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wave number squared
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the k_squared  pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Create prescribed-flux elements and BC elements 
 // from all elements that are adjacent to the boundaries and add them to 
 // Helmholtz_boundary_meshes
 create_flux_elements(0,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 create_flux_elements(1,Bulk_mesh_pt,Helmholtz_inner_boundary_mesh_pt);
 create_outer_bc_elements(2,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);
 create_outer_bc_elements(3,Bulk_mesh_pt,Helmholtz_outer_boundary_mesh_pt);

 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Set pointer to prescribed flux function and DtN mesh
 setup_outer_boundary();
 set_prescribed_incoming_flux_pt(); 
  
 
}// end of actions_after_adapt


//==================start_of_setup_outer_boundary=========================
/// Set pointers for elements on outer boundary
//========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::setup_outer_boundary()
{ 
 // Loop over the flux elements to pass pointer to DtN
 // BC for the outer boundary
 unsigned n_element=Helmholtz_outer_boundary_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Dirichlet to Neumann BC
   if (GlobalParameters::DtN_BC)
    {
     // Upcast from GeneralisedElement to Helmholtz flux element
     HelmholtzDtNBoundaryElement<ELEMENT> *el_pt = 
      dynamic_cast< HelmholtzDtNBoundaryElement<ELEMENT>*>(
       Helmholtz_outer_boundary_mesh_pt->element_pt(e));
     
     // Set pointer to the mesh that contains all the boundary condition
     // elements on this boundary
     el_pt->set_outer_boundary_mesh_pt(Helmholtz_outer_boundary_mesh_pt);
    }
   // ABCs BC
   else
    {
     // Upcast from GeneralisedElement to appropriate type
     HelmholtzAbsorbingBCElement<ELEMENT> *el_pt = 
      dynamic_cast< HelmholtzAbsorbingBCElement<ELEMENT>*>(
       Helmholtz_outer_boundary_mesh_pt->element_pt(e));
     
     // Set pointer to outer radius of artificial boundary
     el_pt->outer_radius_pt()=&GlobalParameters::Outer_radius;
     
     // Set order of absorbing boundary condition
     el_pt->abc_order_pt()=&GlobalParameters::ABC_order;
    }    
  }
}


//==================start_of_set_prescribed_incoming_flux_pt==============
/// Set pointer to prescribed incoming-flux function for all
/// elements in the inner boundary
//========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::set_prescribed_incoming_flux_pt()
{
 // Loop over the flux elements to pass pointer to prescribed 
 // flux function for the inner boundary 
 unsigned n_element=Helmholtz_inner_boundary_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz flux element
   HelmholtzFluxElement<ELEMENT> *el_pt = 
    dynamic_cast<  HelmholtzFluxElement<ELEMENT>*>(
     Helmholtz_inner_boundary_mesh_pt->element_pt(e));
   
   // Set the pointer to the prescribed flux function
   el_pt->flux_fct_pt() =
    &GlobalParameters::prescribed_incoming_flux;
  }

}// end of set prescribed_incoming_flux pt

//=====================start_of_doc=======================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::doc_solution(DocInfo& 
                                              doc_info) 
{ 

 ofstream some_file,some_file2;
 char filename[100];
 
 // Number of plot points
 unsigned npts;
 npts=2; 

 // Total radiated power
 double power=0.0;

 // Compute/output the radiated power
 //----------------------------------
 sprintf(filename,"%s/power%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 
 // Accumulate contribution from elements
 unsigned nn_element=Helmholtz_outer_boundary_mesh_pt->nelement(); 
 for(unsigned e=0;e<nn_element;e++)
  {
   HelmholtzBCElementBase<ELEMENT> *el_pt = 
    dynamic_cast< HelmholtzBCElementBase<ELEMENT>*>(
     Helmholtz_outer_boundary_mesh_pt->element_pt(e)); 
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
 

 // Write power to trace file
 Trace_file << power << std::endl;

 // Do animation of Helmholtz solution
 //-----------------------------------
 // unsigned nstep=40;
 // for (unsigned i=0;i<nstep;i++)
 //  {
 //   sprintf(filename,"%s/helmholtz_animation%i_frame%i.dat",
 //           doc_info.directory().c_str(),
 //           doc_info.number(),i);
 //   some_file.open(filename);
 //   sprintf(filename,"%s/exact_helmholtz_animation%i_frame%i.dat",
 //           doc_info.directory().c_str(),
 //           doc_info.number(),i);
 //   some_file2.open(filename);
 //   double phi=2.0*MathematicalConstants::Pi*double(i)/double(nstep-1);
 //   unsigned nelem=Bulk_mesh_pt->nelement();
 //   for (unsigned e=0;e<nelem;e++)
 //    {
 //     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
 //      Bulk_mesh_pt->element_pt(e));
 //     el_pt->output_real(some_file,phi,npts);    
 //     el_pt->output_real_fct(some_file2,phi,npts,
 //                            GlobalParameters::get_exact_u); 
 //    }
 //   some_file.close();
 //   some_file2.close();
 //  }

} // end of doc

//============start_of_create_flux_elements==================================
/// Create Helmholtz inner Flux Elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object  pointed to by helmholtz_inner_boundary_mesh_pt
//============================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::
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
   HelmholtzFluxElement<ELEMENT>* flux_element_pt = new 
    HelmholtzFluxElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed incoming-flux element to the surface mesh
   helmholtz_inner_boundary_mesh_pt->add_element_pt(flux_element_pt);
   
  } //end of loop over bulk elements adjacent to boundary b
 
} // end of create_flux_elements



//============start_of_create_outer_bc_elements==============================
/// Create outer BC elements on the b-th boundary of 
/// the Mesh object pointed to by bulk_mesh_pt and add the elements 
/// to the Mesh object pointed to by helmholtz_outer_boundary_mesh_pt.
//===========================================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::
create_outer_bc_elements(const unsigned &b, Mesh* const &bulk_mesh_pt,
                         Mesh* const & helmholtz_outer_boundary_mesh_pt)
{
 // Loop over the bulk elements adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b 
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding outer flux element
   
   // Dirichlet to Neumann boundary conditon
   if (GlobalParameters::DtN_BC)
    {
     HelmholtzDtNBoundaryElement<ELEMENT>* flux_element_pt = new 
      HelmholtzDtNBoundaryElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
     helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);
    }
   //  ABCs BC
   else
    {
     HelmholtzAbsorbingBCElement<ELEMENT>* flux_element_pt = new 
      HelmholtzAbsorbingBCElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the flux boundary element to the  helmholtz_outer_boundary_mesh
     helmholtz_outer_boundary_mesh_pt->add_element_pt(flux_element_pt);
    }
  } //end of loop over bulk elements adjacent to boundary b
} // end of create_outer_bc_elements


//============start_of_delete_face_elements================
/// Delete face elements and wipe the boundary mesh
//==========================================================
template<class ELEMENT,class MESH>
void ScatteringProblem<ELEMENT,MESH>::
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

 // Define case to be run
 unsigned i_case=0;
 CommandLineArgs::specify_command_line_flag("--case",&i_case);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Now set flags accordingly
 switch(i_case)
  {
  case 0:
   GlobalParameters::DtN_BC=true;
   break;
   
  case 1:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=1;
   break;

  case 2:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=2;
   break;

  case 3:
   GlobalParameters::DtN_BC=false;
   GlobalParameters::ABC_order=3;
   break;
  }
 
 
 //Set up the problem
 //------------------
 
#ifdef ADAPTIVE
 
 // Typedef the ELEMENT and MESH type
 typedef RefineableQHelmholtzElement<2,3> ELEMENT;
 typedef RefineableQuadFromTriangleMesh<ELEMENT> MESH;
  
#else
 
 // Typedef the ELEMENT and MESH type
 typedef QHelmholtzElement<2,3> ELEMENT;
 typedef QuadFromTriangleMesh<ELEMENT> MESH;
 
#endif
 
 // Set up the problem with 2D nine-node elements from the
 // QHelmholtzElement family.
 ScatteringProblem<ELEMENT,MESH> problem;
 
 // Create label for output
 //------------------------
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");

#ifdef ADAPTIVE

 // Max. number of adaptations
 unsigned max_adapt=1;
 
 // Solve the problem with Newton's method, allowing
 // up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);

#else

 // Solve the problem with Newton's method
 problem.newton_solve();

#endif

 //Output solution
 problem.doc_solution(doc_info);
    
} //end of main


#endif
