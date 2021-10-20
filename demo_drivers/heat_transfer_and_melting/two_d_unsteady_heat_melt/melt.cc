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
//Driver for 2D unsteady heat problem with "melting" constraint
#include <fenv.h> 

//Generic routines
#include "generic.h"

// The unsteady heat equations
#include "unsteady_heat.h"
#include "solid.h"

// Mesh
#include "meshes/triangle_mesh.h"

// The melt equations
#include "heat_transfer_and_melt_elements.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;


/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// 

//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{
 /// Melt-temperature
 double Melt_temperature=0.0;
 
 /// Non-dim density for pseudo-solid
 double Lambda_sq=0.0;

 /// Poisson's ratio for pseudo-solid
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;
 
} // end of ProblemParameters


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//======start_of_ExactSolution========================================
/// Namespace for exact solution 
//====================================================================
namespace ExactSolution
{

 /// Constant/initial temperature
 double U0=0.0;

 /// Growth rate for interface
 double Growth_rate=1.0;
 
 ///  Frequency of co-sinusoidal oscillation of incoming heat flux
 /// (to assess suppression of re-freezing). Set to zero for validation.
 double Omega_cos=0.0;

 /// Exact solution as a Vector
 void get_exact_u_for_unsteady_heat_validation(const double& t, 
                                               const Vector<double>& x, 
                                               Vector<double>& u)
 {
  double X=x[0];
  double Y=x[1];
  u[0]=U0+t*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1);
 }

 /// Exact solution as a scalar
 void get_exact_u_for_unsteady_heat_validation(const double& t, 
                                               const Vector<double>& x, 
                                               double& u)
 {
  Vector<double> u_vect(1);
  get_exact_u_for_unsteady_heat_validation(t,x,u_vect);
  u=u_vect[0];
 }

 /// Source function to make it an exact solution 
 void get_source_for_unsteady_heat_validation(const double& t, 
                                              const Vector<double>& x, 
                                              double& source)
 {
  double X=x[0];
  double Y=x[1];
  source = -2.0*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)-2.0*t*t*t*Y*Y*
Growth_rate*(1.0-cos(2.0*X*0.3141592653589793E1))*cos(2.0*X*
0.3141592653589793E1)+4.0*t*t*t*t*Y*Y*Growth_rate*pow(cos(2.0*X*
0.3141592653589793E1),2.0)*0.3141592653589793E1*0.3141592653589793E1-8.0*t*t*t*
t*Y*Y*Growth_rate*pow(sin(2.0*X*0.3141592653589793E1),2.0)*0.3141592653589793E1
*0.3141592653589793E1-4.0*t*t*Y*Y*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)*0.3141592653589793E1*
0.3141592653589793E1+2.0*t*t*(Y-1.0+Growth_rate*t*t*(1.0-cos(2.0*X*
0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)+4.0*t*t*Y*cos(2.0*X*
0.3141592653589793E1);

 }


 ///  Flux required by the exact solution on a boundary with outer unit
 /// normal n. No dependence on temperature u.
 void prescribed_flux_for_unsteady_heat_validation(const double& t,
                                                   const Vector<double>& x, 
                                                   const Vector<double>& n, 
                                                   const double& u,
                                                   double& flux)
 {
  double X=x[0];
  double Y=x[1];

  //The outer unit normal 
  double Nx =  n[0];
  double Ny =  n[1];

  //The flux in terms of the normal is
  flux=(2.0*t*t*t*t*Y*Y*Growth_rate*sin(2.0*X*0.3141592653589793E1)*
0.3141592653589793E1*cos(2.0*X*0.3141592653589793E1)-2.0*t*t*Y*Y*(Y-1.0+
Growth_rate*t*t*(1.0-cos(2.0*X*0.3141592653589793E1)))*sin(2.0*X*
0.3141592653589793E1)*0.3141592653589793E1)*Nx+(2.0*t*t*Y*(Y-1.0+Growth_rate*t*
t*(1.0-cos(2.0*X*0.3141592653589793E1)))*cos(2.0*X*0.3141592653589793E1)+t*t*Y*
Y*cos(2.0*X*0.3141592653589793E1))*Ny;

  double melt_flux=
2.0*Growth_rate*t*(1.0-cos(2.0*X*0.3141592653589793E1))/sqrt(1.0+4.0
*Growth_rate*Growth_rate*t*t*t*t*pow(sin(2.0*X*0.3141592653589793E1),2.0)*
0.3141592653589793E1*0.3141592653589793E1);

  flux+=melt_flux*cos(Omega_cos*t);
 }

 /// Height of melting surface
 double melting_surface_height(const double& t, const double& x)
 {
  return 1.0-Growth_rate*t*t*(1.0-cos(2.0*x*0.3141592653589793E1));
 }
 
} // end of ExactSolution

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// UnsteadyHeat problem 
//====================================================================
template<class ELEMENT>
class UnsteadyHeatMeltProblem : public Problem
{

public:

 /// Constructor
 UnsteadyHeatMeltProblem();
 
 /// Destructor (empty)
 ~UnsteadyHeatMeltProblem(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 ///  Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
  
 /// Actions before adapt: wipe flux elements
 void actions_before_adapt() 
  {
   // Make backup of surface mesh
   Backed_up_surface_melt_mesh_pt=new BackupMeshForProjection<TElement<1,3> >
    (Surface_melt_mesh_pt,Melt_boundary_id);
      
   // Kill the  elements and wipe surface mesh
   delete_flux_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }

 ///  Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   // Create flux elements
   create_flux_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();

   // Rebuild elements
   complete_problem_setup();

   // Now project from backup of original contact mesh to new one
   Backed_up_surface_melt_mesh_pt->project_onto_new_mesh(
    Surface_melt_mesh_pt);
   
   // Wipe Lagrange multiplier pressure 
   unsigned n_element=Surface_melt_mesh_pt->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     SurfaceMeltElement<ELEMENT>* el_pt = 
      dynamic_cast<SurfaceMeltElement<ELEMENT>*>(
       Surface_melt_mesh_pt->element_pt(e));
     el_pt->set_lagrange_multiplier_pressure_to_zero();
    } 

   // Kill backed up mesh
   delete Backed_up_surface_melt_mesh_pt;
   Backed_up_surface_melt_mesh_pt=0; // hierher add to constructor

  }
 
 /// Doc the solution
 void doc_solution();
 
 /// Dummy global error norm for adaptive time-stepping
 double global_temporal_error_norm(){return 0.0;}

private:

 
 /// Create flux elements
 void create_flux_elements()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Melt_boundary_id; 
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding prescribed-flux element
     SurfaceMeltElement<ELEMENT>* flux_element_pt = new 
      SurfaceMeltElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the prescribed-flux element to the surface mesh
     Surface_melt_mesh_pt->add_element_pt(flux_element_pt);
     
    } //end of loop over bulk elements adjacent to boundary b
  }



 /// Delete flux elements
 void delete_flux_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_melt_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_melt_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_melt_mesh_pt->flush_element_and_node_storage();
  }


 ///  Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup()
  {

   // Set (pseudo-)solid mechanics properties for all elements
   //---------------------------------------------------------
   unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a solid element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     
     // Set the constitutive law
     el_pt->constitutive_law_pt() =
      ProblemParameters::Constitutive_law_pt;
     
     // Set density to zero
     el_pt->lambda_sq_pt()=&ProblemParameters::Lambda_sq;
     
     // Disable inertia
     el_pt->disable_inertia();

     // Set source function
     el_pt->source_fct_pt() = 
      &ExactSolution::get_source_for_unsteady_heat_validation;
    }

   // Apply boundary conditions for solid
   //------------------------------------

   // Bottom: completely pinned
   unsigned b=Bottom_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
    }

   // Sides: Symmetry bcs
   b=Left_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }
   b=Right_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }
   
   // No melting: Pin all solid positions and the Lagrange multipliers
   if (CommandLineArgs::command_line_flag_has_been_set("--disable_melting"))
    {
     unsigned nnod=Bulk_mesh_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       SolidNode* nod_pt=Bulk_mesh_pt->node_pt(j);
       nod_pt->pin_position(0);
       nod_pt->pin_position(1);
       unsigned nval=nod_pt->nvalue();
       if (nval==3)
        {
         nod_pt->pin(1);
        }
      }
    }

   // Assign the Lagrangian coordinates -- sensible
   // because we've completely rebuilt the mesh 
   // and haven't copied across any Lagrange multipliers
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
   n_element=Surface_melt_mesh_pt->nelement();
   for (unsigned e=0;e<n_element;e++)
    {
     SurfaceMeltElement<ELEMENT>* el_pt = 
      dynamic_cast<SurfaceMeltElement<ELEMENT>*>(
       Surface_melt_mesh_pt->element_pt(e));
     el_pt->set_lagrange_multiplier_pressure_to_zero();
    } 

   // Loop over the flux elements to pass pointer to prescribed flux function
   //------------------------------------------------------------------------
   n_element=Surface_melt_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to SurfaceMeltElement flux element
     SurfaceMeltElement<ELEMENT> *el_pt = 
      dynamic_cast<SurfaceMeltElement<ELEMENT>*>(
       Surface_melt_mesh_pt->element_pt(e));
     
     // Pass validation flux
     el_pt->flux_fct_pt()=
      &ExactSolution::prescribed_flux_for_unsteady_heat_validation;

     // Set melt temperature
     el_pt->melt_temperature_pt()=&ProblemParameters::Melt_temperature;

     // Suppress melting?
     if (CommandLineArgs::command_line_flag_has_been_set("--disable_melt_flux"))
      {
       el_pt->disable_melting();
      }
    }
  }
 
 /// Pointers to specific mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_melt_mesh_pt;

 /// ID of melt boundary
 unsigned Melt_boundary_id;

 /// ID of bottom boundary
 unsigned Bottom_boundary_id;

 /// ID of left boundary
 unsigned Left_boundary_id;

 /// ID of right boundary
 unsigned Right_boundary_id;

 /// Trace file
 ofstream Trace_file;

 // Setup labels for output
 DocInfo Doc_info;

 ///  Backup of Surface_melt_mesh_pt so the Lagrange multipliers
 /// and melt rate can be projected across
 BackupMeshForProjection<TElement<1,3> >*  Backed_up_surface_melt_mesh_pt;
 

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT>
UnsteadyHeatMeltProblem<ELEMENT>::UnsteadyHeatMeltProblem()
{ 

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0;

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Allow for crap initial guess
 Problem::Max_residuals=10000.0;
 Problem::Max_newton_iterations=10000;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon
  
 // The boundary is bounded by five distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(4);
 
 // Vertex coordinates on boundary
 Vector<Vector<double> > bound_coords(2);
 
 // Left boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=1.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=0.0;
 bound_coords[1][1]=0.0;

 // Build the boundary polyline
 Left_boundary_id=0;
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,
                                                  Left_boundary_id);

 // Bottom boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=0.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=1.0;
 bound_coords[1][1]=0.0;

 // Build the boundary polyline
 Bottom_boundary_id=1;
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,
                                                  Bottom_boundary_id);
 

 // Right boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=1.0;
 bound_coords[0][1]=0.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=1.0;
 bound_coords[1][1]=1.0;
 
 // Build the boundary polyline
 Right_boundary_id=2;
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_coords,
                                                  Right_boundary_id);

 // Flux boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=1.0;
 bound_coords[0][1]=1.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=0.0;
 bound_coords[1][1]=1.0;
 
 // Build boundary poly line
 Melt_boundary_id=3;
 boundary_polyline_pt[3]=new TriangleMeshPolyLine(bound_coords,
                                                  Melt_boundary_id);
 
 // Create the triangle mesh polygon for outer boundary
 //----------------------------------------------------
 TriangleMeshPolygon *outer_polygon =
  new TriangleMeshPolygon(boundary_polyline_pt);
  
 // Set the pointer
 closed_curve_pt = outer_polygon;
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the maximum area element
 double uniform_element_area=0.002;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 Bulk_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                       time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set element size limits
 Bulk_mesh_pt->max_element_size()=0.2;
 Bulk_mesh_pt->min_element_size()=0.0002; // hierher 0.002; 

 
 // Create the surface mesh as an empty mesh
 Surface_melt_mesh_pt=new Mesh;
 
 // Build 'em 
 create_flux_elements();

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_melt_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // Set the initial conditions
 unsigned nnod = Bulk_mesh_pt->nnode();
 for(unsigned j=0;j<nnod;j++)
  {
   Bulk_mesh_pt->node_pt(j)->set_value(0,ExactSolution::U0); 
  } 

 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnsteadyHeatMeltProblem<ELEMENT>::doc_solution()
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
 
 cout << std::endl;
 cout << "=================================================" << std::endl;
 cout << "Docing solution for t=" << time_pt()->time() << std::endl;
 cout << "=================================================" << std::endl;
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();
 
 // Output flux with melt
 sprintf(filename,"%s/flux_with_melt%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_melt_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<SurfaceMeltElement<ELEMENT>* >(
    Surface_melt_mesh_pt->element_pt(e))->output_melt(some_file);
  }
 some_file.close();
 
 
 // Output Number of Newton iterations in form that can be visualised
 // as vector in paraview
 sprintf(filename,"%s/newton_iter%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 some_file << "0 0 0 " << Nnewton_iter_taken << std::endl;
 some_file.close();
 
 
 // Output exact solution
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(
  some_file,npts,time_pt()->time(),
  ExactSolution::get_exact_u_for_unsteady_heat_validation); 
 some_file.close();
 
 // Output exact position of melting line
 sprintf(filename,"%s/exact_height%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nplot=100;
 for (unsigned j=0;j<nplot;j++)
  {
   double x=double(j)/double(nplot-1);
   some_file << x << " " 
             << ExactSolution::melting_surface_height(time_pt()->time(),x)
             << std::endl;
  }
 some_file.close();

 
 // Write norm of solution to trace file
 double norm=0.0;
 Bulk_mesh_pt->compute_norm(norm); 
 Trace_file  << norm << std::endl;
 
 //Increment counter for solutions 
 Doc_info.number()++;

} // end of doc_solution



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
///  Driver code for unsteady heat equation
//========================================================================
int main(int argc, char* argv[])
{

 FiniteElement::Accept_negative_jacobian=true;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--no_adapt");

 // Suppress melting?
 CommandLineArgs::specify_command_line_flag("--disable_melting");

 // Suppress melt flux?
 CommandLineArgs::specify_command_line_flag("--disable_melt_flux");
  
 // Non-default initial temperature
 CommandLineArgs::specify_command_line_flag("--theta_init",
                                            &ExactSolution::U0);

 //  Frequency of co-sinusoidal oscillation of incoming heat flux
 // (to assess suppression of re-freezing).
 CommandLineArgs::specify_command_line_flag("--omega_cos",
                                            &ExactSolution::Omega_cos);
  
 // Validation
 CommandLineArgs::specify_command_line_flag("--validate");
  
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 ProblemParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&ProblemParameters::Nu);
  
 // Build problem
 UnsteadyHeatMeltProblem<
  ProjectableUnsteadyHeatElement<
  PseudoSolidNodeUpdateElement<TUnsteadyHeatElement<2,3>,TPVDElement<2,3> > > >
 problem;

 
 // Choose simulation interval and timestep
 double t_max=0.5;
 double dt=0.01;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);
 
 // Set impulsive IC
 problem.assign_initial_values_impulsive(dt);
 
 //Output initial condition
 problem.doc_solution();
 
 unsigned max_adapt=1;
 bool first=false;
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }

 // Find number of steps
 unsigned nstep = unsigned(t_max/dt);

 // Validation?
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   nstep=5;
  }

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   cout << " Timestep " << istep << std::endl;
   
   // Dummy double adaptivity (timestep is always accepted because
   // tolerance is set to huge value; mainly used to automatically
   // re-solve with smaller timestep increment after non-convergence
   double epsilon_t=DBL_MAX;
   unsigned suppress_resolve_after_spatial_adapt_flag=0;
   double next_dt=
    problem.doubly_adaptive_unsteady_newton_solve(
     dt,
     epsilon_t,
     max_adapt,
     suppress_resolve_after_spatial_adapt_flag,
     first);
   oomph_info << "Suggested next dt: " << next_dt << std::endl;
   dt = next_dt; 

   //Output solution
   problem.doc_solution();

  }
 
} // end of main
