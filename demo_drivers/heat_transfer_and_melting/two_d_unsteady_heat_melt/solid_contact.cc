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
// Driver for 2D contact problem -- displacement controlled penetrator


//===========================================================
// Jonathan to do:
// - doc projection of contact pressure (and do by default)
// - explore --resolve_after_adapt
//===========================================================

#include <fenv.h> 

//Generic routines
#include "generic.h"

// The solid elements
#include "solid.h"

// Mesh
#include "meshes/triangle_mesh.h"

// Contact stuff
#include "contact_elements.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 

//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{
 
 /// Non-dim density for solid
 double Lambda_sq=0.0;

 /// Poisson's ratio for solid
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Radius of penetrator
 double Radius=0.1;

 /// Initial/max y-position
 double Y_c_max=0.0;

 /// Position of centre of penetrator
 Vector<double> Centre;

 /// Penetrator
 Penetrator* Penetrator_pt=0;

 /// Initial/max element area
 double El_area=0.002;

} // end of ProblemParameters


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// Problem class
//====================================================================
template<class ELEMENT>
class ContactProblem : public Problem
{

public:

 /// Constructor
 ContactProblem();
 
 /// Destructor (empty)
 ~ContactProblem(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
  

 /// Actions before next timestep
 void actions_before_implicit_timestep()
  {
   // Amplitude of oscillation
   double amplitude=0.1; 
   
   // Update position of centre -- amplitude of oscillation increases
   // by 5% per period...
   double time=time_pt()->time();
   ProblemParameters::Centre[1]=ProblemParameters::Y_c_max-amplitude*
    (1.0+0.05*time)*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*time));
   
   oomph_info << "Solving for y_c = " 
              << ProblemParameters::Centre[1] 
              << " for time: " << time << std::endl;
  }

 /// Actions before adapt: wipe contact elements
 void actions_before_adapt() 
  {
   // Make backup of surface mesh
   Backed_up_surface_contact_mesh_pt=
    new BackupMeshForProjection<TElement<1,3> >(
     Surface_contact_mesh_pt,Contact_boundary_id);

   // // Output contact elements
   // ofstream some_file;
   // char filename[100];
   // sprintf(filename,"contact_before.dat");
   // some_file.open(filename);
   // unsigned nel=Surface_contact_mesh_pt->nelement();
   // for (unsigned e=0;e<nel;e++)
   //  {
   //   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
   //    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
   //  }
   // some_file.close();
   
   // // Kill the  elements and wipe surface mesh
   delete_contact_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   // Create contact elements
   create_contact_elements();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
   
   // Rebuild elements
   complete_problem_setup();

   // Now project from backup of original contact mesh to new one
   oomph_info << "Projecting contact pressure.\n";
   Backed_up_surface_contact_mesh_pt->project_onto_new_mesh(
    Surface_contact_mesh_pt);

   // Kill backed up mesh
   delete Backed_up_surface_contact_mesh_pt;
   Backed_up_surface_contact_mesh_pt=0;

   // // Output contact elements
   // ofstream some_file;
   // char filename[100];
   // sprintf(filename,"contact_after.dat");
   // some_file.open(filename);
   // unsigned nel=Surface_contact_mesh_pt->nelement();
   // for (unsigned e=0;e<nel;e++)
   //  {
   //   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
   //    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
   //  }
   // some_file.close();
   // //pause("done");

  }
 
 /// Doc the solution
 void doc_solution();
 
 /// Dummy global error norm for adaptive time-stepping
 double global_temporal_error_norm(){return 0.0;}

private:

 
 /// Create contact elements
 void create_contact_elements()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Contact_boundary_id; 
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding contact element
     NonlinearSurfaceContactElement<ELEMENT>* contact_element_pt = new 
      NonlinearSurfaceContactElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the contact element to the surface mesh
     Surface_contact_mesh_pt->add_element_pt(contact_element_pt);
     
    } //end of loop over bulk elements adjacent to boundary b
  }



 /// Delete contact elements
 void delete_contact_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_contact_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_contact_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_contact_mesh_pt->flush_element_and_node_storage();
  }


 /// Helper function to (re-)set boundary condition
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


   // hierher
   // if (!CommandLineArgs::command_line_flag_has_been_set("--proper_elasticity"))
   //  {
   //   // Assign the Lagrangian coordinates -- sensible
   //   // because we've completely rebuilt the mesh 
   //   // and haven't copied across any Lagrange multipliers
   //   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
   //  }


   // Loop over the contact elements to pass pointer to penetrator
   //-------------------------------------------------------------
   n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     NonlinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
    }
  }
 
 /// Pointer to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_contact_mesh_pt;
 
 /// ID of contact boundary
 unsigned Contact_boundary_id;

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

 /// Backup of Surface_contact_mesh_pt so the Lagrange multipliers
 /// can be projected across
 BackupMeshForProjection<TElement<1,3> >* Backed_up_surface_contact_mesh_pt;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for contact problem in square domain
//========================================================================
template<class ELEMENT>
ContactProblem<ELEMENT>::ContactProblem()
{ 

 // Initialise
 Backed_up_surface_contact_mesh_pt=0;

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0;

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Allow for crap initial guess
 Problem::Max_residuals=10000.0;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);


 double x_ll=0.0;
 double x_ur=1.0;
 double y_ll=0.0;
 double y_ur=1.0;

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
 bound_coords[0][0]=x_ll;
 bound_coords[0][1]=y_ur;

 bound_coords[1].resize(2);
 bound_coords[1][0]=x_ll;
 bound_coords[1][1]=y_ll;

 // Build the boundary polyline
 Left_boundary_id=0;
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,
                                                  Left_boundary_id);

 // Bottom boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=x_ll;
 bound_coords[0][1]=y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=x_ur;
 bound_coords[1][1]=y_ll;

 // Build the boundary polyline
 Bottom_boundary_id=1;
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,
                                                  Bottom_boundary_id);
 

 // Right boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=x_ur;
 bound_coords[0][1]=y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=x_ur;
 bound_coords[1][1]=y_ur;
 
 // Build the boundary polyline
 Right_boundary_id=2;
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_coords,
                                                  Right_boundary_id);


 // Contact boundary
 unsigned npt_contact=4; 
 Vector<Vector<double> > contact_bound_coords(npt_contact);
 contact_bound_coords[0].resize(2);
 contact_bound_coords[0][0]=x_ur;
 contact_bound_coords[0][1]=y_ur;
 for (unsigned j=1;j<npt_contact-1;j++)
  {  
   contact_bound_coords[j].resize(2);
   contact_bound_coords[j][0]=x_ur-(x_ur-x_ll)*double(j)/double(npt_contact-1);
   contact_bound_coords[j][1]=y_ur;
  }
 contact_bound_coords[npt_contact-1].resize(2);
 contact_bound_coords[npt_contact-1][0]=x_ll;
 contact_bound_coords[npt_contact-1][1]=y_ur;

 
 // Build boundary poly line
 Contact_boundary_id=3;
 TriangleMeshPolyLine*contact_boundary_pt=
  new TriangleMeshPolyLine(contact_bound_coords,
                           Contact_boundary_id);
 boundary_polyline_pt[3]=contact_boundary_pt;
 

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
 double uniform_element_area=ProblemParameters::El_area;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 

 // Create the mesh
 Bulk_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                       time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Create the surface mesh as an empty mesh
 Surface_contact_mesh_pt=new Mesh;
 
 // Build 'em 
 create_contact_elements();
 
 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_contact_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();
 
 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ContactProblem<ELEMENT>::doc_solution()
{ 

 oomph_info << "Outputting for step: " << Doc_info.number() << std::endl;
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
 
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
 
 
 // Output contact elements
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
    Surface_contact_mesh_pt->element_pt(e))->output(some_file,20);
  }
 some_file.close();
 

 // Output penetrator
 sprintf(filename,"%s/penetrator%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nplot=500;
 ProblemParameters::Penetrator_pt->output(some_file,nplot);
 some_file.close();
  
 // Write mesh "volume" to trace file
 Trace_file  << time_pt()->time() << " " << Bulk_mesh_pt->total_size() 
             << std::endl;
 
 //Increment counter for solutions 
 Doc_info.number()++;

} // end of doc_solution



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// Driver code
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
    
 // Initial element size
 CommandLineArgs::specify_command_line_flag("--el_area",
                                            &ProblemParameters::El_area);
    
 // Resolve after adaptation
 // hierher CommandLineArgs::specify_command_line_flag("--resolve_after_adapt");

 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--validate");
    
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 ProblemParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&ProblemParameters::Nu);
  
 // Define centre of penetrator
 ProblemParameters::Centre.resize(2);
 ProblemParameters::Centre[0]=0.5;
 ProblemParameters::Centre[1]=1.15;

 // Initial/max y-position: Lift off the surface by a bit...
 ProblemParameters::Y_c_max=1.0+ProblemParameters::Radius+0.01;
   
 // Create penetrator
 ProblemParameters::Penetrator_pt =
  new CircularPenetrator(&ProblemParameters::Centre,
                         ProblemParameters::Radius);
  
 // Build problem
 ContactProblem<ProjectablePVDElement<TPVDElement<2,3> > > problem;

  //Output initial condition
 problem.doc_solution();
 
 unsigned max_adapt=1;
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }


 // Number of parameter increments per period
 unsigned nstep_for_period=40; // 100;

 // Parameter variation
 unsigned nperiod=3;

 // Initial timestep
 double dt=1.0/double(nstep_for_period);

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);

 
 double t_max=double(nperiod);
 if (CommandLineArgs::command_line_flag_has_been_set("--validate"))
  {
   t_max=10.0*dt;
  }
 

 //while (ProblemParameters::Centre[1]>1.08)
 while (problem.time_pt()->time()<t_max)
  {
   // Dummy double adaptivity (timestep is always accepted because
   // tolerance is set to huge value; mainly used to automatically
   // re-solve with smaller timestep increment after non-convergence
   double epsilon_t=DBL_MAX;
   bool first=false;
   unsigned suppress_resolve_after_spatial_adapt_flag=1;
   if (CommandLineArgs::command_line_flag_has_been_set("--resolve_after_adapt"))
    {
     suppress_resolve_after_spatial_adapt_flag=0;
    }

   // hierher
   suppress_resolve_after_spatial_adapt_flag=0;
  
   double next_dt=
    problem.doubly_adaptive_unsteady_newton_solve(
     dt,
     epsilon_t,
     max_adapt,
     suppress_resolve_after_spatial_adapt_flag,
     first);
   dt = next_dt; 
   
   //Output solution
   problem.doc_solution();
  }
  

} // end of main
