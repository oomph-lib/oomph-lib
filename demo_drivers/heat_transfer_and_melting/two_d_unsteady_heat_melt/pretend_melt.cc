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

// Mesh
#include "meshes/triangle_mesh.h"

// The unsteady heat equations
#include "unsteady_heat_flux_pseudo_melt_elements.h"

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
 
 /// Get flux applied along boundary x=0.
 void flux(const double& time, const Vector<double>& x, double& flux)
 {
  flux = 20.0*sin(2.0*4.0*MathematicalConstants::Pi*time);
 }
 
 /// Melt-temperature
 double Melt_temperature=0.0;
 
 
} // end of ProblemParameters

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// UnsteadyHeat problem 
//====================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public Problem
{

public:

 /// Constructor
 UnsteadyHeatProblem();
 
 /// Destructor (empty)
 ~UnsteadyHeatProblem(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 ///  Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
  
 /// Actions before adapt: wipe flux elements
 void actions_before_adapt() 
  {
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
  }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
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
     UnsteadyHeatFluxPseudoMeltElement<ELEMENT>* flux_element_pt = new 
      UnsteadyHeatFluxPseudoMeltElement<ELEMENT>(bulk_elem_pt,face_index);
     
     //Add the prescribed-flux element to the surface mesh
     Surface_melt_mesh_pt->add_element_pt(flux_element_pt);
     
    } //end of loop over bulk elements adjacent to boundary b

   
   // How many bulk elements are adjacent to boundary b?
   b=Flux_boundary_id; 
   n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     if (CommandLineArgs::command_line_flag_has_been_set("--melt_everywhere"))
      {
       // Build the corresponding prescribed-flux element
       UnsteadyHeatFluxPseudoMeltElement<ELEMENT>* flux_element_pt = new 
        UnsteadyHeatFluxPseudoMeltElement<ELEMENT>(bulk_elem_pt,face_index);
       
       //Add the prescribed-flux element to the surface mesh
       Surface_mesh_pt->add_element_pt(flux_element_pt);
      }
     else
      {
       // Build the corresponding prescribed-flux element
       UnsteadyHeatFluxElement<ELEMENT>* flux_element_pt = new 
        UnsteadyHeatFluxElement<ELEMENT>(bulk_elem_pt,face_index);
       
       //Add the prescribed-flux element to the surface mesh
       Surface_mesh_pt->add_element_pt(flux_element_pt);
      }

      

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

   // How many surface elements are in the surface mesh
   n_element = Surface_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_mesh_pt->flush_element_and_node_storage();
  }


 ///  Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup()
  {
   // Loop over the flux elements to pass pointer to prescribed flux function
   unsigned n_element=Surface_melt_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to UnsteadyHeat flux element
     UnsteadyHeatFluxPseudoMeltElement<ELEMENT> *el_pt = 
      dynamic_cast<UnsteadyHeatFluxPseudoMeltElement<ELEMENT>*>(
       Surface_melt_mesh_pt->element_pt(e));
     
     // Set the pointer to the prescribed flux function
     el_pt->flux_fct_pt() = &ProblemParameters::flux;
          
     // Set melt temperature
     el_pt->melt_temperature_pt()=&ProblemParameters::Melt_temperature;

     // Suppress melting?
     if (CommandLineArgs::command_line_flag_has_been_set("--disable_melting"))
      {
       el_pt->disable_melting();
      }
    }
   
   // Loop over the flux elements to pass pointer to prescribed flux function
   n_element=Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     if (CommandLineArgs::command_line_flag_has_been_set("--melt_everywhere"))
      {
       // Upcast from GeneralisedElement to UnsteadyHeat flux element
       UnsteadyHeatFluxPseudoMeltElement<ELEMENT> *el_pt = 
        dynamic_cast<UnsteadyHeatFluxPseudoMeltElement<ELEMENT>*>(
         Surface_mesh_pt->element_pt(e));
       
       // Set the pointer to the prescribed flux function
       el_pt->flux_fct_pt() = &ProblemParameters::flux;
       
       // Set melt temperature
       el_pt->melt_temperature_pt()=&ProblemParameters::Melt_temperature;

       // Suppress melting?
       if (CommandLineArgs::command_line_flag_has_been_set("--disable_melting"))
        {
         el_pt->disable_melting();
        }
      }
     else
      {
       // Upcast from GeneralisedElement to UnsteadyHeat flux element
       UnsteadyHeatFluxElement<ELEMENT> *el_pt = 
        dynamic_cast<UnsteadyHeatFluxElement<ELEMENT>*>(
         Surface_mesh_pt->element_pt(e));
       
       // Set the pointer to the prescribed flux function
       el_pt->flux_fct_pt() = &ProblemParameters::flux;
      }
    }
   
  }
 
 /// Pointers to specific mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_melt_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

 /// ID of flux boundary
 unsigned Flux_boundary_id;

 /// ID of melt boundary
 unsigned Melt_boundary_id;

 /// Trace file
 ofstream Trace_file;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem()
{ 

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
  
 // The boundary is bounded by three distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(3);
 
 // Vertex coordinates on boundary
 Vector<Vector<double> > bound_coords(4);
 
 // First part of the boundary 
 //---------------------------
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=0.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=1.0;
 bound_coords[1][1]=0.0;

 bound_coords[2].resize(2);
 bound_coords[2][0]=1.0;
 bound_coords[2][1]=1.0;

 bound_coords[3].resize(2);
 bound_coords[3][0]=0.0;
 bound_coords[3][1]=1.0;
 
 // Build the 1st boundary polyline
 unsigned boundary_id=0;
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,boundary_id);
 
 // Second part of the boundary
 //----------------------------
 bound_coords.resize(2);
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=1.0;

 bound_coords[1].resize(2);
 bound_coords[1][0]=0.0;
 bound_coords[1][1]=0.5;
 

 // Flux boundary
 Flux_boundary_id=1;

 // Build the 2nd boundary polyline
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,
                                                  Flux_boundary_id);
 

 // Third part of the boundary
 //----------------------------
 bound_coords.resize(2);
 bound_coords[0].resize(2);
 bound_coords[0][0]=0.0;
 bound_coords[0][1]=0.5;

 bound_coords[1].resize(2);
 bound_coords[1][0]=0.0;
 bound_coords[1][1]=0.0;
 

 // Flux boundary
 Melt_boundary_id=2;

 // Build the 3rd boundary polyline
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_coords,
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
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                                  time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set element size limits
 Bulk_mesh_pt->max_element_size()=0.2;
 Bulk_mesh_pt->min_element_size()=0.000002; 
 
 // Create the surface mesh as an empty mesh
 Surface_melt_mesh_pt=new Mesh;
 Surface_mesh_pt=new Mesh;

 // Build 'em 
 create_flux_elements();

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_melt_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // Set the initial conditions
 unsigned nnod = Bulk_mesh_pt->nnode();
 for(unsigned j=0;j<nnod;j++)
  {
   Bulk_mesh_pt->node_pt(j)->set_value(0,-0.5); 
  } 
  
 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
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
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();

 // Output flux with melt
 sprintf(filename,"%s/flux_with_melt%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned nnod=dynamic_cast<UnsteadyHeatFluxPseudoMeltElement<ELEMENT>*>(
  Surface_melt_mesh_pt->element_pt(0))->nnode();
 Surface_melt_mesh_pt->output(some_file,nnod);
 if (CommandLineArgs::command_line_flag_has_been_set("--melt_everywhere"))
  {
   Surface_mesh_pt->output(some_file,nnod);
  }
 some_file.close();


 // Output Number of Newton iterations in form that can be visualised
 // as vector in paraview
 sprintf(filename,"%s/newton_iter%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 some_file << "0 0 0 " << Nnewton_iter_taken << std::endl;
 some_file.close();

 // Write norm of solution to trace file
 double norm=0.0;
 Bulk_mesh_pt->compute_norm(norm); 
 Trace_file  << norm << std::endl;
 
} // end of doc_solution



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
///  Driver code for unsteady heat equation
//========================================================================
int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--no_adapt");

 // Melt everywhere?
 CommandLineArgs::specify_command_line_flag("--melt_everywhere");

 // Suppress melting?
 CommandLineArgs::specify_command_line_flag("--disable_melting");
  
 // Validation
 CommandLineArgs::specify_command_line_flag("--validate");
  
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Build problem
 UnsteadyHeatProblem<ProjectableUnsteadyHeatElement
                     <TUnsteadyHeatElement<2,3> > > problem;
 
 // Setup labels for output
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory("RESLT");
 
 // Output number
 doc_info.number()=0;
 
 // Choose simulation interval and timestep
 double t_max=0.5;
 double dt=0.01;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 problem.initialise_dt(dt);
 
 // Set impulsive IC
 problem.assign_initial_values_impulsive(dt);
 
 //Output initial condition
 problem.doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

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
   
   // Take timestep
   problem.unsteady_newton_solve(dt,max_adapt,first);
   
   //Output solution
   problem.doc_solution(doc_info);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }
 

}; // end of main
