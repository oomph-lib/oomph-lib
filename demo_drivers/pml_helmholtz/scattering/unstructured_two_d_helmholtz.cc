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
// Driver for a specific 2D Helmholtz problem with
// perfectly matched layer treatment for the exterior boundaries

#include<fenv.h>

#include "math.h"
#include <complex>

// Generic routines
#include "generic.h"

// The pml Helmholtz equations
#include "pml_helmholtz.h"

// The meshes needed in the PML constructions
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"

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
 double Wavenumber = sqrt(50.0);

 /// Square of the wavenumber (also known as k^2)
 double K_squared = Wavenumber * Wavenumber;

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
 ~PMLProblem(){}

 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 /// Create PML meshes
 void create_pml_meshes();

 // Apply boundary conditions
 void apply_boundary_conditions();

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
 double a=0.2;
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
 triangle_mesh_parameters.element_area() = 0.1;

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

 // Create the main triangular mesh
 add_sub_mesh(Bulk_mesh_pt);

 // Create PML meshes and add them to the global mesh
 create_pml_meshes();

 // Build the entire mesh from its submeshes
 build_global_mesh();

 // Let's have a look where the boundaries are
 this->mesh_pt()->output("global_mesh.dat");
 this->mesh_pt()->output_boundaries("global_mesh_boundary.dat");

 // Complete the build of all elements so they are fully functional
 unsigned n_element = this->mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to Helmholtz bulk element
   PMLHelmholtzEquations<2> *el_pt =
    dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

   //Set the k_squared double pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Apply boundary conditions
 apply_boundary_conditions();

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

 // Build PML meshes  and add them to the global mesh
 create_pml_meshes();

 // Build the entire mesh from its submeshes
 rebuild_global_mesh();

 // Complete the build of all elements so they are fully functional

 // Loop over the entire mesh elements to set up element-specific
 // things that cannot be handled by constructor
 unsigned n_element = this->mesh_pt()->nelement();

 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to PMLHelmholtz bulk element
   PMLHelmholtzEquations<2> *el_pt =
    dynamic_cast<PMLHelmholtzEquations<2>*>(mesh_pt()->element_pt(e));

   //Set the frequency function pointer
   el_pt->k_squared_pt() = &GlobalParameters::K_squared;
  }

 // Re-apply boundary conditions
 apply_boundary_conditions();

}// end of actions_after_adapt

#endif

//==================start_of_apply_boundary_conditions====================
/// Apply boundary conditions
//========================================================================
template<class ELEMENT>
void PMLProblem<ELEMENT>::apply_boundary_conditions()
{

 // Boundary conditions are set on the surface of the circle
 // as a constant nonzero Dirichlet boundary condition
 unsigned n_bound = Bulk_mesh_pt->nboundary();

 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = Bulk_mesh_pt->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     if ((b==0) || (b==1))
      {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,n);
       nod_pt->pin(0);
       nod_pt->pin(1);

       nod_pt->set_value(0,0.1);
       nod_pt->set_value(1,0.0);
      }
    }
  }

}// end of apply_boundary_conditions


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

 // Output solution
 //-----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output coarse solution
 //-----------------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned npts_coarse=2;
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();


 // Output solution within pml domains
 //-----------------------------------
 sprintf(filename,"%s/pml_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
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




 // Write norm of solution to trace file
 double norm=0.0;
 Bulk_mesh_pt->compute_norm(norm);
 Trace_file  << norm << std::endl;

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

 // PML width in elements for the left layer
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

 // Set up the problem with 2D ten-node elements from the
 // TPMLHelmholtzElement family.
 // PMLProblem<ProjectablePMLHelmholtzElement
 //  <TPMLHelmholtzElement<2,4> > > problem;

#else

 // Set up the problem with 2D six-node elements from the
 // TPMLHelmholtzElement family.
 PMLProblem<TPMLHelmholtzElement<2,3> >  problem;

 // Set up the problem with 2D ten-node elements from the
 // TPMLHelmholtzElement family.
 //   PMLProblem<TPMLHelmholtzElement<2,4> >  problem;

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
